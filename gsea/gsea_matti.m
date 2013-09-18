h
% GSEA     Perform gene set enrichment analysis based on expression data
%
%    GSEA(A, B, GSETS) calculates enrichment scores for the gene sets GSETS
%    based on the differential expression between datasets A and B. Permutation
%    tests are performed by permuting the class labels A and B.
%
%    GSEA(A, [], GSETS) calculates enrichment scores using data in A without
%    the reference B.

% Author: Matti Annala <matti.annala@tut.fi>

function [] = gsea_matti(A, B, gsets, varargin)

if ~isempty(B) && size(A, 1) ~= size(B, 1)
	error 'Matrices for classes A and B must have an equal number of rows.';
end


if isnumeric(gsets)
    %only one gene set, containing probe set ranks
	genes = gsets;
	gsets = struct;
	gsets.Name = { 'foo' };
	gsets.genes = { genes };
end

gset_size = [1 Inf];
permutations = 10000;
significance = 0.001;
early_termination = 1e-9;
whitelist = true(length(gsets.Name), 1);

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'GsetSize')
		gset_size = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'Permutations')
		permutations = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'Significance')
		significance = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'Whitelist')
		if isnumeric(varargin{k+1})
			whitelist = false(length(gsets.Name), 1);
			whitelist(varargin{k+1}) = true;
		%elseif islogical(whitelist)
        elseif islogical(varargin{k+1})
			whitelist = varargin{k+1};
		else
			error 'Whitelist must be given as a logical or index vector.';
		end
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

if 1/permutations >= significance
	error(['Number of permutations is too low to achieve the desired ' ...
	       'significance level.']);
end

if isempty(B)
	ratio_label = 'A';
	diff = mean(A, 2);
	AB = A;
else
	ratio_label = 'A/B';
	diff = mean(A, 2) ./ mean(B, 2);
	AB = [A B];
end

bad = any(isnan(AB), 2);
AB = AB(~bad, :);

labels = [zeros(1, size(A, 2)) ones(1, size(B, 2))];

enriched = struct;
enriched.Name = {};
enriched.Genes = {};
enriched.Pval = [];
enriched.Escore = [];

%shufflings = nchoosek(size(AB, 2), sum(labels == 0));
%if shufflings < permutations * 10
%	fprintf(1, ['WARNING: Number of samples is too low for ' ...
%	            'sample label permutation. Permuting genes...\n']);
%end

progress = Progress;
whitelist = find(whitelist);

for k = 1:length(whitelist)
	idx = whitelist(k);
	progress.update(k / length(whitelist));
	
	g = gsets.genes(idx);
	if length(g) < gset_size(1), continue, end

	logG = false(size(A, 1), 1);
	logG(g) = true;
	g = find(logG(~bad));

	if length(g) < gset_size(1) || length(g) > gset_size(2)
		continue;
	end
	
	[escore, pval] = gsea_single(AB, labels, g, permutations, significance, ...
		early_termination);
		
	if pval <= significance
		enriched.Name{end+1} = gsets.Name{idx};
		enriched.Genes{end+1} = g;
		enriched.Pval(end+1) = pval;
		enriched.Escore(end+1) = escore;
	end
end

[temp, order] = sort(enriched.Pval,'ascend');
for k = 1:length(order)
	idx = order(k);
	fprintf(1, '%s (%d genes):\n', enriched.Name{idx}, ...
		length(enriched.Genes{idx}));
	fprintf(1, 'P-value = %.4f\t\tE-score = %.2f\t\t%s = %f\n', ...
		enriched.Pval(idx), enriched.Escore(idx), ratio_label, ...
		nanmean(diff(enriched.Genes{idx})));
end



function [escore, pval] = gsea_single(AB, labels, S, permutations, ...
	significance, early_termination)

L = rank_features(AB, labels);
escore = compute_enrichment(L, S);
perm_escores = zeros(1, permutations);

terminated_early = false;
perms_higher = 0;

if all(labels == 0)
	M = AB;    % Do this outside the loop to make things faster.
end

for p = 1:permutations
	if all(labels == 0)
		% For experiments with only one class, we permute the genes (rows)
		% within each individual sample.
		for s = 1:size(AB, 2)
			M(:, s) = AB(randperm(size(AB, 1)), s);
		end
		[temp, L] = sort(mean(M, 2), 'descend');
	else
		% For experiments with two classes, the permutation test is
		% done by permuting the class labels A and B.
		plabels = zeros(size(labels));
		perm = randperm(length(labels));
		plabels(perm(labels == 1)) = 1;
		L = rank_features(AB, plabels);
	end
	
	perm_escores(p) = compute_enrichment(L, S);
	
	% For the first iterations of the permutation test, we check if
	% we get a lot of permutation e-scores higher than the observed one.
	% If we do, then we can stop the permutation test early, because
	% very likely the observed e-score will not be significantly high.
	if p <= permutations / 10 && early_termination > 0 && ...
		((escore >= 0 && perm_escores(p) > escore) || ...
		(escore < 0 && perm_escores(p) < escore))
		perms_higher = perms_higher + 1;
		if 1 - binocdf(perms_higher, p, significance) < early_termination
			terminated_early = true;
			break;
		end
	end
end

if terminated_early
	pval = NaN;
else
	if escore >= 0
		pval = sum(perm_escores > escore) / length(perm_escores);
	else
		pval = sum(perm_escores < escore) / length(perm_escores);
	end
	pval = max(pval, 1 / permutations);
end



function L = rank_features(AB, labels)

diff = zeros(size(AB, 1), 1);
A_labels = find(labels == 0);
B_labels = find(labels == 1);
for k = 1:length(A_labels)
	diff = diff + AB(:, A_labels(k)) / length(A_labels);
end
for k = 1:length(B_labels)
	diff = diff - AB(:, B_labels(k)) / length(B_labels);
end
[temp, L] = sort(diff, 'descend');



function escore = compute_enrichment(L, S)

hit = ismember(L, S);
Nh = sum(hit);

% TODO: Steps in the random walk should be weighed by correlation.
Phit = cumsum(hit) / Nh;
Pmis = cumsum(~hit) / (length(L) - Nh);

[temp, idx] = max(abs(Phit - Pmis));
escore = Phit(idx) - Pmis(idx);