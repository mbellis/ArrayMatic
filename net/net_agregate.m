% AGREGATE - Find the clusters
% The current cluster is seeded with the  probe set which has the maximum
% number of correlation higher than a given correlation value.
%INPUT
% Corr: a square matrix with positive correlation values in the range 0-100
%       (correlation between two probe set)
% CorrSel: a logical square matrix, the same size of Corr which indicates which correlation must be used
% CorrInfLimit: if AgregationType='limit', probe sets with corrletations > CorrInfLimit are merged with
%               the current cluster
% CorrStdFactor: if AgregationType='variance', probe sets with corrletations > median - std*CorrStdFactor
%               are merged with the current cluster
% ConnStdFactor: probe sets with a mean connectivity < median - std*ConnStdFactor
%               are not merged to the current cluster
% AgregationType: the type of criterium used to put a probe set in an
%                 existing cluster {'variance','limit','complex'}
% MaxIndex: Index of the probe set that has, in corr, the maximum
%           number of correlation higher than a given correlation value
% CliqueFactor:
% LoadFlag: if = 0 , net_agregate works on a matrix, if = 1 agregate
% receive Corr and CorrSel as a vector and load the nedded columns of the
% matrix.
function [Cluster,varargout]=net_agregate(Corr,Anti,CorrSel,CorrInfLimit,ConnStdFactor,CorrStdFactor,MaxIndex,CliqueFactor,AgregationType,LoadFlag,CliqueFlag,varargin)
global K
%seed the cluster with the ps having the greatest number of correlation higher than CorrSupLimit (parameter MaxIndex)
%find the probe sets that have a link with this seeding ps with a correlation
%higher than CorrSupLimit and agregate them to the seed probe set(=>Merged)


if LoadFlag
    ModelRank=varargin{1};
    NetRank=varargin{2};
    PsNb=varargin{3};
end
if LoadFlag
    Merged=find(CorrSel==1);
else
    Merged=find(CorrSel(:,MaxIndex)==1);
end
MergedNb=length(Merged);
if MergedNb>0
    if MergedNb>1
    %start the counting of Ps in the current cluster
    %each Ps is counted in Cluster as many times it is linked to another Ps in the
    %current cluster (intra cluster node degree)
    Cluster=[ones(length(Merged),1)*MaxIndex;Merged];
    
        %test if the correlation is sufficient
        if isequal(AgregationType,'variance')
            if LoadFlag
                MedianCorr=median(Corr(Merged));
                StdCorr=std(double(Corr(Merged)));
            else
                MedianCorr=median(Corr(MaxIndex,Merged));
                StdCorr=std(double(Corr(MaxIndex,Merged)));
            end
        end
        %process each putative new member of the current cluster
        for i=1:length(Merged)-1;

            SearchedPos=Merged(i+1:end);
            %merge ps that have a correlation > CorrInfLimit
            %and add all the probesets that are correlated with it with a
            %corr higher

            %selection only on Corr
            if isequal(AgregationType,'corr limit')
                if LoadFlag
                    NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',NetRank));
                    Corr=load_data(sprintf('c_m%03u_n%05u.4mat',ModelRank,NetRank),NetDir,PsNb,PsNb,'uint8','ieee-le',1:PsNb,Merged(i));
                    Corr(Merged(i))=0;
                    CurrPos=find(Corr(SearchedPos)>CorrInfLimit);
                else
                    CurrPos=find(Corr(SearchedPos,Merged(i))>CorrInfLimit);
                end
            elseif isequal(AgregationType,'variance')
                if LoadFlag
                    NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',NetRank));
                    Corr=load_data(sprintf('c_m%03u_n%05u.4mat',ModelRank,NetRank),NetDir,PsNb,PsNb,'uint8','ieee-le',1:PsNb,Merged(i));
                    Corr(Merged(i))=0;
                    CurrPos=find(Corr(SearchedPos)>MedianCorr-StdCorr*CorrStdFactor);
                else
                    CurrPos=find(Corr(SearchedPos,Merged(i))>MedianCorr-StdCorr*CorrStdFactor);
                end
            elseif isequal(AgregationType,'complex')
                if LoadFlag
                    CurrPos=COMPLEX_SEL(Merged(i),LoadFlag,ModelRank,NetRank);
                else
                    CurrPos=COMPLEX_SEL(Merged(i),LoadFlag,Corr,Anti);
                end
            end
            if ~isempty(CurrPos)
                Cluster=[Cluster;ones(length(CurrPos),1)*Merged(i);SearchedPos(CurrPos)];
            end
        end
        %eliminate ps that have few correlation with others
        LinkNb=histc(Cluster,unique(Cluster));
        MedianLinkNb=median(LinkNb);
        StdLinkNb=std(LinkNb);
        Cluster=unique(Cluster);
        %at this step Cluster and LinkNb are in correspondance : linkNb(i) is
        %the number of link of Cluster(i)        
        if CliqueFlag
            %find the greatest clique
            CurrLimit=length(Cluster)-1;
            while min(LinkNb)<CurrLimit
                ClearIndex=find(LinkNb==min(LinkNb));
                Cluster(ClearIndex)=[];
                LinkNb(ClearIndex)=[];
                CurrLimit=length(Cluster)-1;
            end
        else
            CurrLimit=MedianLinkNb-ConnStdFactor*StdLinkNb;
            while min(LinkNb)<CurrLimit
                ClearIndex=find(LinkNb==min(LinkNb));
                Cluster(ClearIndex)=[];
                LinkNb(ClearIndex)=[];
                MedianLinkNb=median(LinkNb);
                StdLinkNb=std(LinkNb);
                CurrLimit=MedianLinkNb-ConnStdFactor*StdLinkNb;
            end
        end       
        %normalized deviation from the median
        if ~isempty(Cluster)
            if CliqueFlag==0
                if length(Cluster)>2
                    Deviation=(LinkNb-(length(Cluster)-1))/(length(Cluster)-1);
                else
                    Deviation=zeros(length(Cluster),1);
                end
            end
        else
            Deviation=[];
        end
        if nargout==3
            %find quasi-cliques
            Corr=Corr(Cluster,Cluster);
            CorrSel=uint8(zeros(size(Corr)));
            CorrSel(Corr>CorrInfLimit)=1;
            NodeDegree=sum(CorrSel);
            Clique=find(NodeDegree>min(ceil(length(Corr)*CliqueFactor),length(Corr)-1));
            varargout{1}=Clique;
            varargout{2}=Deviation;
        end
    else %MergeNb=1
        if LoadFlag
            if Corr(Merged(1))>CorrInfLimit
                Cluster=sort([MaxIndex;Merged(1)]);
                Deviation=0;
                if nargout==3
                    varargout{1}=1;
                    varargout{2}=Deviation;
                end
            else
                Cluster=[];
                Deviation=[];
                if nargout==3
                    varargout{1}=0;
                    varargout{2}=Deviation;
                end
            end
        else
            if Corr(Merged(1),MaxIndex)>CorrInfLimit
                Cluster=sort([MaxIndex;Merged(1)]);
                Deviation=0;
                if nargout==3
                    varargout{1}=1;
                    varargout{2}=Deviation;
                end
            else
                Cluster=[];
                Deviation=[];
                if nargout==3
                    varargout{1}=0;
                    varargout{2}=Deviation;
                end
            end
        end                
    end
else
    Cluster=MaxIndex;
    Deviation=0;
    if nargout==3
        varargout{1}=1;
    end
end


function CurrPos=COMPLEX_SEL(PsRank,LoadFlag,varargin)
CORR_LIMIT=0;
ANTI_LIMIT=20;
if LoadFlag
    ModelRank=varargin{1};
    NetRank=varargin{2};
    NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',NetRank));
    Corr=load_data(sprintf('c_m%03u_n%05u.4mat',ModelRank,NetRank),NetDir,PsNb,PsNb,'uint8','ieee-le',1:PsNb,PsRank);
    Corr(PsRank)=0;
    Anti=load_data(sprintf('a_m%03u_n%05u.4mat',ModelRank,NetRank),NetDir,PsNb,PsNb,'uint8','ieee-le',1:PsNb,PsRank);
    PsRank=1;
else
    Corr=varargin{1};
    Anti=varargin{2};
end

CurrPos=find(Corr(:,PsRank)>CORR_LIMIT&Anti(:,PsRank)<ANTI_LIMIT);