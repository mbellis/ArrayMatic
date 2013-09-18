%========================%
% FUNCTION KEGG2ENSEMBL  %
%========================%

% KEGG2ENSEMBL interrogates Kegg database to find correspondance between Kegg gene IDs and
% Ensembl gene IDs

%INPUT PARAMETERS

% 1 KegSpecies: KEGG species (hsa, mmu, rn0 ...)

%OUTPUT FILE

% file %KeggSpecies_kegg containing EnsGeneRank and KeggGeneRank

%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%
%                          c) Michel Bellis                                                %
%                          michel.bellis@crbm.cnrs.fr                                      %
%            Affiliation:  CNRS (Centre National de la Recherche Scientifique - France)    %                               
%  Bioinformatic Project:  ARRAYMATIC => http://code.google.com/p/arraymatic               %
%        Code Repository:  GITHUB => http://github.com/mbellis                             %
%          Personal Page:  http://bns.crbm.cnrs.fr                                         %
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
%  THIS CODE IS DISTRIBUTED UNDER THE CeCILL LICENSE, WHICH IS COMPATIBLE WITH       %
%  THE GNU GENERAL PUBLIC LICENCE AND IN ACCORDANCE WITH THE EUROPEAN LEGISLATION.   %
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%


function kegg2ensembl(KeggSpecies)
global K

cd(K.dir.gene)
if ~ exist('@KEGG')
    url = 'http://soap.genome.jp/KEGG.wsdl';
    createClassFromWsdl(url);
end

Kegg=KEGG;


%recover list of Kegg genes
eval(sprintf('GeneNb=get_number_of_genes_by_organism(Kegg,''%s'')',KeggSpecies));
eval(sprintf('KGeneId=get_genes_by_organism(Kegg,''%s'',1,%u);',KeggSpecies,GeneNb));
EnsGeneRank=zeros(GeneNb,1);
KeggGeneRank=zeros(GeneNb,1);

%load Kegg gene information by series of 100 genes (maximum allowed)
BlocNb=ceil(GeneNb/100);
for BlocL=1:BlocNb
    BlocL
    %construct command
    Offset=(BlocL-1)*100;
    Command='GeneInfo=bget(Kegg,''';    
    for GeneL=1:min(100,GeneNb-Offset)
        Command=[Command,' ',KGeneId{GeneL+Offset}];
    end
    Command=[Command,''');'];
    eval(Command);
    %find EnsGeneRank
    EnsemblPos=findstr('            Ensembl: ENS',GeneInfo);
    for GeneL=1:min(100,GeneNb-Offset)
        GenePos1=findstr(sprintf('ENTRY       %s ',KGeneId{GeneL+Offset}(5:end)),GeneInfo);
        if GeneL<min(100,GeneNb-Offset)
            try
                GenePos2=findstr(sprintf('ENTRY       %s ',KGeneId{GeneL+Offset+1}(5:end)),GeneInfo);
            catch
                GenePos2=length(GeneInfo)+1;
            end
        else
            GenePos2=length(GeneInfo)+1;
        end
        if isempty(GenePos1)
            h=warndlg(sprintf('%s not found',KGeneId{GeneL+Offset}));
            waitfor(h)
        else
            CurrKeggGeneId=regexp(KGeneId{GeneL+Offset},'[0-9]+','match');
            KeggGeneRank(GeneL+Offset)=str2num(CurrKeggGeneId{1});
            Pos=find(EnsemblPos>GenePos1&EnsemblPos<GenePos2);
            if ~isempty(Pos)
                CurrEnsemblGeneId=GeneInfo(EnsemblPos(Pos)+20:EnsemblPos(Pos)+38);
                CurrEnsemblGeneId=regexp(CurrEnsemblGeneId,'[0-9]+','match');                
                EnsGeneRank(GeneL+Offset)=str2num(CurrEnsemblGeneId{1});
            end
        end
    end
end

eval(sprintf('save %s_kegg EnsGeneRank KeggGeneRank',KeggSpecies))

