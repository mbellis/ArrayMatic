%====================%
% FUNCTION KEGGLIST  %
%====================%

% KEGGLIST recover Kegg gene IDs and Ensembl gene IDs corresponding to a list of KEGG
% pathways

%INPUT PARAMETERS

% 1 KegSpecies: KEGG species (hsa, mmu ...)
% 2  PathwayId: list of pathway (numerical part only : [3050,20,4010])
% 3   ListName: name of the saved file

%OUTPUT

% save Pathway under ListName in K.dir.list
% Pathway(PathRank).keggGeneRank: KEGG gene IDS belonging to the pathway
% Pathway(PathRank).endGeneRank: numeric part of Ensembl gene IDS belonging to the pathway
% Pathway(PathRank).name: KEGG pathway name
% 

%EXAMPLE
%kegglist('mmu',[20,3050,4010],'mouse_krebs_proteasome_mapk_ens')
%kegglist('hsa',[20,3050,4010],'mouse_krebs_proteasome_mapk_ens')
%kegglist('hsa',4010,'human_mapk_ens')

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


function kegglist(KeggSpecies,PathwayId,ListName)
global K

PathwayNb=length(PathwayId);

%create eventually KEGG class
cd(K.dir.gene)
if ~ exist('@KEGG')
    url = 'http://soap.genome.jp/KEGG.wsdl';
    createClassFromWsdl(url);
end

%load correspondance Kegg gene id <=> Ensembl gene id
if  exist(sprintf('%s_kegg.mat',KeggSpecies),'file')
    eval(sprintf('load %s_kegg',KeggSpecies))
else
    h=errordlg(sprintf('%s_kegg.mat does not exist',KeggSpecies));
    waitfor(h)
end

%create KEGG object
Kegg=KEGG;

%recover list of pathways
eval(sprintf('Pathways=list_pathways(Kegg,''%s'');',KeggSpecies));

%recover list of kegg gene id
Pathway=[];
for PathL=1:PathwayNb    
    eval(sprintf('CurrKeggGeneId=get_genes_by_pathway(Kegg,''path:%s%05u'');',KeggSpecies,PathwayId(PathL)))        
    Pathway(PathL).keggGeneRank=zeros(length(CurrKeggGeneId),1);
    for GeneL=1:length(CurrKeggGeneId)
        CurrKeggGeneRank=regexp(CurrKeggGeneId{GeneL},'[0-9]+','match');
        Pathway(PathL).keggGeneRank(GeneL)=str2num(CurrKeggGeneRank{1});
    end
    Pathway(PathL).keggGeneRank=sort(Pathway(PathL).keggGeneRank);
    for PathL1=1:length(Pathways)
        if isequal(sprintf('path:%s%05u',KeggSpecies,PathwayId(PathL)),Pathways(PathL1).entry_id)
            Pathway(PathL).name=Pathways(PathL1).definition;
            break
        end
    end
end

%find corresponding Ensembl gene id
for PathL=1:PathwayNb
    CurrKeggGeneRank=Pathway(PathL).keggGeneRank;
    CurrEnsGeneRank=[];
    for GeneL=1:length(CurrKeggGeneRank)
        GenePos=find(KeggGeneRank==CurrKeggGeneRank(GeneL));
        if EnsGeneRank(GenePos)>0
            CurrEnsGeneRank(end+1,1)=EnsGeneRank(GenePos);
        end
    end
    Pathway(PathL).ensGeneRank=sort(CurrEnsGeneRank);
end

cd(K.dir.list)
eval(sprintf('save %s Pathway',ListName))

