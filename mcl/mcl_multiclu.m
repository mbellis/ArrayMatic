% indicate the multiple cluster assignation of probe sets
%=======================%
% FUNCTION MCL_MULTICLU %
%=======================%

% INPUT PARAMETERS
% 1  ChipRank: chip rank
% 2   RefPostFix: file name of MCL results to used as reference
% 3   PostFix: file name of MCL results to be analyzed
% 4   TypePos: type positions to be analyzed
% 5   ListPos: list positionsto be analyzed 





%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%
%                         c) Michel Bellis                                                %
%                         michel.bellis@crbm.cnrs.fr                                      %
%            Affiliation:  CNRS (Centre National de la Recherche Scientifique - France)    %
%  Bioinformatic Project:  ARRAYMATIC => http://code.google.com/p/arraymatic               %
%        Code Repository:  GITHUB => http://github.com/mbellis                             %
%          Personal Page:  http://bns.crbm.cnrs.fr                                         %
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
%  THIS CODE IS DISTRIBUTED UNDER THE CeCILL LICENSE, WHICH IS COMPATIBLE WITH       %
%  THE GNU GENERAL PUBLIC LICENCE AND IN ACCORDANCE WITH THE EUROPEAN LEGISLATION.   %
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

%mcl_multiclu(27,'n88ton108A','n118A',[2,2],[1,1])
%mcl_multiclu(27,'n88ton108A','n164G',[2,2],[1,2])

function mcl_multiclu(ChipRank,PostFix,RefPostFix,TypePos,ListPos)
global K

NetDir=fullfile(K.dir.net,sprintf('m%u',ChipRank),'mcl');
cd(NetDir)

%load MCL clusters
load(sprintf('m%u_mcl_%s.mat',ChipRank,PostFix))
if ~exist('CorrLimit')
    CLimit{1}=Info.limits{TypePos};
else
    CLimit{1}=CorrLimit;    
end
Cluster=Clu{TypePos(1)}{ListPos(1)};

%load reference MCL
load(sprintf('m%u_mcl_%s.mat',ChipRank,RefPostFix))
if ~exist('CorrLimit')
    CLimit{2}=Info.limits{TypePos};
else
    CLimit{2}=CorrLimit;
end
%assume that the first network is used 
CurrClu{2}=Clu{TypePos(2)}{ListPos(2)}{1}(:,1);
MaxClu=max(CurrClu{2});
clear Clu

'stop'

NetNb=length(Cluster);
PsNb=length(Cluster{1});
%construct a new clustering result by gathering the first column of each network and by
%reassigning cluster rank to the reference clusters
NewClu=zeros(PsNb,NetNb);
for NetL=1:NetNb
    CurrClu{1}=Cluster{NetL}(:,1);
    [SimClu,Overlap,TestSize]=mcl_findsim(CurrClu,{[1,max(unique(CurrClu{1}))],[1,MaxClu]},'jaccard',1);
    for CluL=1:max(unique(CurrClu{1}))
        PsPos=find(CurrClu{1}==CluL);
        NewClu(PsPos,NetL)=SimClu(CluL);
    end
end
'stop'
Combinations=unique(NewClu,'rows');
CombNb=zeros(size(Combinations,1),1);
for CombL=1:length(Combinations)
    FoundComb=1:PsNb;
    for NetL=1:NetNb
        Pos=find(NewClu(:,NetL)==Combinations(CombL,NetL));
        FoundComb=intersect(FoundComb,Pos);
        %length(FoundComb)
    end
    CombNb(CombL)=length(FoundComb);
end