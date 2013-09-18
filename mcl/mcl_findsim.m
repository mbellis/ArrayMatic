%find the one to one similarity betweens two list of clusters
%
%======================%
% FUNCTION MCL_FINDSIM %
%======================%
%
% INPUT PARAMETERS
%
%          Clu: two lists of clusterized probe sets (Clu{1} and Clu{2})
%      CluList: two lists of clusters that must be searched for correspondance
%      SimType: indicates which statistics is used to calculate overlap between two clusters
%               -    sqrt: #intersect*100/sqrt(#clu1*#clu2)
%               -     min: #intersect*100/min(#clu1,#clu2)
%               - jaccard: #intersect*100/#union(clu1,clu2)
% OneToOneFlag:
% 
% OUTPUT RESULTS
%
%   SimClu: list of clusters of the second list that are in correspondance with
%           the first list
%  Overlap: overlap between two clusters according to the similarity type selected
% TestSize: size derived from the sizes of the two combined clusters used to evaluate overlap

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

function [SimClu,Overlap,TestSize]=mcl_findsim(Clu,CluList,SimType,OneToOneFlag)


CluNb=CluList{1}(2)-CluList{1}(1)+1;
if OneToOneFlag
    if CluList{2}(2)-CluList{2}(1)+1<CluNb
        h=errordlg(sprintf('if OneToOneFlag=1, number of \nclusters in first list must\n be lower than number\nof clusters in second list'));
        waitofr(h)
        error('process canceled')
    end
end
%overlap between two clusters according to the similarity type selected
Overlap=zeros(1,CluNb);
%list of clusters of the second list that are in correspondance with the first list
SimClu=zeros(1,CluNb);
%size derived from the sizes of the two combined clusters used to evaluate overlap
TestSize=zeros(1,CluNb);
%each cluster of second list is used at most once
UsedClu=[];
for CluL1=CluList{1}(1):CluList{1}(2)
    Ps1=find(Clu{1}==CluL1);
    %find the best overlap
    SearchedClu=0;
    Score=0;
    for CluL2=CluList{2}(1):CluList{2}(2)
        if isempty(find(UsedClu==CluL2))
            Ps2=find(Clu{2}==CluL2);
            switch SimType
                case 'sqrt'
                    CurrScore=round(length(intersect(Ps1,Ps2))*100/sqrt(length(Ps1)*length(Ps2)));
                case 'min'
                    CurrScore=round(length(intersect(Ps1,Ps2))*100/min(length(Ps1),length(Ps2)));
                case 'jaccard'
                    CurrScore=round(length(intersect(Ps1,Ps2))*100/length(union(Ps1,Ps2)));
                    [CluL1 CluL2 length(Ps1) length(Ps2) length(intersect(Ps1,Ps2)) length(union(Ps1,Ps2)) CurrScore];
            end
            if CurrScore>=Score
                Score=CurrScore;
                SearchedClu=CluL2;
            end
        end
    end
    SimClu(CluL1)=SearchedClu;
    Ps2=find(Clu{2}==SearchedClu);
    if OneToOneFlag
        UsedClu=[UsedClu,SearchedClu];
    end
    switch SimType
        case 'sqrt'
            Overlap(CluL1)=round(length(intersect(Ps1,Ps2))*100/sqrt(length(Ps1)*length(Ps2)));
            TestSize(CluL1)=sqrt(length(Ps1)*length(Ps2));
        case 'min'
            Overlap(CluL1)=round(length(intersect(Ps1,Ps2))*100/min(length(Ps1),length(Ps2)));
            TestSize(CluL1)=min(length(Ps1),length(Ps2));
        case 'jaccard'
            Overlap(CluL1)=round(length(intersect(Ps1,Ps2))*100/length(union(Ps1,Ps2)));
            TestSize(CluL1)=length(union(Ps1,Ps2));
    end
end