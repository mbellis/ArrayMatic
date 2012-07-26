%GSEA_ONLIST
%find significative gene sets in a list of results

%INPUT PARAMETERS
% RoundNb: how many random combinations between two groups of biological conditions
% Comparison: either '2vs1' or '1vs2' (indicates which biological group is test (first
%             position) and which is control (snd position) in test versus control comparison
% ModelRank: chip rank
% Values : values to be searched for GSEA
% CluFlag=1 : exist clusters in used gene set (no clusters in cliques
%          ; clusters in Kegg pathways)
% AbsentIndex: probe sets that must not be used in tested gene sets

%EXTERNAL FILES

%OUTPUT PARAMETERS


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

%function [GsIndex,GsClasses]=gsea_onlist(ModelRank,ZVarValues,FdrValues,MeanRanks,CluFlag,AbsentIndex,GeneSetFile,GeneSetDir,BestFile,BestDir,GeneId,Display_flag,PrintFlag)
function gsea_onlist(RoundNb,Comparison,varargin)
%gsea_findgeneset(11,OperonZVar,OperonFdr,OperonMeanRanks,0,AbsentIndex,'m11_cliques','/home/mbellis/array2/mbellis/sosma/arraymatic/amcdata/cliques','m11_cliques_type1_r5489_c2000','/home/mbellis/array2/mbellis/sosma/arraymatic/amcdata/cliques',TairId)
%gsea_findgeneset(11,OperonZVar,OperonFdr,OperonMeanRanks,0,AbsentIndex,'m11_athgsdb11',K.dir.affyMetadata,'athgsdb11_type1_r5489_c2000',K.dir.affyMetadata,TairId)
global K P
PsNb=P.chip.currProbeSetNb;
if nargin==3
    RandSeed=varargin{1};
else
    RandSeed=5489;
end
rand('twister',RandSeed)

[ModelRank,NetRank,NetPos]=select_net('unique',P.chip.chipRank);
if length(K.net{ModelRank}.biolRank{NetPos})==2
    BiolRank{1}=K.net{ModelRank}.biolRank{NetPos}{1};
    BiolRank{2}=K.net{ModelRank}.biolRank{NetPos}{2};
else
    h=errodlg('select a network with two groups of biological conditions');
    waitfor(h)
    error('process canceled')
end
BiolIndex=find(P.net.biolIndex);
BiolPos=cell(1,2);
for GrpL=1:2
    BiolPos{GrpL}=[];
    for BiolL=1:length(BiolRank{GrpL})
        CurrBiolPos=find(BiolIndex==BiolRank{GrpL}(BiolL));
        if ~isempty(CurrBiolPos)
            BiolPos{GrpL}=[BiolPos{GrpL};CurrBiolPos];
        end
    end
end
BiolNb=min(length(BiolPos{1}),length(BiolPos{2}));

%positions in comparison matrix
TABLE_SIZE=1000;
UsedBiolNb=length(find(P.net.biolIndex));
TableFdrPos=uint32(zeros(UsedBiolNb));
Pos=0;
for BiolL1=1:UsedBiolNb-1
    for BiolL2=BiolL1+1:UsedBiolNb
        Pos=Pos+1;
        TableFdrPos(BiolL1,BiolL2)=Pos;
    end
end
TableBiolPos=zeros(1,BiolNb*BiolNb);
FdrPos=reshape(1:BiolNb*BiolNb,BiolNb,BiolNb);
%FdrPos=FdrPos';
Inversion=[];
Pos=0;
for BiolL1=1:BiolNb
    for BiolL2=1:BiolNb
        Pos=Pos+1;
        if BiolPos{1}(BiolL1)<BiolPos{2}(BiolL2)
            TableBiolPos(Pos)=TableFdrPos(BiolPos{1}(BiolL1),BiolPos{2}(BiolL2));
            if isequal(Comparison,'2vs1')
                Inversion=[Inversion,Pos];
            end
        else
            TableBiolPos(Pos)=TableFdrPos(BiolPos{2}(BiolL2),BiolPos{1}(BiolL1));
            if isequal(Comparison,'1vs2')
                Inversion=[Inversion,Pos];
            end
        end
    end
end

%load FDR value 
TableNb=ceil(PsNb/TABLE_SIZE);
%recover fdr data
Fdr=single(ones(PsNb,BiolNb*BiolNb));
cd(P.dir.data)
for TableL=1:TableNb
%for TableL=1
    eval(sprintf('load m%u_%u_1',ModelRank,TableL))
    if TableL<TableNb     
        Fdr(((TableL-1)*TABLE_SIZE)+1:TableL*TABLE_SIZE,:)=Table1(:,TableBiolPos);
    else
        Fdr(((TableL-1)*TABLE_SIZE)+1:end,:)=Table1(:,TableBiolPos);
    end
end
clear Table1
Fdr(:,Inversion)=-Fdr(:,Inversion);


PathGsea=cell(RoundNb,1);
CliGsea=cell(RoundNb,1);
for RoundL=1:RoundNb
    PathGsea{RoundL}=cell(BiolNb,1);
    CliGsea{RoundL}=cell(BiolNb,1);
    RandPos=randperm(BiolNb);
    for CompL=1:BiolNb
        CompL
        [PathGsea{RoundL}{CompL},Temp]=gsea_findgeneset(ModelRank,[],Fdr(:,FdrPos(CompL,RandPos(CompL))),[],0,[],'m11_athgsdb11',K.dir.affyMetadata,'athgsdb11_type1_r5489_c2000',K.dir.affyMetadata,{},0,0);
        [CliGsea{RoundL}{CompL},Temp]=gsea_findgeneset(ModelRank,[],Fdr(:,FdrPos(CompL,RandPos(CompL))),[],0,[],'m11_cliques',K.dir.affyMetadata,'m11_cliques_type1_r5489_c2000',K.dir.affyMetadata,{},0,0);
    end
end
'stop'
