%GSEA_BESTDISTRIUBTION
%Calculate the binned enrichment scores (BES) for a series of control data. Control data may be either a series of variation values (comparisons between biological
%conditions) or a series of rank values (in a series of biological conditions). Gene sets may be either external lists (Kegg, MsigDB ...) or
%Arraymatic defined cliques and modules.
%As the number of control data is in general huge, only a subset of them is used (either randomly or explicitly selected).

%INPUT PARAMETERS
% GsFile = the geneset file to be loaded
% GsDir = the directory of geneset file
% GsName = the name of saved results
% TotalNb = nb of control data present in a Table
% UsedNb = nb of control data used to calculate the BES
% ModelRank = rank of chipset model used
% NetRank = rank of the network used if modules are used
% CluFlag=1 => The geneset has been analysed for the presence of clusters
% ChangeFlag = 0 => calculate BEST for only Increase and Decrease
%            = 1 => calculate BEST for Increas/Decrease abd Change/NotChange
%RandState : state of random generator (twister) at the begining (by default = 5489)
% DataType: 'ZVar' or 'Fdr' or 'Sensitivity'

% In this case exists Gs.clusters{1,2,3} (position in Gs.index, corresponding, respectively,  to the first cluster, the second cluster and the reliquate

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


function gsea_bestdistribution(GsFile,GsDir,GsName,TotalNb,UsedNb,ModelRank,NetRank,CluFlag,ChangeFlag,DataType,RandState,varargin)
global K
tic
%gsea_bestdistribution('m11_cliques.mat',K.dir.cliques,'m11_cliques',186355,2000,11,22,0,1,'Fdr',5489)
%gsea_bestdistribution('m11_cliques.mat',K.dir.cliques,'m11_cliques',186355,2000,11,22,0,1,'Fdr',5490)
%gsea_bestdistribution('m11_cliques.mat',K.dir.cliques,'m11_cliques',186355,2000,11,22,0,1,'Fdr',5491)
%gsea_bestdistribution('m11_athgsdb11.mat',K.dir.affyMetadata,'athgsdb11',186355,2000,11,22,0,1,'Fdr',5489)

% cd(fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',NetRank)))
% eval(sprintf('load m%un%u_cliques_50.mat',ModelRank,NetRank)
% eval(sprintf('load m%un%u_quasi-cliques_clusters_size5_50.mat',ModelRank,NetRank))
% 
% CliqueNb=0;
% RegionNb=length(Cliques);
% for RegionL=1:RegionNb
%     CliqueNb=CliqueNb+length(Cliques{RegionL});
% end
% Gs.name=cell(CliqueNb,1);
% Gs.index=cell(CliqueNb,1);
% Gs.region=zeros(CliqueNb,1);
% Gs.cliqueRank=zeros(CliqueNb,1);
% Gs.cliquePos=zeros(CliqueNb,1);
% GsPos=0;
% for RegionL=1:RegionNb
%     CurrCliques=Cliques{RegionL};
%     CliquePos=0;
%     for CliL=1:length(CurrCliques)
%         CliquePos=CliquePos+1;
%         CliqueRank=CurrCliques(CliL);
%         GsPos=GsPos+1;
%         Gs.name{GsPos}=sprintf('r%02u_cp%03u_cr%04u',RegionL,CliquePos,CliqueRank);
%         Gs.region(GsPos)=RegionL;
%         Gs.cliqueRank(GsPos)=CliqueRank;
%         Gs.cliquePos(GsPos)=CliquePos;
%         Gs.index{GsPos}=find(Clu==CliqueRank);
%     end
% end

DISPLAY_IT=1;


if nargin==11
    %Control data are randomly selected
    RandFlag=1;
    rand('twister',RandState)
    Positions=ceil(rand(UsedNb,1)*TotalNb);
else
    %Control data are selected by list
    RandFlag=0;
    Positions=varargin{1};
    UsedNb=length(Positions);
end

PsNb=K.chip.probeSetNbs{ModelRank}(1);
BlocNb=K.net{ModelRank}.blocNb(NetRank);
BlocSize=K.net{ModelRank}.blocSize(NetRank);

%load gene sets
cd(GsDir)
load(GsFile)
%recover non empty gene sets
GsPos=[];
for i=1:length(Gs.name)
    if ~isempty(Gs.index{i})
        GsPos=[GsPos;i];
    end
end
GsNb=length(GsPos);


Bindex=zeros(PsNb,1);
Bin=ceil(PsNb/100);
BinIndex=[[Bin:Bin:Bin*10],[Bin*10+Bin*4:Bin*4:Bin*10+Bin*24],[PsNb-Bin*10-Bin*24:Bin*4:PsNb-Bin*10-Bin*4],[PsNb-Bin*10:Bin:PsNb-Bin]];
BinNb=length(BinIndex);

%load data
cd(K.dir.tables)
cd(sprintf('m%u',ModelRank))
cd(sprintf('m%u_data',ModelRank))
Data=single(zeros(PsNb,UsedNb));
for BlocL=1:BlocNb
    load(sprintf('m%u_%u_1.mat',ModelRank,BlocL))
    if BlocL<BlocNb
        if size(Table1,1)~=BlocSize|size(Table1,2)~=TotalNb
            h=errordlg(sprintf('Table %u is %ux%u instead of %ux%u',BlocL,size(Table1,1),size(Table1,2),BlocSize,TotalNb));
            waitfor(h)
            error('process canceled')
        end
        Data(((BlocL-1)*BlocSize)+1:BlocL*BlocSize,:)=Table1(:,Positions);
    else
        if size(Table1,1)~=mod(PsNb,BlocSize)|size(Table1,2)~=TotalNb
            h=errordlg(sprintf('Table %u is %ux%u instead of %ux%u',BlocL,size(Table1,1),size(Table1,2),mod(PsNb,BlocSize),TotalNb));
            waitfor(h)
            error('process canceled')
        end

        Data(((BlocL-1)*BlocSize)+1:PsNb,:)=Table1(:,Positions);
    end
end
clear Table1


if ChangeFlag
    RoundNb=2;
else
    RoundNb=1;
end
for RoundL=1:RoundNb
    if CluFlag==1
        TypeNb=4;
    else
        TypeNb=1;
    end
    for TypeL=1:TypeNb
        BinCount{RoundL}=cell(GsNb,1);
        for  GsL=1:GsNb
            BinCount{RoundL}{GsL}=single(zeros(UsedNb,BinNb));
        end

        for GsL=1:GsNb
            RoundL
            GsL
            GsRank=GsPos(GsL);
            switch TypeL
                case 1
                    Index=Gs.index{GsRank};
                case 2
                    Index=Gs.index{GsRank}(Gs.clusters{GsRank}{1});
                case 3
                    Index=Gs.index{GsRank}(Gs.clusters{GsRank}{2});
                case 4
                    Index=Gs.index{GsRank}(Gs.clusters{GsRank}{3});
            end
            for DataL=1:UsedNb
                %matrix of zeros
                Val=Bindex;
                %mark the position where the geneset is present;
                Val(Index)=1;
                %recover a particular result
                CurrData=Data(:,DataL);
                if RoundL==2
                    %FIRST ROUND : INCREASED/DECREASED
                    CurrData=abs(CurrData);
                end

                if isequal(DataType,'ZVar')
                    [Temp SortIndex]=sort(CurrData,'descend');
                else                    
                    NullPos=find(CurrData==0);
                    if ~isempty(NullPos)
                        h=errordlg('Exist null fdr => can not know it they are increased or decreased. Restart and used ZVar values.');
                        waitfor(h)
                        error('process canceled')
                    end
                    if RoundL==2
                        [Temp SortIndex]=sort(CurrData);
                    else
                        Pos=find(CurrData>0);
                        [Temp CurrSortIndex]=sort(CurrData(Pos));
                        SortIndex=Pos(CurrSortIndex);
                        Neg=find(CurrData<0);
                        [Temp CurrSortIndex]=sort(CurrData(Neg));
                        SortIndex=[SortIndex;Neg(CurrSortIndex)];
                    end
                end
                CurrVal=Val(SortIndex);                
                BinCount{RoundL}{GsL}(DataL,:)=single(BINNED_RUNNING_SUM(CurrVal,BinIndex));
            end
        end
    end
end
cd(GsDir)
eval(sprintf('save %s_type%u_r%u_c%u BinCount',GsName,TypeL,RandState,UsedNb))
Mem=BinCount;
for RoundL=1:2
    BinCount{RoundL}=cell(GsNb,1);
    for  GsL=1:GsNb
        BinCount{RoundL}{GsL}=Mem{RoundL}{GsL}(:,[1:32]);
    end
end
    
toc
if DISPLAY_IT
    GENE_SET=[1,2,3,130];
    for RoundL=1:2
        for GsL=1:length(GENE_SET)
            CurrGs=GENE_SET(GsL);
            MidPos=floor(BinNb/2);
            h=figure;            
            if RoundL==1
                set(h,'name',sprintf('BEST for Inc/Dec for gene set %u',CurrGs))
            else
                set(h,'name',sprintf('BEST for Ch/NotCh for gene set %u',CurrGs))
            end
            set(gcf,'color',[1,1,1])
            hold on            
            if mod(BinNb,2)==1
                plot(BinIndex(1:MidPos),BinCount{RoundL}{CurrGs}(:,1:MidPos)')
                plot(BinIndex(MidPos:end-1),BinCount{RoundL}{CurrGs}(:,MidPos:end-1)')
                YLim=get(gca,'ylim');
                line([BinIndex(MidPos),BinIndex(MidPos)],YLim,'color','b','linestyle','-.','linewidth',3)
            else
                plot(BinIndex(1:MidPos),BinCount{RoundL}{CurrGs}(:,1:MidPos)')
                plot(BinIndex(MidPos+1:end-1),BinCount{RoundL}{CurrGs}(:,MidPos+1:end-1)')
                YLim=get(gca,'ylim');
                line([BinIndex(MidPos),BinIndex(MidPos)],YLim,'color','b','linestyle','-.','linewidth',3)
                line([BinIndex(MidPos+1),BinIndex(MidPos+1)],YLim,'color','b','linestyle','-.','linewidth',3)    
            end         
            set(gca,'box','on')
            set(gca,'xtick',BinIndex)
            set(gca,'xticklabel','')            
            xlabel('position from top list')
            ylabel('binned enrichment score from top (BEST)')
            if RoundL==1
                title(sprintf('BEST for Inc/Dec for gene set %u',CurrGs))
            else
                title(sprintf('BEST for Ch/NotCh for gene set %u',CurrGs))
            end

            h=figure;
            set(gcf,'color',[1,1,1])
            if RoundL==1
                set(h,'name',sprintf('BEST HIST for Inc/Dec for gene set %u',CurrGs))
            else
                set(h,'name',sprintf('BEST HIST for Ch/NotCh for gene set %u',CurrGs))
            end
            ES=BinCount{RoundL}{CurrGs}';
            ValNb=size(ES,1)-1;
            ColNb=ceil(sqrt(ValNb));
            RawNb=round(ValNb/ColNb);
            if RawNb*ColNb<ValNb
                ColNb=ColNb+1;
            end
            for ValL=1:ValNb
                if length(unique(ES(ValL,1:end-1)))>1
                subplot(RawNb,ColNb,ValL)
                hist(ES(ValL,1:end-1),100)
                set(gca,'xlim',[min(ES(ValL,:))-0.01,max(ES(ValL,:)+0.01)])
                set(gca,'tickdir','out')
                title(sprintf('%u',BinIndex(ValL)))
                end
            end
        end
    end
end

% %% BINNED RUNNING SUM
% function [ES]=BINNED_RUNNING_SUM(Val,BinIndex)
% 
% DISPLAY_IT=0;
% %one running sum from start to end
% ES=zeros(1,length(BinIndex));
% BinNb=(length(BinIndex)-1)/2;
% BinIndex=[1,BinIndex];
% SumVal=zeros(size(Val));
% HitPos=find(Val);
% HitNb=length(HitPos);
% MissPos=find(Val==0);
% MissNb=length(MissPos);
% SumVal(MissPos)=-1/MissNb;
% SumVal(HitPos)=1/HitNb;
% ESSum=cumsum(SumVal);
% for i=1:BinNb
%     ES(i)=max(ESSum(BinIndex(i):BinIndex(i+1)));
% end
% [t,Pos]=max(abs(ESSum(BinIndex(BinNb+1):BinIndex(BinNb+1+1))));
% ES(BinNb+1)=ESSum(BinIndex(i)+Pos);
% for i=BinNb+2:2*BinNb+1
%     ES(i)=min(ESSum(BinIndex(i):BinIndex(i+1)));
% end
% 
% if DISPLAY_IT
%     BinPos=[];
%     for i=1:2*BinNb+1
%         BinPos(i)=BinIndex(i+1)-(BinIndex(i+1)-BinIndex(i))/2;
%     end
%     for i=BinNb+2:2*BinNb+1
%         BinPos(i)=BinIndex(i+1)-(BinIndex(i+1)-BinIndex(i))/2;
%     end
%     figure
%     hold on
%     plot(ESSum,'k')
%     plot(BinPos,ES,'ro')
% end


%% BINNED RUNNING SUM
function [ES]=BINNED_RUNNING_SUM(Val,BinIndex)

%one running sum from start to end
% BinNb=length(BinIndex);
% ES=zeros(1,BinNb);
HitPos=find(Val);
HitNb=length(HitPos);
MissNb=length(Val)-HitNb;
MissVal=-1/MissNb;
HitVal=1/HitNb;
SumVal=ones(size(Val))*MissVal;
SumVal(HitPos)=HitVal;
ESSum=cumsum(SumVal);
ES=ESSum(BinIndex);
% 
% for BinL=1:BinNb
%     ES(BinL)=length(find(Val(1:BinIndex(BinL))))*HitVal+length(find(Val(1:BinIndex(BinL))==0))*MissVal;
% end
