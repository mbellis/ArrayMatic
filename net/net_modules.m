%OBJECTIVE
%Allows to find modules by using net_agregate function
%which clusterizes probesets around hubs

%PARAMETERS
%Corr
%Anti
%MinCluSize
%LoadFlag
%ModelRAnk
%NetRank


%VERSION
%0-0-1 : uses directly the sum of corr_anti to find hubs
%0-0-0 : select used edges on corr-anti limit value, then
%           works on binary matrix (Corr =< CorrSel

function net_modules(Corr,Anti,MinCluSize,LoadFlag,varargin)
global K
%TEST_MYCLUSTER(8,5,20,2);
if LoadFlag
    ModelRank=varargin{1};
    NetRank=varargin{2};
end

NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',NetRank));


%select the paramaters to be used
Parameters=inputdlg({'CorrSupLimit','CorrInfLimit','AntiLimit','ConnStdFactor','CorrStdFactor','DiffFlag','CliqueFactor'},'',1,{''40','30','5','1','1','1','0.9'});
CorrSupLimit=str2num(Parameters{1});
CorrInfLimit=str2num(Parameters{2});
AntiLimit=str2num(Parameters{3});
ConnStdFactor=str2num(Parameters{4});
CorrStdFactor=str2num(Parameters{5});
DiffFlag=str2num(Parameters{6});
CliqueFactor=str2num(Parameters{7});



Clusters={};
Hubs={};
Cliques={};
Deviations={};
PsRanks=[];
CluSizes=[];
if LoadFlag==0
    for j=1:length(Corr)
        Corr(j,j)=0;
    end
else

end



Continue=1;
RoundL=0;
LastRound=0;
InfLimit=CorrInfLimit;

%find primary clusters
while Continue==1
    RoundL=RoundL+1;
    CorrSel=uint8(zeros(size(Corr)));
    if LastRound==1
        InfLimit=1;
    end

    if DiffFlag==0
        CorrSel(find(Corr>=CorrSupLimit&Anti<=AntiLimit))=1;
        [Clusters{RoundL},Hubs{RoundL},Cliques{RoundL},Deviations{RoundL}]=TWO_CLUSTERS(Corr,CorrSel,InfLimit,0,ConnStdFactor,CorrStdFactor,0,CliqueFactor,ModType);
    else
        CorrSel(find(Corr-Anti>=CorrSupLimit))=1;
        [Clusters{RoundL},Hubs{RoundL},Cliques{RoundL},Deviations{RoundL}]=TWO_CLUSTERS(Corr-Anti,CorrSel,InfLimit,0,ConnStdFactor,CorrStdFactor,0,CliqueFactor,ModType);
    end

    %what is changed if another hub is used first
    Sizes{RoundL}=[];
    NewClusters{RoundL}={};
    NewDeviations{RoundL}={};
    NewCliques{RoundL}={};
    if DiffFlag==0
        for HubL=1:length(Hubs{RoundL})
            [NewClusters{RoundL}{HubL},NewDeviations{RoundL}{HubL},NewCliques{RoundL}{HubL}]=net_agregate(Corr,CorrSel,InfLimit,ConnStdFactor,CorrStdFactor,Hubs{RoundL}(HubL),CliqueFactor,ModType);
            Sizes{RoundL}=[Sizes{RoundL};[length(Clusters{RoundL}{HubL}),length(Cliques{RoundL}{HubL}),length(NewClusters{RoundL}{HubL}),length(NewCliques{RoundL}{HubL})]];
        end
    else
        for HubL=1:length(Hubs{RoundL})
            [NewClusters{RoundL}{HubL},NewDeviations{RoundL}{HubL},NewCliques{RoundL}{HubL}]=net_agregate(Corr-Anti,CorrSel,InfLimit,ConnStdFactor,CorrStdFactor,Hubs{RoundL}(HubL),CliqueFactor,ModType);
            Sizes{RoundL}=[Sizes{RoundL};[length(Clusters{RoundL}{HubL}),length(Cliques{RoundL}{HubL}),length(NewClusters{RoundL}{HubL}),length(NewCliques{RoundL}{HubL})]];
        end
    end


    %how many hubs are not used
    if length(Hubs{RoundL})>1
        HubSizes=sum(CorrSel);
        [HubSizes,CurrHubs]=sort(HubSizes);
        try
            LastHubPos=find(CurrHubs==Hubs{RoundL}(end-1));
        catch
            'stop'
        end
        NotUsedHubs=[];
        for HubL=LastHubPos:length(CurrHubs)
            if isempty(find(Hubs{RoundL}==CurrHubs(HubL)))
                NotUsedHubs=[NotUsedHubs;CurrHubs(HubL)];
            end
        end
        sprintf('%u hubs larger than the last high hub (hub %u with %u ps) are not used as such',length(NotUsedHubs),Hubs{RoundL}(end-1),length(Clusters{RoundL}{end-1}))
    end

    %find ps that are shared and reassign them to the cluster where they have
    %the highest connectivity degree

    PsInfo{RoundL}=zeros(length(PosIndex),8);
    PsInfo{RoundL}(:,7)=ones(length(PosIndex),1)*NaN;
    if length(Hubs{RoundL})>1
        for HubL=1:length(Hubs{RoundL})-1
            for PsL=1:length(Clusters{RoundL}{HubL})
                %fill information for the current round
                PsInfo{RoundL}(Clusters{RoundL}{HubL}(PsL),1)=Hubs{RoundL}(HubL);
                PsInfo{RoundL}(Clusters{RoundL}{HubL}(PsL),2)=HubL;
                PsInfo{RoundL}(Clusters{RoundL}{HubL}(PsL),3)=Deviations{RoundL}{HubL}(PsL);
                if ~isempty(find(Cliques{RoundL}{HubL}==Clusters{RoundL}{HubL}(PsL)))
                    PsInfo{RoundL}(Clusters{RoundL}{HubL}(PsL),4)=1;
                end
            end
        end
        for HubL=1:length(Hubs{RoundL})-1
            %Integrate NewClustering information
            for PsL=1:length(NewClusters{RoundL}{HubL})
                %info not yet written
                SwapIt=0;
                if isnan(PsInfo{RoundL}(NewClusters{RoundL}{HubL}(PsL),7))
                    SwapIt=1;
                elseif NewDeviations{RoundL}{HubL}(PsL)>PsInfo{RoundL}(NewClusters{RoundL}{HubL}(PsL),7)
                    SwapIt=1;
                end
                if SwapIt==1
                    PsInfo{RoundL}(NewClusters{RoundL}{HubL}(PsL),5)=Hubs{RoundL}(HubL);
                    PsInfo{RoundL}(NewClusters{RoundL}{HubL}(PsL),6)=HubL;
                    PsInfo{RoundL}(NewClusters{RoundL}{HubL}(PsL),7)=NewDeviations{RoundL}{HubL}(PsL);
                    if ~isempty(find(NewCliques{RoundL}{HubL}==NewClusters{RoundL}{HubL}(PsL)))
                        PsInfo{RoundL}(NewClusters{RoundL}{HubL}(PsL),8)=1;
                    end
                end
            end
        end



        %find cliques for new assignations
        NextCliques{RoundL}=cell(length(Hubs{RoundL}),1);
        for HubL=1:length(Hubs{RoundL})-1
            CliqueCorr=Corr(find(PsInfo{RoundL}(:,6)==HubL),find(PsInfo{RoundL}(:,6)==HubL));
            CorrSel=uint8(zeros(size(CliqueCorr)));
            CorrSel(find(CliqueCorr>=InfLimit))=1;
            NodeDegree=sum(CorrSel);
            NextCliques{RoundL}{HubL}=find(NodeDegree>=min(ceil(length(CliqueCorr)*CliqueFactor),length(CliqueCorr)-1));
        end

        for HubL=1:length(Hubs{RoundL})-1
            Sizes{RoundL}(HubL,5)=0;
            Sizes{RoundL}(HubL,6)=length(find(PsInfo{RoundL}(:,6)==HubL));
            Sizes{RoundL}(HubL,7)=0;
            Sizes{RoundL}(HubL,8)=0;
        end
        Sizes{RoundL}(end,4)=0;
        Sizes{RoundL}(end,5)=0;

        %some previously assigned ps could miss ?
        MissingPs{RoundL}=[];
        for HubL=1:length(Hubs{RoundL})
            Sizes{RoundL}(HubL,5)=length(setdiff(Clusters{RoundL}{HubL},NewClusters{RoundL}{HubL}));
            Sizes{RoundL}(HubL,7)=length(NextCliques{RoundL}{HubL});
            Sizes{RoundL}(HubL,8)=length(setdiff(Clusters{RoundL}{HubL},find(PsInfo{RoundL}(:,6)==HubL)));
        end
        %             'end'
        %
        b=find(Sizes{RoundL}(:,6)>=MinCluSize);
        d=find(Sizes{RoundL}(:,6)<MinCluSize);
        sprintf('region %u : %u Clu>=%u (%u ps), %u Clu<%u (%u ps)', RegionL,length(b),MinCluSize,sum(Sizes{RoundL}(b,6)),length(d),MinCluSize,sum(Sizes{RoundL}(d,6)))
    end
    if LastRound==1|length(Hubs{RoundL})==1
        Continue=0;
        if length(Hubs{RoundL})>1
            SupPos=find(Sizes{RoundL}(:,6)>0);
            [temp SortIndex]=sort(Sizes{RoundL}(SupPos,6));
            SortIndex=flipud(SortIndex);
            SupPos=SupPos(SortIndex);
            if ~isempty(SupPos)
                for SupL=1:length(SupPos)
                    PsRanks=[PsRanks;RankIndex(find(PsInfo{RoundL}(:,6)==SupPos(SupL)))];
                    CluSizes=[CluSizes;length(find(PsInfo{RoundL}(:,6)==SupPos(SupL)))];
                end
            end
        else
            try
                PsRanks=[PsRanks;RankIndex];
                CluSizes=[CluSizes;length(RankIndex)];
            catch
                'stop'
            end
        end
    else
        if sum(Sizes{RoundL}(b,6))/CurrNb<0.10
            LastRound=1;
        end
        %write right ps Rank
        SupPos=find(Sizes{RoundL}(:,6)>=MinCluSize);
        [temp SortIndex]=sort(Sizes{RoundL}(SupPos,6));
        SortIndex=flipud(SortIndex);
        SupPos=SupPos(SortIndex);
        if ~isempty(SupPos)
            for SupL=1:length(SupPos)
                PsRanks=[PsRanks;RankIndex(find(PsInfo{RoundL}(:,6)==SupPos(SupL)))];
                CluSizes=[CluSizes;length(find(PsInfo{RoundL}(:,6)==SupPos(SupL)))];
            end
        end
        %start a new round with non clusterized ps
        %recover all non clusterized ps plus ps in clusters<MinCluSize
        InfPos=find(Sizes{RoundL}(:,6)<MinCluSize);
        PosIndex=[];
        if ~isempty(InfPos)
            for InfL=1:length(InfPos)
                PosIndex=[PosIndex;find(PsInfo{RoundL}(:,6)==InfPos(InfL))];
            end
        end
        PosIndex=[PosIndex;find(isnan(PsInfo{RoundL}(:,7)))];
        PosIndex=unique(PosIndex);
        CurrNb=length(PosIndex)
        'end'
    end
end %of while Continue==1
%construct definitive index of ps
CluSizes(end+1)=0;

for PsL=1:length(MemRankIndex);
    if isempty(find(PsRanks==MemRankIndex(PsL)))
        PsRanks=[PsRanks;MemRankIndex(PsL)];
        CluSizes(end)=CluSizes(end)+1;
    end
end
Region{RegionL}.ranks=PsRanks;
Region{RegionL}.sizes=CluSizes;

else %of if ~isempty(RankIndex)
    Region{RegionL}.ranks=[];
    Region{RegionL}.sizes=[];
end
end %of RegionL
'end'
% Ps=[];
% Sizes=[];
% for i=1:length(Region)
%     CurrSize=Region{i}.sizes;
%     CurrSize(find(CurrSize==0))=[];
%     Ps=[Ps;Region{i}.ranks];
%     Sizes=[Sizes;CurrSize];
% end
%
% figure
% subplot(2,1,1)
% plot(Sizes)
% subplot(2,1,2)
% plot(sort(Sizes))
% a=find(Sizes>=MinCluSize);
% b=find(Sizes<MinCluSize);
% length(a)
% length(b)
% sum(Sizes(a))
% sum(Sizes(b))
%Complete the Region Info
if RegionNb<max(Clu)
    for i=RegionNb+1:max(Clu)
        Region{i}.ranks=find(Clu==i);
        Region{i}.sizes=length(find(Clu==i));
    end
end
Region{max(Clu)+1}.ranks=find(Clu==0);
Region{max(Clu)+1}.sizes=length(find(Clu==0));



eval(sprintf('save %s_cs%u_ci%u_a%u_nstd%u_rstd%u_d%u_ct%u.mat Region -MAT',FileRoot,CorrSupLimit,CorrInfLimit,AntiLimit,strrep(ConnStdFactor,'.',''),strrep(ConnStdFactor,'.',''),DiffFlag,ModType))

