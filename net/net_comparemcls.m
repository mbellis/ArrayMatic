%==========================%
% FUNCTION NET_COMPAREMCLS %
%==========================%

%NET_COMPAREMCLS: compares several results of MCL clustering
%INPUT PARAMETERS
% 1  ModelRank: list of two chip ranks
% 2   FileName: list of two MCL results
% 3   CompList: intra chip probe set list comparisons (one comparison per line)
% 4  ClassName: list legend
% 5 PsListFlag: indicates for each chip if the probe set lists must be loaded directly (1)
%               or reconstituted from Clustering results
% 6    SimType: how is calculated similarity of two clusters
%               either : 'min' (#intersect/min(clusters))
%                   or : 'geom' (#intersect/geometric mean)
%                   or : 'jaccard' (#intersect/union(clusters) (Jaccard Index)) 



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


function net_comparemcls(ChipRank,FileName,CompList,ClassName,PsListFlag,SimType)
%net_comparemcls([8,27],[4,3],[1,2;1,3;1,4;1,5;1,6;2,5;3,5;4,5],{'all','ss','ms','cx','s','rest','allc','ssc','msc','cxc','sc','restc'};)
%net_comparemcls([8,27],[1,1])
%net_comparemcls([27,27],[1,2])
%net_comparemcls([27],[2])
%net_comparemcls([27,27],{'n165','n166'})
%net_comparemcls([27,27],{'n165','n164'})
%1prc-10prc
%net_comparemcls([8,27],{'n229','n165'})
%1prc
%net_comparemcls([8,27],{'n227','n166'})
%1prm
% limit 10:10:60
%net_comparemcls([8,27],{'n228','n164'})
%limit 20:2:30
%net_comparemcls([8,27],{'n228A','n164A'})
%limit 12:2:18
%net_comparemcls([8,27],{'n228B','n164B'},[1,2],{'all','s'})
%net_comparemcls([27,27],{'n169A','n170A'},[1,2],{'all','s'},[0,0],'min')
global K

if length(ChipRank)>2
    h=errordlg('net_comparemcls need at most two chips');
    waitfor(h)
    error('process canceled')
end


%% RECOVER INFORMATION AND RESULTS OF MCL CLUSTERING
PsNb=[];
Information={};
ChipNb=length(ChipRank);
DiffChipFlag=0;
if ChipNb==2 & ChipRank(1)~=ChipRank(2)
    DiffChipFlag=1;
end
for ChipL=1:ChipNb
    ChipPos(ChipL)=find(K.chip.rank==ChipRank(ChipL));
    CmlDir=fullfile(K.dir.net,sprintf('m%u',ChipRank(ChipL)),'mcl');
    cd(CmlDir)
    eval(sprintf('load m%u_mcl_%s.mat',ChipRank(ChipL),FileName{ChipL}))
    Information{ChipL}=Info;
    Cluster{ChipL}=Clu;
    ClusterNb{ChipL}=CluNb;
    PsNb(ChipL)=K.chip.probesetNb(ChipPos(ChipL));
end
clear Info Clu CluNb

if ChipNb==1
    ListNb(1)=Information{1}.listNb;
    if ListNb(1)==1
        PlotPos{1}=1;
        RowNb=1;
        ColNb=1;
    elseif ListNb<=6
        PlotPos{1}=[1:ListNb(1)];
        RowNb=2;
        ColNb=ceil(ListNb(1)/2);
    else
        PlotPos{1}=[1,3,5,7,9,11,2,4,6,8,10,12];
        RowNb=6;
        ColNb=2;
    end
else
    ListNb(1)=Information{1}.listNb;
    ListNb(2)=Information{2}.listNb;
    if ListNb(1)==1
        PlotPos{1}=1;
        PlotPos{2}=2;
        RowNb=1;
        ColNb=2;
    elseif ListNb(1)<=6
        PlotPos{1}=[1:2:(ListNb(1)*2)-1];
        PlotPos{2}=[2:2:(ListNb(1)*2)];
        RowNb=ListNb(1);
        ColNb=2;
    else
        PlotPos{1}=[1,4,7,10,13,16,2,5,8,11,14,17];
        PlotPos{2}=[3,6,9,12,15,18];
        RowNb=6;
        ColNb=3;
    end
end

%process types of units used to measure correlation
TypeName={'C','C-A','rawC','raw(C-A)'};
Types=zeros(ChipNb,4);
MaxType=0;
for ChipL=1:ChipNb
    for TypeL=1:4
        if Information{ChipL}.type(TypeL)
            Types(ChipL,TypeL)=1;
            MaxType=max(MaxType,TypeL);
        end
    end
end

%recover ps lists
for ChipL=1:ChipNb
    if PsListFlag(ChipL)
        cd(fullfile(K.dir.net,sprintf('m%u',ChipRank(ChipL)),'list'))
        for ListL=1:ListNb(ChipL)
            fid=fopen(sprintf('m%u_pslist%u.u32',ChipRank(ChipL),Information{ChipL}.listRank(ListL)),'r','ieee-le');
            PsRankList{ChipL}{ListL}=fread(fid,inf,'uint32');
            ListSize{ChipL}(ListL)=length(PsRankList{ChipL}{ListL});
            fclose(fid);
        end
    else
        NetNb(ChipL)=length(Information{ChipL}.netRanks);
        ListSize{ChipL}=zeros(ListNb(ChipL),1);
        PsRankList{ChipL}=cell(ListNb(ChipL),1);
        for TypeL=1:4
            if length(Cluster{ChipL})>=TypeL
                if ~isempty(Cluster{ChipL}{TypeL})
                    for ListL=1:length(Cluster{ChipL}{TypeL})
                        if ~isempty(Cluster{ChipL}{TypeL}{ListL})
                            for NetL=1:NetNb(ChipL)
                                if ~isempty(Cluster{ChipL}{TypeL}{ListL}{NetL})
                                    for LimitL=1:size(Cluster{ChipL}{TypeL}{ListL}{NetL},2)
                                        if ListSize{ChipL}(ListL)<length(find(Cluster{ChipL}{TypeL}{ListL}{NetL}(:,LimitL)))
                                            ListSize{ChipL}(ListL)=length(find(Cluster{ChipL}{TypeL}{ListL}{NetL}(:,LimitL)));
                                            PsRankList{ChipL}{ListL}=find(Cluster{ChipL}{TypeL}{ListL}{NetL}(:,LimitL));
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% INTRA CHIP INTER PS LIST REPRODUCIBILITY

Legend={};
for CompL=1:size(CompList,1)
    Legend{end+1}=sprintf('%s-%s',ClassName{CompList(CompL,1)},ClassName{CompList(CompL,2)});
end
IntraChipCluOrder=cell(1,2);
for TypeL=1:MaxType
    DoIt=0;
    if Information{1}.type(TypeL)==1
        DoIt=1;
        if ChipNb==2
            if Information{2}.type(TypeL)==0
                DoIt=0;
            end
        end
    end
    if DoIt
        IntraChipCluOrder{TypeL}=cell(1,2);
        h=figure;
        set(h,'name',sprintf('%s: REPRODUCIBILITY Ps List vs Ps List',TypeName{TypeL}))        
        set(gcf,'color',[1,1,1])
        Colors=colors(colormap,length(CompList));
        MaxLimitNb=max(length(Information{1}.limits{TypeL}),length(Information{2}.limits{TypeL}));
        WeightedScore=cell(1,ChipNb);
        for ChipL=1:ChipNb
            IntraChipCluOrder{TypeL}{ChipL}=cell(size(CompList,1),1);
            LimitNb=length(Information{ChipL}.limits{TypeL});            
            if ~isempty(Types(ChipL,TypeL))
                WeightedScore{ChipL}{TypeL}=zeros(LimitNb,size(CompList,1));
                for CompL=1:size(CompList,1)
                    IntraChipCluOrder{TypeL}{ChipL}{CompL}=zeros((LimitNb+1)*4,10);
                    CurrComp=CompList(CompL,:);
                    CReprod{ChipL}=zeros(10*LimitNb,4*LimitNb);
                    for LimitL=1:LimitNb
                        UsedClu=[];                        
                        CRep=zeros(1,10);
                        SClu=zeros(1,10);
                        CluSize1=zeros(1,10);
                        CluSize2=zeros(1,10);
                        CurrSize=zeros(1,10);
                        for CluL1=1:10
                            a=find(Cluster{ChipL}{TypeL}{CurrComp(1)}{1}(:,LimitL)==CluL1);
                            CluSize1(CluL1)=length(a);
                            %find the best overlap
                            SearchedClu=0;
                            Score=0;
                            for CluL2=1:10
                                if isempty(find(UsedClu==CluL2))                                   
                                    b=find(Cluster{ChipL}{TypeL}{CurrComp(2)}{1}(:,LimitL)==CluL2);
                                    CluSize2(CluL2)=length(b);
                                    switch SimType
                                        case 'sqrt'
                                            CurrScore=round(length(intersect(a,b))*100/sqrt(length(a)*length(b)));
                                        case 'min'
                                            CurrScore=round(length(intersect(a,b))*100/min(length(a),length(b)));
                                        case 'jaccard'
                                            CurrScore=round(length(intersect(a,b))*100/length(union(a,b)));
                                    end
                                    if CurrScore>=Score
                                        Score=CurrScore;
                                        SearchedClu=CluL2;
                                    end
                                end
                            end
                            SClu(CluL1)=SearchedClu;                          
                            b=find(Cluster{ChipL}{TypeL}{CurrComp(2)}{1}(:,LimitL)==SearchedClu);
                            switch SimType
                                case 'sqrt'
                                    CRep(CluL1)=round(length(intersect(a,b))*100/sqrt(length(a)*length(b)));
                                    CurrSize(CluL1)=sqrt(length(a)*length(b));
                                case 'min'
                                    CRep(CluL1)=round(length(intersect(a,b))*100/min(length(a),length(b)));
                                    CurrSize(CluL1)=min(length(a),length(b));
                                case 'jaccard'
                                    CRep(CluL1)=round(length(intersect(a,b))*100/length(union(a,b)));
                                    CurrSize(CluL1)=length(union(a,b));
                            end
                        end
                        IntraChipCluOrder{TypeL}{ChipL}{CompL}(LimitL,:)=SClu;
                        IntraChipCluOrder{TypeL}{ChipL}{CompL}(LimitL+LimitNb+1,:)=CluSize1;
                        IntraChipCluOrder{TypeL}{ChipL}{CompL}(LimitL+LimitNb*2+2,:)=CluSize2(SClu);
                        IntraChipCluOrder{TypeL}{ChipL}{CompL}(LimitL+LimitNb*3+3,:)=CRep;
                        figure(h)
                        if ChipL==1    
                            subplot(MaxLimitNb,2,LimitL*2-1)
                        else
                            subplot(MaxLimitNb,2,LimitL*2)
                        end
                        hold on
                        plot([1:10],CRep,'color',Colors(CompL,:))
                        if CompL==1
                            if LimitL==1
                                title(sprintf('m%u-%s Limit = %u',ChipRank(ChipL),FileName{ChipL},Information{ChipL}.limits{TypeL}(LimitL)))
                            else
                                title(sprintf('Limit = %u',Information{ChipL}.limits{TypeL}(LimitL)))
                            end
                            %xlabel(sprintf('C-A m%u cluster rank',ChipRank(ChipL)))
                            set(gca,'xtick',[1:10])
                            set(gca,'xlim',[1,10])
                            set(gca,'ylim',[0,100])
                            set(gca,'xticklabel',SClu)
                            set(gca,'box','on')
                        end
                        %fill WeightedScore
                        WeightedScore{ChipL}{TypeL}(LimitL,CompL)=sum(CurrSize(1:7).*CRep(1:7))/sum(CurrSize(1:7));
                    end
                end
                legend(Legend)
            end
        end
    end
end

%% INTRA CHIP REPRODUCIBILITY BETWEEN C AND C-A
for ChipL=1:ChipNb
    if Types(ChipL,1)&Types(ChipL,2)
        h1=figure;
        set(h1,'name',sprintf('Chip m%u: REPRODUCIBILITY C & C-A: Cluster sizes',ChipRank(ChipL)))
        set(gcf,'color',[1,1,1])
        h2=figure;
        set(h2,'name',sprintf('Chip m%u: REPRODUCIBILITY C & C-A',ChipRank(ChipL)))
        set(gcf,'color',[1,1,1])

        CurrListNb=length(Information{ChipL}.listRank);
        LimitNb=length(Information{ChipL}.limits{1});
        Colors=colors(colormap,CurrListNb);
        CurrRowNb=ceil(sqrt(LimitNb));
        CurrColNb=round(LimitNb/CurrRowNb);
        if CurrRowNb*CurrColNb<LimitNb
            CurrColNb=CurrColNb+1;
        end
        CReprod{ChipL}=zeros(10*LimitNb,4*LimitNb);

        for LimitL=1:LimitNb
            for ListL=1:CurrListNb
                LimitPos=find(Information{ChipL}.limits{2}==Information{ChipL}.limits{1}(LimitL));
                if ~isempty(LimitPos)
                    UsedClu=[];
                    Size{1}=zeros(1,10);
                    Size{2}=zeros(1,10);
                    CRep=zeros(1,10);
                    SClu=zeros(1,10);
                    for CluL1=1:10
                        a=find(Cluster{ChipL}{1}{ListL}{1}(:,LimitL)==CluL1);
                        %find the best overlap
                        SearchedClu=0;
                        Score=0;
                        for CluL2=1:10
                            if isempty(find(UsedClu==CluL2))
                                b=find(Cluster{ChipL}{2}{ListL}{1}(:,LimitPos)==CluL2);
                                CurrScore=round(length(intersect(a,b))*100/max(length(a),length(b)));
                                if CurrScore>=Score
                                    Score=CurrScore;
                                    SearchedClu=CluL2;
                                end
                            end
                        end
                        SClu(CluL1)=SearchedClu;
                        UsedClu=[UsedClu,SearchedClu];
                        b=find(Cluster{ChipL}{2}{ListL}{1}(:,LimitPos)==SearchedClu);
                        CReprod{ChipL}(10*(ListL-1)+CluL1,4*(LimitL-1)+1)=length(a);
                        Size{1}(CluL1)=length(a);
                        CReprod{ChipL}(10*(ListL-1)+CluL1,4*(LimitL-1)+2)=SearchedClu;
                        CReprod{ChipL}(10*(ListL-1)+CluL1,4*(LimitL-1)+3)=length(b);
                        Size{2}(CluL1)=length(b);
                        CReprod{ChipL}(10*(ListL-1)+CluL1,4*(LimitL-1)+4)=round(length(intersect(a,b))*100/min(length(a),length(b)));
                        CRep(CluL1)=round(length(intersect(a,b))*100/min(length(a),length(b)));
                    end
                    figure(h1)
                    subplot(CurrRowNb,CurrColNb,LimitL)
                    hold on
                    plot([1:10],Size{1},'-','color',Colors(ListL,:))
                    plot([1:10],Size{2},':','color',Colors(ListL,:))
                    if ListL==CurrListNb
                        title(sprintf('Limit = %u',Information{ChipL}.limits{1}(LimitL)))
                        xlabel(sprintf('C-A m%u cluster rank',ChipRank(ChipL)))
                        set(gca,'xtick',[1:10])
                        set(gca,'xlim',[1,10])
                        set(gca,'xticklabel',SClu)
                        set(gca,'box','on')
                        if LimitL==LimitNb
                            try
                            legend(ClassName(Information{ChipL}.listRank))
                            catch
                                'stop'
                            end
                        end
                    end
                    figure(h2)
                    subplot(CurrRowNb,CurrColNb,LimitL)
                    hold on
                    plot([1:10],CRep,'color',Colors(ListL,:))
                    if ListL==CurrListNb
                        title(sprintf('Limit = %u',Information{ChipL}.limits{1}(LimitL)))
                        xlabel(sprintf('C-A m%u cluster rank',ChipRank(ChipL)))
                        set(gca,'xtick',[1:10])
                        set(gca,'xlim',[1,10])
                        set(gca,'xticklabel',SClu)
                        set(gca,'box','on')
                        if LimitL==LimitNb
                            legend(ClassName(Information{ChipL}.listRank))
                        end
                    end
                end
            end
        end
    end
end


%% LOAD CORESPONDANCE BETWEEN CHIPS IF NECESSARY

if DiffChipFlag
    cd(K.dir.chip)
    for ChipL=1:2
        eval(sprintf('load m%u_combinedps_corr60;',ChipRank(ChipL)))
        MatchPsRank{ChipL}=zeros(PsNb(ChipL));
        for PsL=1:length(PsRank)
            MatchPsRank{ChipL}(PsRank{PsL})=PsL;
        end
        CommonPsRank{ChipL}=PsRank;
    end
    clear PsRank
else
    for ChipL=1:2
        MatchPsRank{ChipL}=[1:PsNb(ChipL)]';
    end
end

% if  ChipRank(1)~=ChipRank(2)
%     cd(K.dir.affyMetadata)
%     FileName=sprintf('m%u_m%u_commonps.mat',min(ChipRank),max(ChipRank));
%     if exist(FileName,'file')
%         load(FileName)
%         if ChipRank(1)>ChipRank(2)
%             Temp=ComPsRank;
%             ComPsRank(:,1)=ComPsRank(:,2);
%             ComPsRank(:,2)=Temp(:,1);
%             clear Temp
%         end
%         if ~isempty(ChipRank==8)&~isempty(ChipRank==27)
%             %if min(Information{1}.listRank)>=7
%             %   MatchPsRank=ComPsRank(:,2);
%             %else
%             MatchPsRank=zeros(max(ComPsRank(:,1)),1);
%             MatchPsRank(ComPsRank(:,1))=ComPsRank(:,2);
%             %end
%         else
%             MatchPsRank=zeros(max(ComPsRank(:,1)),1);
%             MatchPsRank(ComPsRank(:,1))=ComPsRank(:,2);
%         end
%     else
%         h=errordlg(sprintf('no correspondance file between m%y and m%u',ChipRank(1),ChipRank(2)));
%         waitfor(h)
%         error('process canceled')
%     end
% else
%     ComPsRank(:,1)=[1:PsNb(1)]';
%     ComPsRank(:,2)=[1:PsNb(1)]';
%     MatchPsRank=ComPsRank;
% end

% STATISTICS ACCORDING TO TYPE, PS LIST, LIMITS ...
%% NB OF PROCESSED NETWORKS ACCORDING TO TYPE
h1=figure;
if ChipNb==1
    set(h1,'name',sprintf('NB OF PROCESSED NETWORKS IN m%u BY TYPE',ChipRank(1)))
else
    set(h1,'name',sprintf('NB OF PROCESSED NETWORKS IN m%u and m%u BY TYPE',ChipRank(1),ChipRank(2)))
end
Legend={};
for ChipL=1:ChipNb
    CluNbs={};
    hs=zeros(1,4);
    Limits=[];
    for TypeL=1:MaxType
        Limits=union(Limits,Information{ChipL}.limits{TypeL});
    end
    for TypeL=1:MaxType
        if Types(ChipL,TypeL)
            if isempty(strmatch(TypeName{TypeL},Legend,'exact'))
                Legend{end+1}=TypeName{TypeL};
            end
            set(gcf,'color',[1,1,1])
            Colors=colors(colormap,4);
            LimitNb=length(Information{ChipL}.limits{TypeL});
            NetNbs{TypeL}=[];
            if ~isempty(ClusterNb{ChipL}{TypeL})
                for ListL=1:ListNb(ChipL)
                    %try
                        if ~isempty(ClusterNb{ChipL}{TypeL}{ListL})
                            subplot(RowNb,ColNb,PlotPos{ChipL}(ListL))
                            hold on
                            CluNb=zeros(1,LimitNb);
                            CurrNetNb=zeros(LimitNb,1);;
                            for NetL=1:NetNb(ChipL)
                                for LimitL=1:LimitNb
                                    try
                                        if length(ClusterNb{ChipL}{TypeL}{ListL}{NetL}{LimitL})>0
                                            CurrNetNb(LimitL)=CurrNetNb(LimitL)+1;
                                        end
                                    catch
                                    end
                                end
                            end
                            h=plot(Information{ChipL}.limits{TypeL},CurrNetNb,'color',Colors(TypeL,:));
                            if ListL==ListNb(ChipL)
                                hs(TypeL)=h;
                            end
                            plot(Information{ChipL}.limits{TypeL},CurrNetNb,'+','color',Colors(TypeL,:))
                            NetNbs{TypeL}(ListL,:)=CurrNetNb;
                            if ListL==1
                                title(sprintf('m%u-%s, %s',ChipRank(ChipL),FileName{ChipL},ClassName{ListL}))
                            else
                                title(sprintf('%s',ClassName{ListL}))
                            end
                            xlabel('limit')
                            ylabel('net nb')
                            set(gca,'box','on')
                            set(gca,'xlim',[0,Limits(end)])
                            set(gca,'xgrid','on')
                            set(gca,'ygrid','on')
                            set(gca,'ylim',[0,NetNb(ChipL)])
                        end
                    %catch
                    %end
                end
            end
        end
    end
end
hs(find(hs==0))=[];
legend(hs,Legend)



%% NB OF CLUSTERS ACCORDING TO TYPE OF CORRELATION
h=figure;
if ChipNb==1
    set(h,'name',sprintf('NB OF CLUSTERS IN m%u BY TYPE',ChipRank(1)))
else
    set(h,'name',sprintf('NB OF CLUSTERS IN m%u and m%u BY TYPE',ChipRank(1),ChipRank(2)))
end
for ChipL=1:ChipNb
    CluNbs={};
    hs=zeros(1,4);
    for TypeL=1:MaxType
        if Types(ChipL,TypeL)
            set(gcf,'color',[1,1,1])
            Colors=colors(colormap,4);
            LimitNb=length(Information{ChipL}.limits{TypeL});
            CluNbs{TypeL}=[];
            if ~isempty(ClusterNb{ChipL}{TypeL})
                for ListL=1:ListNb(ChipL)
                    %try
                    if ~isempty(ClusterNb{ChipL}{TypeL}{ListL})
                        subplot(RowNb,ColNb,PlotPos{ChipL}(ListL))
                        hold on
                        CluNb=zeros(1,LimitNb);
                        CurrNetNb=zeros(1,LimitNb);
                        for NetL=1:NetNb(ChipL)
                            for LimitL=1:LimitNb
                                try
                                    CluNb(LimitL)=CluNb(LimitL)+length(ClusterNb{ChipL}{TypeL}{ListL}{NetL}{LimitL});
                                    if length(ClusterNb{ChipL}{TypeL}{ListL}{NetL}{LimitL})>0
                                        CurrNetNb(LimitL)=CurrNetNb(LimitL)+1;
                                    end
                                catch
                                end
                            end
                        end
                        h=plot(Information{ChipL}.limits{TypeL},CluNb./CurrNetNb,'color',Colors(TypeL,:));
                        if ListL==ListNb(ChipL)
                            hs(TypeL)=h;
                        end
                        plot(Information{ChipL}.limits{TypeL},CluNb./CurrNetNb,'+','color',Colors(TypeL,:))
                        %CluNbs{TypeL}(ListL,:)=round(CluNb./CurrNetNb);
                        CluNbs{TypeL}(ListL,:)=CluNb./CurrNetNb;
                        if ListL==1
                            title(sprintf('m%u-%s, %s',ChipRank(ChipL),FileName{ChipL},ClassName{ListL}))
                        else
                            title(sprintf('%s',ClassName{ListL}))
                        end                        
                        xlabel('limit')
                        ylabel('cluster nb')
                        set(gca,'box','on')
                        set(gca,'xgrid','on')
                        set(gca,'ygrid','on')
                    end
                    %catch
                    %end
                end
            end
        end
    end
    MaxCluNb=0;
    for TypeL=1:MaxType
        if Types(ChipL,TypeL)
            for ListL=1:ListNb(ChipL)
                MaxCluNb=max(MaxCluNb,max(CluNbs{TypeL}(ListL,:)));
            end
            subplot(RowNb,ColNb,PlotPos{ChipL}(ListL))
            set(gca,'ylim',[0,min(1500,MaxCluNb)])
        end
    end
end
hs(find(hs==0))=[];
legend(hs,Legend)

%% NB OF CLUSTERS ACCORDING TO PS LIST
for ChipL=1:ChipNb
    h=figure;
    set(h,'name',sprintf('NB OF CLUSTERS IN m%u-%s BY PS LIST',ChipRank(ChipL),FileName{ChipL}))
    CurrPlotPos=0;
    CurrTypeNb=sum(Types(ChipL,:));
    for TypeL=1:MaxType
        if Types(ChipL,TypeL)
            set(gcf,'color',[1,1,1])
            Colors=colors(colormap,ListNb(ChipL));
            %Colors=colors(colormap,6);
            LimitNb=length(Information{ChipL}.limits{TypeL});
            if TypeL==MaxType
                hs=zeros(1,ListNb(ChipL));
            end
            CurrPlotPos=CurrPlotPos+1;
            switch CurrTypeNb
                case 1
                    subplot(1,1,CurrPlotPos)
                case 2
                    subplot(1,2,CurrPlotPos)
                otherwise
                    subplot(2,2,CurrPlotPos)
            end                                                                    
            hs=zeros(1,ListNb(ChipL));
            hold on
            if ~isempty(ClusterNb{ChipL}{TypeL})               
                for ListL=1:ListNb(ChipL)
                    try
                        if ~isempty(ClusterNb{ChipL}{TypeL}{ListL})
                            CluNb=zeros(1,LimitNb);
                            CurrNetNb=zeros(1,LimitNb);
                            for NetL=1:NetNb(ChipL)
                                for LimitL=1:LimitNb
                                    try
                                        CluNb(LimitL)=CluNb(LimitL)+length(ClusterNb{ChipL}{TypeL}{ListL}{NetL}{LimitL});
                                        if length(ClusterNb{ChipL}{TypeL}{ListL}{NetL}{LimitL})>0
                                            CurrNetNb(LimitL)=CurrNetNb(LimitL)+1;
                                        end
                                    catch
                                    end
                                end
                            end
                            h=plot(Information{ChipL}.limits{TypeL},CluNb./CurrNetNb,'color',Colors(ListL,:));
                            if TypeL==MaxType
                                hs(ListL)=h;
                            end
                            plot(Information{ChipL}.limits{TypeL},CluNb./CurrNetNb,'+','color',Colors(ListL,:));
                        end
                    catch
                    end
                end
            end
            if ListL==1
                title(sprintf('m%u-%s, %s',ChipRank(ChipL),FileName{ChipL},TypeName{ListL}))
            else
                title(sprintf('%s',TypeName{ListL}))
            end
            xlabel('limit')
            ylabel('cluster nb')
            set(gca,'box','on')
            set(gca,'xgrid','on')
            set(gca,'ygrid','on')
        end
    end
    %legend(hs,ClassName)
    legend(hs,ClassName(1:ListNb(ChipL)))
    CurrPlotPos=0;
    if Types(ChipL,TypeL)
        CurrPlotPos=CurrPlotPos+1;
        switch CurrTypeNb
            case 1
                subplot(1,1,CurrPlotPos)
            case 2
                subplot(1,2,CurrPlotPos)
            otherwise
                subplot(2,2,CurrPlotPos)
        end
        set(gca,'ylim',[0,min(1500,max(max(CluNbs{TypeL})))])
    end
end


%% CLUSTER SIZE

% plot each network (1) or mean of all networks (0)
AllNetFlag=0;
% plot each cluster (0) of cumulative sum (1)
CumulFlag=1;
% plot cluster size (0) of percentage of processed probe sets (1)
PercFlag=1;
% CluRange (clusters to be processed)
CluRange=[1:50];
%merge all limits
AllLimits=[];
for ChipL=1:ChipNb
    for TypeL=1:MaxType
        AllLimits=union(AllLimits,Information{ChipL}.limits{TypeL});
    end
end
LimitNb=length(AllLimits);
Colors=colors(colormap,LimitNb);
LimitNames=cell(1,LimitNb);
for LimitL=1:LimitNb
    LimitNames{LimitL}=num2str(AllLimits(LimitL));
end


for TypeL=1:MaxType
    if ~isempty(find(Types(:,TypeL)))
        h=figure;
        if ChipRank(1)==ChipRank(2)
            set(h,'name',sprintf('CLUSTER SIZE IN m%u-%s - TYPE %s',ChipRank(ChipL),FileName{ChipL},TypeName{TypeL}))
        else
            set(h,'name',sprintf('CLUSTER SIZE IN m%u-%s and m%u-%s - TYPE %s',ChipRank(1),FileName{1},ChipRank(2),FileName{2},TypeName{TypeL}))
        end
        set(gcf,'color',[1,1,1])

        for ChipL=1:ChipNb
            LimitNb=length(Information{ChipL}.limits{TypeL});
            hs=zeros(1,LimitNb);
            if ~isempty(ClusterNb{ChipL}{TypeL})
                for ListL=1:ListNb(ChipL)
                    %for ListL=7
                    
                        if ~isempty(ClusterNb{ChipL}{TypeL}{ListL})
                            subplot(RowNb,ColNb,PlotPos{ChipL}(ListL))
                            hold on
                            for LimitL=1:LimitNb
                                LimitPos=find(AllLimits==Information{ChipL}.limits{TypeL}(LimitL));
                                Values=[];
                                for NetL=1:NetNb(ChipL)
                                    try
                                        if isempty(Values)
                                            if ~isempty(CluRange)
                                                Values=ClusterNb{ChipL}{TypeL}{ListL}{NetL}{LimitL}(1:min(length(ClusterNb{ChipL}{TypeL}{ListL}{NetL}{LimitL}),length(CluRange)))';
                                            else
                                                Values=ClusterNb{ChipL}{TypeL}{ListL}{NetL}{LimitL}';
                                            end
                                        else
                                            if isempty(CluRange)
                                                if length(ClusterNb{ChipL}{TypeL}{ListL}{NetL}{LimitL})>size(Values,2)
                                                    Values=[Values,zeros(size(Values,1),length(ClusterNb{ChipL}{TypeL}{ListL}{NetL}{LimitL})-size(Values,2))];
                                                    Values(end+1,:)=ClusterNb{ChipL}{TypeL}{ListL}{NetL}{LimitL}';
                                                elseif length(ClusterNb{ChipL}{TypeL}{ListL}{NetL}{LimitL})<size(Values,2)
                                                    Values(end+1,:)=[ClusterNb{ChipL}{TypeL}{ListL}{NetL}{LimitL}',zeros(1,size(Values,2)-length(ClusterNb{ChipL}{TypeL}{ListL}{NetL}{LimitL}))];
                                                else
                                                    Values(end+1,:)=ClusterNb{ChipL}{TypeL}{ListL}{NetL}{LimitL}';
                                                end

                                            else
                                                if length(ClusterNb{ChipL}{TypeL}{ListL}{NetL}{LimitL})>=length(CluRange)
                                                    Values(end+1,:)=ClusterNb{ChipL}{TypeL}{ListL}{NetL}{LimitL}(1:length(CluRange))';
                                                else
                                                    Values(end+1,1:length(ClusterNb{ChipL}{TypeL}{ListL}{NetL}{LimitL}))=ClusterNb{ChipL}{TypeL}{ListL}{NetL}{LimitL}';
                                                end
                                            end
                                        end
                                    catch
                                    end
                                end
                                if CumulFlag
                                    if size(Values,1)>1
                                        Values=Values';
                                        Values=cumsum(Values);
                                        Values=Values';
                                    else
                                        Values=cumsum(Values);
                                    end
                                end
                                if ~isempty(Values)
                                    if PercFlag
                                        Values=Values*100/ListSize{ChipL}(ListL);
                                    end
                                    if size(Values,1)>1
                                        if AllNetFlag==0
                                            h=plot([1:size(Values,2)],mean(Values),'color',Colors(LimitPos,:));
                                            if CumulFlag==0|~isempty(CluRange)
                                                plot([1:size(Values,2)],mean(Values),'+','color',Colors(LimitPos,:));
                                            end
                                        else
                                            for NetL=1:size(Values,1)
                                                h=plot([1:size(Values,2)],Values(NetL,:),'color',Colors(LimitPos,:));
                                            end
                                        end
                                    else
                                        h=plot([1:size(Values,2)],Values,'color',Colors(LimitPos,:));
                                        if CumulFlag==0|~isempty(CluRange)
                                            plot([1:size(Values,2)],Values,'+','color',Colors(LimitPos,:));
                                        end
                                    end
                                    if hs(LimitL)==0
                                        hs(LimitL)=h;
                                    end
                                    if size(Values,1)>1 & AllNetFlag==0
                                        plot([1:size(Values,2)],mean(Values)+std(Values),':','color',Colors(LimitPos,:));
                                        plot([1:size(Values,2)],mean(Values)-std(Values),':','color',Colors(LimitPos,:));
                                    end
                                end
                            end
                            if ListL==1
                                title(sprintf('m%u-%s, %s(%u ps)',ChipRank(ChipL),FileName{ChipL},ClassName{ListL},ListSize{ChipL}(ListL)))
                            else
                                title(sprintf('%s(%u ps)',ClassName{ListL},ListSize{ChipL}(ListL)))
                            end
                            xlabel('cluster rank')
                            if PercFlag
                                ylabel('probe set %')
                            else
                                ylabel('probe set nb')
                            end
                            set(gca,'box','on')
                            set(gca,'xgrid','on')
                            set(gca,'ygrid','on')
                            if ChipL==1
                                legend(hs,LimitNames);
                            end
                            if PercFlag & CumulFlag==1
                                set(gca,'ylim',[0,100])
                            end
                        end
                    
                end
            end
        end
    end
end

%% INTRA CHIP REPRODUCIBILITY ACCORDING TO TYPE
%ChipRange indicates the chips that have their MCL clustering made on several networks
ChipRange=[];
if NetNb(1)>1 & NetNb(2)>1
    ChipRange=[1,2];
elseif NetNb(1)>1
    ChipRange=1;
elseif NetNb(2)>1
    ChipRange=2;
end

if ~isempty(ChipRange)
    h=figure;
    CLU=1;
    set(h,'name',sprintf('INTRA m%u-%s CHIP REPRODUCIBILITY OF CLU %u ACCORDING TO CORR TYPE',ChipRank(ChipL),FileName{ChipL},CLU))
    Colors=colors(colormap,4);
    for ChipL=ChipRange
        hs=zeros(1,4);
        Limits=[];
        for TypeL=1:MaxType
            Limits=union(Limits,Information{ChipL}.limits{TypeL});
        end
        NetFreq{ChipL}=cell(1,4);
        for TypeL=1:MaxType
            if Types(ChipL,TypeL)
                %for TypeL=2
                CurrLimits=Information{ChipL}.limits{TypeL};
                set(gcf,'color',[1,1,1])
                LimitNb=length(Information{ChipL}.limits{TypeL});
                NetFreq{ChipL}{TypeL}=cell(1,ListNb(ChipL));
                if ~isempty(ClusterNb{ChipL}{TypeL})
                    for ListL=1:ListNb(ChipL)
                        %for ListL=1
                        NetFreq{ChipL}{TypeL}{ListL}=zeros(LimitNb,NetNb(ChipL));
                        if ~isempty(ClusterNb{ChipL}{TypeL}{ListL})
                            subplot(RowNb,ColNb,PlotPos{ChipL}(ListL))
                            hold on
                            Reprod=zeros(NetNb(ChipL)*(NetNb(ChipL)-1)/2,LimitNb);
                            for LimitL=1:LimitNb
                                ResPos=0;
                                PsMat=zeros(PsNb(ChipL),NetNb(ChipL));
                                CluPos=zeros(1,NetNb(ChipL));
                                %find correspondance (can be different in different limit; explain
                                %that some plot are not smooth (local v or ^ shape between
                                %three successive points)
                                FoundNet=0;
                                for NetL=1:NetNb(ChipL)
                                    try
                                        Ps1=find(Cluster{ChipL}{TypeL}{ListL}{NetL}(:,LimitL)==CLU);
                                        if length(Ps1)>0
                                            CluPos(NetL)=CLU;
                                            FoundNet=1;
                                            break
                                        end
                                    catch
                                    end
                                end
                                if FoundNet
                                    for NetL=1:NetNb(ChipL)
                                        if CluPos(NetL)==0
                                            if size(Cluster{ChipL}{TypeL}{ListL}{NetL},2)>=LimitL
                                                if ~isempty(find(Cluster{ChipL}{TypeL}{ListL}{NetL}(:,LimitL))>0)
                                                    InterSize=0;
                                                    for CluL=1:10
                                                        %assumes that CLU=1 to 5 or 6 and that the
                                                        %corresponding cluster in less than 10
                                                        Ps2=find(Cluster{ChipL}{TypeL}{ListL}{NetL}(:,LimitL)==CluL);
                                                        if InterSize<length(intersect(Ps1,Ps2))
                                                            InterSize=length(intersect(Ps1,Ps2));
                                                            CluPos(NetL)=CluL;
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                                for NetL1=1:NetNb(ChipL)-1
                                    Ps1=[];
                                    try
                                        if CluPos(NetL1)>0
                                            Ps1=find(Cluster{ChipL}{TypeL}{ListL}{NetL1}(:,LimitL)==CluPos(NetL1));
                                            PsMat(Ps1,NetL1)=1;
                                        end
                                    catch
                                    end
                                    for NetL2=NetL1+1:NetNb(ChipL)
                                        ResPos=ResPos+1;
                                        Ps2=[];
                                        try
                                            if CluPos(NetL2)>0
                                                Ps2=find(Cluster{ChipL}{TypeL}{ListL}{NetL2}(:,LimitL)==CluPos(NetL2));
                                                if length(Ps1)>0&length(Ps2>0)
                                                    Reprod(ResPos,LimitL)=length(intersect(Ps1,Ps2))*100/min(length(Ps1),length(Ps2));
                                                end
                                                if NetL2==NetNb(ChipL)&length(Ps2>0)
                                                    PsMat(Ps2,NetL2)=1;
                                                end
                                            end
                                        catch
                                        end
                                    end
                                end
                                PsMat=sum(PsMat,2);
                                PsMat(find(PsMat==0))=[];
                                if ~isempty(PsMat)
                                    NetFreq{ChipL}{TypeL}{ListL}(LimitL,:)=histc(PsMat,[1:NetNb(ChipL)])';
                                end
                            end

                            MeanReprod=zeros(1,LimitNb);
                            StdReprod=zeros(1,LimitNb);
                            for LimitL=1:LimitNb
                                MeanReprod(LimitL)=mean(Reprod(find(Reprod(:,LimitL)>0),LimitL));
                                StdReprod(LimitL)=std(Reprod(find(Reprod(:,LimitL)>0),LimitL));
                                if isnan(MeanReprod(LimitL))
                                    MeanReprod(LimitL)=0;
                                    StdReprod(LimitL)=0;
                                end
                            end
                            h=plot(CurrLimits,MeanReprod,'-+','color',Colors(TypeL,:));
                            plot(CurrLimits,MeanReprod+StdReprod,':','color',Colors(TypeL,:))
                            plot(CurrLimits,MeanReprod-StdReprod,':','color',Colors(TypeL,:))
                            if hs(TypeL)==0
                                hs(TypeL)=h;
                            end
                            if ListL==1
                                title(sprintf('m%u-%s, %s',ChipRank(ChipL),FileName{ChipL},ClassName{ListL}))
                            else
                                title(sprintf('%s',ClassName{ListL}))
                            end                        
                            xlabel('limit')
                            ylabel('reproducibility')
                            set(gca,'box','on')
                            set(gca,'xlim',[0,Limits(end)])
                            set(gca,'xtick',Limits)
                            set(gca,'ylim',[0,100])
                            set(gca,'xgrid','on')
                            set(gca,'ygrid','on')
                        end
                    end
                end
            end
        end
        legend(hs,TypeName)
    end

    for TypeL=1:MaxType
        if ~isempty(find(Types(:,TypeL)))
            h=figure;
            set(gcf,'color',[1,1,1])
            set(h,'name',sprintf('%s: INTRA m%u-%s CHIP NETWORK REPRODUCIBILITY OF CLU %u ACCORDING TO TYPE %s',TypeName{TypeL},ChipRank(ChipL),FileName{ChipL},CLU,TypeName{TypeL}))
            for ChipL=ChipRange
                if Types(ChipL,TypeL)
                    hs=zeros(1,length(Limits));
                    LimitNb=length(Information{ChipL}.limits{TypeL});
                    LimitNames=cell(1,LimitNb);
                    for LimitL=1:LimitNb
                        LimitNames{LimitL}=num2str(Information{ChipL}.limits{TypeL}(LimitL));
                    end
                    Colors=colors(colormap,LimitNb);
                    for ListL=1:ListNb(ChipL)
                        subplot(RowNb,ColNb,PlotPos{ChipL}(ListL))
                        hold on
                        for LimitL=1:LimitNb
                            try
                                h=plot(1:NetNb(ChipL),cumsum(NetFreq{ChipL}{TypeL}{ListL}(LimitL,:)*100/sum(NetFreq{ChipL}{TypeL}{ListL}(LimitL,:))),'-+','color',Colors(LimitL,:));
                            catch
                            end
                        end
                        if ListL==1
                            title(sprintf('m%u-%s, %s',ChipRank(ChipL),FileName{ChipL},ClassName{ListL}))
                        else
                            title(sprintf('%s',ClassName{ListL}))
                        end                                               
                        xlabel('net nb')
                        ylabel('frequency')
                        set(gca,'box','on')
                        set(gca,'ylim',[0,100])
                        set(gca,'xlim',[1,NetNb(ChipL)])
                        set(gca,'xtick',[1:NetNb(ChipL)])
                        %set(gca,'xgrid','on')
                        %set(gca,'ygrid','on')
                    end
                    legend(LimitNames)
                end
            end
        end
    end



    LIMIT=2;
    h=figure;
    set(gcf,'color',[1,1,1])
    set(h,'name',sprintf('INTRA m%u CHIP NETWORK REPRODUCIBILITY OF CLU %u ACCORDING TO PS LIST TYPE FOR LIMIT %u',ChipRank(ChipL),CLU,Information{ChipL}.limits{TypeL}(LIMIT)))

    ChipL=1;
    Colors=colors(colormap,6);
    for TypeL=1:MaxType
        if ~isemtpy(Types(1,TypeL))
            subplot(ChipNb,MaxType,TypeL)
            hold on

            for ListL=1:ListNb(ChipL)

                try
                    plot(1:NetNb(ChipL),cumsum(NetFreq{ChipL}{TypeL}{ListL}(LIMIT,:)*100/sum(NetFreq{ChipL}{TypeL}{ListL}(LIMIT,:))),'-+','color',Colors(ListL,:));
                catch
                end

                title(TypeName{TypeL})
                xlabel('net nb')
                ylabel('frequency')
                set(gca,'box','on')
                set(gca,'ylim',[0,100])
                set(gca,'xlim',[1,NetNb(ChipL)])
                set(gca,'xtick',[1:NetNb(ChipL)])
            end
        end

        %set(gca,'xgrid','on')
        %set(gca,'ygrid','on')
    end
    legend(ClassName(1:ListNb(ChipL)))

    % for TypeL=1:MaxType
    %     subplot(3,4,4+TypeL)
    %     hold on
    %
    %     for ListL=7:12
    %
    %         try
    %             plot(1:NetNb(ChipL),cumsum(NetFreq{ChipL}{TypeL}{ListL}(LIMIT,:)*100/sum(NetFreq{ChipL}{TypeL}{ListL}(LIMIT,:))),'-+','color',Colors(ListL-6,:));
    %         catch
    %         end
    %
    %         title(TypeName{TypeL})
    %         xlabel('net nb')
    %         ylabel('frequency')
    %         set(gca,'box','on')
    %         set(gca,'ylim',[0,100])
    %         set(gca,'xlim',[1,NetNb(ChipL)])
    %         set(gca,'xtick',[1:NetNb(ChipL)])
    %     end
    %
    %     %set(gca,'xgrid','on')
    %     %set(gca,'ygrid','on')
    % end
    % legend(ClassName(7:12))
    if ChipNb==2
        ChipL=2;
        for TypeL=1:MaxType
            if ~isempty(Types(2,ChipL))
                subplot(2,MaxType,8+TypeL)
                hold on

                for ListL=1:ListNb(ChipL)

                    try
                        plot(1:NetNb(ChipL),cumsum(NetFreq{ChipL}{TypeL}{ListL}(LIMIT,:)*100/sum(NetFreq{ChipL}{TypeL}{ListL}(LIMIT,:))),'-+','color',Colors(ListL,:));
                    catch
                    end

                    title(TypeName{TypeL})
                    xlabel('net nb')
                    ylabel('frequency')
                    set(gca,'box','on')
                    set(gca,'ylim',[0,100])
                    set(gca,'xlim',[1,NetNb(ChipL)])
                    set(gca,'xtick',[1:NetNb(ChipL)])
                end
            end

            %set(gca,'xgrid','on')
            %set(gca,'ygrid','on')
        end
        legend(ClassName(1:ListNb(ChipL)))
    end



    %% REPRODUCIBILITY MATRIX

    for ChipL=1:ChipNb
        PlotNb=NetNb(ChipL)*(NetNb(ChipL)-1)/2;
        ColNbMat=round(sqrt(PlotNb/2));
        RowNbMat=round(PlotNb/ColNbMat);
        if PlotNb>RowNbMat*ColNbMat
            RowNbMat=RowNbMat+1;
        end
        for TypeL=1:MaxType
            if Types(ChipL,TypeL)
                for ListL=1:ListNb(ChipL)
                    h=figure;
                    set(h,'name',sprintf('INTRA m%u CHIP REPRODUCIBILITY MATRIX ACCORDING TO TYPE %s AND PS LIST %s',ChipRank(ChipL),TypeName{TypeL},ClassName{ListL}))
                    set(gcf,'color',[1,1,1])
                    if ~isempty(ClusterNb{ChipL}{TypeL})
                        if ~isempty(ClusterNb{ChipL}{TypeL}{ListL})
                            LimitL=1;
                            NetPos=0;
                            for NetL1=1:NetNb(ChipL)-1
                                for NetL2=NetL1+1:NetNb(ChipL)
                                    NetPos=NetPos+1;
                                    try
                                        subplot(RowNbMat,ColNbMat,NetPos)
                                        SimMcl=zeros(11,11);
                                        for CluL1=1:10
                                            Ps1=find(Cluster{ChipL}{TypeL}{ListL}{NetL1}(:,LimitL)==CluL1);
                                            for CluL2=1:10
                                                Ps2=find(Cluster{ChipL}{TypeL}{ListL}{NetL2}(:,LimitL)==CluL2);
                                                SimMcl(CluL1,CluL2)=length(intersect(Ps1,Ps2))*100/length(union(Ps1,Ps2));
                                            end
                                        end
                                        h=pcolor(SimMcl);
                                        set(h,'linestyle','none')
                                        set(gca,'xticklabel','')
                                        set(gca,'yticklabel','')
                                        xlabel(sprintf('n%u',Information{ChipL}.netRanks(NetL1)))
                                        ylabel(sprintf('n%u',Information{ChipL}.netRanks(NetL2)))
                                    catch
                                    end
                                end
                            end
                        end
                    end
                end
            end %of ListL
        end %of TypeL
    end





%% INTRA CHIP REPRODUCIBILITY ACCORDING TO PS LIST

    for ChipL=1:ChipNb
        h=figure;
        CLU=1;
        set(h,'name',sprintf('INTRA m%u CHIP REPRODUCIBILITY OF CLU %u ACCORDING TO PS LIST',ChipRank(ChipL),CLU))
        Colors=colors(colormap,6);
        hs=zeros(1,ListNb(ChipL));
        Limits=[];
        for TypeL=1:MaxType
            if Types(ChipL,TypeL)
                Limits=union(Limits,Information{ChipL}.limits{TypeL});
            end
        end
        for TypeL=1:MaxType
            if Types(ChipL,TypeL)
                CurrLimits=Information{ChipL}.limits{TypeL};
                set(gcf,'color',[1,1,1])
                LimitNb=length(Information{ChipL}.limits{TypeL});
                if MaxType==2
                    subplot(1,2,TypeL)
                elseif MaxType>2    
                    subplot(2,2,TypeL)
                end
                hold on
                if ~isempty(ClusterNb{ChipL}{TypeL})
                    for ListL=1:ListNb(ChipL)
                        if ~isempty(ClusterNb{ChipL}{TypeL}{ListL})
                            Reprod=zeros(NetNb(ChipL)*(NetNb(ChipL)-1)/2,LimitNb);
                            for LimitL=1:LimitNb
                                ResPos=0;
                                FoundNet=0;
                                CluPos=zeros(1,NetNb(ChipL));
                                for NetL=1:NetNb(ChipL)
                                    Ps1=find(Cluster{ChipL}{TypeL}{ListL}{NetL}(:,LimitL)==CLU);
                                    if length(Ps1)>0
                                        CluPos(NetL)=CLU;
                                        FoundNet=1;
                                        break
                                    end
                                end
                                if FoundNet
                                    for NetL=1:NetNb(ChipL)
                                        if CluPos(NetL)==0
                                            if size(Cluster{ChipL}{TypeL}{ListL}{NetL},2)>=LimitL
                                                if ~isempty(find(Cluster{ChipL}{TypeL}{ListL}{NetL}(:,LimitL))>0)
                                                    InterSize=0;
                                                    for CluL=1:10
                                                        %assumes that CLU=1 to 5 or 6 and that the
                                                        %corresponding cluster in less than 10
                                                        Ps2=find(Cluster{ChipL}{TypeL}{ListL}{NetL}(:,LimitL)==CluL);
                                                        if InterSize<length(intersect(Ps1,Ps2))
                                                            InterSize=length(intersect(Ps1,Ps2));
                                                            CluPos(NetL)=CluL;
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                                for NetL1=1:NetNb(ChipL)-1
                                    Ps1=[];
                                    try
                                        if CluPos(NetL1)>0
                                            Ps1=find(Cluster{ChipL}{TypeL}{ListL}{NetL1}(:,LimitL)==CluPos(NetL1));
                                        end
                                    catch
                                    end
                                    for NetL2=NetL1+1:NetNb(ChipL)
                                        ResPos=ResPos+1;
                                        Ps2=[];
                                        try
                                            if CluPos(NetL2)>0
                                                Ps2=find(Cluster{ChipL}{TypeL}{ListL}{NetL2}(:,LimitL)==CluPos(NetL2));
                                                if length(Ps1)>0&length(Ps2>0)
                                                    Reprod(ResPos,LimitL)=length(intersect(Ps1,Ps2))*100/min(length(Ps1),length(Ps2));
                                                end
                                            end
                                        catch
                                        end
                                    end
                                end
                            end
                            %eliminate empty lines
                            MeanReprod=zeros(1,LimitNb);
                            StdReprod=zeros(1,LimitNb);
                            for LimitL=1:LimitNb
                                MeanReprod(LimitL)=mean(Reprod(find(Reprod(:,LimitL)>0),LimitL));
                                StdReprod(LimitL)=std(Reprod(find(Reprod(:,LimitL)>0),LimitL));
                                if isnan(MeanReprod(LimitL))
                                    MeanReprod(LimitL)=0;
                                    StdReprod(LimitL)=0;
                                end
                            end
                            h=plot(CurrLimits,MeanReprod,'-+','color',Colors(ListL,:));
                            if hs(ListL)==0
                                hs(ListL)=h;
                            end
                            plot(CurrLimits,MeanReprod+StdReprod,':','color',Colors(ListL,:))
                            plot(CurrLimits,MeanReprod-StdReprod,':','color',Colors(ListL,:))

                        end
                    end
                end
                title(sprintf('%s',TypeName{TypeL}))
                xlabel('limit')
                ylabel('reproducibility')
                set(gca,'box','on')
                set(gca,'xlim',[0,Limits(end)])
                set(gca,'xtick',Limits)
                set(gca,'ylim',[0,100])

                legend(hs,ClassName(1:ListNb(ChipL)))
            end
        end %TypeL

        for ListL=1:ListNb(ChipL)
            h=figure;
            set(gcf,'color',[1,1,1])
            set(h,'name',sprintf('INTRA m%u CHIP NETWORK REPRODUCIBILITY OF CLU %u ACCORDING TO PS LIST %s',ChipRank(ChipL),CLU,ClassName{ListL}))
            for TypeL=1:MaxType
                if Types(ChipL,TypeL)
                    LimitNb=length(Information{ChipL}.limits{TypeL});
                    LimitNames=cell(1,LimitNb);
                    for LimitL=1:LimitNb
                        LimitNames{LimitL}=num2str(Information{ChipL}.limits{TypeL}(LimitL));
                    end
                    Colors=colors(colormap,LimitNb);

                    if MaxType==2
                        subplot(1,2,TypeL)
                    elseif MaxType>2
                        subplot(2,2,TypeL)
                    end
                    hold on
                    for LimitL=1:LimitNb
                        try
                            h=plot(1:NetNb(ChipL),cumsum(NetFreq{ChipL}{TypeL}{ListL}(LimitL,:)*100/sum(NetFreq{ChipL}{TypeL}{ListL}(LimitL,:))),'-+','color',Colors(LimitL,:));
                        catch
                        end
                    end

                    title(TypeName{TypeL})
                    xlabel('net nb')
                    ylabel('frequency')
                    set(gca,'box','on')
                    set(gca,'ylim',[0,100])
                    set(gca,'xtick',[1:NetNb(ChipL)])
                    set(gca,'xgrid','on')
                    set(gca,'ygrid','on')
                    legend(LimitNames)
                end
            end
        end

    end %ChipL


%% CONSENSUS REGIONS
    %constructs union of the same region in all the networks at different level of
    %reproducibility


    for ChipL=1:ChipNb
        Region{ChipL}={};
        for TypeL=1:MaxType
            if Types(ChipL,TypeL)
                Region{ChipL}{TypeL}={};
                for ListL=1:ListNb(ChipL)
                    Region{ChipL}{TypeL}{ListL}={};
                    for RegL=1:ListNb(ChipL)
                        for LimitL=1:2
                            Region{ChipL}{TypeL}{ListL}{RegL}{LimitL}.ps=cell(1,10);
                            if ListL==1 & LimitL==1
                                CurrReg=RegL;
                            else
                                %search corresponding cluster
                                Ps1=Region{ChipL}{TypeL}{1}{RegL}{1}.ps{1};
                                CurrReg=0;
                                for NetL=1:NetNb(ChipL)
                                    if size(Cluster{ChipL}{TypeL}{ListL}{NetL},2)>=LimitL
                                        if ~isempty(find(Cluster{ChipL}{TypeL}{ListL}{NetL}(:,LimitL))>0)
                                            InterSize=0;
                                            for CluL=1:10
                                                %assumes that CLU=1 to 5 or 6 and that the
                                                %corresponding cluster in less than 10
                                                Ps2=find(Cluster{ChipL}{TypeL}{ListL}{NetL}(:,LimitL)==CluL);
                                                if InterSize<length(intersect(Ps1,Ps2))
                                                    InterSize=length(intersect(Ps1,Ps2));
                                                    CurrReg=CluL;
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            Region{ChipL}{TypeL}{ListL}{RegL}{LimitL}.clu=CurrReg;
                            %search the first network
                            FoundNet=0;
                            CluPos=zeros(1,NetNb(ChipL));
                            for NetL=1:NetNb(ChipL)
                                Ps1=find(Cluster{ChipL}{TypeL}{ListL}{NetL}(:,LimitL)==CurrReg);
                                if length(Ps1)>0
                                    CluPos(NetL)=CurrReg;
                                    FoundNet=1;
                                    break
                                end
                            end
                            if FoundNet
                                for NetL=1:NetNb(ChipL)
                                    if CluPos(NetL)==0
                                        if size(Cluster{ChipL}{TypeL}{ListL}{NetL},2)>=LimitL
                                            if ~isempty(find(Cluster{ChipL}{TypeL}{ListL}{NetL}(:,LimitL))>0)
                                                InterSize=0;
                                                for CluL=1:10
                                                    %assumes that CLU=1 to 5 or 6 and that the
                                                    %corresponding cluster in less than 10
                                                    Ps2=find(Cluster{ChipL}{TypeL}{ListL}{NetL}(:,LimitL)==CluL);
                                                    if InterSize<length(intersect(Ps1,Ps2))
                                                        InterSize=length(intersect(Ps1,Ps2));
                                                        CluPos(NetL)=CluL;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            PsMat=zeros(PsNb(ChipL),NetNb(ChipL));
                            for NetL1=1:NetNb(ChipL)-1
                                Ps1=[];
                                try
                                    if CluPos(NetL1)>0
                                        Ps1=find(Cluster{ChipL}{TypeL}{ListL}{NetL1}(:,LimitL)==CluPos(NetL1));
                                        PsMat(Ps1,NetL1)=1;
                                    end
                                catch
                                end
                                for NetL2=NetL1+1:NetNb(ChipL)
                                    Ps2=[];
                                    try
                                        if CluPos(NetL2)>0
                                            Ps2=find(Cluster{ChipL}{TypeL}{ListL}{NetL2}(:,LimitL)==CluPos(NetL2));
                                            if length(Ps1)>0&length(Ps2>0)
                                            end
                                            if NetL2==NetNb(ChipL)&length(Ps2>0)
                                                PsMat(Ps2,NetL2)=1;
                                            end
                                        end
                                    catch
                                    end
                                end
                            end
                            PsMat=sum(PsMat,2);
                            for RepL=1:10
                                Region{ChipL}{TypeL}{ListL}{RegL}{LimitL}.ps{RepL}=find(PsMat>=RepL);
                            end

                        end
                    end
                end
            end
        end
    end

    % REGIONS 1 TO 6
    for RegL=1:ListNb(ChipL)
        %REGION=1;
        LIMIT=1;
        MAXREP=10;
        Colors=colors(colormap,4);
        h=figure;
        set(h,'name',sprintf('SIZE OF REGION %u ACCORDING TO UNION THRESHOLD FOR LIMIT %u',RegL,LIMIT));
        set(gcf,'color',[1,1,1])
        TypeName1={'A all','C-A all','rawC all','raw(C-A) all','A core','C-A core','rawC core','raw(C-A) core'}
        hs=zeros(1,8);
        for TypeL=1:MaxType
            for ChipL=1:ChipNb
                if Types(ChipL,TypeL)
                    for ListL=1:ListNb(ChipL)
                        subplot(RowNb,ColNb,PlotPos{ChipL}(ListL))
                        hold on
                        RegSize=zeros(1,MAXREP);
                        RegCoreSize=zeros(1,MAXREP);
                        for RepL=1:MAXREP
                            RegSize(RepL)=length(Region{ChipL}{TypeL}{ListL}{RegL}{LIMIT}.ps{RepL});
                            RegCoreSize(RepL)=length(Region{ChipL}{TypeL}{ListL}{RegL}{LIMIT}.core{RepL});
                        end
                        plot([1:MAXREP],RegSize,'color',Colors(TypeL,:))
                        h=plot([1:MAXREP],RegSize,'+','color',Colors(TypeL,:));
                        if hs(TypeL)==0
                            hs(TypeL)=h;
                        end
                        plot([1:MAXREP],RegCoreSize,':','color',Colors(TypeL,:))
                        h=plot([1:MAXREP],RegCoreSize,'o','color',Colors(TypeL,:));
                        if hs(TypeL+4)==0
                            hs(TypeL+4)=h;
                        end
                        title(sprintf('%s',ClassName{ListL}))
                        ylabel('#ps')
                        xlabel('>= #networks')
                        set(gca,'box','on')
                        set(gca,'xgrid','on')
                        set(gca,'ygrid','on')
                    end
                end
            end
        end
        legend(hs,TypeName1)
    end

    % RELIQUAT
    LIMIT=1;
    MAXREP=10;
    Colors=colors(colormap,4);
    h=figure;
    set(h,'name',sprintf('SIZE OF RELIQUAT ACCORDING TO UNION THRESHOLD FOR LIMIT %u',RegL,LIMIT));
    set(gcf,'color',[1,1,1])
    TypeName1={'A all','C-A all','rawC all','raw(C-A) all','A core','C-A core','rawC core','raw(C-A) core'};
    hs=zeros(1,8);
    for TypeL=1:MaxType
        for ChipL=1:ChipNb
            if Types(ChipL,TypeL)
                for ListL=1:ListNb(ChipL)
                    subplot(RowNb,ColNb,PlotPos{ChipL}(ListL))
                    hold on
                    ReliquatSize=zeros(1,MAXREP);
                    ReliquatCoreSize=zeros(1,MAXREP);
                    for RepL=1:MAXREP
                        RegPs=[];
                        RegCorePs=[];
                        for RegL=1:ListNb(ChipL)
                            RegPs=union(RegPs,Region{ChipL}{TypeL}{ListL}{RegL}{LIMIT}.ps{RepL});
                            RegCorePs=union(RegCorePs,Region{ChipL}{TypeL}{ListL}{RegL}{LIMIT}.core{RepL});
                        end
                        if ChipL==1
                            if ListL<7
                                ReliquatPs=setdiff(PsRankList{ChipL}{1},RegPs);
                                ReliquatCorePs=setdiff(PsRankList{ChipL}{1},RegCorePs);
                            else
                                ReliquatPs=setdiff(PsRankList{ChipL}{7},RegPs);
                                ReliquatCorePs=setdiff(PsRankList{ChipL}{7},RegCorePs);
                            end
                        else
                            ReliquatPs=setdiff(PsRankList{ChipL}{1},RegPs);
                            ReliquatCorePs=setdiff(PsRankList{ChipL}{1},RegCorePs);
                        end
                        Region{ChipL}{TypeL}{ListL}{7}{LIMIT}.ps{RepL}=ReliquatPs;
                        Region{ChipL}{TypeL}{ListL}{7}{LIMIT}.core{RepL}=ReliquatCorePs;
                        ReliquatSize(RepL)=length(ReliquatPs);
                        ReliquatCoreSize(RepL)=length(ReliquatCorePs);
                    end
                    plot([1:MAXREP],ReliquatSize,'color',Colors(TypeL,:))
                    h=plot([1:MAXREP],ReliquatSize,'+','color',Colors(TypeL,:));
                    if hs(TypeL)==0
                        hs(TypeL)=h;
                    end
                    plot([1:MAXREP],ReliquatCoreSize,':','color',Colors(TypeL,:))
                    h=plot([1:MAXREP],ReliquatCoreSize,'o','color',Colors(TypeL,:));
                    if hs(TypeL+4)==0
                        hs(TypeL+4)=h;
                    end
                    title(sprintf('%s',ClassName{ListL}))
                    ylabel('#ps')
                    xlabel('>= #networks')
                    set(gca,'box','on')
                    set(gca,'xgrid','on')
                    set(gca,'ygrid','on')
                end
            end
        end
    end
    legend(hs,TypeName1)




%% OVERLAP BETWEEN CONSENSUS REGIONS
    for ChipL=1:ChipNb
        for TypeL=1:MaxType
            if Types(ChipL,TypeL)
                for ListL=1:ListNb(ChipL)
                    for RegL1=1:ListNb(ChipL)
                        for LimitL=1:2
                            Region{ChipL}{TypeL}{ListL}{RegL1}{LimitL}.core=cell(1,10);
                            for RepL=1:10
                                CurrPs=Region{ChipL}{TypeL}{ListL}{RegL1}{LimitL}.ps{RepL};
                                %find intersection
                                InterReg=[];
                                for RegL2=1:ListNb(ChipL)
                                    if RegL2~=RegL1
                                        try
                                            InterReg=union(InterReg,intersect(CurrPs,Region{ChipL}{TypeL}{ListL}{RegL2}{LimitL}.ps{RepL}));
                                        catch
                                        end
                                    end
                                end
                                Region{ChipL}{TypeL}{ListL}{RegL1}{LimitL}.core{RepL}=setdiff(Region{ChipL}{TypeL}{ListL}{RegL1}{LimitL}.ps{RepL},InterReg);
                            end
                        end
                    end
                end
            end
        end
    end

    REGION=4;
    LIMIT=1;
    MAXREP=10;
    Colors=colors(colormap,4);
    h=figure;
    set(h,'name',sprintf('OVERLAP BETWEEN REGION %u AND OTHER REGIONS FOR LIMIT %u',REGION,LIMIT));
    set(gcf,'color',[1,1,1])
    hs=zeros(1,4);
    for TypeL=1:MaxType
        for ChipL=1:ChipNb
            if Types(ChipL,TypeL)
                for ListL=1:ListNb(ChipL)
                    subplot(RowNb,ColNb,PlotPos{ChipL}(ListL))
                    hold on
                    Overlap=zeros(1,MAXREP);
                    for RepL=1:MAXREP
                        Overlap(RepL)=100-(length(Region{ChipL}{TypeL}{ListL}{REGION}{LIMIT}.core{RepL})*100/length(Region{ChipL}{TypeL}{ListL}{REGION}{LIMIT}.ps{RepL}));
                    end
                    h=plot([1:MAXREP],Overlap,'color',Colors(TypeL,:));
                    if hs(TypeL)==0
                        hs(TypeL)=h;
                    end
                    plot([1:MAXREP],Overlap,'+','color',Colors(TypeL,:))
                    title(sprintf('%s',ClassName{ListL}))
                    ylabel('Overlap')
                    xlabel('>= #networks')
                    set(gca,'box','on')
                    set(gca,'xgrid','on')
                    set(gca,'ygrid','on')
                end
            end
        end
    end
    legend(hs,TypeName)



%% INTRA CHIP REPRODUCIBILITY

    REGION=1;
    LIMIT=1;
    MAXREP=10;
    OVERLAPFLAG=1;
    Colors=colors(colormap,4);
    h=figure;
    set(h,'name',sprintf('INTRA REPRODUCIBILITY OF REGION %u ACCORDING TO UNION THRESHOLD FOR LIMIT %u',REGION,LIMIT));
    set(gcf,'color',[1,1,1])
    hs=zeros(1,4);
    for ChipL=1:ChipNb
        for TypeL=1:MaxType
            if Types(ChipL,TypeL)
                LoopNb=1;
                ColNb=2;
                if ListNb(1)==12
                    ColNb=3;
                    if ChipL==1
                        LoopNb=2;
                    end
                end
                for LoopL=1:LoopNb
                    for ListL=2:6
                        if  ChipL==1
                            if LoopL==1
                                RepPlotPos=[0,1,4,7,10,13];
                                Offset=0;
                            else
                                RepPlotPos=[0,2,5,8,11,14];
                                Offset=6;
                            end
                        else
                            if ListNb(1)==12
                                RepPlotPos=[0,3,6,9,12,15];
                            else
                                RepPlotPos=[0,2,5,8,11,14];
                            end
                            Offset=0;
                        end
                        subplot(8,ColNb,RepPlotPos(ListL))
                        hold on
                        RegRep=zeros(1,MAXREP);
                        for RepL=1:MAXREP
                            if OVERLAPFLAG
                                RegRep(RepL)=length(intersect(Region{ChipL}{TypeL}{ListL+Offset}{REGION}{LIMIT}.core{RepL},Region{ChipL}{TypeL}{1}{REGION}{LIMIT}.core{RepL}))*100/min(length(Region{ChipL}{TypeL}{ListL+Offset}{REGION}{LIMIT}.core{RepL}),length(Region{ChipL}{TypeL}{1}{REGION}{LIMIT}.core{RepL}));
                            else
                                RegRep(RepL)=length(intersect(Region{ChipL}{TypeL}{ListL+Offset}{REGION}{LIMIT}.ps{RepL},Region{ChipL}{TypeL}{1}{REGION}{LIMIT}.ps{RepL}))*100/min(length(Region{ChipL}{TypeL}{ListL+Offset}{REGION}{LIMIT}.ps{RepL}),length(Region{ChipL}{TypeL}{1}{REGION}{LIMIT}.ps{RepL}));
                            end
                        end
                        try
                            h=plot([1:MAXREP],RegRep,'color',Colors(TypeL,:));
                            if hs(TypeL)==0
                                hs(TypeL)=h;
                            end
                            plot([1:MAXREP],RegRep,'+','color',Colors(TypeL,:))
                        catch
                        end
                        title(sprintf('%s vs %s',ClassName{ListL+Offset},ClassName{1}))
                        ylabel('reproducibility')
                        xlabel('>= #networks')
                        set(gca,'box','on')
                        set(gca,'xgrid','on')
                        set(gca,'ygrid','on')
                    end
                end
            end

            for LoopL=1:LoopNb
                for ListL=2:4
                    if  ChipL==1
                        if LoopL==1
                            RepPlotPos=[0,16,19,22];
                            Offset=0;
                        else
                            RepPlotPos=[0,17,20,23];
                            Offset=6;
                        end
                    else
                        if ListNb(1)==12
                            RepPlotPos=[0,18,21,24];
                        else
                            RepPlotPos=[0,17,20,23];
                        end
                        Offset=0;
                    end
                    subplot(8,ColNb,RepPlotPos(ListL))
                    hold on
                    RegRep=zeros(1,MAXREP);
                    for RepL=1:MAXREP
                        RegRep(RepL)=length(intersect(Region{ChipL}{TypeL}{ListL+Offset}{REGION}{LIMIT}.ps{RepL},Region{ChipL}{TypeL}{5}{REGION}{LIMIT}.ps{RepL}))*100/min(length(Region{ChipL}{TypeL}{ListL+Offset}{REGION}{LIMIT}.ps{RepL}),length(Region{ChipL}{TypeL}{5}{REGION}{LIMIT}.ps{RepL}));
                    end
                    try
                        h=plot([1:MAXREP],RegRep,'color',Colors(TypeL,:));
                        if hs(TypeL)==0
                            hs(TypeL)=h;
                        end
                        plot([1:MAXREP],RegRep,'+','color',Colors(TypeL,:))
                    catch
                    end
                    title(sprintf('%s vs %s',ClassName{ListL+Offset},ClassName{5}))
                    ylabel('reproducibility')
                    xlabel('>= #networks')
                    set(gca,'box','on')
                    set(gca,'xgrid','on')
                    set(gca,'ygrid','on')
                end
            end
        end
    end
    legend(hs,TypeName)





%% INTER CHIP REPRODUCIBILITY    (EXIST OVERLAP BETWEENS REGIONS CONSTRUCTED BY UNION OF SEVERAL NETWORKS)

    REGION1=1;
    REGION2=1;

    LIMIT=1;
    MAXREP=10;
    SIMREGFLAG=1;
    OVERLAPFLAG=1;

    for REGION1=7
        REGION2=REGION1;


        Colors=colors(colormap,4);
        h=figure;
        set(h,'name',sprintf('INTER CHIP REPRODUCIBILITY OF REGION %u ACCORDING TO UNION THRESHOLD FOR LIMIT %u',REGION1,LIMIT));
        set(gcf,'color',[1,1,1])
        hs=zeros(1,4);
        for TypeL=1:MaxType
            if ~isempty(Types(1,TypeL))&~isempty(Types(2,TypeL))
                LoopNb=1;
                ColNb=1;
                if ListNb(1)==12
                    ColNb=2;
                    LoopNb=2;
                end
                for LoopL=1:LoopNb
                    for ListL=1:ListNb(ChipL)
                        if LoopL==1
                            RepPlotPos=[1,3,5,7,9,11];
                            Offset=0;
                        else
                            RepPlotPos=[2,4,6,8,10,12];
                            Offset=6;
                        end
                        subplot(6,ColNb,RepPlotPos(ListL))
                        hold on
                        RegRep=zeros(1,MAXREP);
                        if SIMREGFLAG
                            for RepL=1:MAXREP
                                MaxRep=0;
                                CurrSimReg=0;
                                for RegL=1:ListNb(ChipL)
                                    if OVERLAPFLAG
                                        Ps1=MatchPsRank{1}(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.core{RepL});
                                        Ps1(find(Ps1==0))=[];
                                        Ps2=MatchPsRank{2}(Region{2}{TypeL}{ListL}{RegL}{LIMIT}.core{RepL});
                                        Ps2(find(Ps2==0))=[];                                        
                                    else
                                        Ps1=MatchPsRank{1}(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.ps{RepL});
                                        Ps1(find(Ps1==0))=[];
                                        Ps2=MatchPsRank{2}(Region{2}{TypeL}{ListL}{RegL}{LIMIT}.ps{RepL});
                                        Ps2(find(Ps2==0))=[];
                                    end
                                    CurrRegRep=length(intersect(Ps1,Ps2))*100/min(length(Ps1),length(Ps2));
                                    if CurrRegRep>MaxRep
                                        CurrSimReg=RegL;
                                        MaxRep=CurrRegRep;
                                    end
                                end
                                %                     if CurrSimReg==0
                                %                         CurrSimReg=REGION2;
                                %                     end
                                if CurrSimReg>0
                                    if OVERLAPFLAG
                                        Ps1=MatchPsRank{1}(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.core{RepL});
                                        Ps1(find(Ps1==0))=[];
                                        Ps2=MatchPsRank{2}(Region{2}{TypeL}{ListL}{CurrSimReg}{LIMIT}.core{RepL});
                                        Ps2(find(Ps2==0))=[];                                                                                         
                                    else
                                        Ps1=MatchPsRank{1}(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.ps{RepL});
                                        Ps1(find(Ps1==0))=[];
                                        Ps2=MatchPsRank{2}(Region{2}{TypeL}{ListL}{CurrSimReg}{LIMIT}.ps{RepL});
                                        Ps2(find(Ps2==0))=[];                                            
                                    end
                                    RegRep(RepL)=length(intersect(Ps1,Ps2))*100/min(length(Ps1),length(Ps2));
                                end
                            end
                        else
                            for RepL=1:MAXREP
                                if OVERLAPFLAG
                                    Ps1=MatchPsRank{1}(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.core{RepL});
                                    Ps1(find(Ps1==0))=[];
                                    Ps2=MatchPsRank{2}(Region{2}{TypeL}{ListL}{REGION2}{LIMIT}.core{RepL});
                                    Ps2(find(Ps2==0))=[];                                                                                                                                 
                                else
                                    Ps1=MatchPsRank{1}(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.ps{RepL});
                                    Ps1(find(Ps1==0))=[];
                                    Ps2=MatchPsRank{2}(Region{2}{TypeL}{ListL}{REGION2}{LIMIT}.ps{RepL});
                                    Ps2(find(Ps2==0))=[];                                            
                                end
                                RegRep(RepL)=length(intersect(Ps1,Ps2))*100/min(length(Ps1),length(Ps2));
                            end
                        end
                        try
                            h=plot([1:MAXREP],RegRep,'color',Colors(TypeL,:));
                            if hs(TypeL)==0
                                hs(TypeL)=h;
                            end
                            plot([1:MAXREP],RegRep,'+','color',Colors(TypeL,:))
                        catch
                        end
                        title(sprintf('%s',ClassName{ListL+Offset}))
                        ylabel('reproducibility')
                        xlabel('>= #networks')
                        set(gca,'box','on')
                        set(gca,'xlim',[1,MAXREP])
                        set(gca,'ylim',[0,100])
                        set(gca,'xgrid','on')
                        set(gca,'ygrid','on')
                    end
                end
            end

        end

        legend(hs,TypeName)
    end




%% RANK OF THE MOST SIMILAR REGION


    LIMIT=1;
    REGION1=1;
    MAXREP=10;
    Colors=colors(colormap,4);
    h=figure;
    set(h,'name',sprintf('MOST SIMILAR REGION FOR REGION %u FOR LIMIT %u',REGION,LIMIT));
    set(gcf,'color',[1,1,1])
    hs=zeros(1,4);
    for TypeL=1:MaxType
        if ~isempty(Types(1,TypeL))
            LoopNb=1;
            ColNb=1;
            if ListNb(1)==12
                ColNb=2;
                LoopNb=2;
            end
            for LoopL=1:LoopNb
                for ListL=1:ListNb(ChipL)
                    if LoopL==1
                        RepPlotPos=[1,3,5,7,9,11];
                        Offset=0;
                    else
                        RepPlotPos=[2,4,6,8,10,12];
                        Offset=6;
                    end
                    subplot(8,ColNb,RepPlotPos(ListL))
                    hold on
                    SimReg=zeros(1,MAXREP);
                    for RepL=1:MAXREP
                        try
                            MaxRep=0;
                            CurrSimReg=0;
                            for RegL=1:ListNb(ChipL)
                                Ps1=MatchPsRank{1}(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.ps{RepL});
                                Ps1(find(Ps1==0))=[];
                                Ps2=MatchPsRank{2}(Region{2}{TypeL}{ListL}{RegL}{LIMIT}.ps{RepL});
                                Ps2(find(Ps2==0))=[];
                                CurrRegRep=length(intersect(Ps1,Ps2))*100/min(length(Ps1),length(Ps2));
                                if CurrRegRep>MaxRep
                                    CurrSimReg=RegL;
                                    MaxRep=CurrRegRep;
                                end
                            end                         
                            SimReg(RepL)=CurrSimReg;
                        catch
                        end
                    end
                    try
                        h=plot([1:MAXREP],SimReg,'color',Colors(TypeL,:));
                        if hs(TypeL)==0
                            hs(TypeL)=h;
                        end
                        plot([1:MAXREP],SimReg,'+','color',Colors(TypeL,:))
                    catch
                    end
                    title(sprintf('%s',ClassName{ListL+Offset}))
                    ylabel('most similar')
                    xlabel('>= #networks')
                    set(gca,'box','on')
                    set(gca,'xgrid','on')
                    set(gca,'ygrid','on')
                end
            end
        end

    end

    legend(hs,TypeName)




%% RANDOM MODULES

    TYPE=1;
    REGION=5;
    LIMIT=1;
    REPROD=5;
    %PS CLASS=5;
    CHIP=1;
    NET=1;
    cd('/home/mbellis/net/clum8')
    NetRanks{1}=[24:38];
    tic
    %for NetL=1:NetNb(CHIP)


    C=load_data(sprintf('c_m%u_n%u.4mat',ChipRank(CHIP),NetRanks{CHIP}(NET)),'./',PsNb(CHIP),PsNb(CHIP),'uint8','ieee-le',Region{CHIP}{TYPE}{LIST}{REGION}{LIMIT}.ps{REPROD},Region{CHIP}{TYPE}{LIST}{REGION}{LIMIT}.ps{REPROD});
    if TYPE==2 |TYPE==4
        A=load_data(sprintf('a_m%u_n%u.4mat',ChipRank(CHIP),NetRanks{CHIP}(NET)),'./',PsNb(CHIP),PsNb(CHIP),'uint8','ieee-le',Region{CHIP}{TYPE}{LIST}{REGION}{LIMIT}.ps{REPROD},Region{CHIP}{TYPE}{LIST}{REGION}{LIMIT}.ps{REPROD});
    end
    if TYPE==3 | TYPE==4
        F=load_data(sprintf('f_m%u_n%u.4mat',ChipRank(CHIP),NetRanks{CHIP}(NET)),'./',PsNb(CHIP),PsNb(CHIP),'uint8','ieee-le',Region{CHIP}{TYPE}{LIST}{REGION}{LIMIT}.ps{REPROD},Region{CHIP}{TYPE}{LIST}{REGION}{LIMIT}.ps{REPROD});
    end
    Mod3=[];
    Mod5=[];
    Mod10=[];
    RegSize=length(C);

    RandNb=1000;
    MaxMod=30;
    for i=1:MaxMod
        Pos=uint16(zeros(10,RandNb));
    end
    for PsL=1:MaxMod
        Pos(PsL,:)=uint16(ceil(rand(RandNb,1)*RegSize));
    end
    %modules of size 3
    C=uint16(C);
    Density=uint8(zeros(4,RandNb));
    CurrC=uint16(zeros(1,RandNb));
    for k=1:RandNb
        CurrC(k)=uint16(sum(sum(C(Pos(1:3,k),Pos(1:3,k))))-300);
    end
    Density(1,:)=uint8(round(CurrC/6));

    %modules of size 5
    CurrC=uint16(zeros(1,RandNb));
    for k=1:RandNb
        CurrC(k)=uint16(sum(sum(C(Pos(1:5,k),Pos(1:5,k))))-500);
    end
    Density(2,:)=uint8(round(CurrC/20));

    %modules of size 10
    CurrC=uint16(zeros(1,RandNb));
    for k=1:RandNb
        CurrC(k)=uint16(sum(sum(C(Pos(1:10,k),Pos(1:10,k))))-1000);
    end
    Density(3,:)=uint8(round(CurrC/90));

    %modules of size 30
    CurrC=uint32(zeros(1,RandNb));
    for k=1:RandNb
        CurrC(k)=uint32(sum(sum(C(Pos(1:30,k),Pos(1:30,k))))-3000);
    end
    Density(4,:)=uint8(round(CurrC/870));

    figure
    plot([0:100],histc(Density(1,:),[0:100]),'b')
    hold on
    plot([0:100],histc(Density(2,:),[0:100]),'g')
    plot([0:100],histc(Density(3,:),[0:100]),'r')
    plot([0:100],histc(Density(4,:),[0:100]),'k')

    toc
    figure
    Pos=find(C>0&C<100);
    plot([0:100],histc(C(Pos),[0:100]),'c')
end

%% INTERCHIP COMPARISONS

%only one network for each chip
if ChipNb==2
    CurrTypeNb=length(find(sum(Types)));
    if isempty(ChipRange)
        AllLimits=[];
        for ChipL=1:ChipNb
            for TypeL=1:MaxType
                AllLimits=union(AllLimits,Information{ChipL}.limits{TypeL});
            end
        end
    end
    LimitNb=length(AllLimits);
    LimitNames=cell(1,LimitNb);
    for LimitL=1:LimitNb
        LimitNames{LimitL}=num2str(AllLimits(LimitL));
    end

    %only one network for each chip
    InterChipCluOrder=cell(1,2);
    InterChipCluOrder{1}=cell(1,6);
    InterChipCluOrder{2}=cell(1,6);

    MinSize=10;
    for ListL=1:ListNb(ChipL)
        h=figure;
        set(h,'color',[1,1,1])
        set(h,'name',sprintf('Liste %u (%s)',ListL,ClassName{ListL}))
        CurrTypePos=0;
        %used to resize individual plots
        XSizes=cell(MaxType);
        YSizes=cell(MaxType);
        TypePos=0;
        for TypeL=1:MaxType
            if Types(1,TypeL)==1&Types(2,TypeL)==1
                TypePos=TypePos+1;
                XSizes{TypeL}=zeros(LimitNb,1);
                YSizes{TypeL}=zeros(LimitNb,1);
                InterChipCluOrder{TypeL}{ListL}=zeros((LimitNb+1)*4,20);
                for LimitL=1:LimitNb
                    %number of clusters that have a size >=MinSize
                    MaxCluPos=[0,0];
                    CluSize1=zeros(1,20);
                    CluSize2=zeros(1,20);
                    for ChipL=1:2
                        CurrLimitPos=find(Information{ChipL}.limits{TypeL}==AllLimits(LimitL));
                        if ~isempty(CurrLimitPos)
                            LimitPos(ChipL)=CurrLimitPos;
                            try
                                CurrClusterNb=ClusterNb{ChipL}{TypeL}{ListL}{1}{LimitPos(ChipL)};
                                if ~isempty(CurrClusterNb)
                                    CurrMaxCluPos=find(CurrClusterNb>=MinSize);
                                    MaxCluPos(ChipL)=CurrMaxCluPos(end);
                                    if ChipL==1
                                        YSizes{TypeL}(LimitL)=CurrMaxCluPos(end);
                                    else
                                        XSizes{TypeL}(LimitL)=CurrMaxCluPos(end);
                                    end
                                end
                            catch
                            end
                        end
                    end
                    if MaxCluPos(1)>0&MaxCluPos(2)>0
                        Reprod=zeros(MaxCluPos(1)+1,MaxCluPos(2)+1);
                        %nb of common ps between any two clusters
                        Used=zeros(MaxCluPos(1)+1,MaxCluPos(2)+1);
                        %nb of ps in selected (>=MinSize) clusters
                        Total1=zeros(MaxCluPos(1)+1,1);
                        Total2=zeros(1,MaxCluPos(2)+1);
                        %reorder Reprod columns
                        Sort=[];
                        for CluL1=1:MaxCluPos(1)
                            Ps1=MatchPsRank{1}(find(Cluster{1}{TypeL}{ListL}{1}(:,LimitPos(1))==CluL1));
                            Ps1(find(Ps1==0))=[];
                            CluSize1(CluL1)=length(Ps1);
                            Total1(CluL1)=length(Ps1);
                            MaxPos=0;
                            MaxVal=0;
                            %find the cluster with the maximal overlap
                            for CluL2=1:MaxCluPos(2)
                                Ps2=MatchPsRank{2}(find(Cluster{2}{TypeL}{ListL}{1}(:,LimitPos(2))==CluL2));
                                Ps2(find(Ps2==0))=[];
                                CluSize2(CluL2)=length(Ps2);
                                MinReprod=round(length(intersect(Ps1,Ps2))*100/max(length(Ps1),length(Ps2)));
                                %Reprod(CluL1,CluL2)=round(length(intersect(Ps1,Ps2))*100/min(length(Ps1),length(Ps2)));
                                %Reprod(CluL1,CluL2)=round(length(intersect(Ps1,Ps2))*100/sqrt(length(Ps1)*length(Ps2)));
                                Reprod(CluL1,CluL2)=round(length(intersect(Ps1,Ps2))*100/length(union(Ps1,Ps2)));
                                Used(CluL1,CluL2)=length(intersect(Ps1,Ps2));
                                if  CluL1==1
                                    Total2(CluL2)=length(Ps2);
                                end
                                if MinReprod>MaxVal
                                    %don'use twice the same column
                                    if isempty(find(Sort==CluL2))
                                        Sort(CluL1)=CluL2;
                                        MaxVal=MinReprod;
                                    end
                                end
                            end
                        end
                        NullPos=find(Sort==0);
                        if MaxCluPos(1)>MaxCluPos(2)
                            Sort(NullPos)=[];
                        else
                            if ~isempty(NullPos)
                                if isempty(intersect(NullPos,Sort))
                                    Sort(NullPos)=NullPos;
                                else
                                    NotUsed=setdiff([1:MaxCluPos(2)],Sort);
                                    UsedRank=0;
                                    for PosL=1:length(NullPos)
                                        if isempty(find(Sort==NullPos(PosL)))
                                            Sort(NullPos(PosL))=NullPos(PosL);
                                            NotUsed=setdiff(NotUsed,NullPos(PosL));
                                        else
                                            UsedRank=UsedRank+1;
                                            Sort(NullPos(PosL))=NotUsed(UsedRank);
                                        end
                                    end
                                end
                            end
                        end
                        Reprod(:,1:length(Sort))=Reprod(:,Sort);
                        EndPos=min(20,length(Sort));
                        CluSize2=CluSize2(Sort);
                        InterChipCluOrder{TypeL}{ListL}(LimitL,1:EndPos)=Sort(1:EndPos);
                        InterChipCluOrder{TypeL}{ListL}(LimitL+LimitNb+1,1:EndPos)=CluSize1(1:EndPos);
                        InterChipCluOrder{TypeL}{ListL}(LimitL+LimitNb*2+2,1:EndPos)=CluSize2(1:EndPos);
                        InterChipCluOrder{TypeL}{ListL}(LimitL+LimitNb*3+3,1:EndPos)=diag(Reprod(1:EndPos,1:EndPos))';

                        Used(:,1:length(Sort))=Used(:,Sort);
                        AllUsed=sum(sum(Used));
                        AllTotal1=sum(Total1);
                        Total2(:,1:length(Sort))=Total2(:,Sort);
                        AllTotal2=sum(Total2);
                        subplot(LimitNb,CurrTypeNb,CurrTypeNb*(LimitL-1)+TypePos)                        
                        h=pcolor(Reprod);
                        if MaxCluPos(1)*MaxCluPos(2)>1000
                            set(h,'linestyle','none')
                        end
                        %print Limit
                        % percentage of all ps common to corresponding clusters (>=MinSize) (i.e. main diagonal of Used) relatively to all the ps contained in the (>=MinSize) clusters of one chip
                        % nb of ps contained in the (>=minSize) clusters of one chip
                        % percentage of ps contained in the (>=minSize) clusters of one chiprelatively to the number of ps contained in the current list
                        xlabel(sprintf('%u (%u(%u) %u(%u))',Information{2}.limits{TypeL}(LimitPos(2)),round(sum(diag(Used))*100/AllTotal2),AllTotal2,round(AllTotal2*100/ListSize{2}(ListL)),ListSize{2}(ListL)))
                        ylabel(sprintf('%u (%u(%u) %u(%u))',Information{1}.limits{TypeL}(LimitPos(1)),round(sum(diag(Used))*100/AllTotal1),AllTotal1,round(AllTotal1*100/ListSize{1}(ListL)),ListSize{1}(ListL)))
%                         if MaxType==2
%                             %process the first 20 clusters that are on the main diagonal
%                             UsedA=Used(1:min(20,size(Reprod,1)),1:min(20,size(Reprod,2)));
%                             DiagUsedA=sum(diag(UsedA));
%                             AllUsedA=sum(sum(UsedA));
%                             TotalA1=Total1(1:min(20,size(Reprod,1)));
%                             TotalA2=Total2(1:min(20,size(Reprod,2)));
%                             AllTotalA1=sum(TotalA1);
%                             AllTotalA2=sum(TotalA2);
%                             %process all other clusters (find the best overlap for each
%                             MaxClu1=size(Reprod,1);
%                             MaxClu2=size(Reprod,2);
%                             if MaxClu1>20|MaxClu2>20
%                                 AllTotalB1=AllTotal1-AllTotalA1;
%                                 AllTotalB2=AllTotal2-AllTotalA2;
%                                 BestUsedB1=0;
%                                 BestUsedB2=0;
%                                 if MaxClu1>20
%                                     for CluL=21:MaxClu1
%                                         BestUsedB1=BestUsedB1+max(Reprod(CluL,21:end));
%                                     end
%                                 end
%                                 if MaxClu2>20
%                                     for CluL=21:MaxClu2
%                                         BestUsedB2=BestUsedB2+max(Reprod(21:end,CluL));
%                                     end
%                                 end
%                                 ReprodScoreB1=round(BestUsedB1*100/AllTotalB1);
%                                 ReprodScoreB2=round(BestUsedB2*100/AllTotalB2);
%                             else
%                                 ReprodScoreB1=0;
%                                 ReprodScoreB2=0;
%                                 AllTotalB1=0;
%                                 AllTotalB2=0;
%                             end
%                             subplot(LimitNb,4,(LimitL-1)*4+TypeL+2)
%                             h=pcolor(Reprod(1:min(20,size(Reprod,1)),1:min(20,size(Reprod,2))));
%                             xlabel(sprintf('%u (%u %u %u(%u) %u(%u))',Information{2}.limits{TypeL}(LimitPos(2)),round(AllTotalA2*100/AllTotal2),round(AllUsedA*100/AllTotalA2),round(DiagUsedA*100/AllTotalA2),AllTotalA2,ReprodScoreB2,AllTotalB2))
%                             ylabel(sprintf('%u (%u %u %u(%u) %u(%u))',Information{1}.limits{TypeL}(LimitPos(1)),round(AllTotalA1*100/AllTotal1),round(AllUsedA*100/AllTotalA1),round(DiagUsedA*100/AllTotalA1),AllTotalA1,ReprodScoreB1,AllTotalB1))
%                         end
                    end
                end
            end
            %resize plots
        end
    end
end


%% WRITE LISTS
WriteList=questdlg('Do you want to write a list','','yes','no','no');
while isequal(WriteList,'yes')
    %write correspondance between chips
    cd(fullfile(K.dir.net,sprintf('m%u',ChipRank(1)),'list'));
    fid=fopen(sprintf('m%u_m%u_cluster_correspondance.txt',ChipRank(1),ChipRank(2)),'a');
    CluSel=inputdlg({'Type','LimitRank','ListRank','CluRank','I,U,C'});
    CurrType=str2num(CluSel{1});
    CurrLimit1=str2num(CluSel{2});
    CurrLimit2=find(Information{2}.limits{CurrType}==Information{2}.limits{CurrType}(CurrLimit1));
    CurrList=str2num(CluSel{3});
    CurrClu1=str2num(CluSel{4});
    CurrSel=CluSel{5};
    CurrClu2=InterChipCluOrder{CurrType}{CurrList}(CurrLimit1,CurrClu1);
    PsRank1=MatchPsRank{1}(find(Cluster{1}{CurrType}{CurrList}{1}(:,CurrLimit1)==CurrClu1));
    PsRank1(find(PsRank1==0))=[];
    PsRank1in2=MatchPsRank{2}(find(Cluster{1}{CurrType}{CurrList}{1}(:,CurrLimit1)==CurrClu1));
    PsRank1in2(find(PsRank1in2==0))=[];
    if ~isempty(CurrClu2)&~isempty(CurrLimit2)
        PsRank2=find(Cluster{2}{CurrType}{CurrList}{1}(:,CurrLimit2)==CurrClu2);
    else
        PsRank2=[];
    end

    if isequal(CurrSel,'I')
        if ~isempty(PsRank2)
            [Inter,Index1,Index2]=intersect(PsRank1in2,PsRank2);
            ListRank=sprintf('11%u%u%u',CurrType,CurrList,CurrClu1);
            ListRank1=str2num(ListRank);
            ListRank=sprintf('11%u%u%u',CurrType,CurrList,CurrClu2);
            ListRank2=str2num(ListRank);
        else
            h=warndlg('no ps in second chip cluster');
            waitfor(h)
        end
    elseif isequal(CurrSel,'U')
        if ~isempty(PsRank2)
            Index2=union(PsRank1in2,PsRank2);
            %calculate Index1  to be done
            ListRank=sprintf('12%u%u%u',CurrType,CurrList,CurrClu1);
            ListRank1=str2num(ListRank);
            ListRank=sprintf('12%u%u%u',CurrType,CurrList,CurrClu2);
            ListRank2=str2num(ListRank);
        else
            h=warndlg('no ps in second chip cluster');
            waitfor(h)
        end
    elseif isequal(CurrSel,'C')
        Index1=PsRank1;
        if ~isempty(PsRank2)
            Index2=PsRank2;
        else
            Index2=[];
        end
    else
        h=warndlg(sprintf('%s not allowed',CurrSel));
        waitfor(h)
    end
    cd(fullfile(K.dir.net,sprintf('m%u',ChipRank(1)),'list'));
    save_data(Index1,sprintf('m%u_pslist%u.u32',ChipRank(1),ListRank1),'./','w','uint32','ieee-le');
    if ~isempty(PsRank2)
        fprintf(fid,'%u\t%u\n',CurrClu1,CurrClu2)
        cd(fullfile(K.dir.net,sprintf('m%u',ChipRank(2)),'list'));
        save_data(Index2,sprintf('m%u_pslist%u.u32',ChipRank(2),ListRank2),'./','w','uint32','ieee-le');
    end
    WriteList=questdlg('Do you want to write another list','','yes','no','no');
    if ~isequal(WriteList,'yes')
        fclose(fid)
    end
end


%% DISPLAY CVM
DisplayCVM=questdlg('Do you want to display CVM','','yes','no','no');
%Cluster{ChipL}{TypeL}{ListL}{NetL}(:,LimitL)
while isequal(DisplayCVM,'yes')
    %select in list 1 (All)
    CluSel=inputdlg({'Type','LimitRank','ListRank','CluRanks'});
    CurrType=str2num(CluSel{1});
    CurrLimit=str2num(CluSel{2});
    CurrList=str2num(CluSel{3});
    CurrClu=str2num(CluSel{4});
    OutName{1}=sprintf('m%u_t%u_l%u_L1_c%s',ChipRank(1),CurrType,Information{1}.limits{CurrType}(CurrLimit),strrep(CluSel{4},',','-'));
    OutName{2}=sprintf('m%u_m%u_t%u_l%u_L1_c%s',ChipRank(2),ChipRank(1),CurrType,Information{1}.limits{CurrType}(CurrLimit),strrep(CluSel{4},',','-'));
    %recover different Ps classes contained in All clusters
    %Limit={};
    Limit={};
    Ps{1}={};
    Ps{2}={};
    PsPos={};
    PsIndex{2}={};
    PsIndex{1}={};
    for GrpL=[2,3,4,6]            
        PsIndex{1}{GrpL}=[];
        PsIndex{2}{GrpL}=[];
        PsPos{GrpL}=[0,0];
    end
    %Class2to4=union(union(PsRankList{1}{2},PsRankList{1}{3}),PsRankList{1}{4});
    Class{2}=PsRankList{1}{2};
    Class{3}=PsRankList{1}{3};
    Class{4}=PsRankList{1}{4};
    Class{6}=PsRankList{1}{6};
    for CluL=1:length(CurrClu)        
        Ps{1}{CluL}={};
        Ps{2}{CluL}={};
        AllPs=find(Cluster{1}{CurrType}{1}{1}(:,CurrLimit)==CurrClu(CluL));
        for GrpL=[2,3,4,6]            
            CurrClass=Class{GrpL};
            if DiffChipFlag
                Ps1=MatchPsRank{1}(intersect(AllPs,CurrClass));
                Ps1(find(Ps1==0))=[];            
                Ps{1}{CluL}{GrpL}=zeros(1,length(Ps1));
                Ps{2}{CluL}{GrpL}=zeros(1,length(Ps1));
                %keep only one ps for each common position
                for PsL=1:length(Ps1)
                    Ps{1}{CluL}{GrpL}(PsL)=CommonPsRank{1}{Ps1(PsL)}(1);
                    Ps{2}{CluL}{GrpL}(PsL)=CommonPsRank{2}{Ps1(PsL)}(1);
                end
            else
                Ps{1}{CluL}{GrpL}(PsL)=intersect(AllPs,CurrClass);
                Ps{2}{CluL}{GrpL}=Ps{1}{CluL}{GrpL};
            end
            PsIndex{1}{GrpL}=[PsIndex{1}{GrpL},Ps{1}{CluL}{GrpL}'];
            PsIndex{2}{GrpL}=[PsIndex{2}{GrpL},Ps{2}{CluL}{GrpL}'];
            PsPos{GrpL}(1)=PsPos{GrpL}(2)+1;
            PsPos{GrpL}(2)=PsPos{GrpL}(1)+length(Ps{1}{CluL}{GrpL})-1;
            Limit{CluL}{GrpL}=PsPos{GrpL};
        end
            
%         for ClassL=[2,3,4,6]
%             Ps{1}{CluL}{ClassL}=intersect(AllPs,PsRankList{1}{ClassL});
%             PsIndex{1}=[PsIndex{1},Ps{1}{CluL}{ClassL}'];
%             PsPos(1)=PsPos(2)+1;
%             PsPos(2)=PsPos(1)+length(Ps{1}{CluL}{ClassL})-1;
%             Limit{1}{CluL}{ClassL}=PsPos;
%             %recover corresponding ps in second chip in the same order
%             Ps{2}{CluL}{ClassL}=MatchPsRank(Ps{1}{CluL}{ClassL});
%         end
    end
%     for GrpL=[2,3,4,6]            
%         PsIndex{2}{GrpL}=MatchPsRank(PsIndex{1}{GrpL});
%     end




    %calculate ps order in each cluster by running anneal clust in first chip
    ChipL=1;
    cd(fullfile(K.dir.net,sprintf('m%u',ChipRank(ChipL)),sprintf('n%u',Information{ChipL}.netRanks(ChipL))))
    FileNameList={};
    PsOrder={};
    tic
    for CluL=1:length(CurrClu)
        for GrpL=[2,3,4,6]        
            c=load_data(sprintf('c_m%u_n%u.4mat',ChipRank(ChipL),Information{ChipL}.netRanks(1)),'.',...
                K.chip.probesetNb(ChipPos(ChipL)),K.chip.probesetNb(ChipPos(ChipL)),'uint8','ieee-le',...
                Ps{1}{CluL}{GrpL},Ps{1}{CluL}{GrpL});
            if CurrType==1
                [PsOrder{CluL}{GrpL},temp,temp]=annealclust(single(c)/100);
            else
                a=load_data(sprintf('a_m%u_n%u.4mat',ChipRank(ChipL),Information{ChipL}.netRanks(1)),'.',...
                    K.chip.probesetNb(ChipPos(ChipL)),K.chip.probesetNb(ChipPos(ChipL)),'uint8','ieee-le',...
                    Ps{1}{CluL}{GrpL},Ps{1}{CluL}{GrpL});
                [PsOrder{CluL}{GrpL},temp,temp]=annealclust(single(c-a)/100);
            end
        end
    end
    toc


    
    %distribution anti/corr
    CurrNet=[111,110];
    for ChipL=1:2
        if CurrNet(ChipL)>0
            h=figure;
            set(h,'name',sprintf('m%u n%u: rank distribution',ChipRank(ChipL),CurrNet(ChipL)))
            set(gcf,'color',[1,1,1])
            cd(fullfile(K.dir.net,sprintf('m%u',ChipRank(ChipL)),sprintf('n%u',CurrNet(ChipL))))
            GrpPos=0;
            for GrpL=[2,3,4,6]
                GrpPos=GrpPos+1;
                c=load_data(sprintf('c_m%u_n%u.4mat',ChipRank(ChipL),CurrNet(ChipL)),'.',...
                    K.chip.probesetNb(ChipPos(ChipL)),K.chip.probesetNb(ChipPos(ChipL)),'uint8','ieee-le',...
                    PsIndex{ChipL}{GrpL},Ps{ChipL}{1}{GrpL});
                a=load_data(sprintf('a_m%u_n%u.4mat',ChipRank(ChipL),CurrNet(ChipL)),'.',... 
                    K.chip.probesetNb(ChipPos(ChipL)),K.chip.probesetNb(ChipPos(ChipL)),'uint8','ieee-le',...
                    PsIndex{ChipL}{GrpL},Ps{ChipL}{1}{GrpL});               
                %statistics on c and a
%                 Pos=find(c);
%                 Neg=find(a);
%                 subplot(1,4,GrpPos);
%                 Val=histc(c(Pos),[1:100]);
%                 CMedianVal=median(c(Pos));
%                 SVal=cumsum(Val);
%                 C75th=find(SVal<=0.75*sum(Val));
%                 C75th=C75th(end);
%                 CPerc=round(SVal(70)*100/sum(Val));
%                 plot([1:100],Val,'r-')
%                 hold on
%                 line([CMedianVal,CMedianVal],[0,Val(CMedianVal)],'color','r')
%                 line([C75th,C75th],[0,Val(C75th)],'color','r','linestyle',':')
%                 Val=histc(a(Neg),[1:100]);
%                 AMedianVal=median(a(Pos));
%                 SVal=cumsum(Val);
%                 A75th=find(SVal<=0.75*sum(Val));
%                 A75th=A75th(end);
%                 APerc=round(SVal(60)*100/sum(Val));
%                 plot([1:100],Val,'b-')
%                 line([AMedianVal,AMedianVal],[0,Val(AMedianVal)-1],'color','b')
%                 line([A75th,A75th],[0,Val(A75th)],'color','b','linestyle',':')
%                 xlabel(sprintf('c-a:%u-%u & %u-%u  (c:%u%% a:%u%%)',CMedianVal,C75th,AMedianVal,A75th,round(length(Pos)*100/(size(c,1)*size(c,2))),round(length(Neg)*100/(size(a,1)*size(a,2)))));
%                 title(sprintf('m%u grp%u (70=c%uth 60=a%uth)',ChipRank(ChipL),GrpL,CPerc,APerc));
                %statistics on c-a
                %eliminate 100 values (on diagonal where corr is between the same probeset)
                ValNb=size(c,1)*size(c,2);
                Range=[-100:100];
                Pos=find(c==100);
                c(Pos)=0;
                Pos=find(c);
                c=c(Pos);
                a=a(Pos);
                c=double(c)-double(a);
                subplot(1,4,GrpPos);
                Val=histc(c,Range);
                CMedianVal=median(c);
                SVal=cumsum(Val);
                C50thPos=find(SVal<=0.50*sum(Val));
                C50thPos=C50thPos(end);
                C50th=Range(C50thPos);
                C75thPos=find(SVal<=0.75*sum(Val));
                C75thPos=C75thPos(end);
                C75th=Range(C75thPos);
                CPerc=round(SVal(131)*100/sum(Val));
                plot([-100:100],Val,'k-')
                hold on
                line([C50th,C50th],[0,Val(C50thPos)],'color','k')
                line([C75th,C75th],[0,Val(C75thPos)],'color','k','linestyle',':')                
                %set(gca,'xlim',[-50,50])
                xlabel(sprintf('med:%u 75th:%u  (c:%u%%)',CMedianVal,C75th,round(length(Pos)*100/ValNb)));
                title(sprintf('m%u grp%u (30=%uth)',ChipRank(ChipL),GrpL,CPerc));
            end
        end
    end


    %write individual bmp files (Clu vs Clu)
    MaxVal=30;    
    Colors=colors(colormap,MaxVal);                        
    for ChipL=1:2
        tic
        cd(fullfile(K.dir.net,sprintf('m%u',ChipRank(ChipL)),sprintf('n%u',Information{ChipL}.netRanks(1))))
        for CluL1=1:length(CurrClu)
            CluL1
            GrpPos=0;
            for GrpL=[2,3,4,6]
                GrpPos=GrpPos+1;
                c=load_data(sprintf('c_m%u_n%u.4mat',ChipRank(ChipL),Information{ChipL}.netRanks(1)),'.',...
                    K.chip.probesetNb(ChipPos(ChipL)),K.chip.probesetNb(ChipPos(ChipL)),'uint8','ieee-le',...
                    PsIndex{ChipL}{GrpL},Ps{ChipL}{CluL1}{GrpL});
                a=load_data(sprintf('a_m%u_n%u.4mat',ChipRank(ChipL),Information{ChipL}.netRanks(1)),'.',...                K.chip.probesetNb(ChipPos(ChipL)),K.chip.probesetNb(ChipPos(ChipL)),'uint8','ieee-le',...
                    K.chip.probesetNb(ChipPos(ChipL)),K.chip.probesetNb(ChipPos(ChipL)),'uint8','ieee-le',...
                    PsIndex{ChipL}{GrpL},Ps{ChipL}{CluL1}{GrpL});
                %make a and c horizontal
                c=c';
                a=a';
                %order ps lines
                c=c(PsOrder{CluL1}{GrpL},:);
                a=a(PsOrder{CluL1}{GrpL},:);
                %order ps column
                for CluL2=1:length(CurrClu)                
                    ColOrder=[Limit{CluL2}{GrpL}(1):Limit{CluL2}{GrpL}(2)];
                    c(:,Limit{CluL2}{GrpL}(1):Limit{CluL2}{GrpL}(2))=c(:,ColOrder(PsOrder{CluL2}{GrpL}));
                    a(:,Limit{CluL2}{GrpL}(1):Limit{CluL2}{GrpL}(2))=a(:,ColOrder(PsOrder{CluL2}{GrpL}));
                end
                
                %keep c horizontal                
                %make a vertical
                a=a';
                %process all clu and save bmp files
                for CluL2=CluL1:length(CurrClu)
                     CurrC=c(:,Limit{CluL2}{GrpL}(1):Limit{CluL2}{GrpL}(2));
                     CurrA=a(Limit{CluL2}{GrpL}(1):Limit{CluL2}{GrpL}(2),:);
                     if CluL2==CluL1                           
                        %make composit bloc (C and A values)
                        CurrC=triu(CurrC,1)+tril(CurrA,-1);                      
                    end                    
                    FileName=sprintf('%s_clu%u_vs_clu%u_grp%u.bmp',OutName{ChipL},CluL1,CluL2,GrpL);
                    imwrite(CurrC,Colors,FileName,'bmp')
                    if CluL2~=CluL1
                        FileName=sprintf('%s_clu%u_vs_clu%u_grp%u.bmp',OutName{ChipL},CluL2,CluL1,GrpL);
                        imwrite(CurrA,Colors,FileName,'bmp')
                    end
                end
            end
        end
        toc
    end
    
    %write bmp files for qlimit file
    MaxVal=30;
    CurrNet=[211,0];
    for ChipL=1:2
        if CurrNet(ChipL)>0
        tic
        cd(fullfile(K.dir.net,sprintf('m%u',ChipRank(ChipL)),sprintf('n%u',CurrNet(ChipL))))
        for CluL1=1:length(CurrClu)
            CluL1
            for GrpL=[2,3,4,6]
                c=load_data(sprintf('c_m%u_n%u.4mat',ChipRank(ChipL),CurrNet(ChipL)),'.',...
                    K.chip.probesetNb(ChipPos(ChipL)),K.chip.probesetNb(ChipPos(ChipL)),'uint8','ieee-le',...
                    PsIndex{ChipL}{GrpL},Ps{ChipL}{CluL1}{GrpL});
                a=load_data(sprintf('a_m%u_n%u.4mat',ChipRank(ChipL),CurrNet(ChipL)),'.',...
                    K.chip.probesetNb(ChipPos(ChipL)),K.chip.probesetNb(ChipPos(ChipL)),'uint8','ieee-le',...
                    PsIndex{ChipL}{GrpL},Ps{ChipL}{CluL1}{GrpL});
                %make a and c horizontal
                c=c';
                a=a';
                %order ps lines
                c=c(PsOrder{CluL1}{GrpL},:);
                a=a(PsOrder{CluL1}{GrpL},:);
                %order ps column
                for CluL2=1:length(CurrClu)
                    ColOrder=[Limit{CluL2}{GrpL}(1):Limit{CluL2}{GrpL}(2)];
                    c(:,Limit{CluL2}{GrpL}(1):Limit{CluL2}{GrpL}(2))=c(:,ColOrder(PsOrder{CluL2}{GrpL}));
                    a(:,Limit{CluL2}{GrpL}(1):Limit{CluL2}{GrpL}(2))=a(:,ColOrder(PsOrder{CluL2}{GrpL}));
                end
                %keep c horizontal
                %make a vertical
                a=a';
                %process all clu and save bmp files
                for CluL2=CluL1:length(CurrClu)
                    CurrC=c(:,Limit{CluL2}{GrpL}(1):Limit{CluL2}{GrpL}(2));
                    CurrA=a(Limit{CluL2}{GrpL}(1):Limit{CluL2}{GrpL}(2),:);
                    if CluL2==CluL1
                        %make composit bloc (C and A values)
                        CurrC=triu(CurrC,1)+tril(CurrA,-1);
                        if CluL1==1&CluL2==1
                            if MaxVal==0
                                MaxVal=max(max(CurrC));
                            end
                            Colors=colors(colormap,MaxVal+1);
                        end
                    end                    
                    FileName=sprintf('%s_clu%u_vs_clu%u_grp%u.bmp',OutName{ChipL},CluL1,CluL2,GrpL);
                    imwrite(CurrC,Colors,FileName,'bmp')
                    if CluL2~=CluL1
                        FileName=sprintf('%s_clu%u_vs_clu%u_grp%u.bmp',OutName{ChipL},CluL2,CluL1,GrpL);
                        imwrite(CurrA,Colors,FileName,'bmp')
                    end
                end
            end
        end
        toc
        end
    end
    

    for GrpL=[2,3,4,6]
        for ChipL=1:2
            FileNameList{ChipL}{GrpL}={};
            for CluL1=1:length(CurrClu)
                for CluL2=CluL1:length(CurrClu)
                    FileName=sprintf('%s_clu%u_vs_clu%u_grp%u.bmp',OutName{ChipL},CluL1,CluL2,GrpL);
                    FileNameList{ChipL}{GrpL}{end+1,1}=FileName;
                    if CluL2~=CluL1
                        FileName=sprintf('%s_clu%u_vs_clu%u_grp%u.bmp',OutName{ChipL},CluL2,CluL1,GrpL);
                        FileNameList{ChipL}{GrpL}{end+1,1}=FileName;
                    end
                end
            end
        end
    end
    

    FileOrder=zeros(length(CurrClu),length(CurrClu));
    Rank=0;
    for CluL1=1:length(CurrClu)
        for CluL2=CluL1:length(CurrClu)
            Rank=Rank+1;
            FileOrder(CluL1,CluL2)=Rank;
            if CluL2~=CluL1
                Rank=Rank+1;
                FileOrder(CluL2,CluL1)=Rank;
            end
        end
    end
    FileOrder=FileOrder';

    
    BlocNb=length(CurrClu);
    %BlocNb=2;
    for GrpL=[2,3,4,6]
        for ChipL=1:2
            cd(fullfile(K.dir.net,sprintf('m%u',ChipRank(ChipL)),sprintf('n%u',Information{ChipL}.netRanks(1))))
            fid=fopen(sprintf('%s_imagelist_grp%u.txt',OutName{ChipL},GrpL),'w');
            CurrFileOrder=FileOrder(1:BlocNb,1:BlocNb);
            CurrFileNameList=FileNameList{ChipL}{GrpL}(CurrFileOrder(:));
            for NameL=1:length(CurrFileNameList)
                fprintf(fid,'%s\n',CurrFileNameList{NameL});
            end
            fclose(fid)
            fid=fopen(sprintf('%s_grp%u.bat',OutName{ChipL},GrpL),'w');
            fprintf(fid,'montage @%s_imagelist_grp%u.txt -tile %ux%u -geometry "10x10<" -mode Concatenate  %s_grp%u.bmp',OutName{ChipL},GrpL,BlocNb,BlocNb,OutName{ChipL},GrpL)
            fclose(fid)
        end
    end
    DisplayCVM=questdlg('Do you want to display CVM','','yes','no','no');
end


%     for ChipL=1:2
%         tic       
%         cd(fullfile(K.dir.net,sprintf('m%u',ChipRank(ChipL)),sprintf('n%u',Information{ChipL}.netRanks(1))))
%         for CluL1=1:length(CurrClu)
%             CluL1
%             for ClassL1=[2,3,4,6]
%                 c=load_data(sprintf('c_m%u_n%u.4mat',ChipRank(ChipL),Information{ChipL}.netRanks(1)),'.',...
%                     K.chip.probesetNb(ChipPos(ChipL)),K.chip.probesetNb(ChipPos(ChipL)),'uint8','ieee-le',...
%                     PsIndex{ChipL},Ps{ChipL}{CluL1}{ClassL1});
%                 a=load_data(sprintf('a_m%u_n%u.4mat',ChipRank(ChipL),Information{ChipL}.netRanks(1)),'.',...                K.chip.probesetNb(ChipPos(ChipL)),K.chip.probesetNb(ChipPos(ChipL)),'uint8','ieee-le',...
%                     K.chip.probesetNb(ChipPos(ChipL)),K.chip.probesetNb(ChipPos(ChipL)),'uint8','ieee-le',...
%                     PsIndex{ChipL},Ps{ChipL}{CluL1}{ClassL1});
%                 %make a and c horizontal
%                 c=c';
%                 a=a';
%                 %order ps lines
%                 c=c(PsOrder{CluL1}{ClassL1},:);
%                 a=a(PsOrder{CluL1}{ClassL1},:);
%                 %order ps column
%                 for CluL2=1:length(CurrClu)
%                     for ClassL2=[2,3,4,6]
%                         ColOrder=[Limit{1}{CluL2}{ClassL2}(1):Limit{1}{CluL2}{ClassL2}(2)];
%                         c(:,Limit{1}{CluL2}{:,ClassL2}(1):Limit{1}{CluL2}{ClassL2}(2))=c(:,ColOrder(PsOrder{CluL2}{ClassL2}));
%                         a(:,Limit{1}{CluL2}{:,ClassL2}(1):Limit{1}{CluL2}{ClassL2}(2))=a(:,ColOrder(PsOrder{CluL2}{ClassL2}));
%                     end
%                 end
%                 %keep c horizontal
%                 %make a vertical
%                 a=a';
%                 %process all clu and save bmp files
%                 for CluL2=CluL1:length(CurrClu)
%                     for ClassL2=[2,3,4,6]
%                         CurrC=c(:,Limit{1}{CluL2}{ClassL2}(1):Limit{1}{CluL2}{ClassL2}(2));
%                         CurrA=a(Limit{1}{CluL2}{ClassL2}(1):Limit{1}{CluL2}{ClassL2}(2),:);
%                         if CluL2==CluL1&ClassL2==ClassL1
%                             %make composit bloc (C and A values)
%                             CurrA=triu(CurrC,1)+tril(CurrA,-1);
%                         end
%                         MaxVal=max(max(CurrA));
%                         Colors=colors(colormap,MaxVal+1);
%                         FileName=sprintf('%s_clu%u_class%u_vs_clu%u_class%u.bmp',OutName{ChipL},CluL1,ClassL1,CluL2,ClassL2);
%                         imwrite(CurrA,Colors,FileName,'bmp')           
%                         if CluL2~=CluL1|ClassL1~=ClassL2
%                             MaxVal=max(max(CurrC));
%                             FileName=sprintf('%s_clu%u_class%u_vs_clu%u_class%u.bmp',OutName{ChipL},CluL2,ClassL2,CluL1,ClassL1);
%                             imwrite(CurrC,Colors,FileName,'bmp')
%                         end
%                     end
%                 end
%             end
%         end
%         toc
%     end
% 
% 
% 
%     for ChipL=1:2
%         FileNameList{ChipL}={};
%         for CluL1=1:length(CurrClu)
%             for ClassL1=[2,3,4,6]
%                 for CluL2=CluL1:length(CurrClu)
%                     for ClassL2=[2,3,4,6]
%                         FileName=sprintf('%s_clu%u_class%u_vs_clu%u_class%u.bmp',OutName{ChipL},CluL1,ClassL1,CluL2,ClassL2);
%                         FileNameList{ChipL}{end+1,1}=FileName;
%                         if CluL2~=CluL1|ClassL1~=ClassL2
%                             FileName=sprintf('%s_clu%u_class%u_vs_clu%u_class%u.bmp',OutName{ChipL},CluL2,ClassL2,CluL1,ClassL1);
%                             FileNameList{ChipL}{end+1,1}=FileName;
%                         end
%                     end
%                 end
%             end
%         end
%     end
% 
% 
%     FileOrder=zeros(length(CurrClu)*4,length(CurrClu)*4);
%     Rank=0;
%     for CluL1=1:length(CurrClu)
%         for ClassL1=1:4
%             for CluL2=CluL1:length(CurrClu)
%                 for ClassL2=1:4
%                     Rank=Rank+1;
%                     FileOrder((CluL1-1)*4+ClassL1,(CluL2-1)*4+ClassL2)=Rank;K
%                     if CluL2~=CluL1|ClassL1~=ClassL2
%                         Rank=Rank+1;
%                         FileOrder((CluL2-1)*4+ClassL2,(CluL1-1)*4+ClassL1)=Rank;
%                     end
%                 end
%             end
%         end
%     end
%     FileOrder=FileOrder';
% 
% 
%     BlocNb=length(CurrClu)*4;
%     %BlocNb=5;
%     for ChipL=1:2
%         cd(fullfile(K.dir.net,sprintf('m%u',ChipRank(ChipL)),sprintf('n%u',Information{ChipL}.netRanks(1))))
%         fid=fopen(sprintf('%s_imagelist.txt',OutName{ChipL}),'w');
%         CurrFileOrder=FileOrder(1:BlocNb,1:BlocNb);
%         CurrFileNameList=FileNameList{ChipL}(CurrFileOrder(:));
%         for NameL=1:length(CurrFileNameList)
%             fprintf(fid,'%s\n',CurrFileNameList{NameL});
%         end
%         fclose(fid)
%         fid=fopen(sprintf('%s.bat',OutName{ChipL}),'w');
%         fprintf(fid,'montage @%s_imagelist.txt -tile %ux%u -geometry "10x10<" -mode Concatenate  %s.bmp',OutName{ChipL},BlocNb,BlocNb,OutName{ChipL})
%         fclose(fid)
%     end
%     DisplayCVM=questdlg('Do you want to display CVM','','yes','no','no');









