%==========================%
% FUNCTION NET_COMPAREMCLS %
%==========================%

%NET_COMPAREMCLS: compares several results of MCL clustering
%INPUT PARAMETERS
% 1 ModelRanks: list of two chip ranks

% net_comparemcls(8)
 
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


function net_comparemcls(ChipRank,OutRank)
%net_comparemcls([8,27],[4,3])
%net_comparemcls([8,27],[1,1])
%net_comparemcls([27,27],[1,2])
%net_comparemcls([27],[2])
global K

if length(ChipRank)>2
    h=errordlg('net_comparemcls need at most two chips');
    waitfor(h)
    error('process canceled')
end

PsNb=[];
Information={};
ChipNb=length(ChipRank);
for ChipL=1:ChipNb
    CmlDir=fullfile(K.dir.net,sprintf('m%u',ChipRank(ChipL)),'mcl');
    cd(CmlDir)
    eval(sprintf('load m%u_mcl%u.mat',ChipRank(ChipL),OutRank(ChipL)))
    Information{ChipL}=Info;
    Cluster{ChipL}=Clu;
    ClusterNb{ChipL}=CluNb;
    if ChipRank(ChipL)==8
        if min(Information{ChipL}.listRank)>=7
            PsNb(ChipL)=22690;
        else
            PsNb(ChipL)=K.chip.probesetNb(ChipRank(ChipL));
        end
    else
        PsNb(ChipL)=K.chip.probesetNb(ChipRank(ChipL));
    end
end
clear Info Clu CluNb

if ~isempty(find(ChipRank==8))
    ListName={'all','ss','ms','cx','s','rest','allc','ssc','msc','cxc','sc','restc'};
else
    ListName={'all','ss','ms','cx','s','rest'};
end
ListNb=length(ListName);
TypeName={'C','C-A','rawC','raw(C-A)'};




if ChipNb==1
    ListNb(1)=Information{1}.listNb;
    if ListNb==6
        PlotPos{1}=[1:6];
        RowNb=2;
        ColNb=3;
    else
        PlotPos{1}=[1,3,5,7,9,11,2,4,6,8,10,12];
        RowNb=6;
        ColNb=2;
    end
else
    ListNb(1)=Information{1}.listNb;
    ListNb(2)=Information{2}.listNb;
    if ListNb(1)==6
        PlotPos{1}=[1,3,5,7,9,11];
        PlotPos{2}=[2,4,6,8,10,12];
        RowNb=6;
        ColNb=2;
    else
        PlotPos{1}=[1,4,7,10,13,16,2,5,8,11,14,17];
        PlotPos{2}=[3,6,9,12,15,18];
        RowNb=6;
        ColNb=3;
    end
end

%load ps lists
for ChipL=1:ChipNb
    NetNb(ChipL)=length(Information{ChipL}.netRanks);
    ListSize{ChipL}=zeros(ListNb(ChipL),1);
    PsRankList{ChipL}=cell(ListNb(ChipL),1);
    cd(K.dir.chip)
    for ListL=1:ListNb(ChipL)
        fid=fopen(sprintf('m%u_pslist%u.u32',ChipRank(ChipL),Information{ChipL}.listRank(ListL)),'r','ieee-le');
        PsRankList{ChipL}{ListL}=fread(fid,inf,'uint32');
        ListSize{ChipL}(ListL)=length(PsRankList{ChipL}{ListL});
        fclose(fid);        
    end    
end


cd(K.dir.affyMetadata)
FileName=sprintf('m%u_m%u_commonps.mat',min(ChipRank),max(ChipRank));
if exist(FileName,'file')
    load(FileName)
    if ChipRank(1)>ChipRank(2)
        Temp=ComPsRank;
        ComPsRank(:,1)=ComPsRank(:,2);
        ComPsRank(:,2)=Temp(:,1);
        clear Temp
    end
    if ~isempty(ChipRank==8)&~isempty(ChipRank==27)
        if min(Information{1}.listRank)>=7
            MatchPsRank=ComPsRank(:,2);
        else
            MatchPsRank=zeros(max(ComPsRank(:,1)),1);
            MatchPsRank(ComPsRank(:,1))=ComPsRank(:,2);
        end
    else
        MatchPsRank=zeros(max(ComPsRank(:,1)),1);
        MatchPsRank(ComPsRank(:,1))=ComPsRank(:,2);
    end
else
    h=errordlg(sprintf('no correspondance file between m%y and m%u',ChipRank(1),ChipRank(2)));
    waitfor(h)
    error('process canceled')
end

%% STATISTICS ACCORDING TO TYPE, LIST, LIMITS ...

%% NB OF PROCESSED NETWORKS ACCORDING TO TYPE
h=figure;
if ChipNb==1
    set(h,'name',sprintf('NB OF PROCESSED NETWORKS IN m%u BY TYPE',ChipRank(1)))
else
    set(h,'name',sprintf('NB OF PROCESSED NETWORKS IN m%u and m%uBY TYPE',ChipRank(1),ChipRank(2)))
end
for ChipL=1:ChipNb
    CluNbs={};
    hs=zeros(1,4);
        Limits=[];
    for TypeL=1:4
        Limits=union(Limits,Information{ChipL}.limits{TypeL});
    end
    for TypeL=1:4                
        set(gcf,'color',[1,1,1])
        Colors=colors(colormap,4);
        LimitNb=length(Information{ChipL}.limits{TypeL});
        NetNbs{TypeL}=[];
        if ~isempty(ClusterNb{ChipL}{TypeL})
            for ListL=1:ListNb(ChipL)
                try
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
                        title(sprintf('%s',ListName{ListL}))
                        xlabel('limit')
                        ylabel('net nb')
                        set(gca,'box','on') 
                        set(gca,'xlim',[0,Limits(end)])
                        set(gca,'xgrid','on')
                        set(gca,'ygrid','on')
                        set(gca,'ylim',[0,NetNb(ChipL)])
                    end
                catch
                end
            end
        end
    end        
end
legend(hs,TypeName)



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
    for TypeL=1:4
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
                    title(sprintf('%s',ListName{ListL}))
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
    for TypeL=1:4
        for ListL=1:ListNb(ChipL)
            subplot(RowNb,ColNb,PlotPos{ChipL}(ListL))
            set(gca,'ylim',[0,min(1500,max([CluNbs{1}(ListL,:),CluNbs{2}(ListL,:),CluNbs{3}(ListL,:),CluNbs{4}(ListL,:)]))])
        end
    end

end
legend(hs,TypeName)

%% NB OF CLUSTERS ACCORDING TO LIST
for ChipL=1:ChipNb
    h=figure;
    set(h,'name',sprintf('NB OF CLUSTERS IN m%u BY LIST',ChipRank(ChipL)))
    for TypeL=1:4
        set(gcf,'color',[1,1,1])
        Colors=colors(colormap,ListNb(ChipL));
        %Colors=colors(colormap,6);
        LimitNb=length(Information{ChipL}.limits{TypeL});
        if TypeL==4
            hs=zeros(1,ListNb(ChipL));
        end
        subplot(2,2,TypeL)
        hs=zeros(1,ListNb(ChipL));
        hold on
        if ~isempty(ClusterNb{ChipL}{TypeL})
            %for ListL=1:ListNb(ChipL)
            for ListL=1:6
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
                        if TypeL==4
                            hs(ListL)=h;
                        end
                        plot(Information{ChipL}.limits{TypeL},CluNb./CurrNetNb,'+','color',Colors(ListL,:));
                    end
                catch
                end
            end
        end
        title(sprintf('%s',TypeName{TypeL}))
        xlabel('limit')
        ylabel('cluster nb')
        set(gca,'box','on')
        set(gca,'xgrid','on')
        set(gca,'ygrid','on')
    end   
    %legend(hs,ListName)
    legend(hs,ListName(1:6))
    for TypeL=1:4
        subplot(2,2,TypeL)
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
CluRange=[1:10];
%merge all limits
AllLimits=[];
for ChipL=1:ChipNb
    for TypeL=1:4
        AllLimits=union(AllLimits,Information{ChipL}.limits{TypeL});
    end
end
LimitNb=length(AllLimits);
Colors=colors(colormap,LimitNb);
LimitNames=cell(1,LimitNb);
for LimitL=1:LimitNb
    LimitNames{LimitL}=num2str(AllLimits(LimitL));
end


for TypeL=1:4
    h=figure;
    set(h,'name',sprintf('CLUSTER SIZE IN m%u - TYPE %s',ChipRank(ChipL),TypeName{TypeL}))
    set(gcf,'color',[1,1,1])

    for ChipL=1:ChipNb
        LimitNb=length(Information{ChipL}.limits{TypeL});
        hs=zeros(1,LimitNb);        
        if ~isempty(ClusterNb{ChipL}{TypeL})
            for ListL=1:ListNb(ChipL)
                %for ListL=7
                try
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
                        title(sprintf('%s (%u ps)',ListName{ListL},ListSize{ChipL}(ListL)))
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
                catch
                end
            end
        end
    end
end

%% INTRA CHIP REPRODUCIBILITY ACCORDING TO TYPE
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
    set(h,'name',sprintf('INTRA m%u CHIP REPRODUCIBILITY OF CLU %u ACCORDING TO CORR TYPE',ChipRank(ChipL),CLU))
    Colors=colors(colormap,4);
    for ChipL=ChipRange
        hs=zeros(1,4);
        Limits=[];
        for TypeL=1:4
            Limits=union(Limits,Information{ChipL}.limits{TypeL});
        end
        NetFreq{ChipL}=cell(1,4);
        for TypeL=1:4
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
                        title(sprintf('%s',ListName{ListL}))
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
        legend(hs,TypeName)
    end

    for TypeL=1:4
        h=figure;
        set(gcf,'color',[1,1,1])
        set(h,'name',sprintf('%s: INTRA m%u CHIP NETWORK REPRODUCIBILITY OF CLU %u ACCORDING TO TYPE %s',TypeName{TypeL},ChipRank(ChipL),CLU,TypeName{TypeL}))
        for ChipL=ChipRange
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
                title(ListName{ListL})
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



    LIMIT=2;
    h=figure;
    set(gcf,'color',[1,1,1])
    set(h,'name',sprintf('INTRA m%u CHIP NETWORK REPRODUCIBILITY OF CLU %u ACCORDING TO LIST TYPE FOR LIMIT %u',ChipRank(ChipL),CLU,Information{ChipL}.limits{TypeL}(LIMIT)))

    ChipL=1;
    Colors=colors(colormap,6);
    for TypeL=1:4
        subplot(ChipNb,4,TypeL)
        hold on

        for ListL=1:6

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

        %set(gca,'xgrid','on')
        %set(gca,'ygrid','on')
    end
    legend(ListName(1:6))

    % for TypeL=1:4
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
    % legend(ListName(7:12))
    if ChipNb==2
        ChipL=2;
        for TypeL=1:4
            subplot(2,4,8+TypeL)
            hold on

            for ListL=1:6

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

            %set(gca,'xgrid','on')
            %set(gca,'ygrid','on')
        end
        legend(ListName(1:6))
    end



    %% REPRODUCIBILITY MATRIX

    for ChipL=1:ChipNb
        PlotNb=NetNb(ChipL)*(NetNb(ChipL)-1)/2;
        ColNbMat=round(sqrt(PlotNb/2));
        RowNbMat=round(PlotNb/ColNbMat);
        if PlotNb>RowNbMat*ColNbMat
            RowNbMat=RowNbMat+1;
        end
        for TypeL=1:4
            for ListL=1:ListNb(ChipL)
                h=figure;
                set(h,'name',sprintf('INTRA m%u CHIP REPRODUCIBILITY MATRIX ACCORDING TO TYPE %s AND LIST %s',ChipRank(ChipL),TypeName{TypeL},ListName{ListL}))
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
            end %of ListL
        end %of TypeL
    end





    %% INTRA CHIP REPRODUCIBILITY ACCORDING TO LIST

    for ChipL=1:ChipNb
        h=figure;
        CLU=1;
        set(h,'name',sprintf('INTRA m%u CHIP REPRODUCIBILITY OF CLU %u ACCORDING TO LIST',ChipRank(ChipL),CLU))
        Colors=colors(colormap,6);
        hs=zeros(1,ListNb(ChipL));
        Limits=[];
        for TypeL=1:4
            Limits=union(Limits,Information{ChipL}.limits{TypeL});
        end
        for TypeL=1:4
            CurrLimits=Information{ChipL}.limits{TypeL};
            set(gcf,'color',[1,1,1])
            LimitNb=length(Information{ChipL}.limits{TypeL});
            subplot(2,2,TypeL)
            hold on
            if ~isempty(ClusterNb{ChipL}{TypeL})
                for ListL=1:6
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

            legend(hs,ListName(1:6))
        end %TypeL

        for ListL=1:ListNb(ChipL)
            h=figure;
            set(gcf,'color',[1,1,1])
            set(h,'name',sprintf('INTRA m%u CHIP NETWORK REPRODUCIBILITY OF CLU %u ACCORDING TO LIST %s',ChipRank(ChipL),CLU,ListName{ListL}))
            for TypeL=1:4
                LimitNb=length(Information{ChipL}.limits{TypeL});
                LimitNames=cell(1,LimitNb);
                for LimitL=1:LimitNb
                    LimitNames{LimitL}=num2str(Information{ChipL}.limits{TypeL}(LimitL));
                end
                Colors=colors(colormap,LimitNb);

                subplot(2,2,TypeL)
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

    end %ChipL


    %% CONSENSUS REGIONS
    %constructs union of the same region in all the networks at different level of
    %reproducibility


    for ChipL=1:ChipNb
        Region{ChipL}={};
        for TypeL=1:4
            Region{ChipL}{TypeL}={};
            for ListL=1:ListNb(ChipL)
                Region{ChipL}{TypeL}{ListL}={};
                for RegL=1:6
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

    % REGIONS 1 TO 6
    for RegL=1:6
        %REGION=1;
        LIMIT=1;
        MAXREP=10;
        Colors=colors(colormap,4);
        h=figure;
        set(h,'name',sprintf('SIZE OF REGION %u ACCORDING TO UNION THRESHOLD FOR LIMIT %u',RegL,LIMIT));
        set(gcf,'color',[1,1,1])
        TypeName1={'A all','C-A all','rawC all','raw(C-A) all','A core','C-A core','rawC core','raw(C-A) core'}
        hs=zeros(1,8);
        for TypeL=1:4
            for ChipL=1:ChipNb
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
                    title(sprintf('%s',ListName{ListL}))
                    ylabel('#ps')
                    xlabel('>= #networks')
                    set(gca,'box','on')
                    set(gca,'xgrid','on')
                    set(gca,'ygrid','on')
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
    for TypeL=1:4
        for ChipL=1:ChipNb
            for ListL=1:ListNb(ChipL)
                subplot(RowNb,ColNb,PlotPos{ChipL}(ListL))
                hold on
                ReliquatSize=zeros(1,MAXREP);
                ReliquatCoreSize=zeros(1,MAXREP);
                for RepL=1:MAXREP
                    RegPs=[];
                    RegCorePs=[];
                    for RegL=1:6
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
                title(sprintf('%s',ListName{ListL}))
                ylabel('#ps')
                xlabel('>= #networks')
                set(gca,'box','on')
                set(gca,'xgrid','on')
                set(gca,'ygrid','on')
            end
        end
    end
    legend(hs,TypeName1)




    %% OVERLAP BETWEEN CONSENSUS REGIONS
    for ChipL=1:ChipNb
        for TypeL=1:4
            for ListL=1:ListNb(ChipL)
                for RegL1=1:6
                    for LimitL=1:2
                        Region{ChipL}{TypeL}{ListL}{RegL1}{LimitL}.core=cell(1,10);
                        for RepL=1:10
                            CurrPs=Region{ChipL}{TypeL}{ListL}{RegL1}{LimitL}.ps{RepL};
                            %find intersection
                            InterReg=[];
                            for RegL2=1:6
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

    REGION=4;
    LIMIT=1;
    MAXREP=10;
    Colors=colors(colormap,4);
    h=figure;
    set(h,'name',sprintf('OVERLAP BETWEEN REGION %u AND OTHER REGIONS FOR LIMIT %u',REGION,LIMIT));
    set(gcf,'color',[1,1,1])
    hs=zeros(1,4);
    for TypeL=1:4
        for ChipL=1:ChipNb
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
                title(sprintf('%s',ListName{ListL}))
                ylabel('Overlap')
                xlabel('>= #networks')
                set(gca,'box','on')
                set(gca,'xgrid','on')
                set(gca,'ygrid','on')
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
        for TypeL=1:4
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
                    title(sprintf('%s vs %s',ListName{ListL+Offset},ListName{1}))
                    ylabel('reproducibility')
                    xlabel('>= #networks')
                    set(gca,'box','on')
                    set(gca,'xgrid','on')
                    set(gca,'ygrid','on')
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
                    title(sprintf('%s vs %s',ListName{ListL+Offset},ListName{5}))
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





    %% INTER CHIP REPRODUCIBILITY    (EXIST OVERLAP BETWEENS REGIONS CONSTRUCTED BY UNION OF
    % SEVERAL NETWORKS)

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
        for TypeL=1:4
            LoopNb=1;
            ColNb=1;
            if ListNb(1)==12
                ColNb=2;
                LoopNb=2;
            end
            for LoopL=1:LoopNb
                for ListL=1:6
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
                            for RegL=1:6
                                if OVERLAPFLAG
                                    CurrRegRep=length(intersect(MatchPsRank(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.core{RepL}),Region{2}{TypeL}{ListL}{RegL}{LIMIT}.core{RepL}))*100/min(length(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.core{RepL}),length(Region{2}{TypeL}{ListL}{RegL}{LIMIT}.core{RepL}));
                                else
                                    CurrRegRep=length(intersect(MatchPsRank(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.ps{RepL}),Region{2}{TypeL}{ListL}{RegL}{LIMIT}.ps{RepL}))*100/min(length(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.ps{RepL}),length(Region{2}{TypeL}{ListL}{RegL}{LIMIT}.ps{RepL}));
                                end
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
                                    RegRep(RepL)=length(intersect(MatchPsRank(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.core{RepL}),Region{2}{TypeL}{ListL}{CurrSimReg}{LIMIT}.core{RepL}))*100/min(length(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.core{RepL}),length(Region{2}{TypeL}{ListL}{CurrSimReg}{LIMIT}.core{RepL}));
                                else
                                    RegRep(RepL)=length(intersect(MatchPsRank(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.ps{RepL}),Region{2}{TypeL}{ListL}{CurrSimReg}{LIMIT}.ps{RepL}))*100/min(length(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.ps{RepL}),length(Region{2}{TypeL}{ListL}{CurrSimReg}{LIMIT}.ps{RepL}));
                                end
                            end
                        end
                    else
                        for RepL=1:MAXREP
                            if OVERLAPFLAG
                                RegRep(RepL)=length(intersect(MatchPsRank(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.core{RepL}),Region{2}{TypeL}{ListL}{REGION2}{LIMIT}.core{RepL}))*100/min(length(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.core{RepL}),length(Region{2}{TypeL}{ListL}{REGION2}{LIMIT}.core{RepL}));
                            else
                                RegRep(RepL)=length(intersect(MatchPsRank(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.ps{RepL}),Region{2}{TypeL}{ListL}{REGION2}{LIMIT}.ps{RepL}))*100/min(length(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.ps{RepL}),length(Region{2}{TypeL}{ListL}{REGION2}{LIMIT}.ps{RepL}));
                            end
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
                    title(sprintf('%s',ListName{ListL+Offset}))
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
    for TypeL=1:4
        LoopNb=1;
        ColNb=1;
        if ListNb(1)==12
            ColNb=2;
            LoopNb=2;
        end
        for LoopL=1:LoopNb
            for ListL=1:6
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
                        for RegL=1:6
                            CurrRegRep=length(intersect(MatchPsRank(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.ps{RepL}),Region{2}{TypeL}{ListL}{RegL}{LIMIT}.ps{RepL}))*100/min(length(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.ps{RepL}),length(Region{2}{TypeL}{ListL}{RegL}{LIMIT}.ps{RepL}));
                            if CurrRegRep>MaxRep
                                CurrSimReg=RegL;
                                MaxRep=CurrRegRep;
                            end
                        end
                        %RegRep(1,end+1)=length(intersect(MatchPsRank(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.ps{RepL}),Region{2}{TypeL}{ListL}{REGION2}{LIMIT}.ps{RepL}))*100/min(length(Region{1}{TypeL}{ListL+Offset}{REGION1}{LIMIT}.ps{RepL}),length(Region{2}{TypeL}{ListL}{REGION2}{LIMIT}.ps{RepL}));
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
                title(sprintf('%s',ListName{ListL+Offset}))
                ylabel('most similar')
                xlabel('>= #networks')
                set(gca,'box','on')
                set(gca,'xgrid','on')
                set(gca,'ygrid','on')
            end
        end

    end

    legend(hs,TypeName)




    %% RANDOM MODULES

    TYPE=1;
    REGION=5;
    LIMIT=1;
    REPROD=5;
    LIST=5;
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

if isempty(ChipRange)    
    AllLimits=[];
for ChipL=1:ChipNb
    for TypeL=1:4
        AllLimits=union(AllLimits,Information{ChipL}.limits{TypeL});
    end
end
LimitNb=length(AllLimits);
Colors=colors(colormap,LimitNb);
LimitNames=cell(1,LimitNb);
for LimitL=1:LimitNb
    LimitNames{LimitL}=num2str(AllLimits(LimitL));
end
   %only one network for each chip
   MinSize=10;
   for ListL=1:6
       h=figure;
       set(h,'color',[1,1,1])
       set(h,'name',sprintf('Liste %u (%s)',ListL,ListName{ListL}))
       for TypeL=1:2
           for LimitL=1:LimitNb
               MaxCluPos=[0,0];
               for ChipL=1:2
                   CurrLimitPos=find(Information{ChipL}.limits{TypeL}==AllLimits(LimitL));
                   if ~isempty(CurrLimitPos)
                       LimitPos(ChipL)=CurrLimitPos;
                       try
                           CurrClusterNb=ClusterNb{ChipL}{TypeL}{ListL}{1}{LimitPos(ChipL)};
                           if ~isempty(CurrClusterNb)
                               CurrMaxCluPos=find(CurrClusterNb>=MinSize);
                               MaxCluPos(ChipL)=CurrMaxCluPos(end);
                           end
                       catch
                       end
                   end
               end
               if MaxCluPos(1)>0&MaxCluPos(2)>0
                   Reprod=zeros(MaxCluPos(1)+1,MaxCluPos(2)+1);
                   for CluL1=1:MaxCluPos(1)
                       Ps1=MatchPsRank(find(Cluster{1}{TypeL}{ListL}{1}(:,LimitPos(1))==CluL1));
                       for CluL2=1:MaxCluPos(2)
                           Ps2=find(Cluster{2}{TypeL}{ListL}{1}(:,LimitPos(2))==CluL2);
                           Reprod(CluL1,CluL2)=round(length(intersect(Ps1,Ps2))*100/min(length(Ps1),length(Ps2)));
                       end
                   end

                   subplot(LimitNb,2,(LimitL-1)*2+TypeL)
                   h=pcolor(Reprod);
                   if MaxCluPos(1)*MaxCluPos(2)>1000
                       set(h,'linestyle','none')
                   end
                   xlabel(sprintf('%u',Information{2}.limits{TypeL}(LimitPos(2))))
                   ylabel(sprintf('%u',Information{1}.limits{TypeL}(LimitPos(1))))
               end
           end
       end
   end
end
