%recover CORR/ANTI values for networks constructed on transcription factors on m8
%(networks 56 to 109) and for random networks
%plot correlation of all probe sets with the TF probe set either in random networks or in FT
%network

ModelRank=8;
RandNetRanks=[7:21,24:38,39:53];
TfNetRanks=[56:109];
PsNb=K.chip.probeSetNbs{ModelRank};

TfNetNb=length(TfNetRanks);
RandNetNb=length(RandNetRanks);

CORR=cell(TfNetNb,1);
ANTI=cell(TfNetNb,1);

%RECOVER CORR/ANTI VALUES
TfNames={};
TfPsRanks=[];
NetDir=fullfile('/home/mbellis/array1/sosma/net',sprintf('m%03u',ModelRank));
%for TfNetL=1:TfNetNb
for TfNetL=53:TfNetNb
    TfNetL
    tic
    TfNetRank=TfNetRanks(TfNetL);
    TfNetPos=find(K.net{ModelRank}.rank==TfNetRank);
    NetName=K.net{ModelRank}.name{TfNetPos};
    PsPos=findstr(NetName,' ps');
    TfPsRank=str2num(NetName(PsPos+3:end));
    TfName=NetName(1:PsPos-1);
    TfNames{end+1,1}=TfName;
    TfPsRanks(end+1,1)=TfPsRank;
    CORR{TfNetL}=zeros(RandNetNb+1,PsNb);
    ANTI{TfNetL}=zeros(RandNetNb+1,PsNb);
    %load data of random networks
    for RandNetL=1:RandNetNb      
        RandNetRank=RandNetRanks(RandNetL);
        cd(fullfile(NetDir,sprintf('n%05u',RandNetRank)))
        CORR{TfNetL}(RandNetL,:)=load_data(sprintf('c_m%u_n%u.4mat',ModelRank,RandNetRank),'./',PsNb,PsNb,'uint8','ieee-le',[1:PsNb],TfPsRank);
        ANTI{TfNetL}(RandNetL,:)=load_data(sprintf('a_m%u_n%u.4mat',ModelRank,RandNetRank),'./',PsNb,PsNb,'uint8','ieee-le',[1:PsNb],TfPsRank);
    end
    %load data of tf network
    try
        cd(fullfile(NetDir,sprintf('n%05u',TfNetRank)))
        CORR{TfNetL}(RandNetNb+1,:)=load_data(sprintf('c_m%u_n%u.4mat',ModelRank,TfNetRank),'./',PsNb,PsNb,'uint8','ieee-le',[1:PsNb],TfPsRank);
        ANTI{TfNetL}(RandNetNb+1,:)=load_data(sprintf('a_m%u_n%u.4mat',ModelRank,TfNetRank),'./',PsNb,PsNb,'uint8','ieee-le',[1:PsNb],TfPsRank);
    catch
    end
    toc
end

[TfNames,SortOrder]=sort(TfNames);
TfPsRanks=TfPsRanks(SortOrder);
CORR=CORR(SortOrder);
ANTI=ANTI(SortOrder);

%LOAD PsMatrix
NetRanks=[7:21];
NetNb=15;
cd('/home/mbellis/sosma/data/psawn/mldata/mouse/m8_1p')
eval(sprintf('load m%u_n%u_netnb%u_probenb1_newps_stat_netprc100_pvcorr1',ModelRank,NetRanks(1),NetNb));
%LOAD PROBE SET IDS
cd('/home/mbellis/sosma/data/psawn/rawdata');
PsNames=textread('m8_probeset.txt','%s','delimiter','\t');


%LOAD TF INFORMATION
cd('/home/mbellis/sosma/data/tf/rawdata')
[UsedTF,UsedGeneName,UsedPsRank]=textread(sprintf('TF_m%u_transcript_ids.txt',ModelRank),'%s%s%u','delimiter','\t');
%[TF,EnsId,GeneName]=textread(sprintf('TF_Gene_interactions_090925_EnsIds.txt',ModelRank),'%s%s%s','delimiter','\t');
[TF,EnsId]=textread(sprintf('TF_Gene_interactions.txt',ModelRank),'%s%s','delimiter','\t');


%FIND TARGETS
TfTargPsRanks=cell(length(TfNames),1);
PsRanksClass=cell(length(TfNames),1);
for TfL=1:length(TfNames)
    %find correspondant transcript id
    TfPos=find(UsedPsRank==TfPsRanks(TfL));
    if ~isempty(TfPos)
        %use the corresponding transcript ID
        CurrTrsId=UsedTF{TfPos(1)};
        %to find the targeted genes
        TfPos=strmatch(CurrTrsId,TF,'exact');
        length(TfPos)
        if ~isempty(TfPos)
            for PosL=1:length(TfPos)
                %verify if the curent gene exit on the chip
                GenePos=strmatch(EnsId{TfPos(PosL)},Genes.name,'exact');
                if ~isempty(GenePos)
                    %find the probe set that have been assigned to this gene
                    CurrPsRanks=find(PsMatrix(:,1)==GenePos);
                    for PsL=1:length(CurrPsRanks)
                        TfTargPsRanks{TfL}(end+1,:)=[GenePos,CurrPsRanks(PsL)];           
                        PsRanksClass{TfL}(end+1)=PsMatrix(CurrPsRanks(PsL),17);
                    end
                end
            end
            TfTargPsRanks{TfL}=unique(TfTargPsRanks{TfL},'rows');
        end
    end
end




Colors=colors(colormap,15);
h=figure;
set(h,'name','CORR')
set(gcf,'color',[1,1,1])
for p=29:54
    subplot(7,4,p-28)
    if ~isempty(CORR{p})
        [CorrVal,SortOrder]=sort(CORR{p}(46,:));
        hold on
        NullCorrPerc=zeros(15,1);
        PosCorrPerc=zeros(15,1);
        NullAntiPerc=zeros(15,1);
        PosAntiPerc=zeros(15,1);
        CorrPerc=zeros(15,1);
        AntiPerc=zeros(15,1);
        for i=1:15
            plot(CORR{p}(i,SortOrder),'marker','.','linestyle','none','color',Colors(i,:),'markersize',3)
            NullCorrPerc(i)=length(find(CORR{p}(i,:)>ANTI{p}(i,:)&CORR{p}(i,:)>0&CORR{p}(46,:)==0))*100/length(find(CORR{p}(i,:)>ANTI{p}(i,:)&CORR{p}(i,:)>0));
            NullAntiPerc(i)=length(find(CORR{p}(i,:)<ANTI{p}(i,:)&CORR{p}(i,:)>0&CORR{p}(46,:)==0))*100/length(find(CORR{p}(i,:)<ANTI{p}(i,:)&CORR{p}(i,:)>0));
            PosCorrPerc(i)=length(find(CORR{p}(i,:)==0&CORR{p}(46,:)>0&CORR{p}(46,:)>ANTI{p}(46,:)))*100/length(find(CORR{p}(46,:)>0&CORR{p}(46,:)>ANTI{p}(46,:)));
            PosAntiPerc(i)=length(find(CORR{p}(i,:)==0&CORR{p}(46,:)>0&CORR{p}(46,:)<ANTI{p}(46,:)))*100/length(find(CORR{p}(46,:)>0&CORR{p}(46,:)<ANTI{p}(46,:)));
            CorrPerc(i)=length(find(CORR{p}(i,:)>ANTI{p}(i,:)))*100/PsNb;
            AntiPerc(i)=length(find(CORR{p}(i,:)<ANTI{p}(i,:)))*100/PsNb;
        end
        plot(CorrVal,'k-')
        %plot(CorrVal,'k.')

        CurrCorr=CORR{p}([1:15,46],SortOrder);
        %find position corresponding to targets
        CurrPs=zeros(1,size(TfTargPsRanks{p}),1);
        for PsL=1:size(TfTargPsRanks{p},1)
            CurrPs(PsL)=find(SortOrder==TfTargPsRanks{p}(PsL,2));
        end
        MaxCorr=max(CurrCorr(1:15,CurrPs));
        plot(CurrPs,MaxCorr,'r+')
        plot(CurrPs,CurrCorr(16,CurrPs),'ro')
        set(gca,'box','on')
        set(gca,'xticklabel','')
        set(gca,'yticklabel','')
        title(sprintf('%s ps%u (%u)',TfNames{p},TfPsRanks(p),size(TfTargPsRanks{p},1)))
        xlabel(sprintf('%u%%=>0=>%u%% - %u%%=>0=>%u%%',round(mean(NullCorrPerc)),round(mean(PosCorrPerc)),round(mean(NullAntiPerc)),round(mean(PosAntiPerc))))
        ylabel(sprintf('%u%% %u%% - %u%% %u%%',round(mean(CorrPerc)),round(mean(AntiPerc)),round(length(find(CORR{p}(46,:)>ANTI{p}(46,:)))*100/PsNb),round(length(find(CORR{p}(46,:)<ANTI{p}(46,:)))*100/PsNb)));
    end
end

h=figure;
set(h,'name','ANTI')
set(gcf,'color',[1,1,1])
for p=1:54
    subplot(9,6,p)
    if ~isempty(ANTI{p})
        [AntiVal,SortOrder]=sort(ANTI{p}(46,:));
        hold on
        NullPerc=zeros(15,1);
        PosPerc=zeros(15,1);
        for i=1:15
            plot(ANTI{p}(i,SortOrder),'marker','.','linestyle','none','color',Colors(i,:),'markersize',3)
            NullPerc(i)=length(find(ANTI{p}(i,:)>0&ANTI{p}(46,:)==0))*100/length(find(ANTI{p}(i,:)>0));
            PosPerc(i)=length(find(ANTI{p}(i,:)==0&ANTI{p}(46,:)>0))*100/length(find(ANTI{p}(46,:)>0));            
        end
        plot(AntiVal,'k-')
        %plot(AntiVal,'k.')
        CurrAnti=ANTI{p}([1:15,46],SortOrder);
        %find position corresponding to targets
        CurrPs=zeros(1,size(TfTargPsRanks{p}),1);
        for PsL=1:size(TfTargPsRanks{p},1)
            CurrPs(PsL)=find(SortOrder==TfTargPsRanks{p}(PsL,2));
        end
        MaxAnti=max(CurrAnti(1:15,CurrPs));
        plot(CurrPs,MaxAnti,'b+')
        plot(CurrPs,CurrAnti(16,CurrPs),'bo')                
        set(gca,'box','on')
        set(gca,'xticklabel','')
        set(gca,'yticklabel','')
        title(sprintf('%s ps%u (%u)',TfNames{p},TfPsRanks(p),size(TfTargPsRanks{p},1)))
        xlabel(sprintf('%u%%=>0=>%u%%',round(mean(NullPerc)),round(mean(PosPerc))))
    end
end


h=figure;
set(h,'name','TFvsRand CORR')
set(gcf,'color',[1,1,1])
for p=1:54
    subplot(9,6,p)
    if ~isempty(CORR{p})        
        hold on
        CorrNb=zeros(1,15);
        AntiNb=zeros(1,15);
        for i=1:15
            plot(CORR{p}(i,:),CORR{p}(46,:),'marker','.','linestyle','none','color',Colors(i,:),'markersize',6)
            CorrNb(i)=length(find(CORR{p}(i,:)>ANTI{p}(i,:)));
            AntiNb(i)=length(find(CORR{p}(i,:)<ANTI{p}(i,:)));
        end
        plot(CORR{p}(1,TfTargPsRanks{p}),CORR{p}(46,TfTargPsRanks{p}),'ro')
        set(gca,'box','on')
        set(gca,'xticklabel','')
        set(gca,'yticklabel','')
        set(gca,'ylim',[0,100])
        title(sprintf('%s ps%u (%u)',TfNames{p},TfPsRanks(p),size(TfTargPsRanks{p},1)))
        xlabel(sprintf('%u c %u a',round(mean(CorrNb)),round(mean(AntiNb))))
        ylabel(sprintf('%u c %u a',length(find(CORR{p}(46,:)>ANTI{p}(46,:))),length(find(ANTI{p}(46,:)>CORR{p}(46,:)))))
    end
end

h=figure;
set(h,'name','EXAMPLE CORR')
set(gcf,'color',[1,1,1])
PlotRank=0;
%PlotRanks=[1,5,9,13,17,2,6,10,14,18,3,7,11,15,19,4,8,12,16,20];
PlotRanks=[1,2,3,4];
%for p=[19,5,48,2]
for p=[16,10,14,5]
%for p=[2,16,27,44,49,1,4,6,7,8,13,14,25,31,48,5,19,26,38]
    PlotRank=PlotRank+1;
    %for TypeL=1:2
    for TypeL=1
        %CORR & ANTI
        if TypeL==1
            VAL=CORR;
            %subplot(4,4,PlotRank);
            subplot(2,2,PlotRank);
            %subplot(5,4,PlotRanks(PlotRank));
        else
            VAL=ANTI;
            subplot(4,4,8+PlotRank);
        end
        [Val,SortOrder]=sort(VAL{p}(46,:));
        hold on
        NullPerc=zeros(15,1);
        PosPerc=zeros(15,1);
        Perc=zeros(15,1);
        for i=1:15
            %plot(VAL{p}(i,SortOrder),'marker','.','linestyle','none','color',Colors(i,:),'markersize',3)
            plot(VAL{p}(i,SortOrder),'marker','.','linestyle','none','color',Colors(i,:))
            NullPerc(i)=length(find(VAL{p}(i,:)>0&VAL{p}(46,:)==0))*100/length(find(VAL{p}(i,:)>0));
            PosPerc(i)=length(find(VAL{p}(i,:)==0&VAL{p}(46,:)>0))*100/length(find(VAL{p}(46,:)>0));            
            Perc(i)=length(find(VAL{p}(i,:)>0))*100/PsNb;
        end
        plot(Val,'k-')
        plot(Val,'k.')

%         CurrVal=VAL{p}([1:15,46],SortOrder);
%         %find position corresponding to targets
%         CurrPs=zeros(1,size(TfTargPsRanks{p}),1);
%         for PsL=1:size(TfTargPsRanks{p},1)
%             CurrPs(PsL)=find(SortOrder==TfTargPsRanks{p}(PsL,2));
%         end
%         MaxCorr=max(CurrVal(1:15,CurrPs));
%         if TypeL==1
%             plot(CurrPs,MaxCorr,'mv','markersize',12,'linewidth',2)
%             plot(CurrPs,CurrVal(16,CurrPs),'mo','markersize',12,'linewidth',2)
%         else
%             plot(CurrPs,MaxCorr,'bx','markersize',10,'linewidth',2)
%             plot(CurrPs,CurrVal(16,CurrPs),'bo','markersize',10,'linewidth',2)
%         end
        set(gca,'box','on')
        set(gca,'xticklabel','')
        set(gca,'yticklabel','')
        %title(sprintf('Corr with TF %s (ps%u - %u targets)',TfNames{p},TfPsRanks(p),size(TfTargPsRanks{p},1)))
        title(sprintf('TF %s (ps%u)',TfNames{p},TfPsRanks(p)))
        xlabel(sprintf('%u%% p->0 - %u%% 0->p',round(mean(NullPerc)),round(mean(PosPerc))),'fontsize',12)
        ylabel(sprintf('%u%% in random networks\n%u%% in TF network (black)',round(mean(Perc)),round(length(find(Val))*100/PsNb)))
        
        %TF vs RAND
%         if TypeL==1        
%             subplot(4,4,4+PlotRank);
%         else
%             subplot(4,4,12+PlotRank);
%         end
% 
%         hold on
%         for i=1:15
%             %plot(VAL{p}(i,:),VAL{p}(46,:),'marker','.','linestyle','none','color',Colors(i,:),'markersize',3)
%             plot(VAL{p}(i,:),VAL{p}(46,:),'marker','.','linestyle','none','color',Colors(i,:))
%         end
%         plot(VAL{p}(1,TfTargPsRanks{p}),VAL{p}(46,TfTargPsRanks{p}),'ro')
%         set(gca,'box','on')
%         set(gca,'xticklabel','')
%         set(gca,'yticklabel','')
%         title(sprintf('%s ps%u (%u)',TfNames{p},TfPsRanks(p),size(TfTargPsRanks{p},1)))
    end
end
