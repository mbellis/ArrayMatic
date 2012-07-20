function net_importmcl(ModelRank,NetRank)
global K

PsNb=K.chip.probeSetNbs{ModelRank}(1);
CmlDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',NetRank),'mcl');
cd(CmlDir)
Dir=dir;
PlotRank=0;
Legend={};
CorrLimits=[];
Clu=[];
CluNb={};
Pos=0;
for FileL=3:length(Dir)
    CurrFile=Dir(FileL).name;
    if ~isempty(findstr('out',CurrFile))
        Pos=Pos+1;
        CorrLimit=regexp(CurrFile,'c[0-9]+','match');
        CorrLimits=[CorrLimits;str2num(CorrLimit{1}(2:end))];
        PlotRank=PlotRank+1;
        Legend=[Legend,CorrLimit{1}];
        CurrClu=zeros(PsNb,1);
        CurrCluNb=[];
        Rank=0;
        InFid=fopen(CurrFile,'r');
        while 1
            Rank=Rank+1;
            CurrPsRank=fgetl(InFid);
            if ~ischar(CurrPsRank), break, end
            CurrPsRank=str2num(CurrPsRank);
            CurrClu(CurrPsRank)=Rank;
            CurrCluNb(Rank,1)=length(CurrPsRank);
        end
        fclose(InFid);
        Clu(:,Pos)=CurrClu;
        CluNb{Pos,1}=CurrCluNb;
    end
end

h=figure;
set(h,'name',sprintf('FIG1 - DISTRIBUTION OF MCL CLUSTER SIZES FOR M%uN%u',ModelRank,NetRank))
set(gcf,'color',[1,1,1])
MclNb=length(CluNb);
Colors=colors(colormap,MclNb);
hs=[];
for MclL=1:MclNb
    h=plot(1:length(CluNb{MclL}),CluNb{MclL},'color',Colors(MclL,:));
    hs=[hs;h];
    hold on
    plot(1:length(CluNb{MclL}),CluNb{MclL},'marker','+','color',Colors(MclL,:))
end
set(gca,'xlim',[0,10])
set(gca,'ylim',[0,max(CluNb{1})+round(max(CluNb{1})/10)])
legend(hs,Legend)
title(sprintf('distribution of MCL cluster sizes for m%un%u',ModelRank,NetRank))

MatFile=sprintf('m%un%u_mcl.mat',ModelRank,NetRank);
eval(sprintf('save %s Clu CluNb CorrLimits',MatFile))