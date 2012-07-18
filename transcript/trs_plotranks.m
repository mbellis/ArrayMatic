function trs_plotranks()
global P DataRanks

if P.flag.loadData   
    DataRanks=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le');
end

if P.flag.testAlgo
    AlgoList=unique(P.point.algo);
    AlgoNb=length(AlgoList);
    for AlgoL=1:AlgoNb
        h(AlgoL)=figure;
        set(h(AlgoL),'name',sprintf('RANK RANK PLOTS - %s',AlgoList{AlgoL}));
        set(h(AlgoL),'color',[1,1,1])
    end
else
    %assumes that all points have been analyzed by the same algorithm
    AlgoList{1}=P.point.algo{1};
    AlgoNb=1;
    h(1)=figure;
    set(h(1),'name',sprintf('RANK RANK PLOTS - %s',AlgoList{1}));
    set(h(1),'color',[1,1,1])
end

%calculate the number of plots
PlotNb=0;
Bindex=zeros(P.point.nb,1);
Bindex(strmatch(AlgoList{1},P.point.algo,'exact'))=1;
for BiolL=1:P.biol.nb
    PointNb=length(find(Bindex&P.point.biolRank==BiolL));
    if PointNb>1
        PlotNb=PlotNb+(PointNb-1)*PointNb/2;
    end
end
ColNb=ceil(sqrt(PlotNb));
RawNb=floor(PlotNb/ColNb);
if PlotNb>RawNb*ColNb
    ColNb=ColNb+1;
end
% initialize corrCoeff field
P.biol.corrCoeff=cell(P.biol.nb,1);
%draw plots
Colors=colors(colormap,P.biol.nb);
cd(P.dir.resProp)
for AlgoL=1:AlgoNb
    Bindex=zeros(P.point.nb,1);
    Bindex(strmatch(AlgoList{AlgoL},P.point.algo,'exact'))=1;
    figure(h(AlgoL))    
    PlotRank=0;
    for BiolL=1:P.biol.nb
        PointPos=find(Bindex&P.point.biolRank==BiolL);
        %BiolL;
        PointNb=length(PointPos);
        if PointNb>1
            for PointL1=1:PointNb-1
                for PointL2=PointL1+1:PointNb               
                    PlotRank=PlotRank+1;
                    subplot(RawNb,ColNb,PlotRank)                    
                    plot(DataRanks(:,PointPos(PointL1)),DataRanks(:,PointPos(PointL2)),'.','color',Colors(BiolL,:),'markersize',3)
                    RankCorr=corrcoef(DataRanks(:,PointPos(PointL1)),DataRanks(:,PointPos(PointL2)));
                    RankCorr=RankCorr(1,2);
                    P.biol.corrCoeff{BiolL}=[P.biol.corrCoeff{BiolL},RankCorr];
                    set(gca,'box','on')
                    xlabel(P.point.name{PointPos(PointL1)})
                    ylabel(P.point.name{PointPos(PointL2)})            
                    title(sprintf('%.2f corr',RankCorr));
                    set(gca,'xlim',[0,100])
                    set(gca,'ylim',[0,100])
                    title(P.biol.name{BiolL})
                    xlabel(P.point.name{PointPos(PointL1)})
                    ylabel(P.point.name{PointPos(PointL2)})
                end
            end
        end
    end
end
for AlgoL=1:AlgoNb
    figure(h(AlgoL))
    set_figsize('1910px')
    saveas(h(AlgoL),sprintf('rankplot_%s_%s_%s',P.project.name,AlgoList{AlgoL},date),'png');
    delete(h(AlgoL))
end

'stop'

    
    