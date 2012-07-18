function trs_plotalgos()

global P M DataRanks


REPLICATE=1;
POINT_NAME=P.point.name{1};
BindexN=zeros(P.point.nb,1);
BindexN(strmatch(POINT_NAME,P.point.name,'exact'))=1;
PsNb=P.chip.currProbeSetNb;




if P.flag.testAlgo
    
    cd(P.dir.data)
    load Comp
    AlgoList=unique(P.point.algo);
    AlgoNb=length(AlgoList);
    Colors=colors(colormap,AlgoNb);
    %find comparison names
    CompName={};
    for CompL=1:length(M{1}.compName)
        AlgoPos=findstr(M{1}.compName{CompL},['_',AlgoList{1}]);
        if ~isempty(AlgoPos)
            CompName{end+1,1}=M{1}.compName{CompL}(1:AlgoPos-1);
        end
    end
    Comparison=listdlg('liststring',CompName,'selectionmode','single','promptstring','select a comparison');
                    
    %CompNb=length(CompName);
    
    
    
    %calculate the number of plots
    PlotNb=AlgoNb;
    ColNb=ceil(sqrt(PlotNb));
    RawNb=floor(PlotNb/ColNb);
    if PlotNb>RawNb*ColNb
        ColNb=ColNb+1;
    end

    
    %COMPARE CURVES
    %for CompL=1:CompNb
        PlotRank=0;
        h=figure;
        set(h,'name',sprintf('%s: FDR,PV,S vs ZVAR CURVES COMPARISON BETWEEN ALGOS',CompName{Comparison}))

        for AlgoL=1:AlgoNb
            Fdr=load_data('Fdr_02.float32le',P.dir.data,PsNb,PsNb,'single','ieee-le',1:PsNb,(Comparison-1)*AlgoNb+AlgoL);
            ZVar=load_data('ZVar_02.float32le',P.dir.data,PsNb,PsNb,'single','ieee-le',1:PsNb,(Comparison-1)*AlgoNb+AlgoL);
            Sensitivity=load_data('Sensitivity_02.float32le',P.dir.data,PsNb,PsNb,'single','ieee-le',1:PsNb,(Comparison-1)*AlgoNb+AlgoL);
            Pv=load_data('Pv_02.float32le',P.dir.data,PsNb,PsNb,'single','ieee-le',1:PsNb,(Comparison-1)*AlgoNb+AlgoL);
            [ZVar,SortOrder]=sort(ZVar);
            Fdr=Fdr(SortOrder);
            Sensitivity=Sensitivity(SortOrder);
            Pv=Pv(SortOrder);
            Neg=find(ZVar<0);
            Pos=find(ZVar>=0);

            PlotRank=PlotRank+1;
            subplot(RawNb,ColNb,PlotRank)
            plot(ZVar(Neg),Fdr(Neg),'r')
            hold on
            plot(ZVar(Pos),Fdr(Pos),'r')
            plot(ZVar(Neg),Pv(Neg),'g')
            plot(ZVar(Pos),Pv(Pos),'g')
            plot(ZVar(Neg),Sensitivity(Neg),'m')
            plot(ZVar(Pos),Sensitivity(Pos),'m')
            x=title(sprintf('%s - %u Inc - %u Dec',AlgoList{AlgoL},length(find(Fdr<=0.10&Fdr>0)),length(find(Fdr>=-0.10&Fdr<0))));            
            set(x,'FontSize',[16])
            set(gca,'xlim',[-100,100])
            set(gca,'ylim',[-1,1])
            
            if AlgoL==1
                text(-50,1.2,sprintf('%s: FDR,PV,S vs ZVAR CURVES Comparison BETWEEN ALGOS',strrep(CompName{Comparison},'_',' ')));
            end
        end
        set(gcf,'color',[1,1,1])
    %end



    %calculate the number of plots
    PlotNb=(AlgoNb-1)*AlgoNb/2;
    ColNb=ceil(sqrt(PlotNb));
    RawNb=floor(PlotNb/ColNb);
    if PlotNb>RawNb*ColNb
        ColNb=ColNb+1;
    end

    % COMPARE RANKS
    %draw plots
    
    cd(P.dir.resCalib)
    PlotRank=0;
    h=figure;
    set(h,'name',sprintf('RANK COMPARISON BETWEEN ALGOS FOR %s',POINT_NAME))
    for AlgoL1=1:AlgoNb-1
        Bindex1=zeros(P.point.nb,1);
        Bindex1(strmatch(AlgoList{AlgoL1},P.point.algo,'exact'))=1;
        PointPos1=find(BindexN&Bindex1&P.point.replicateRank==REPLICATE)
        DataRanks1=load_data('DataRanks.float32le',P.dir.data,PsNb,PsNb,'single','ieee-le',1:PsNb,PointPos1);
        for AlgoL2=AlgoL1+1:AlgoNb
            if AlgoL1==1&AlgoL2==2
                title(sprintf('RANK COMPARISON BETWEEN ALGOS FOR %s',POINT_NAME))
            end
            Bindex2=zeros(P.point.nb,1);
            Bindex2(strmatch(AlgoList{AlgoL2},P.point.algo,'exact'))=1;
            PointPos2=find(BindexN&Bindex2&P.point.replicateRank==REPLICATE)
            PlotRank=PlotRank+1;
            subplot(RawNb,ColNb,PlotRank)            
            try
            DataRanks2=load_data('DataRanks.float32le',P.dir.data,PsNb,PsNb,'single','ieee-le',1:PsNb,PointPos2);
            catch
                'stop'
            end
            plot(DataRanks1,DataRanks2,'.','color',Colors(AlgoL1,:),'markersize',3)
            if AlgoL==1
                title(sprintf('RANK COMPARISON BETWEEN ALGOS FOR %s',POINT_NAME))
            end
            x=xlabel(AlgoList{AlgoL1});
            y=ylabel(AlgoList{AlgoL2});
            set(x,'color',Colors(AlgoL1,:))
            set(x,'FontSize',[16])
            set(y,'color',Colors(AlgoL2,:))
            set(y,'FontSize',[16])
            set(gca,'xlim',[0,100])
            set(gca,'ylim',[0,100])
        end
    end
    set(gcf,'color',[1,1,1])
    

    

    %COMPARE COMPARISONS
    h1=figure;
    set(h1,'name','FDR');
    h2=figure;
    set(h2,'name','ZVAR');
    h3=figure;
    set(h3,'name','FC');
    PlotRank=0;
    for AlgoL1=1:AlgoNb-1
        Fdr1=load_data('Fdr_02.float32le',P.dir.data,PsNb,PsNb,'single','ieee-le',1:PsNb,(Comparison-1)*AlgoNb+AlgoL1);
        ZVar1=load_data('ZVar_02.float32le',P.dir.data,PsNb,PsNb,'single','ieee-le',1:PsNb,(Comparison-1)*AlgoNb+AlgoL1);
        FC1=load_data('Fc_02.float32le',P.dir.data,PsNb,PsNb,'single','ieee-le',1:PsNb,(Comparison-1)*AlgoNb+AlgoL1);
        for AlgoL2=AlgoL1+1:AlgoNb
            Fdr2=load_data('Fdr_02.float32le',P.dir.data,PsNb,PsNb,'single','ieee-le',1:PsNb,(Comparison-1)*AlgoNb+AlgoL2);
             ZVar2=load_data('ZVar_02.float32le',P.dir.data,PsNb,PsNb,'single','ieee-le',1:PsNb,(Comparison-1)*AlgoNb+AlgoL2);
            FC2=load_data('Fc_02.float32le',P.dir.data,PsNb,PsNb,'single','ieee-le',1:PsNb,(Comparison-1)*AlgoNb+AlgoL2);
            PlotRank=PlotRank+1;
            
            figure(h1)            
            subplot(RawNb,ColNb,PlotRank)
            plot(Fdr1,Fdr2,'.','color',Colors(AlgoL1,:),'markersize',3)
            if AlgoL1==1&AlgoL2==2
                title(sprintf('FDR COMPARISON BETWEEN ALGOS FOR %s',strrep(CompName{Comparison},'_',' ')))
            end
            x=xlabel(AlgoList{AlgoL1});
            y=ylabel(AlgoList{AlgoL2});
            set(x,'color',Colors(AlgoL1,:))
            set(x,'FontSize',[16])
            set(y,'color',Colors(AlgoL2,:))
            set(y,'FontSize',[16])
            set(gca,'xlim',[-0.1,0.1])
            set(gca,'ylim',[-0.1,0.1])
            
            figure(h2)            
            subplot(RawNb,ColNb,PlotRank)
            plot(ZVar1,ZVar2,'.','color',Colors(AlgoL1,:),'markersize',3)
            if AlgoL1==1&AlgoL2==2
                title(sprintf('ZVAR COMPARISON BETWEEN ALGOS FOR %s',strrep(CompName{Comparison},'_',' ')))
            end
            x=xlabel(AlgoList{AlgoL1});
            y=ylabel(AlgoList{AlgoL2});
            set(x,'color',Colors(AlgoL1,:))
            set(x,'FontSize',[16])
            set(y,'color',Colors(AlgoL2,:))
            set(y,'FontSize',[16])
            set(gca,'xlim',[-100,100])
            set(gca,'ylim',[-100,100])
            
            figure(h3)
            subplot(RawNb,ColNb,PlotRank)
            plot(FC1,FC2,'.','color',Colors(AlgoL1,:),'markersize',3)
            if AlgoL1==1&AlgoL2==2
                title(sprintf('FC COMPARISON BETWEEN ALGOS FOR %s',strrep(CompName{Comparison},'_',' ')))
            end
            x=xlabel(AlgoList{AlgoL1});
            y=ylabel(AlgoList{AlgoL2});
            set(x,'color',Colors(AlgoL1,:))
            set(x,'FontSize',[16])
            set(y,'color',Colors(AlgoL2,:))
            set(y,'FontSize',[16])          
            set(gca,'xlim',[-10,10])
            set(gca,'ylim',[-10,10])
        end
    end
    set(h1,'color',[1,1,1])
    set(h2,'color',[1,1,1])
    set(h3,'color',[1,1,1])
end


%FC vs ZVAR

PlotNb=AlgoNb;
ColNb=ceil(sqrt(PlotNb));
RawNb=floor(PlotNb/ColNb);
if PlotNb>RawNb*ColNb
    ColNb=ColNb+1;
end

h=figure;
set(h,'name','FC VS ZVAR');
PlotRank=0;
for AlgoL=1:AlgoNb
    PlotRank=PlotRank+1;
    ZVar=load_data('ZVar_02.float32le',P.dir.data,PsNb,PsNb,'single','ieee-le',1:PsNb,(Comparison-1)*AlgoNb+AlgoL);
    FC=load_data('Fc_02.float32le',P.dir.data,PsNb,PsNb,'single','ieee-le',1:PsNb,(Comparison-1)*AlgoNb+AlgoL);

    subplot(RawNb,ColNb,PlotRank)
    plot(ZVar,FC,'.','color',Colors(AlgoL,:),'markersize',3)
    if AlgoL==1
        text(-50,1.2,sprintf('FDR COMPARISON BETWEEN ALGOS FOR %s',strrep(CompName{Comparison},'_',' ')))
    end
    title(AlgoList{AlgoL});
    xlabel('ZVar')
    ylabel('FC')
    set(gca,'xlim',[-100,100])
    set(gca,'ylim',[-10,10])

end
set(h,'color',[1,1,1])


