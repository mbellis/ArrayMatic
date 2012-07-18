% TRS_DISPLAYPROPERTIES - Display point properties

% c) Michel Bellis
% arraymatic@gmail.com

% NO INPUT

function trs_displayproperties
global K P

RunLimit=300;
FigLimit=500;

%process only points that have been loaded at the initiation step of
%the project
UsedPoint=find(P.point.used);

h(1)=figure;
set(gcf,'color',[1,1,1])
%set(h1,'name',sprintf('Properties of dataset %',DataSetName))


% number of nan
if P.point.nb<FigLimit
    subplot(2,2,1)
end
hold on
plot(P.point.nanNb,'r')
plot(UsedPoint,P.point.nanNb(UsedPoint),'k+','markersize',3)
title('# absent S')
ylabel('# probesets')
%xlabel('experimental point rank')
set(gca,'ylim',[0,500])
set(gca,'box','on')

% number of diff(signal)==0 et longest run size
if P.point.nb<FigLimit
    subplot(2,2,2)
else
    h(2)=figure;
    set(gcf,'color',[1,1,1])
end
hold on
plot(P.point.runNb,'r')
if isfield(P.point,'diversity')
    plot(P.point.diversity,'b')
    plot(UsedPoint,P.point.diversity(UsedPoint),'k+','markersize',3)
end
plot(UsedPoint,P.point.runNb(UsedPoint),'k+','markersize',3)
%plot(P.point.delta(:),'b')
title('longest run size (r) & diversity (b) & ')
ylabel('# probesets')
set(gca,'box','on')

% min/max
if P.point.nb<FigLimit
    subplot(2,2,3)
else
    h(3)=figure;
    set(gcf,'color',[1,1,1])
end
hold on
MaxS=P.point.maxSignal;
MaxS(find(MaxS>38000))=38000;
MinS=P.point.minSignal;
MinS(find(MinS<-8000))=-8000;
plot(MaxS,'r')
plot(-MinS,'b')
plot(UsedPoint,MaxS(UsedPoint),'k+','markersize',3)
plot(UsedPoint,-MinS(UsedPoint),'k+','markersize',3)
title('max(S) (r) & min(S) (b)')
ylabel('min - max')
set(gca,'box','on')

% number of values=run
if P.point.nb<FigLimit
    subplot(2,2,4)
else
    h(4)=figure;
    set(gcf,'color',[1,1,1])
end
hold on
Value=P.point.negNb;
plot(Value,'r')
Value=P.point.nullNb;
plot(Value,'g')
Value=P.point.negNb+P.point.nullNb;
Index1=find(P.point.runNb<RunLimit|P.point.runVal==0);
%Index1=find(Value>0);

Bindex=zeros(P.point.nb,1);
Bindex1=Bindex;
Bindex1(UsedPoint)=1;
Bindex2=Bindex;
Bindex2(Index1)=1;
Index2=find(Bindex1&Bindex2);
plot(Index2,Value(Index2),'k+','markersize',3)

Index1=find(P.point.runNb<RunLimit|P.point.runVal==0);
Value=P.point.runNb+P.point.negRunNb;
Value(Index1)=0;
plot(Value,'b')
Index1=find(P.point.runNb>=RunLimit&P.point.runVal~=0);
Bindex=zeros(P.point.nb,1);
Bindex1=Bindex;
Bindex1(UsedPoint)=1;
Bindex2=Bindex;
Bindex2(Index1)=1;
Index2=find(Bindex1&Bindex2);
plot(Index2,Value(Index2),'k+','markersize',3)
title(sprintf('# S<=0 (r) & # S<= value(run|size>%u) (b)',RunLimit))
ylabel('# probesets')
set(gca,'box','on')

DivFactor=50;


for SubL=1:4
    if P.point.nb<FigLimit
        subplot(2,2,SubL)
    else
        figure(h(SubL))
    end
    YLim=get(gca,'ylim');
    if YLim(2)<500
        set(gca,'ylim',[YLim(2),500])
    end
end
for SubL=1:4
    if P.point.nb<FigLimit
        subplot(2,2,SubL)
    else
        figure(h(SubL))
    end
    YLim=get(gca,'ylim');
    YDel=-YLim(2)/DivFactor;
    set(gca,'ylim',[YDel*6,YLim(2)])
    xlabel('experiment rank')
end

% delimitation of experiments
Pos=1;

FirstPos=1;
Color='g';
XTick=[];
XTickLabel={};
XTickPos=0;
ExpRank=P.point.expRank;
CurrExp=ExpRank(1);
PointNb=P.point.nb;
while Pos<PointNb
    while ExpRank(Pos)==CurrExp & Pos<PointNb
        Pos=Pos+1;
    end
    if Pos<PointNb
        EndPos=Pos-1;
    else
        if ExpRank(Pos)==CurrExp
            EndPos=Pos;
        else
            EndPos=Pos-1;
        end
    end
    XTick=[XTick,FirstPos+(EndPos-FirstPos)/2];
    XTickPos=XTickPos+1;
    if mod(XTickPos,5)==0
        XTickLabel{XTickPos,1}=sprintf('%u',CurrExp);
    else
        XTickLabel{XTickPos,1}='';
    end

    for SubL=[1:4]
        if P.point.nb<FigLimit
            subplot(2,2,SubL)
        else
            figure(h(SubL))
        end
        YLim=get(gca,'YLim');
        YDel=-YLim(2)/DivFactor;
        hold on
        %plot([FirstPos,EndPos],[100,100],'color',Color,'Marker','.')
        %plot([FirstPos:EndPos],ones(1,EndPos-FirstPos+1)*100,'color',Color,'Marker','.')
        line([[FirstPos:EndPos];[FirstPos:EndPos]],[ones(1,EndPos-FirstPos+1)*round(YDel*5);ones(1,EndPos-FirstPos+1)*round(YDel*6)],'color',Color)
        %patch([FirstPos,FirstPos,EndPos,EndPos],[round(YDel*5),round(YDel*6),round(YDel*6),round(YDel*5)],Color)
        %plot([FirstPos,EndPos],[100,100],'color',Color)
    end
    FirstPos=EndPos+1;
    CurrExp=ExpRank(Pos);
    if Color=='g'
        Color='y';
    else
        Color='g';
    end
end




% delimitation of biol conditions
Answer=questdlg(sprintf('Do you want to show the %u biological condition (may be long)',P.biol.nb),'','yes','no','yes');
if isequal(Answer,'yes')
    BiolCond=P.point.biolRank;
    Pos=1;
    CurrBiol=P.point.biolRank(1);
    FirstPos=1;
    Color='k';
    while Pos<PointNb
        while BiolCond(Pos)==CurrBiol & Pos<PointNb
            Pos=Pos+1;
        end
        if Pos<PointNb
            EndPos=Pos-1;
        else
            if BiolCond(Pos)==CurrBiol
                EndPos=Pos;
            else
                EndPos=Pos-1;
            end
        end
        for SubL=[1:4]
            if P.point.nb<FigLimit
                subplot(2,2,SubL)
            else
                figure(h(SubL))
            end

            YLim=get(gca,'YLim');
            YDel=-YLim(2)/DivFactor;
            hold on
            %plot([FirstPos:EndPos],ones(1,EndPos-FirstPos+1)*200,'color',Color,'Marker','.')
            line([[FirstPos:EndPos];[FirstPos:EndPos]],[ones(1,EndPos-FirstPos+1)*YDel*4;ones(1,EndPos-FirstPos+1)*YDel*5],'color',Color)
            %plot([FirstPos,EndPos],[200,200],'color',Color,'Marker','.')
            %plot([FirstPos,EndPos],[200,200],'color',Color)
        end
        FirstPos=EndPos+1;
        CurrBiol=BiolCond(Pos);
        if Color=='w'
            Color='k';
        else
            Color='w';
        end
    end
end

Answer=questdlg('Select the algorithm determination method (by name = MAS4, MAS5,... / by number={1=Mas4 & 2=non Mas4)','','by name','by number','by number');
if isequal(Answer,'by name')
    AlgoNames=unique(P.point.algo);
    Colors='rygcbmk';
    ColorNb=length(Colors);
    if ~isempty(AlgoNames)
        for SubL=[1:4]
            if P.point.nb<FigLimit
                subplot(2,2,SubL)
            else
                figure(h(SubL))
            end
            YLim=get(gca,'YLim');
            YDel=-YLim(2)/DivFactor;
            hold on
            for AlgoL=1:length(AlgoNames)
                AlgoPos=strmatch(AlgoNames{AlgoL},P.point.algo,'exact');
                line([AlgoPos;AlgoPos],[ones(size(AlgoPos))*YDel*3;ones(size(AlgoPos))*YDel*4],'color',Colors(max(1,mod(AlgoL,ColorNb))))
            end
            YLim=get(gca,'ylim');
            set(gca,'ylim',[YDel*6,YLim(2)])
            set(gca,'xtick',XTick)
            set(gca,'xticklabel',XTickLabel)
            set(gca,'tickdir','out')
        end
    end
else
    Mas4=find(P.point.algoGrp==1);
    Mas5=find(P.point.algoGrp==2);

    for SubL=[1:4]
        if P.point.nb<FigLimit
            subplot(2,2,SubL)
        else
            figure(h(SubL))
        end

        YLim=get(gca,'YLim');
        YDel=-YLim(2)/DivFactor;
        if ~isempty(Mas5)
            line([Mas5;Mas5],[ones(size(Mas5))*YDel*4;ones(size(Mas5))*YDel*5],'color','m')
        end
        hold on
        if ~isempty(Mas4)
            line([Mas4;Mas4],[ones(size(Mas4))*YDel*3;ones(size(Mas4))*YDel*4],'color','c')
        end
        YLim=get(gca,'ylim');
        set(gca,'ylim',[YDel*6,YLim(2)])
        set(gca,'xtick',XTick)
        set(gca,'xticklabel',XTickLabel)
        set(gca,'tickdir','out')
    end
end

if P.point.nb<FigLimit
    set_figsize('1024px')
end

