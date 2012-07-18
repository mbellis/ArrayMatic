%TRS_DUPPLLICATECORR calculate the correlation existing between pairs of points used in biological conditions

%INPUT PARAMETERS
%DisplayFlag: display figures

%EXTERNAL FILES
%DataRanks.float32le containing ranks

%OUTPUT PARAMETERS
%none


%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%
%                          c) Michel Bellis                                                %
%                          michel.bellis@crbm.cnrs.fr                                      %
%            Affiliation:  CNRS (Centre National de la Recherche Scientifique - France)    %                               
%  Bioinformatic Project:  ARRAYMATIC => http://code.google.com/p/arraymatic               %
%        Code Repository:  GITHUB =>http://github.com/mbellis                              %
%          Personal Page:  http://bns.crbm.cnrs.fr                                         %
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
%  THIS CODE IS DISTRIBUTED UNDER THE CeCILL LICENSE, WHICH IS COMPATIBLE WITH       %
%  THE GNU GENERAL PUBLIC LICENCE AND IN ACCORDANCE WITH THE EUROPEAN LEGISLATION.   %
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

function trs_duplicatecorr(DisplayFlag)
global DataRanks P

BiolNb=P.biol.nb;
PsNb=P.chip.currProbeSetNb;

% (re)initialize corrCoeff field
P.biol.corrCoeff=cell(P.biol.nb,1);
if isequal(P.par.analType,'network')
    for BiolL=1:BiolNb
        if length(P.biol.pairs)>=BiolL
            if ~isempty(P.biol.pairs{BiolL})
                if length(P.biol.pairs{BiolL})==2
                    PointRank1=P.biol.pairs{BiolL}(1);
                    PointRank2=P.biol.pairs{BiolL}(2);
                    Ranks1=load_data('DataRanks.float32le',P.dir.data,PsNb,P.point.nb,'single','ieee-le',1:PsNb,PointRank1);
                    Ranks2=load_data('DataRanks.float32le',P.dir.data,PsNb,P.point.nb,'single','ieee-le',1:PsNb,PointRank2);
                    RankCorr=corrcoef(Ranks1,Ranks2);
                    RankCorr=RankCorr(1,2);
                    P.biol.corrCoeff{BiolL}=RankCorr;
                    if DisplayFlag
                        hfig=figure;
                        title(sprintf('P%03.f vs P%03.f (rank)',PointRank2,PointRank1))
                        plot(Rank1.rank,Rank2.rank,'b.')
                        set(gca,'xlim',[0,100])
                        set(gca,'ylim',[0,100])
                        xlabel(sprintf('Point %03u - %s',PointRank1,strrep(P.point.name{PointRank1},'_',' ')))
                        ylabel(sprintf('Point %03u - %s',PointRank2,strrep(P.point.name{PointRank2},'_',' ')))
                        cd(P.dir.resProp)
                        h1=get(gca,'parent');
                        set(h1,'Color',[1 1 1])
                        set(h1,'name',['Reproducibility : ',P.biol.name{BiolL}])
                        SavedName=sprintf('biol%03u_point%04u_vs_point%04u_corr%u',BiolL,PointRank2,PointRank1,round(P.biol.corrCoeff{BiolL}*100));
                        saveas(h1,SavedName,'png')
                        delete(hfig)
                    end
                end
            end
        end
    end    
else
    for BiolL=1:BiolNb
        BiolPos=find(P.point.biolRank==BiolL);
        if ~isempty(BiolPos)
            PointPos=[];
            PointNb=0;
            for BiolPosL=1:length(BiolPos)
                if ~isempty(find(P.chip.pointIndex==BiolPos(BiolPosL)))
                    PointPos=[PointPos;BiolPos(BiolPosL)];
                    PointNb=PointNb+1;
                end
            end
            if PointNb>1
                for PointL1=1:PointNb-1
                    PointRank1=PointPos(PointL1);
                    if P.flag.loadData
                        Ranks1=load_data('DataRanks.float32le',P.dir.data,PsNb,P.point.nb,'single','ieee-le',1:PsNb,PointRank1);
                    else
                        Ranks1=DataRanks{:,PointRank1};
                    end
                    for PointL2=PointL1+1:PointNb
                        PointRank2=PointPos(PointL2);
                        if P.flag.loadData
                            Ranks2=load_data('DataRanks.float32le',P.dir.data,PsNb,P.point.nb,'single','ieee-le',1:PsNb,PointRank2);
                        else
                            Ranks2=DataRanks{PointRank2};
                        end
                        RankCorr=corrcoef(Ranks1,Ranks2);
                        RankCorr=RankCorr(1,2);
                        P.biol.corrCoeff{BiolL}=[P.biol.corrCoeff{BiolL},RankCorr];
                    end
                end
            end
        end
    end
end

CorrCoeff=[];
BiolPos=[];
for BiolL=1:P.biol.nb
    if ~isempty(P.biol.corrCoeff{BiolL})
    BiolPos=[BiolPos;repmat(BiolL,length(P.biol.corrCoeff{BiolL}),1)];
    CorrCoeff=[CorrCoeff,P.biol.corrCoeff{BiolL}];
    end
end

h=figure;
set(gcf,'color',[1,1,1])
set(h,'color',[1,1,1])
set(h,'name',strrep(sprintf('PAIRS CORRELATION %s %s',P.project.name,date),'_',' '))

subplot(1,3,1)
plot(sort(CorrCoeff),1/length(CorrCoeff):1/(length(CorrCoeff)+1):1,'b')
%bar([0.95:0.001:1],BinVal,'histc')
title(strrep(sprintf('pairs correlation %s',P.project.name),'_',' '))
xlabel('correlation coefficient')
ylabel('cdf')
set(gca,'xlim',[0.90,1])

subplot(1,3,2)
BinVal=histc(CorrCoeff,[floor(min(CorrCoeff)*100)/100:0.005:1]);
plot([floor(min(CorrCoeff)*100)/100:0.005:1],BinVal,'b')
hold on
plot([floor(min(CorrCoeff)*100)/100:0.005:1],BinVal,'g.')
%bar([0.95:0.001:1],BinVal,'histc')
title(strrep(sprintf('pairs correlation %s',P.project.name),'_',' '))
xlabel('correlation coefficient')
ylabel('frequency')

subplot(1,3,3)
BinVal=histc(CorrCoeff,[0.90:0.001:1]);
plot([0.90:0.001:1],BinVal,'b')
hold on
plot([0.90:0.001:1],BinVal,'g.')
%bar([0.95:0.001:1],BinVal,'histc')
title(strrep(sprintf('pairs correlation %s',P.project.name),'_',' '))
xlabel('correlation coefficient')
ylabel('frequency')

cd(P.dir.resCalib)
eval(sprintf('saveas(h,''pairs_correlation_%s_%s.png'',''png'')',P.project.name,date))

Limits=[0.8,0.9,0.92,0.94,0.96,0.98,1];
[BinVal,Pos]=histc(CorrCoeff,Limits);
for BinL=0:length(Limits)-1
    CurrPos=find(Pos==BinL);    
    CurrNb=length(CurrPos);
    if length(CurrPos)>30
        CurrPos=CurrPos(1:30);
    end    
    h=figure;
    set(h,'color',[1,1,1])
    if BinL==0
        set(h,'name',sprintf('rank plot for corr<=%.2f',Limits(1)))
        FigName=sprintf('rank_plot_maxcorr%.2f_%s.png',Limits(1),date);
        Title=sprintf('%u pairs with corr <=%.2f',CurrNb,Limits(1));
    else
        set(h,'name',sprintf('rank plot for corr>%.2f and corr<=%.2f',Limits(BinL),Limits(BinL+1)))
        FigName=sprintf('rank_plot_mincorr%.2f_maxcorr%.2f_%s.png',Limits(BinL),Limits(BinL+1),date);
        Title=sprintf('%u pairs with corr >%.2f and corr<=%.2f',CurrNb,Limits(BinL),Limits(BinL+1));
    end
    set(gcf,'color',[1,1,1])
    ColNb=ceil(sqrt(length(CurrPos)));
    LineNb=floor(length(CurrPos)/ColNb);
    if ColNb*LineNb<length(CurrPos)
        LineNb=LineNb+1;
    end
    for PosL=1:length(CurrPos)
        subplot(LineNb,ColNb,PosL)
        Ranks=P.biol.pairs{BiolPos(CurrPos(PosL))};
        Ranks1=load_data('DataRanks.float32le',P.dir.data,PsNb,P.point.nb,'single','ieee-le',1:PsNb,Ranks(1));
        Ranks2=load_data('DataRanks.float32le',P.dir.data,PsNb,P.point.nb,'single','ieee-le',1:PsNb,Ranks(2));
        plot(Ranks1,Ranks2,'b.','markersize',3);
        title(sprintf('biol cond %u - %.4f',BiolPos(CurrPos(PosL)),CorrCoeff(CurrPos(PosL))))
        xlabel(sprintf('point %u',Ranks(1)))
        ylabel(sprintf('point %u',Ranks(2))) 
        if PosL==1
            text(1,130,Title,'FontSize',20);
        end
    end
    cd(P.dir.resCalib)
    if  length(CurrPos)==30
        set_figsize('1280px')
    elseif length(CurrPos)>10
        set_figsize('1024px')
    end
    saveas(h,FigName,'png')
end

CorrLimits=inputdlg({'inferior correlation (>=)','superior correlation (<=)'},'',1,{'0.9','1'});
P.biol.corrLimits=[str2num(CorrLimits{1}),str2num(CorrLimits{2})];
cd(P.dir.project)
eval(sprintf('save %s P',P.project.name));