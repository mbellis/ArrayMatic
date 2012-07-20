%=======================%
% FUNCTION NET_COMPCORR %
%=======================%

%NET_COMPCORR compare CORR and ANTI values between several networks

%INPUT PARAMETERS


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
%calculate mean corr, anti and connectivity

%net_compcorr(27,[88:108],1250)
%net_compcorr(27,[25:45],1250)
function net_comprank(ChipRank,NetRanks,PsRank)
global K
PsNb=K.chip.probesetNb(ChipRank);
NetNb=length(NetRanks);
Corr=zeros(NetNb,PsNb);
Anti=zeros(NetNb,PsNb);
NormFlag=1;
TypeNb=1;
for NetL=1:NetNb
    CurrNet=NetRanks(NetL);
    CmlDir=fullfile(K.dir.net,sprintf('m%u',ChipRank),sprintf('n%u',CurrNet));
    cd(CmlDir)
    Corr(NetL,:)=load_data(sprintf('c_m%u_n%u.4mat',ChipRank,CurrNet),'./',PsNb,PsNb,'int8','ieee-le',[1:PsNb],PsRank);
    Anti(NetL,:)=load_data(sprintf('a_m%u_n%u.4mat',ChipRank,CurrNet),'./',PsNb,PsNb,'int8','ieee-le',[1:PsNb],PsRank);
    if NormFlag
        try
            Norm(NetL,:)=load_data(sprintf('f_m%u_n%u.4mat',ChipRank,CurrNet),'./',PsNb,PsNb,'int8','ieee-le',[1:PsNb],PsRank);
        catch
            NormFlag=0;
            TypeNb=1;
        end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    end
end



for TypeL=1:TypeNb
    CompNb=NetNb*(NetNb-1)/2;
    Colors=colors(colormap,CompNb);
    h=figure;
    set(gcf,'color',[1,1,1])
    if TypeL==1
        set(h,'name','CORR and ANTI PLOTS');
        CurrCorr=Corr;
        CurrAnti=Anti;
    else
        set(h,'name','RAW(CORR) and RAW(ANTI) PLOTS');
        CurrCorr=double(Corr).*double(Norm)/100;
        CurrAnti=double(Anti).*double(Norm)/100;
    end
    CurrComp=0;
    for NetL1=1:NetNb-1
        for NetL2=NetL1+1:NetNb
            CurrComp=CurrComp+1;
            subplot(3,1,1)
            hold on
            plot(CurrCorr(NetL1,:),CurrCorr(NetL2,:),'.','color',Colors(CurrComp,:),'markersize',3)
            subplot(3,1,2)
            hold on
            plot(CurrAnti(NetL1,:),CurrAnti(NetL2,:),'.','color',Colors(CurrComp,:),'markersize',3)
            subplot(3,1,3)
            hold on
            plot(CurrCorr(NetL1,:)-CurrAnti(NetL1,:),CurrCorr(NetL2,:)-CurrAnti(NetL2,:),'.','color',Colors(CurrComp,:),'markersize',3)
        end
    end
    for PlotL=1:3
        subplot(3,1,PlotL)
        set(gca,'box','on')
        if PlotL==1
            title(sprintf('CORR plots for chip m%u Ps %u and nets %u to %u',ChipRank,PsRank,NetRanks(1),NetRanks(end)))
        elseif PlotL==2
            title(sprintf('ANTI plots for chip m%u Ps %u and nets %u to %u',ChipRank,PsRank,NetRanks(1),NetRanks(end)))
        else
            title(sprintf('CORR-ANTI plots for chip m%u Ps %u and nets %u to %u',ChipRank,PsRank,NetRanks(1),NetRanks(end)))
        end

    end
end

for TypeL=1:TypeNb
    h=figure;
    set(gcf,'color',[1,1,1])
    if TypeL==1
        set(h,'name','CORR and ANTI DISTRIBUTIONS')
        CurrCorr=Corr;
        CurrAnti=Anti;
        Range1=0:10:100;
        Range2=-100:10:100;
    else
        set(h,'name','RAW(CORR) and RAW(ANTI) DISTRIBUTIONS')
        CurrCorr=round(double(Corr).*double(Norm)/100);
        CurrAnti=round(double(Anti).*double(Norm)/100);
        Range1=0:5:50;
        Range2=-50:5:50;
    end

    for NetL1=1:NetNb-1
        for NetL2=NetL1+1:NetNb
            Colors=colors(colormap,11);
            %for Val=100:-10:0
            ColorPos=0;
            for Val=Range1
                ColorPos=ColorPos+1;
                ValPos=find(CurrCorr(NetL1,:)==Val);
                %ValPos=find(CurrCorr(NetL1,:)>=Range1(ValL)&CurrCorr(NetL1,:)<Range1(ValL+1));
                if ~isempty(ValPos)
                    Values=CurrCorr(NetL2,ValPos);
                    UValues=unique(Values);
                    Freq=histc(Values,UValues);
                    subplot(3,1,1)
                    hold on
                    plot(UValues,Freq./length(Values)+0.05*Val/10,'color',Colors(ColorPos,:))
                end
                ValPos=find(CurrAnti(NetL1,:)==Val);
                %ValPos=find(CurrAnti(NetL1,:)>=Range1(ValL)&CurrAnti(NetL1,:)<Range1(ValL+1));
                if ~isempty(ValPos)
                    Values=CurrAnti(NetL2,ValPos);
                    UValues=unique(Values);
                    Freq=histc(Values,UValues);
                    subplot(3,1,2)
                    hold on
                    plot(UValues,Freq./length(Values)+0.05*Val/10,'color',Colors(ColorPos,:))
                end
            end
            Colors=colors(colormap,21);
            ColorPos=0;
            %for Val=100:-10:-100
            for Val=Range2
                ColorPos=ColorPos+1;
                ValPos=find(CurrCorr(NetL1,:)-CurrAnti(NetL1,:)==Val);
                %ValPos=find(CurrCorr(NetL1,:)-CurrAnti(NetL1,:)>=Range1(ValL)&CurrCorr(NetL1,:)-CurrAnti(NetL1,:)<Range1(ValL+1));
                if ~isempty(ValPos)
                    Values=CurrCorr(NetL2,ValPos)-CurrAnti(NetL2,ValPos);
                    UValues=unique(Values);
                    Freq=histc(Values,UValues);
                    subplot(3,1,3)
                    hold on
                    plot(UValues,Freq./length(Values)+0.05*(ColorPos-1),'color',Colors(ColorPos,:))
                end
            end
        end
    end
    for PlotL=1:3
        subplot(3,1,PlotL)
        set(gca,'box','on')
        if PlotL==1
            title(sprintf('CORR distributions for chip m%u Ps %u and nets %u to %u',ChipRank,PsRank,NetRanks(1),NetRanks(end)))
        elseif PlotL==2
            title(sprintf('ANTI distributions for chip m%u Ps %u and nets %u to %u',ChipRank,PsRank,NetRanks(1),NetRanks(end)))
        else
            title(sprintf('CORR-ANTI distributions for chip m%u Ps %u and nets %u to %u',ChipRank,PsRank,NetRanks(1),NetRanks(end)))
        end
    end
end

for TypeL=1:TypeNb
    h=figure;
    if TypeL==1
        set(h,'name','CORR and ANTI HEATHMAPS')
        CurrCorr=Corr;
        CurrAnti=Anti;
    else
        set(h,'name','RAW(CORR) and RAW(ANTI) HEATHMAPS')
        CurrCorr=double(Corr).*double(Norm)/100;
        CurrAnti=double(Anti).*double(Norm)/100;
    end
    set(gcf,'color',[1,1,1])
    subplot(3,1,1)
    h=pcolor(CurrCorr);
    set(h,'linestyle','none')
    if TypeL==1
        title(sprintf('CORR heathmap for chip m%u Ps %u and nets %u to %u',ChipRank,PsRank,NetRanks(1),NetRanks(end)))
    else
        title(sprintf('RAW(CORR) heathmap for chip m%u Ps %u and nets %u to %u',ChipRank,PsRank,NetRanks(1),NetRanks(end)))
    end
    subplot(3,1,2)
    h=pcolor(CurrAnti);
    set(h,'linestyle','none')
    if TypeL==1
        title(sprintf('ANTI heathmap for chip m%u Ps %u and nets %u to %u',ChipRank,PsRank,NetRanks(1),NetRanks(end)))
    else
        title(sprintf('RAW(ANTI) heathmap for chip m%u Ps %u and nets %u to %u',ChipRank,PsRank,NetRanks(1),NetRanks(end)))
    end
    subplot(3,1,3)
    h=pcolor(CurrCorr-CurrAnti);
    set(h,'linestyle','none')
    if TypeL==1
        title(sprintf('CORR-ANTI heathmap for chip m%u Ps %u and nets %u to %u',ChipRank,PsRank,NetRanks(1),NetRanks(end)))
    else
        title(sprintf('RAW(CORR)-RAW(ANTI) heathmap for chip m%u Ps %u and nets %u to %u',ChipRank,PsRank,NetRanks(1),NetRanks(end)))
    end
end

