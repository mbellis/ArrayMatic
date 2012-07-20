%GSEA_FINDGENESET
%find significative gene sets in a result (ordered list of fdr, zvar,...)

%INPUT PARAMETERS
% ModelRank: chip rank
% Values : values to be searched for GSEA
% CluFlag=1 : exist clusters in used gene set (no clusters in cliques
%          ; clusters in Kegg pathways)
% AbsentIndex: probe sets that must not be used in tested gene sets

%EXTERNAL FILES

%OUTPUT PARAMETERS


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

function [GsIndex,GsClasses,GsName]=gsea_findgeneset(ModelRank,ZVarValues,FdrValues,MeanRanks,CluFlag,AbsentIndex,GeneSetFile,GeneSetDir,BestFile,BestDir,GeneId,DisplayFlag,PrintFlag,PrintName)
%gsea_findgeneset(11,OperonZVar,OperonFdr,OperonMeanRanks,0,AbsentIndex,'m11_cliques','/home/mbellis/array2/mbellis/sosma/arraymatic/amcdata/cliques','m11_cliques_type1_r5489_c2000','/home/mbellis/array2/mbellis/sosma/arraymatic/amcdata/cliques',TairId)
%gsea_findgeneset(11,OperonZVar,OperonFdr,OperonMeanRanks,0,AbsentIndex,'m11_athgsdb11',K.dir.affyMetadata,'athgsdb11_type1_r5489_c2000',K.dir.affyMetadata,TairId)
global K P
RandSeed=regexp(BestFile,'r[0-9]+','match');
cd(GeneSetDir)
eval(sprintf('load %s', GeneSetFile))
cd(BestDir)
eval(sprintf('load %s',BestFile))
Best{1}=BinCount;
clear BinCount
PsNb=K.chip.probeSetNbs{ModelRank}(1);

%recover non empty gene sets
GsPos=[];
for i=1:length(Gs.name)
    if ~isempty(Gs.index{i})
        GsPos=[GsPos;i];
    end
end
GsNb=length(GsPos);



Bin=ceil(PsNb/100);
BinIndex=[[Bin:Bin:Bin*10],[Bin*10+Bin*4:Bin*4:Bin*10+Bin*24],[PsNb-Bin*10-Bin*24:Bin*4:PsNb-Bin*10-Bin*4],[PsNb-Bin*10:Bin:PsNb-Bin]];
BinNb=length(BinIndex);



ResNb=length(ZVarValues);

if CluFlag==0
    TypeNb=1;
    UseAll=1;
else
    TypeNb=4;
    UseAll=0;
end

%merge selected results
if ~isempty(ZVarValues)
    [IdZVarValues IdSortIndex]=sort(ZVarValues,'descend');
    [Temp IdRevIndex]=sort(IdSortIndex);
    [Temp ChSortIndex]=sort(abs(ZVarValues),'descend');
    [Temp ChRevIndex]=sort(ChSortIndex);
    [Temp ContSortIndex]=sort(MeanRanks(:,1),'descend');
    [Temp ContRevIndex]=sort(ContSortIndex);
    [Temp TestSortIndex]=sort(MeanRanks(:,2),'descend');
    [Temp TestRevIndex]=sort(TestSortIndex);
    IdFdrValues=FdrValues(IdSortIndex);

else
    NullFdr=find(FdrValues==0);
    if ~isempty(NullFdr)
        h=errordlg('Exist null fdr => can not know it they are increased or decreased. Restart and used ZVar values.');
        waitfor(h)
        error('process canceled')
    end
    PosFdr=find(FdrValues>0);
    [Temp SortIndex]=sort(FdrValues(PosFdr));
    IdSortIndex=PosFdr(SortIndex);
    NegFdr=find(FdrValues<0);
    [Temp SortIndex]=sort(FdrValues(NegFdr));
    IdSortIndex=[IdSortIndex;NegFdr(SortIndex)];
    IdFdrValues=FdrValues(IdSortIndex);
    [Temp ChSortIndex]=sort(abs(FdrValues));
    [Temp ChRevIndex]=sort(ChSortIndex);
end
PsRank=[1:PsNb]';
IdPsRank=PsRank(IdSortIndex);

Bindex=zeros(PsNb,1);
if DisplayFlag
    h=figure;
    set(gcf,'color',[1,1,1])
    plot(ZVarValues,FdrValues,'b.')



    Bindex=zeros(PsNb,1);
    h1=figure;
    set(gcf,'color',[1,1,1])
    for TypeL=1:4
        if TypeL==1
            Pos=find(ZVarValues>0&MeanRanks(:,2)>MeanRanks(:,1));
        elseif TypeL==2
            Pos=find(ZVarValues<0&MeanRanks(:,2)<MeanRanks(:,1));
        elseif TypeL==3
            Pos=find(ZVarValues>0&MeanRanks(:,2)<MeanRanks(:,1));
        else
            Pos=find(ZVarValues<0&MeanRanks(:,2)>MeanRanks(:,1));
        end
        PosNb=length(Pos);

        IdPos=Bindex';
        IdPos(Pos)=IdRevIndex(Pos);
        IdPos=IdPos(find(IdPos));

        ContPos=Bindex';
        ContPos(Pos)=ContRevIndex(Pos);
        ContPos=ContPos(find(ContPos));

        TestPos=Bindex';
        TestPos(Pos)=TestRevIndex(Pos);
        TestPos=TestPos(find(TestPos));

        figure(h1)
        subplot(4,1,TypeL)
        line([IdPos;TestPos;ContPos],[ones(1,PosNb)*3;ones(1,PosNb)*2;ones(1,PosNb)]);
        set(gca,'box','on')
        set(gca,'ytick',[1,2,3])
        set(gca,'yticklabel',{'WT','MU','ZVar'})
        set(gca,'xlim',[1,PsNb])
        title(sprintf('%u ps',PosNb))

    end
end




if DisplayFlag
    CurrFdrPos=find(IdFdrValues>=0.10);
    FdrPos=CurrFdrPos(1);
    CurrFdrPos=find(IdFdrValues<=-0.10);
    FdrPos=[FdrPos,CurrFdrPos(end)];
    CurrFdrPos=find(IdFdrValues>=0.30);
    FdrPos=[FdrPos,CurrFdrPos(1)];
    CurrFdrPos=find(IdFdrValues<=-0.30);
    FdrPos=[FdrPos,CurrFdrPos(end)];
end

IdEs=cell(TypeNb);
ChEs=IdEs;


for TypeL=1:TypeNb
    IdEs{TypeL}=zeros(GsNb,BinNb);
    ChEs{TypeL}=zeros(GsNb,BinNb);
    Bindex=zeros(PsNb,1);
    for GsL=1:GsNb
        GsRank=GsPos(GsL);
        switch TypeL
            case 1
                Index=Gs.index{GsRank};
            case 2
                %Clu1
                Index=Gs.index{GsRank}(Gs.clusters{GsRank}{1});
            case 3
                %Clu2
                Index=Gs.index{GsRank}(Gs.clusters{GsRank}{2});
            case 4
                %Reliquat
                Index=Gs.index{GsRank}(Gs.clusters{GsRank}{3});
        end
        CurrValues=Bindex;
        CurrValues(Index)=1;
        %don't use probe in AbsentIndex
        CurrValues(AbsentIndex)=0;
        if ~isempty(find(CurrValues))
            CurrValues=CurrValues(IdSortIndex);
            %FIRST ROUND : INCREASED/DECREASED
            IdEs{TypeL}(GsL,:)=BINNED_RUNNING_SUM(CurrValues,BinIndex);
            %SND ROUND : CHANGED/NOT CHANGED
            CurrValues=CurrValues(ChSortIndex);
            ChEs{TypeL}(GsL,:)=BINNED_RUNNING_SUM(CurrValues,BinIndex);
        end
    end
end


if DisplayFlag
    %find the limit of the regions
    if isfield(Gs,'region')
        CurrR=Gs.region(1);
        Mark=0.3;
        RLim=Mark;
        GsL=1;
        Continue=1;
        while CurrR<max(Gs.region)
            GsL=GsL+1;
            if CurrR~=Gs.region(GsL);
                CurrR=Gs.region(GsL);
                if Mark==0.3
                    Mark=0.6;
                else
                    Mark=0.3;
                end
            end
            if CurrR<max(Gs.region)
                RLim=[RLim,Mark];
            else
                for i=GsL:GsNb
                    RLim=[RLim,0.45];
                end
            end
        end
        RLim(end)=0;
        RLim(end-1)=1;
        RLim=repmat(RLim,10,1);
        SubY=4;
        RLimFlag=1;
        FirstPlot=2;
        SndPlot=4;
    else
        SubY=2;
        RLimFlag=0;
        FirstPlot=1;
        SndPlot=2;
    end

    GsName=cell(length(Gs.name),1);
    if isfield(Gs,'region')
        for GsL=1:length(Gs.name)
            GsName{GsL}=sprintf('%s_%u_%u',Gs.name{GsL}(1:9),length(Gs.index{GsL}),length(setdiff(Gs.index{GsL},AbsentIndex)));
        end
    else
        for GsL=1:length(Gs.name)
            GsName{GsL}=sprintf('%s_%u_%u',Gs.name{GsL},length(Gs.index{GsL}),length(setdiff(Gs.index{GsL},AbsentIndex)));
        end
    end


    h=figure;
    set(gcf,'color',[1,1,1])
    set(h,'name','heatmap of ES for all gene sets')
    if RLimFlag
        subplot(SubY,1,1)
        h=pcolor(RLim);
        set(h,'linestyle','none')
        set(gca,'position',[0.05,0.95,0.90,0.03])
        set(gca,'ytick',[])
        set(gca,'xtick',[])
    end
    subplot(SubY,1,FirstPlot)
    image(abs(IdEs{1}*255)');
    %h=pcolor(abs(IdEs{1})');
    set(gca,'position',[0.05,0.60,0.90,0.34])
    %set(h,'linestyle','none')
    title('INCREASED/DECREASED')
    set(gca,'clim',[-0.5,0.5])
    if RLimFlag
        %colorbar
        subplot(SubY,1,3)
        h=pcolor(RLim);
        set(h,'linestyle','none')
        set(gca,'position',[0.05,0.45,0.90,0.03])
        set(gca,'ytick',[])
        set(gca,'xtick',[])
    end
    subplot(SubY,1,SndPlot)
    image(abs(ChEs{1}*255)');
    %h=pcolor(abs(ChEs{1})');
    set(gca,'position',[0.05,0.10,0.90,0.34])
    %set(h,'linestyle','none')
    title('CHANGED/NOT CHANGED')
    set(gca,'clim',[-0.5,0.5])
end




%calculate Pv at all bin position
IdEsPv=cell(TypeNb);
ChEsPv=IdEsPv;



for TypeL=1:TypeNb
    PointNb=size(Best{TypeL}{1}{1},1);
    IdEsPv{TypeL}=ones(GsNb,BinNb);
    ChEsPv{TypeL}=ones(GsNb,BinNb);

    for GsL=1:GsNb
        if sum(IdEs{TypeL}(GsL,:)~=0)
            %INC / CHANGED FOR FIRST INTERVALS
            for BinL=1:BinNb/2
                %INC
                CurrEs=single(IdEs{TypeL}(GsL,BinL));
                EsDist=sort(Best{TypeL}{1}{GsL}(:,BinL));
                %if CurrEs>min(EsDist)
                Pos=find(round(EsDist*1000)>=round(CurrEs*1000));
                if isempty(Pos)
                    IdEsPv{TypeL}(GsL,BinL)=1/PointNb;
                else
                    IdEsPv{TypeL}(GsL,BinL)=(PointNb-Pos(1)+1)/PointNb;
                end
                %end

                %CHANGED
                CurrEs=single(ChEs{TypeL}(GsL,BinL));
                EsDist=sort(Best{TypeL}{2}{GsL}(:,BinL));
                %if CurrEs>min(EsDist)
                Pos=find(round(EsDist*1000)>=round(CurrEs*1000));
                if isempty(Pos)
                    ChEsPv{TypeL}(GsL,BinL)=1/PointNb;
                else
                    ChEsPv{TypeL}(GsL,BinL)=(PointNb-Pos(1)+1)/PointNb;
                end
                %end
            end
            %         %MIDDLE INTERVAL
            %         BinL=((BinNb-1)/2)+1;
            %         %LIMIT INC/DEC
            %         CurrEs=abs(IdEs{TypeL}(GsL,BinL));
            %         EsDist=sort(abs(Best{TypeL}{1}{GsL}(:,BinL)));
            %         if CurrEs>min(EsDist)
            %             Pos=find(EsDist>=CurrEs);
            %             if isempty(Pos)
            %                 IdEsPv{TypeL}(GsL,BinL)=1/PointNb;
            %             else
            %                 IdEsPv{TypeL}(GsL,BinL)=(PointNb-Pos(1)+1)/PointNb;
            %             end
            %         end
            %         %LIMIT CHANGED/NOT CHANGED
            %         CurrEs=abs(ChEs{TypeL}(GsL,BinL));
            %         EsDist=sort(abs(Best{TypeL}{2}{GsL}(:,BinL)));
            %         if CurrEs>min(EsDist)
            %             Pos=find(EsDist>=CurrEs);
            %             if isempty(Pos)
            %                 ChEsPv{TypeL}(GsL,BinL)=1/PointNb;
            %             else
            %                 ChEsPv{TypeL}(GsL,BinL)=(PointNb-Pos(1)+1)/PointNb;
            %             end
            %         end
            %DEC / NOT CHANGED
            for BinL=BinNb/2+1:BinNb
                CurrEs=single(IdEs{TypeL}(GsL,BinL));
                %DEC
                EsDist=sort(Best{TypeL}{1}{GsL}(:,BinL));
                %if CurrEs<max(EsDist)
                Pos=find(round(EsDist*1000)<=round(CurrEs*1000));
                if isempty(Pos)
                    IdEsPv{TypeL}(GsL,BinL)=1/PointNb;
                else
                    IdEsPv{TypeL}(GsL,BinL)=Pos(end)/PointNb;
                end
                %end
                %NOT CHANGED
                CurrEs=single(ChEs{TypeL}(GsL,BinL));
                EsDist=sort(Best{TypeL}{2}{GsL}(:,BinL));
                %if CurrEs<max(EsDist)
                Pos=find(round(EsDist*1000)<=round(CurrEs*1000));
                if isempty(Pos)
                    ChEsPv{TypeL}(GsL,BinL)=1/PointNb;
                else
                    ChEsPv{TypeL}(GsL,BinL)=Pos(end)/PointNb;
                end
                %end
            end
        end
    end
end

%clear Best
%heat map of p-values
if DisplayFlag
    %ALL P-VALUES
    %HISTOGRAMES

    h=figure;
    set(gcf,'color',[1,1,1])
    set(h,'name','histogramme of all and min ES p-values for all gene sets');
    subplot(2,2,1)
    hist(-log10(IdEsPv{1}(:)),100)
    title('-log10(all P-values) INCREASED/DECREASED')
    subplot(2,2,2)
    hist(-log10(ChEsPv{1}(:)),100)
    title('-log10(all P-values)  CHANGED/NOT CHANGED')
    subplot(2,2,3)
    plot(sort(-log10(min(IdEsPv{1},[],2))),'m.')
    title('-log10(min(P-values)) INCREASED/DECREASED')
    subplot(2,2,4)
    plot(sort(-log10(min(IdEsPv{1},[],2))),'g.')
    title('-log10(min(P-values))  CHANGED/NOT CHANGED')


    %HEATMAP

    h=figure;
    set(gcf,'color',[1,1,1])
    set(h,'name','heatmap of ES p-values for all gene sets');
    if RLimFlag
        subplot(SubY,1,1)
        h=pcolor(RLim);
        set(h,'linestyle','none')
        set(gca,'position',[0.05,0.95,0.90,0.03])
        set(gca,'ytick',[])
        set(gca,'xtick',[])
    end

    subplot(SubY,1,FirstPlot)
    h=pcolor(flipud(-log10(IdEsPv{1}')));
    set(gca,'position',[0.05,0.60,0.90,0.34])
    set(h,'linestyle','none')
    set(gca,'yticklabel',flipud(get(gca,'yticklabel')))
    title('INCREASED/DECREASED')
    %set(gca,'clim',[-0.5,0.5])
    %colorbar

    if RLimFlag
        subplot(SubY,1,3)
        h=pcolor(RLim);
        set(h,'linestyle','none')
        set(gca,'position',[0.05,0.45,0.90,0.03])
        set(gca,'ytick',[])
        set(gca,'xtick',[])
    end

    subplot(SubY,1,SndPlot)
    h=pcolor(flipud(-log10(ChEsPv{1}')));
    set(gca,'position',[0.05,0.10,0.90,0.34])
    set(h,'linestyle','none')
    set(gca,'yticklabel',flipud(get(gca,'yticklabel')))
    title('CHANGED/NOT CHANGED')


    %P-VALUES OF INCREASED
    %HISTOGRAMES

    h=figure;
    set(gcf,'color',[1,1,1])
    set(h,'name','histogramme of all and min ES p-values for INC for all gene sets');
    subplot(2,2,1)
    a=IdEsPv{1}(:,1:floor(BinNb/2));
    hist(-log10(a(:)),100)
    title('-log10(all p-values) of INC')
    subplot(2,2,2)
    a=ChEsPv{1}(:,1:floor(BinNb/2));
    hist(-log10(a(:)),100)
    clear a
    title('-log10(all p-values) of CHANGED')
    subplot(2,2,3)
    plot(sort(-log10(min(IdEsPv{1}(:,1:floor(BinNb/2)),[],2))),'m.')
    title('-log10(min(P-values) of INCD')
    subplot(2,2,4)
    plot(sort(-log10(min(IdEsPv{1}(:,1:floor(BinNb/2)),[],2))),'g.')
    title('-log10(min(P-values) of CHANGED')


    %P-VALUES OF DECREASED
    %HISTOGRAMES

    h=figure;
    set(gcf,'color',[1,1,1])
    set(h,'name','histogramme of all and min ES p-values for DEC for all gene sets');
    subplot(2,2,1)
    a=IdEsPv{1}(:,floor(BinNb/2)+1:end);
    hist(-log10(a(:)),100)
    title('-log10(all P-values) of DECREASED')
    subplot(2,2,2)
    a=ChEsPv{1}(:,floor(BinNb/2)+1:end);
    hist(-log10(a(:)),100)
    clear a
    title('-log10(all P-values) of NOT CHANGED')
    subplot(2,2,3)
    plot(sort(-log10(min(IdEsPv{1}(:,floor(BinNb/2)+1:end),[],2))),'m.')
    title('-log10(min(P-values) of DECREASED')
    subplot(2,2,4)
    plot(sort(-log10(min(ChEsPv{1}(:,floor(BinNb/2)+1:end),[],2))),'g.')
    title('-log10(min(P-values) of NOT CHANGED')



    for PvL=1:2
        if PvL==1
            PvLimit=3;
        else
            PvLimit=2;
        end
        for VarL=1:2
            h=figure;
            set(gcf,'color',[1,1,1])
            %find the most significatives for all
            if VarL==1
                LowPv=find(-log10(min(IdEsPv{1},[],2))>=PvLimit);
            else
                LowPv=find(-log10(min(ChEsPv{1},[],2))>=PvLimit);
            end
            if PvL==1
                set(h,'name','heatmap of ES for pv<=0.001');
            else
                set(h,'name','heatmap of ES for pv<=0.01');
            end
            if RLimFlag
                subplot(SubY,1,1)
                h=pcolor(RLim);
                set(h,'linestyle','none')
                set(gca,'position',[0.05,0.95,0.90,0.03])
                set(gca,'ytick',[])
                set(gca,'xtick',[])
            end
            %ES
            subplot(SubY,1,FirstPlot)
            if VarL==1
                h=pcolor(flipud(abs(IdEs{1}(LowPv,:))'));
                title('ES FOR SIGNIFICATIVE INC/DEC')
            else
                h=pcolor(flipud(abs(ChEs{1}(LowPv,:))'));
                title('ES FOR SIGNIFICATIVE CHANGED/NOT CHANGED')
            end
            set(gca,'position',[0.05,0.60,0.90,0.34])
            set(h,'linestyle','none')
            %set(gca,'xticklabel',num2str(LowPv(round(str2num(get(gca,'xticklabel'))))))
            set(gca,'xticklabel','')
            set(gca,'yticklabel',flipud(get(gca,'yticklabel')))
            set(gca,'clim',[-0.5,0.5])
            if RLimFlag
                subplot(SubY,1,3)
                h=pcolor(RLim);
                set(h,'linestyle','none')
                set(gca,'position',[0.05,0.45,0.90,0.03])
                set(gca,'ytick',[])
                set(gca,'xtick',[])
            end
            subplot(SubY,1,SndPlot)
            if VarL==1
                h=pcolor(flipud(-log10(IdEsPv{1}(LowPv,:)')));
                title('PV FOR SIGNIFICATIVE INC/DEC OF INC')
            else
                h=pcolor(flipud(-log10(ChEsPv{1}(LowPv,:))'));
                title('PV FOR SIGNIFICATIVE CHANGED/NOT CHANGED')
            end
            set(h,'linestyle','none')
            set(gca,'position',[0.05,0.10,0.90,0.34])
            %set(gca,'xticklabel',num2str(LowPv(str2num(get(gca,'xticklabel')))))
            set(gca,'xticklabel','')
            set(gca,'yticklabel',flipud(get(gca,'yticklabel')))

        end
    end
end
%definition of different classes of results
PvLimit=2;
%inc/dec
IdIndex=find(-log10(min(IdEsPv{1},[],2))>=PvLimit);
%inc & dec &incdec
IncIndex=[];
DecIndex=[];
IncDecIndex=[];
for GsL=1:length(IdIndex)
    IncPos=find(-log10(IdEsPv{1}(IdIndex(GsL),1:BinNb/2))>=PvLimit);
    DecPos=find(-log10(IdEsPv{1}(IdIndex(GsL),(BinNb/2)+1:BinNb))>=PvLimit);
    if ~isempty(IncPos)&~isempty(DecPos)
        IncDecIndex=[IncDecIndex,IdIndex(GsL)];
    elseif ~isempty(IncPos)
        IncIndex=[IncIndex,IdIndex(GsL)];
    else
        DecIndex=[DecIndex,IdIndex(GsL)];
    end
end
%Ch/NotCh
CncIndex=find(-log10(min(ChEsPv{1},[],2))>=PvLimit);
%inc & dec &incdec
ChIndex=[];
NotChIndex=[];
ChNotChIndex=[];
for GsL=1:length(CncIndex)
    ChPos=find(-log10(ChEsPv{1}(CncIndex(GsL),1:BinNb/2))>=PvLimit);
    NotChPos=find(-log10(ChEsPv{1}(CncIndex(GsL),(BinNb/2)+1:BinNb))>=PvLimit);
    if ~isempty(ChPos)&~isempty(NotChPos)
        ChNotChIndex=[ChNotChIndex,CncIndex(GsL)];
    elseif ~isempty(ChPos)
        ChIndex=[ChIndex,CncIndex(GsL)];
    elseif  ~isempty(NotChPos)
        NotChIndex=[NotChIndex,CncIndex(GsL)];
    end
end

GsIndex{1}=setdiff(IncIndex,union(NotChIndex,union(ChIndex,ChNotChIndex)));
GsIndex{2}=intersect(IncIndex,ChIndex);
GsIndex{3}=intersect(IncIndex,NotChIndex);
GsIndex{4}=intersect(IncIndex,ChNotChIndex);
GsIndex{5}=setdiff(DecIndex,union(NotChIndex,union(ChIndex,ChNotChIndex)));
GsIndex{6}=intersect(DecIndex,ChIndex);
GsIndex{7}=intersect(DecIndex,NotChIndex);
GsIndex{8}=intersect(DecIndex,ChNotChIndex);
GsIndex{9}=setdiff(IncDecIndex,union(NotChIndex,union(ChIndex,ChNotChIndex)));
GsIndex{10}=intersect(IncDecIndex,ChIndex);
GsIndex{11}=intersect(IncDecIndex,NotChIndex);
GsIndex{12}=intersect(IncDecIndex,ChNotChIndex);
GsIndex{13}=setdiff(ChIndex,union(IncIndex,union(DecIndex,IncDecIndex)));
GsIndex{14}=setdiff(ChNotChIndex,union(IncIndex,union(DecIndex,IncDecIndex)));
GsIndex{15}=setdiff(NotChIndex,union(IncIndex,union(DecIndex,IncDecIndex)));
GsIndex{16}=union(IncIndex,union(DecIndex,union(IncDecIndex,union(ChIndex,ChNotChIndex))));
GsClasses={'Inc','Inc & Ch','Inc & NotCh','Inc & Ch  & NotCh',...
    'Dec','Dec & Ch','Dec & NotCh','Dec & ChNotCh',...
    'IncDec','IncDec & Ch','IncDec & NotCh','IncDec & ChNotCh',...
    'Ch','ChNotCh','NotCh','AllCh'};

if DisplayFlag
    %relation between pv and es
    PvLimit1=2;
    PvLimit2=3;
    PlotSize=20;
    PlotRawNb=4;
    PlotColNb=5;
    for IndexL=1:16
        CurrGsIndex=GsIndex{IndexL};
        if ~isempty(CurrGsIndex)
            if IndexL<16
                if IndexL<13 | IndexL==16
                    Es=IdEs{1};
                    EsPv=IdEsPv{1};
                    RevIndex=IdRevIndex;
                else
                    Es=ChEs{1};
                    EsPv=ChEsPv{1};
                    RevIndex=ChRevIndex;
                end
                PlotNb=length(CurrGsIndex);
                BlocNb=ceil(PlotNb/PlotSize);
                for BlocL=1:BlocNb
                    if BlocL==BlocNb
                        CurrPlotNb=mod(PlotNb,PlotSize);
                        if CurrPlotNb==0
                            CurrPlotNb=PlotSize;
                            RawNb=PlotRawNb;
                            ColNb=PlotColNb;
                        else
                            ColNb=ceil(sqrt(CurrPlotNb));
                            RawNb=round(CurrPlotNb/ColNb);
                            if ColNb*RawNb<CurrPlotNb
                                ColNb=ColNb+1;
                            end
                        end
                    else
                        CurrPlotNb=PlotSize;
                        RawNb=PlotRawNb;
                        ColNb=PlotColNb;
                    end
                    h=figure;
                    set(gcf,'color',[1,1,1])
                    set(h,'name',sprintf('pv vs ES for %s -  bloc %u',GsClasses{IndexL},BlocL));
                    %for j=1:length(CurrGsIndex)
                    for PlotL=1+(BlocL-1)*PlotSize:CurrPlotNb+(BlocL-1)*PlotSize
                        VarPos=Bindex';
                        VarPos(Gs.index{CurrGsIndex(PlotL)})=RevIndex(Gs.index{CurrGsIndex(PlotL)});
                        %don't use probe in AbsentIndex
                        VarPos(AbsentIndex)=0;
                        VarPos=VarPos(find(VarPos));

                        IdPos=Bindex';
                        IdPos(Gs.index{CurrGsIndex(PlotL)})=RevIndex(Gs.index{CurrGsIndex(PlotL)});
                        %don't use probe in AbsentIndex
                        IdPos(AbsentIndex)=0;
                        IdPos=IdPos(find(IdPos));



                        ContPos=Bindex';
                        ContPos(Gs.index{CurrGsIndex(PlotL)})=ContRevIndex(Gs.index{CurrGsIndex(PlotL)});
                        %don't use probe in AbsentIndex
                        ContPos(AbsentIndex)=0;
                        ContPos=ContPos(find(ContPos));

                        TestPos=Bindex';
                        TestPos(Gs.index{CurrGsIndex(PlotL)})=TestRevIndex(Gs.index{CurrGsIndex(PlotL)});
                        %don't use probe in AbsentIndex
                        TestPos(AbsentIndex)=0;
                        TestPos=TestPos(find(TestPos));
                        GsSize=length(TestPos);


                        subplot(RawNb,ColNb,PlotL-(BlocL-1)*PlotSize)
                        Factor=3/max(abs(Es(CurrGsIndex(PlotL),:)));
                        plot(BinIndex,Es(CurrGsIndex(PlotL),:)*Factor,'b')
                        hold on
                        plot(BinIndex,Es(CurrGsIndex(PlotL),:)*Factor,'b.')
                        [t a]=max(Es(CurrGsIndex(PlotL),1:16));
                        MaxEs=Es(CurrGsIndex(PlotL),a)*Factor;
                        plot(BinIndex(a),MaxEs,'ko')
                        [t a]=min(Es(CurrGsIndex(PlotL),17:end));
                        MinEs=Es(CurrGsIndex(PlotL),a+16)*Factor;
                        plot(BinIndex(a+16),MinEs,'co')
                        MaxPv=max(-log10(EsPv(CurrGsIndex(PlotL),:)));
                        plot(BinIndex,-log10(EsPv(CurrGsIndex(PlotL),:)),'g')
                        plot(BinIndex,-log10(EsPv(CurrGsIndex(PlotL),:)),'g.')
                        a=find(-log10(EsPv(CurrGsIndex(PlotL),:))>=PvLimit1);
                        if ~isempty(a)
                            plot(BinIndex(a),-log10(EsPv(CurrGsIndex(PlotL),a)),'ro')
                        end
                        a=find(-log10(EsPv(CurrGsIndex(PlotL),:))>=PvLimit2);
                        if ~isempty(a)
                            plot(BinIndex(a),-log10(EsPv(CurrGsIndex(PlotL),a)),'r+')
                        end
                        title(strrep(GsName{CurrGsIndex(PlotL)},'_',' '))
                        line([PsNb/2,PsNb/2],[MinEs-1,MaxEs],'color','k','linestyle',':')
                        line([0,PsNb],[0,0],'color','k','linestyle',':')
                        line([0,PsNb/2],[MinEs-1,MinEs-1],'color','r','linestyle','-')
                        line([PsNb/2,PsNb],[MinEs-1,MinEs-1],'color','b','linestyle','-')
                        line([0,PsNb],[MinEs-2,MinEs-2],'color','m','linestyle','-')
                        line([0,PsNb],[MinEs-3,MinEs-3],'color','c','linestyle','-')
                        XPos=[IdPos;TestPos;ContPos];
                        YPos=[ones(1,GsSize)*(MinEs-1);ones(1,GsSize)*(MinEs-2);ones(1,GsSize)*(MinEs-3)];
                        line(XPos,YPos,'color','k')
                        Pos=find(IdPos<=FdrPos(3));
                        if ~isempty(Pos)
                            line(XPos(:,Pos),YPos(:,Pos),'color','m')
                        end
                        Pos=find(IdPos<=FdrPos(1));
                        if ~isempty(Pos)
                            line(XPos(:,Pos),YPos(:,Pos),'color','r')
                        end
                        Pos=find(IdPos>=FdrPos(2));
                        if ~isempty(Pos)
                            line(XPos(:,Pos),YPos(:,Pos),'color','c')
                        end
                        Pos=find(IdPos>=FdrPos(4));
                        if ~isempty(Pos)
                            line(XPos(:,Pos),YPos(:,Pos),'color','b')
                        end


                        set(gca,'xlim',[0,PsNb])
                        set(gca,'ylim',[MinEs-3.5,max(MaxPv,MaxEs)+0.5])
                        set(gca,'xtick',[])
                        set(gca,'ytick',[])
                        set(gca,'xticklabel','')
                        set(gca,'yticklabel','')
                    end
                end
            end



            Filled=zeros(PsNb,1);
            for TypeL=1:TypeNb
                Es=IdEs{TypeL};
                EsPv=IdEsPv{TypeL};
                SortIndex=IdSortIndex;
                LowPv1=CurrGsIndex;
                %length of separation between contiguous gene set
                LowSize=0;
                IncPerc=zeros(length(LowPv1),1);
                MeanPos=zeros(length(LowPv1),1);
                %status = 1 => Inc/Ch
                %status = 0 => Dec/NotCh
                Status=zeros(length(LowPv1),1);
                for GsL=1:length(LowPv1)
                    CurrGs=LowPv1(GsL);
                    CurrIndex=setdiff(Gs.index{CurrGs},AbsentIndex);
                    if ~isempty(CurrIndex)
                        LowSize=LowSize+length(CurrIndex);
                        ValPos=[];
                        for PsL=1:length(CurrIndex)
                            ValPos=[ValPos,find(IdSortIndex==CurrIndex(PsL))];
                        end
                        IncPerc(GsL)=round(length(find(ValPos<=PsNb/2))*100/length(ValPos));
                        PvPos=find(-log10(EsPv(CurrGs,:))>PvLimit1);
                        if length(find(PvPos<=16))>=length(find(PvPos>=17))
                            Status(GsL)=1;
                        end
                        %if IncPerc(GsL)>=50
                        if Status(GsL)
                            MainPos=find(ValPos<=PsNb/2);
                        else
                            MainPos=find(ValPos>PsNb/2);
                        end
                        MeanPos(GsL)=mean(ValPos(MainPos));
                    end
                end
                SepSize=floor((PsNb-LowSize)/GsL);
                LowPv2=find(-log10(min(EsPv,[],2))>=PvLimit2);
                MemLowPv1=LowPv1;
                LowPv1=setdiff(LowPv1,LowPv2);
                h=figure;
                hold on
                set(h,'color',[1,1,1])
                %separate Inc and Dec in LowPv2
                IncList=[];
                DecList=[];
                IncListPos=[];
                DecListPos=[];
                for GsL=1:length(LowPv2)
                    CurrGs=LowPv2(GsL);
                    CurrPos=find(MemLowPv1==CurrGs);
                    if Status(CurrPos)
                        %if IncPerc(CurrPos)>=50
                        IncList=[IncList,CurrGs];
                        IncListPos=[IncListPos,CurrPos];
                    else
                        DecList=[DecList,CurrGs];
                        DecListPos=[DecListPos,CurrPos];
                    end
                end
                MemDecList=DecList;
                MemDecListPos=DecListPos;

                GsPos=0;

                %process Inc in LowPv2
                CurrPerc=IncPerc(IncListPos);
                CurrMeanPos=MeanPos(IncListPos);
                [Temp,SortIndex]=sort(CurrMeanPos,'descend');
                CurrPerc=CurrPerc(SortIndex);
                CurrList=IncList(SortIndex);
                [Temp,SortIndex]=sort(CurrPerc,'descend');
                CurrList=CurrList(SortIndex);


                for GsL=1:length(CurrList)
                    CurrGs=CurrList(GsL);
                    CurrIndex=setdiff(Gs.index{CurrGs},AbsentIndex);
                    for PsL=1:length(CurrIndex)
                        GsPos=GsPos+1;
                        ValPos=find(IdSortIndex==CurrIndex(PsL));
                        Filled(ValPos)=1;
                        if ValPos<=PsNb/2
                            color='r';
                        else
                            color='m';
                        end
                        line([GsPos,ValPos],[1,0],'color',color)
                    end
                    GsPos=GsPos+SepSize;
                    %             if color=='r'
                    %                 color='m';
                    %             else
                    %                 color='r';
                    %             end
                end


                %separate Inc and Dec in LowPv1
                IncList=[];
                DecList=[];
                IncListPos=[];
                DecListPos=[];
                for GsL=1:length(LowPv1)
                    CurrGs=LowPv1(GsL);
                    CurrPos=find(MemLowPv1==CurrGs);
                    %if IncPerc(CurrPos)>=50
                    if Status(CurrPos)
                        IncList=[IncList,CurrGs];
                        IncListPos=[IncListPos,CurrPos];
                    else
                        DecList=[DecList,CurrGs];
                        DecListPos=[DecListPos,CurrPos];
                    end
                end

                %process Inc in LowPv1
                CurrPerc=IncPerc(IncListPos);
                CurrMeanPos=MeanPos(IncListPos);
                [Temp,SortIndex]=sort(CurrMeanPos,'descend');
                CurrPerc=CurrPerc(SortIndex);
                CurrList=IncList(SortIndex);
                [Temp,SortIndex]=sort(CurrPerc,'descend');
                CurrList=CurrList(SortIndex);

                for GsL=1:length(CurrList)
                    CurrGs=CurrList(GsL);
                    CurrIndex=setdiff(Gs.index{CurrGs},AbsentIndex);
                    for PsL=1:length(CurrIndex)
                        GsPos=GsPos+1;
                        ValPos=find(IdSortIndex==CurrIndex(PsL));
                        Filled(ValPos)=1;
                        if ValPos<=PsNb/2
                            color='g';
                        else
                            color='y';
                        end
                        line([GsPos,ValPos],[1,0],'color',color)
                    end
                    GsPos=GsPos+SepSize;
                    %             if color=='g'
                    %                 color='y';
                    %             else
                    %                 color='g';
                    %             end
                end


                %process Dec in LowPv1
                CurrPerc=IncPerc(DecListPos);
                CurrMeanPos=MeanPos(DecListPos);
                [Temp,SortIndex]=sort(CurrMeanPos,'descend');
                CurrPerc=CurrPerc(SortIndex);
                CurrList=DecList(SortIndex);
                [Temp,SortIndex]=sort(CurrPerc,'descend');
                CurrList=CurrList(SortIndex);

                for GsL=1:length(CurrList)
                    CurrGs=CurrList(GsL);
                    CurrIndex=setdiff(Gs.index{CurrGs},AbsentIndex);
                    for PsL=1:length(CurrIndex)
                        GsPos=GsPos+1;
                        ValPos=find(IdSortIndex==CurrIndex(PsL));
                        Filled(ValPos)=1;
                        if ValPos<=PsNb/2
                            color='c';
                        else
                            color='b';
                        end

                        line([GsPos,ValPos],[1,0],'color',color)
                    end
                    GsPos=GsPos+SepSize;
                    %             if color=='b'
                    %                 color='c';
                    %             else
                    %                 color='b';
                    %             end
                end




                %process dec in LowPv2
                DecList=MemDecList;
                DecListPos=MemDecListPos;
                CurrPerc=IncPerc(DecListPos);
                CurrMeanPos=MeanPos(DecListPos);
                [Temp,SortIndex]=sort(CurrMeanPos,'descend');
                CurrPerc=CurrPerc(SortIndex);
                CurrList=DecList(SortIndex);
                [Temp,SortIndex]=sort(CurrPerc,'descend');
                CurrList=CurrList(SortIndex);

                for GsL=1:length(CurrList)
                    CurrGs=CurrList(GsL);
                    CurrIndex=setdiff(Gs.index{CurrGs},AbsentIndex);
                    for PsL=1:length(CurrIndex)
                        GsPos=GsPos+1;
                        ValPos=find(IdSortIndex==CurrIndex(PsL));
                        Filled(ValPos)=1;
                        if ValPos<=PsNb/2
                            color='r';
                        else
                            color='m';
                        end

                        line([GsPos,ValPos],[1,0],'color',color)
                    end
                    GsPos=GsPos+SepSize;
                    %             if color=='r'
                    %                 color='m';
                    %             else
                    %                 color='r';
                    %             end
                end

                Inc=cumsum(Filled)./[1:PsNb]';
                Dec=flipud(cumsum(flipud(Filled)))./[PsNb:-1:1]';
                plot(Inc,'k-','linewidth',3)
                plot(Dec,'b-','linewidth',3)
                Inc=cumsum(Filled)/length(find(Filled));
                Dec=flipud(cumsum(flipud(Filled)))./length(find(Filled));
                plot(Inc,'k-','linewidth',3)
                plot(Dec,'b-','linewidth',3)
                set(gca,'xlim',[1,PsNb])
                set(gca,'box','on')
                title(sprintf('%s CLIQUES (%u  GENES)',GsClasses{IndexL},LowSize))
                xlabel('GENES ORDERED ON FDR (INC => DEC)')
            end
        end
    end
end


if PrintFlag
    %calculate the number of significative gene sets
    CurrGsNb=0;
    for IndexL=1:15
        CurrGsIndexes=GsIndex{IndexL};
        CurrGsNb=CurrGsNb+length(CurrGsIndexes);
    end
    %fill Result and construct headers
    Result=uint8(zeros(PsNb,CurrGsNb));
    ResPos=0;
    Header1=sprintf('PsRank\tPsId\tTairId\tZVar\tFdr\tIsGs\tPrcInGs');
    Header2=sprintf('\t\t\t\t\t\t');
    for IndexL=1:15
        CurrGsIndexes=GsIndex{IndexL};
        if ~isempty(CurrGsIndexes)
            Header1=sprintf('%s\t%s',Header1,GsClasses{IndexL});
            for GsL=1:length(CurrGsIndexes)
                if GsL>1
                    Header1=sprintf('%s\t',Header1);
                end
                Header2=sprintf('%s\t%s',Header2,Gs.name{CurrGsIndexes(GsL)});
                ResPos=ResPos+1;
                Result(Gs.index{CurrGsIndexes(GsL)},ResPos)=1;
            end
        end
    end
    Result=Result(IdSortIndex,:);
    ExistGs=sum(Result,2);
    ExistGs(find(ExistGs))=1;
    IncGs=cumsum(ExistGs)*100./[1:PsNb]';
    DecGs=flipud(cumsum(flipud(ExistGs)))*100./[PsNb:-1:1]';
    PercGs=max([IncGs,DecGs],[],2);
    fid=fopen(sprintf('gsea_%s_%s_%s.txt',PrintName,RandSeed{1},date),'w');
    fprintf(fid,'%s\n',Header1);
    fprintf(fid,'%s\n',Header2);
    GeneId=GeneId(IdSortIndex);
    PsId=P.chip.probeSetIds(IdSortIndex);
    for PsL=1:PsNb
        fprintf(fid,'%u\t%s\t%s\t%.0f\t%.2f\t%u\t%u',IdPsRank(PsL),PsId{PsL},GeneId{PsL},IdZVarValues(PsL),IdFdrValues(PsL),ExistGs(PsL),round(PercGs(PsL)));
        for GsL=1:CurrGsNb
            fprintf(fid,'\t%u',Result(PsL,GsL));
        end
        fprintf(fid,'\n');
    end
    fclose(fid)
end

GsName=Gs.name;

% %PRINT THE RESULTS
% OutName=inputdlg({'name of output fil'},'',1,{''});
% OutName=OutName{1};
% cd(K.dir.ures)
% for TypeL=1:TypeNb
%     IncCode=-log10(min(IdEsPv{1}{TypeL}(:,1:(BinNb-1)/2),[],2))>=2;
%     if ResNb>1
%         for ResL=2:ResNb
%             IncCode=[IncCode+(2^(ResL-1))*(-log10(min(IdEsPv{TypeL}(:,1:(BinNb-1)/2),[],2))>=2)];
%         end
%     end
%     IncPos=find(IncCode);
%     %IncCode=IncCode(IncPos);
%     DecCode=-log10(min(IdEsPv{1}{TypeL}(:,((BinNb-1)/2)+2:BinNb),[],2))>=2;
%     if ResNb>1
%         for ResL=2:ResNb
%             DecCode=[DecCode+(2^(ResL-1))*(-log10(min(IdEsPv{TypeL}(:,((BinNb-1)/2)+2:BinNb),[],2))>=2)];
%         end
%     end
%     DecPos=find(DecCode);
%     %DecCode=DecCode(IncPos);
%
%     Pos=unique([IncPos;DecPos]);
%
%
%     for BranchL=1:3
%         switch BranchL
%             case 1
%                 fid=fopen(sprintf('bp_%s_t%u.csv',OutName,TypeL),'w');
%             case 2
%                 fid=fopen(sprintf('mf_%s_t%u.csv',OutName,TypeL),'w')
%             case 3
%                 fid=fopen(sprintf('cc_%s_t%u.csv',OutName,TypeL),'w');
%         end
%         %print header INC
%         fprintf(fid,'INC\n')
%         %print specificity code
%         fprintf(fid,'\tSPECIFICITY');
%         for PosL=1:length(Pos)
%             fprintf(fid,'\t%u',IncCode(Pos(PosL)));
%         end
%         fprintf(fid,'\n');
%         fprintf(fid,'ES POSITION / P-VALUE\n');
%
%         for ResL=1:ResNb
%             %print Bin position
%             fprintf(fid,'%s\tPosition',CurrRes.name);
%             for PosL=1:length(Pos)
%                 if ~isempty(find(IncPos==Pos(PosL)))
%                     [MinPv BinIndex]=min(IdEsPv{TypeL}(Pos(PosL),1:(BinNb-1)/2));
%                     if MinPv<=0.01
%                         fprintf(fid,'\t%u',PsPos(BinIndex));
%                     else
%                         fprintf(fid,'\t%u',0);
%                     end
%                 else
%                     fprintf(fid,'\t0')
%                 end
%             end
%             fprintf(fid,'\n');
%             %print Bin p-value
%             fprintf(fid,'\tp-value');
%             for PosL=1:length(Pos)
%                 if ~isempty(find(IncPos==Pos(PosL)))
%                     [MinPv Temp]=min(IdEsPv{TypeL}(Pos(PosL),1:(BinNb-1)/2));
%                     if MinPv<=0.01
%                         fprintf(fid,'\t%.2f',-log10(MinPv));
%                     else
%                         fprintf(fid,'\t%u',0);
%                     end
%
%
%                 else
%                     fprintf(fid,'\t0')
%                 end
%             end
%             fprintf(fid,'\n');
%         end
%
%         %print header DEC
%         fprintf(fid,'DEC\n')
%         %print specificity code
%         fprintf(fid,'\tSPECIFICITY');
%         for PosL=1:length(Pos)
%             fprintf(fid,'\t%u',DecCode(Pos(PosL)));
%         end
%         fprintf(fid,'\n');
%         fprintf(fid,'ES POSITION / P-VALUE\n');
%         %print Bin position
%         for ResL=1:ResNb
%             %print Bin position
%             fprintf(fid,'%s\tPosition',CurrRes.name);
%             for PosL=1:length(Pos)
%                 if ~isempty(find(DecPos==Pos(PosL)))
%                     [MinPv BinIndex]=min(IdEsPv{TypeL}(Pos(PosL),((BinNb-1)/2)+1:end));
%                     if MinPv<=0.01
%                         fprintf(fid,'\t%u',PsNb-PsPos(BinIndex+(BinNb-1)/2));
%                     else
%                         fprintf(fid,'\t0');
%                     end
%
%                 else
%                     fprintf(fid,'\t0')
%                 end
%             end
%             fprintf(fid,'\n');
%             %print Bin p-value
%             fprintf(fid,'\tp-value');
%             for PosL=1:length(Pos)
%                 if ~isempty(find(DecPos==Pos(PosL)))
%                     [MinPv Temp]=min(IdEsPv{TypeL}(Pos(PosL),((BinNb-1)/2)+1:end));
%                     if MinPv<=0.01
%                         fprintf(fid,'\t%.2f',-log10(MinPv));
%                     else
%                         fprintf(fid,'\t0')
%                     end
%                 else
%                     fprintf(fid,'\t0')
%                 end
%             end
%             fprintf(fid,'\n');
%
%         end
%
%         fprintf(fid,'\n');
%         %print split information
%         fprintf(fid,'SPLITTED MODULES\n')
%         for ResL=1:ResNb
%             fprintf(fid,'%s\t',CurrRes.name);
%             for PosL=1:length(Pos)
%                 if ~isempty(find(IncPos==Pos(PosL)))|~isempty(find(DecPos==Pos(PosL)))
%                     if min(IdEsPv{TypeL}(Pos(PosL),1:(BinNb-1)/2))<=0.01&min(IdEsPv{TypeL}(Pos(PosL),((BinNb-1)/2)+1:end))<=0.01
%                         fprintf(fid,'\t3');
%                     elseif min(IdEsPv{TypeL}(Pos(PosL),1:(BinNb-1)/2))<=0.01
%                         fprintf(fid,'\t1');
%                     elseif min(IdEsPv{TypeL}(Pos(PosL),((BinNb-1)/2)+1:end))<=0.01
%                         fprintf(fid,'\t2');
%                     else
%                         fprintf(fid,'\t0');
%                     end
%                 end
%             end
%             fprintf(fid,'\n');
%         end
%
%
%         fprintf(fid,'\n');
%         %print region,module,size
%         fprintf(fid,'REGION\t');
%         for PosL=1:size(Pos)
%             fprintf(fid,'\t%u',Gs.region(Pos(PosL)));
%         end
%         fprintf(fid,'\n');
%         %print region,module,size
%         fprintf(fid,'MODULE\t');
%         for PosL=1:size(Pos)
%             fprintf(fid,'\t%u',Gs.module(Pos(PosL)));
%         end
%         fprintf(fid,'\n');
%         %print region,module,size
%         fprintf(fid,'SIZE\t');
%         for PosL=1:size(Pos)
%             fprintf(fid,'\t%u',length(Gs.index{Pos(PosL)}));
%         end
%         fprintf(fid,'\n');
%
%
%         %print Exist GO information
%         fprintf(fid,'\n')
%         fprintf(fid,'EXIST GO INFO\t')
%         for PosL=1:length(Pos)
%             CurrRegion=Gs.region(Pos(PosL));
%             CurrModule=Gs.module(Pos(PosL));
%             TablePos=find(Tables{BranchL}(1,:)==CurrRegion&Tables{BranchL}(2,:)==CurrModule);
%             if ~isempty(TablePos)
%                 if ~isempty(find(Tables{BranchL}(2:end,TablePos)<=-2))
%                     %if sum(Tables{BranchL}(2:end,TablePos))>0
%                     fprintf(fid,'\t1')
%                 else
%                     fprintf(fid,'\t0')
%                 end
%             else
%                 fprintf(fid,'\t0')
%             end
%         end
%         fprintf(fid,'\n')
%
%         %print GO information
%         TablePos=[];
%         for PosL=1:length(Pos)
%             CurrRegion=Gs.region(Pos(PosL));
%             CurrModule=Gs.module(Pos(PosL));
%             CurrTablePos=find(Tables{BranchL}(1,:)==CurrRegion&Tables{BranchL}(2,:)==CurrModule);
%             if ~isempty(CurrTablePos)
%                 TablePos=[TablePos,CurrTablePos];
%             end
%         end
%         fprintf(fid,'\n');
%         fprintf(fid,'GO ID\tGO TERM\n')
%         for GoL=1:length(GoIDs{BranchL})
%             if ~isempty(find(Tables{BranchL}(GoL+2,TablePos)<=-2))
%                 %if sum(Tables{BranchL}(GoL+2,TablePos))>0
%                 fprintf(fid,'%u',GoIDs{BranchL}(GoL));
%                 fprintf(fid,'\t%s',GoTerms{BranchL}{GoL});
%                 for PosL=1:length(Pos)
%                     CurrRegion=Gs.region(Pos(PosL));
%                     CurrModule=Gs.module(Pos(PosL));
%                     CurrTablePos=find(Tables{BranchL}(1,:)==CurrRegion&Tables{BranchL}(2,:)==CurrModule);
%                     if ~isempty(CurrTablePos)
%                         fprintf(fid,'\t%u',Tables{BranchL}(GoL+2,CurrTablePos));
%                     else
%                         fprintf(fid,'\t0')
%                     end
%                 end
%                 fprintf(fid,'\n');
%             end
%         end
%         fclose(fid)
%     end
% end

%% BINNED RUNNING SUM
function [ES]=BINNED_RUNNING_SUM(Val,BinIndex)

%one running sum from start to end
% BinNb=length(BinIndex);
% ES=zeros(1,BinNb);
HitPos=find(Val);
HitNb=length(HitPos);
MissNb=length(Val)-HitNb;
MissVal=-1/MissNb;
HitVal=1/HitNb;
SumVal=ones(size(Val))*MissVal;
SumVal(HitPos)=HitVal;
ESSum=cumsum(SumVal);
ES=ESSum(BinIndex);

% %% BINNED RUNNING SUM
% function [ES]=BINNED_RUNNING_SUM(Val,BinIndex)
%
% DISPLAY_IT=0;
% %one running sum from start to end
% ES=zeros(1,length(BinIndex));
% BinNb=(length(BinIndex)-1)/2;
% BinIndex=[1,BinIndex];
% SumVal=zeros(size(Val));
% HitPos=find(Val);
% HitNb=length(HitPos);
% MissPos=find(Val==0);
% MissNb=length(MissPos);
% SumVal(MissPos)=-1/MissNb;
% SumVal(HitPos)=1/HitNb;
% ESSum=cumsum(SumVal);
% for i=1:BinNb
%     ES(i)=max(ESSum(BinIndex(i):BinIndex(i+1)));
% end
% [t,Pos]=max(abs(ESSum(BinIndex(BinNb+1):BinIndex(BinNb+1+1))));
% ES(BinNb+1)=ESSum(BinIndex(i)+Pos);
% for i=BinNb+2:2*BinNb+1
%     ES(i)=min(ESSum(BinIndex(i):BinIndex(i+1)));
% end
%
% if DISPLAY_IT
%     BinIndex=[];
%     for i=1:2*BinNb+1
%         BinIndex(i)=BinIndex(i+1)-(BinIndex(i+1)-BinIndex(i))/2;
%     end
%     for i=BinNb+2:2*BinNb+1
%         BinIndex(i)=BinIndex(i+1)-(BinIndex(i+1)-BinIndex(i))/2;
%     end
%     figure
%     hold on
%     plot(ESSum,'k')
%     plot(BinIndex,ES,'ro')
% end

%% ES RUNNING SUM
function [Es]=ES_RUNNING_SUM(Values)
%one running sum from start to end
Es=zeros(size(Values));
HitPos=find(Values);
HitNb=length(HitPos);
MissPos=find(Values==0);
MissNb=length(MissPos);
if MissNb==0
    Es(:)=-1/HitNb;
elseif HitNb==0
    Es(:)=-1/MissNb;
else
    Es(MissPos)=-1/MissNb;
    Es(HitPos)=1/HitNb;
end
Es=cumsum(Es);
