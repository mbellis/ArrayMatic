function pca(Data)
global P
CurrColor=colormap;
ColorNb=size(CurrColor,1);

MarkerList='.ox+*sdv^<>ph';
MarkerNb=length(MarkerList);
Markers='';
for BiolL=1:P.biol.nb
    if mod(BiolL,MarkerNb)
        Markers(BiolL)=MarkerList(mod(BiolL,MarkerNb));
    else
        Markers(BiolL)=MarkerList(end);
    end
end

if P.flag.testAlgo
    LoopNb=length(unique(P.point.algo))+1;
    AlgoNames=unique(P.point.algo);
else
    LoopNb=1;
end
for LoopL=1:LoopNb
    if LoopL==1
        CurrPointNb=P.point.nb;
        CurrPointRanks=1:P.point.nb;
    else
        %P.flag.testAlgo==1 and the same set of points is analyzed
        %by
        %different algorithms
        CurrPointNb=length(strmatch(AlgoNames{1},P.point.algo,'exact'));
        CurrPointRanks=strmatch(AlgoNames{LoopL-1},P.point.algo,'exact');
    end
    CurrDataSignals=Data(:,CurrPointRanks);
    icolor=zeros(CurrPointNb,3);
    imarker=repmat('.',CurrPointNb,1);
    legendcellstr2={cell(P.biol.nb,1)};
    mlegendcolor=zeros(P.biol.nb,3);
    ColorStep=floor(ColorNb/P.biol.nb);
    for BiolL=1:P.biol.nb
        if mod((BiolL-1)*ColorStep+1,ColorNb)
            mlegendcolor(BiolL,:)=CurrColor(mod((BiolL-1)*ColorStep+1,ColorNb),:);
        else
            mlegendcolor(BiolL,:)=CurrColor(end,:);
        end
    end
    for BiolL=1:P.biol.nb
        legendcellstr2{1}{BiolL}=sprintf('%c %s',Markers(BiolL),P.biol.name{BiolL});
        BiolPos=find(P.point.biolRank(CurrPointRanks)==BiolL);
        try
            icolor(BiolPos,:)=repmat(mlegendcolor(BiolL,:),length(BiolPos),1);
        catch
            'stop'
        end
        imarker(BiolPos)=Markers(BiolL);
    end
    paramstruct = struct('viout',3, ...
        'vipcplot',1:2, ...
        'icolor',icolor, ...
        'imarker',imarker,...
        'mlegendcolor',mlegendcolor, ...
        'legendcellstr1',legendcellstr2, ...
        'legendcellstr2',legendcellstr2, ...
        'iscreenwrite',0);

    % Verify if exist negative values
    NegPos=find(CurrDataSignals<=0);
    if ~isempty(NegPos)
        CurrDataSignals(NegPos)=eps;
    end
    CurrDataSignals=log2(CurrDataSignals);
    h=figure;
    if LoopL==1
        set(h,'name','all points - biological conditions');
    else
        set(h,'name',sprintf('points analyzed by %s - biological conditions',AlgoNames{LoopL-1}));
    end
    %curvdatSM(DataLog2Signal,paramstruct)
    %replace NaN by 0 because some operation in curvdatSM are not
    %possible with Nan values
    NanPos=find(isnan(CurrDataSignals));
    if ~isempty(NanPos)
        CurrDataSignals(NanPos)=0;
    end
    curvdatSM_mb(CurrDataSignals,paramstruct)
    cd(P.dir.resTree)
    set_figsize('1024px')

    if LoopL==1
        saveas(h,sprintf('all_biolcond_%s_%s',P.project.name,date),'png');
    else
        saveas(h,sprintf('%s_all_biolcond_%s',AlgoNames{LoopL-1},date),'png');        
    end
    close(h)
    if LoopL==1&length(unique(P.point.algo))>1&P.flag.testAlgo
        AlgoNames=unique(P.point.algo);
        AlgoNb=length(AlgoNames);
        mlegendcolor=zeros(AlgoNb,3);
        ColorStep=floor(ColorNb/AlgoNb);
        for AlgoL=1:AlgoNb
            if mod((AlgoL-1)*ColorStep+1,ColorNb)
                mlegendcolor(AlgoL,:)=CurrColor(mod((AlgoL-1)*ColorStep+1,ColorNb),:);
            else
                mlegendcolor(AlgoL,:)=CurrColor(end,:);
            end
        end
        icolor=zeros(CurrPointNb,3);
        legendcellstr={cell(AlgoNb,1)};
        for AlgoL=1:AlgoNb
            legendcellstr{1}{AlgoL}=strrep(AlgoNames{AlgoL},'_',' ');
            AlgoPos=strmatch(AlgoNames{AlgoL},P.point.algo,'exact');
            try
                icolor(AlgoPos,:)=repmat(mlegendcolor(AlgoL,:),length(AlgoPos),1);
            catch
                'stop'
            end
        end
        paramstruct = struct('viout',3, ...
            'vipcplot',1:4, ...
            'icolor',icolor, ...
            'imarker',imarker,...
            'mlegendcolor',mlegendcolor, ...
            'legendcellstr1',legendcellstr, ...
            'legendcellstr2',legendcellstr2, ...
            'iscreenwrite',0);

        % Verify if exist negative values
        NegPos=find(CurrDataSignals<=0);
        if ~isempty(NegPos)
            CurrDataSignals(NegPos)=eps;
        end
        CurrDataSignals=log2(CurrDataSignals);
        h=figure;
        set(h,'name','all points - algorithms');
        %curvdatSM(DataLog2Signal,paramstruct)
        %replace NaN by 0 because some operation in curvdatSM are not
        %possible with Nan values
        NanPos=find(isnan(CurrDataSignals));
        if ~isempty(NanPos)
            CurrDataSignals(NanPos)=0;
        end
        curvdatSM_mb(CurrDataSignals,paramstruct)
        cd(P.dir.resTree)
        set_figsize('1024px')       
        saveas(h,sprintf('all_allgorithms_%s_%s',P.project.name,date),'png');        
        close(h)
    end
    if isfield(P.point,'factorValues')
        for FactorL=1:length(P.point.factorNames)
            Val=unique(P.point.factorValues{FactorL});
            ValNb=length(Val);

            mlegendcolor=zeros(ValNb,3);
            ColorStep=floor(ColorNb/ValNb);
            for ValL=1:ValNb
                if mod((ValL-1)*ColorStep+1,ColorNb)
                    mlegendcolor(ValL,:)=CurrColor(mod((ValL-1)*ColorStep+1,ColorNb),:);
                else
                    mlegendcolor(ValL,:)=CurrColor(end,:);
                end
            end
            icolor=zeros(CurrPointNb,3);
            legendcellstr={cell(ValNb,1)};
            for ValL=1:ValNb
                if isequal(P.point.factorTypes{FactorL},'str')
                    legendcellstr{1}{ValL}=strrep(Val{ValL},'_',' ');
                    ValPos=strmatch(Val{ValL},P.point.factorValues{FactorL}(CurrPointRanks),'exact');
                else
                    ValPos=find(P.point.factorValues{FactorL}(CurrPointRanks)==Val(ValL));
                    legendcellstr{1}{ValL}=Val(ValL);
                end
                icolor(ValPos,:)=repmat(mlegendcolor(ValL,:),length(ValPos),1);
            end
            paramstruct = struct('viout',3, ...
                'vipcplot',1:2, ...
                'icolor',icolor, ...
                'imarker',imarker,...
                'mlegendcolor',mlegendcolor, ...
                'legendcellstr1',legendcellstr, ...
                'legendcellstr2',legendcellstr2, ...
                'iscreenwrite',0);
            h=figure;
            if LoopL==1
                set(h,'name',sprintf('all points - %s',P.point.factorNames{FactorL}))
            else
                set(h,'name',sprintf('points analyzed by %s - %s',AlgoNames{LoopL-1},P.point.factorNames{FactorL}))
            end
            curvdatSM_mb(CurrDataSignals,paramstruct)
            cd(P.dir.resTree)
            set_figsize('1024px')
            if LoopL==1
                saveas(h,sprintf('%s_%s',P.point.factorNames{FactorL},date),'png');
            else
                saveas(h,sprintf('%s_%s_%s',AlgoNames{LoopL-1},P.point.factorNames{FactorL},date),'png');                
            end
            close(h)
        end
    end
end

