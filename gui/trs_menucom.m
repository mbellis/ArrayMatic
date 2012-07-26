%TRS_MENUCOM
% menu programme
% function amc_menu_com triggered by uicontrols of TRANSCRIPTOME menu in amc_menu.fig

function varargout=amc_menu_com(action,varargin)

global K C D F M N O P R S T X DataRanks

%% action
gco_tag=get(gco,'tag')
Properties=get(gco);
if isfield(Properties,'Callback')
    gco_callback=get(gco,'callback');
end
% empty root tag
set(0,'tag','')
action

switch action
    case 'new project'
        trs_newproject

    case 'load project'
        trs_loadproject

    case 'quit project'
%% quit project
        DataRanks=[];
        P=[];
        R=[];
        S=[];
        T=[];
        set(get(F.gh.amcmenu.TRANSCRIPTOME,'children'),'enable','off')
        set(get(F.gh.amcmenu.TRANSCRIPTOME,'children'),'checked','off')
        TrsChildren=get(F.gh.amcmenu.TRANSCRIPTOME,'children');
        for ChildL=1:length(TrsChildren)
            set(get(TrsChildren(ChildL),'children'),'enable','off')
            set(get(TrsChildren(ChildL),'children'),'checked','off')
        end
        set(F.gh.amcmenu.LOAD_A_PROJECT,'enable','on')
        set(F.gh.amcmenu.NEW_PROJECT,'enable','on')
        set(F.gh.amcmenu.PROJECT_MANAGEMENT,'enable','on')
        K.par.menu='';
        P.tmp.menu={};
        Children1=get(F.gh.amcmenu.TRANSCRIPTOME,'children');
        for ChildL1=1:length(Children1)
            if ~isempty(get(Children1(ChildL1),'TAG'))
                P.tmp.menu{end+1,1}=sprintf('set(F.gh.amcmenu.%s,''Enable'',''%s'')',get(Children1(ChildL1),'TAG'),get(Children1(ChildL1),'Enable'));
                P.tmp.menu{end+1,1}=sprintf('set(F.gh.amcmenu.%s,''Checked'',''%s'')',get(Children1(ChildL1),'TAG'),get(Children1(ChildL1),'Checked'));
            end
            Children2=get(Children1(ChildL1),'children');
            for ChildL2=1:length(Children2)
                if ~isempty(get(Children2(ChildL2),'TAG'))
                    P.tmp.menu{end+1,1}=sprintf('set(F.gh.amcmenu.%s,''Enable'',''%s'')',get(Children2(ChildL2),'TAG'),get(Children2(ChildL2),'Enable'));
                    P.tmp.menu{end+1,1}=sprintf('set(F.gh.amcmenu.%s,''Checked'',''%s'')',get(Children2(ChildL2),'TAG'),get(Children2(ChildL2),'Checked'));
                end
            end
        end


    case 'keep data on disk'
%% keep data on disk
        if P.flag.loadData==0
            DataRanks=[];
            P.flag.loadData=1;
            cd(P.dir.project)
            eval(sprintf('save %s P',P.project.name))
            set(F.gh.amcmenu.keep_data_on_disk,'checked','on')
            set(F.gh.amcmenu.keep_data_on_disk,'enable','off')
            set(F.gh.amcmenu.keep_data_in_memory,'checked','off')
            set(F.gh.amcmenu.keep_data_in_memory,'enable','on')
        end


    case 'keep data in memory'
%% keep data in memory
        if P.flag.loadData
            DataRanks=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le');
            P.flag.loadData=0;
            cd(P.dir.project)
            eval(sprintf('save %s P',P.project.name))
            set(F.gh.amcmenu.keep_data_on_disk,'checked','off')
            set(F.gh.amcmenu.keep_data_on_disk,'enable','on')
            set(F.gh.amcmenu.keep_data_in_memory,'checked','on')
            set(F.gh.amcmenu.keep_data_in_memory,'enable','off')
        end

    case 'same algorithm'
%% same algorithm
        set(F.gh.amcmenu.same_algorithm,'checked','on')
        set(F.gh.amcmenu.same_algorithm,'enable','off')
        set(F.gh.amcmenu.several_algorithms,'checked','off')
        set(F.gh.amcmenu.several_algorithms,'enable','on')

    case 'several algorithms'
%% severalalgorithms
        set(F.gh.amcmenu.several_algorithms,'checked','on')
        set(F.gh.amcmenu.several_algorithms,'enable','off')
        set(F.gh.amcmenu.same_algorithm,'checked','off')
        set(F.gh.amcmenu.same_algorithm,'enable','on')


    case 'select used points'
%% select used points
        %keep in memory the number of used points in biol cond marked as
        %not used
        NotUsedBiol=find(P.biol.used==0);
        UsedPointNbs=zeros(length(NotUsedBiol),1);
        if ~isempty(NotUsedBiol)
            for BiolL=1:length(NotUsedBiol)
                UsedPointNb=0;
                for PointL=1:length(P.biol.pointIndex{NotUsedBiol(BiolL)})
                    if P.point.used(P.biol.pointIndex{NotUsedBiol(BiolL)}(PointL))
                        UsedPointNb=UsedPointNb+1;
                    end
                end
                UsedPointNbs(BiolL)=UsedPointNb;
            end
        end
        %keep in memory the current values which will be eventually changed
        MemUsedPoint=P.point.used;
        %% select not used points
        SelPar=inputdlg({'missing nb >=';'largest run >='},'select the parameters to �liminate wrong points',1,{'500';'500'});
        if ~isempty(SelPar)
            MissLimit=str2num(SelPar{1});
            RunLimit=str2num(SelPar{2});
            if isempty(MissLimit)|isempty(RunLimit)
                h=errordlg('Enter numbers');
                waitfor(h)
                error('process canceled')
            end

            for PointL=1:P.point.nb
                if P.point.nanNb(PointL)>=MissLimit|P.point.runNb(PointL)>=RunLimit
                    P.point.used(PointL)=0;
                end
            end
        end
        %% manual selection
        Sel=questdlg('Do you want to eliminate other points by manual selection','','yes','no','yes');
        if isequal(Sel,'yes')
            SelType=questdlg('How do you want to select points','','points','biological conditions','experiments','points');
            switch SelType
                case 'points'
                    [SelIndex,Ok]=select_points('multiple','SELECT THE POINT TO BE MARKED AS NOT USED');
                    if Ok
                        P.point.used(SelIndex)=0;
                    end
                    %update used biol conditions
                    %not used => used
                    NotUsedBiol=find(P.biol.used==0);
                    if ~isempty(NotUsedBiol)
                        for BiolL=1:length(NotUsedBiol)
                            UsedPointNb=0;
                            for PointL=1:length(P.biol.pointIndex{NotUsedBiol(BiolL)})
                                if P.point.used(P.biol.pointIndex{NotUsedBiol(BiolL)}(PointL))
                                    UsedPointNb=UsedPointNb+1;
                                end
                            end
                            if UsedPointNb>0 && UsedPointNb~=UsedPointNbs(BiolL)
                                UseFlag=questdlg(sprintf('%s was marked as not used with %u used points but has now %u used points',P.biol.name{NotUsedBiol(BiolL)},UsedPointNbs(BiolL),UsedPointNb),'mark it as used','keep as not used','keep as not used')
                                if isequal(UseFlag,'mark it as used')
                                    P.biol.used(NotUsedBiol(BiolL))=1;
                                end
                            end
                        end
                    end
                    %used => not used
                    NotUsedPoint=find(P.point.used==0);
                    if ~isempty(NotUsedPoint)
                        BiolRanks=P.point.biolRank(NotUsedPoint);
                        BiolRank=unique(BiolRanks);
                        for BiolL=1:length(BiolRank)
                            BiolPos=find(BiolRanks==BiolRank(BiolL));
                            CurrPoints=NotUsedPoint(BiolPos);
                            KeptPoint=setdiff(P.biol.pointIndex{BiolRank(BiolL)},CurrPoints);
                            if isempty(KeptPoint)
                                P.biol.used(BiolRank(BiolL))=0;
                            end
                        end
                    end
                    cd(P.dir.project)
                    eval(sprintf('save %s P',P.project.name))

                case 'biological conditions'
                    [SelIndex,Ok]=select_biolconditions('multiple','SELECT THE POINT TO BE MARKED AS NOT USED');
                    if Ok
                        P.biol.used(SelIndex)=0;
                        for BiolL=1:length(SelIndex)
                            P.point.used(P.biol.pointIndex{SelIndex(BiolL)})=0;
                        end
                    end

                case 'experiments'
                    [SelIndex,Ok]=select_experiments('multiple','SELECT THE POINT TO BE MARKED AS NOT USED');
                    if Ok
                        P.exp.used(SelIndex)=0;
                        for ExpL=1:length(SelIndex)
                            BiolIndex=unique(P.point.biolRank(P.exp.pointIndex{SelIndex(ExpL)}));
                            for BiolL=1:length(BiolIndex)
                                P.biol.used(BiolIndex(BiolL))=0;
                                P.point.used(P.biol.pointIndex{BiolIndex(BiolL)})=0;
                            end
                        end
                    end
            end
        end
        cd(P.dir.project)
        eval(sprintf('save %s P',P.project.name))
       



    case 'plot signal vs rank'
%% plot signal vs rank


        %recover color map

        %figure with all curves rank vs log2(signal)
        h=figure;
        set(h,'color',[1,1,1])
        set(h,'name','all signal vs rank curves')
        subplot(1,2,1)
        hold on
        subplot(1,2,2)
        hold on
        [Legend,Colors]=fig_com('legend',h,'EastOutside');

        for PointL=1:P.point.nb
            if P.flag.loadData
                CurrRanks=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,PointL,'single','ieee-le',[],PointL);
            else
                CurrRanks=DataRanks(:,PointL);
            end
            CurrSignals=load_data('DataSignals.float32le',P.dir.data,P.chip.currProbeSetNb,PointL,'single','ieee-le',[],PointL);

            %shift Signals not to have negative values
            CurrSignals=CurrSignals+abs(min(CurrSignals))+1;

            [CurrRanks SortIndex]=sort(CurrRanks);
            CurrSignals=CurrSignals(SortIndex);

            subplot(1,2,1)
            if P.flag.testAlgo
                ColorPos=strmatch(P.point.algo{PointL},Legend,'exact');
            else
                ColorPos=strmatch(strrep(P.biol.name{P.point.biolRank(PointL)},'_',' '),Legend,'exact');
            end
            plot(CurrRanks,log2(CurrSignals),'color',Colors(ColorPos,:))
            if PointL==1
                title('log2(signal) vs rank')
                ylabel('log2(Signal+min(Signal)+1)')
                xlabel('rank')
            end


            CurrSignals=log2(CurrSignals);
            BinVal=histc(CurrSignals,[0:16/100:16]);
            figure(h)
            subplot(1,2,2)
            plot([0:16/100:16],BinVal,'color',Colors(ColorPos,:));
            set(gca,'box','on')
            if PointL==1
                title('log2(signal) distribution')
                xlabel('log2(CurrSignals)')
                ylabel('frequency')
            end

        end



        figure(h)
        set(h,'units','normalized')
        subplot(1,2,1)
        set(gca,'box','on')
        set(gca,'xlim',[0,100])
        set(gca,'ylim',[0,16])
        set(gca,'box','on')
        subplot(1,2,2)
        set(gca,'box','on')
        set(gca,'xlim',[0,16])

        set_figsize('1024px')
        set(gcf,'Color',[1 1 1])
        cd(P.dir.resCalib)

        h=msgbox(sprintf('resize and save as \nsignals_vs_ranks_%s_%s.png',P.project.name,date));
        waitfor(h)
        plot2svg(sprintf('%s_%s.svg',P.project.name,date))


    case 'calculate cv'
%% CALCULATE CV

        BlOC_SIZE=1000;
        BlocNb=ceil(P.chip.probeSetNb/BlOC_SIZE);
        cd(P.dir.data)
        P.chip.cv=zeros(P.chip.currProbeSetNb,1);
        if BlocNb==1
            CurrRanks=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le');
            MeanRank=mean(CurrRanks,2);
            StdRank=std(CurrRanks,1,2);
            P.chip.cv(1:end)=MeanRank./StdRank;
        else
            for BlocL=1:BlocNb-1
                CurrRanks=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le',((BlocL-1)*BlOC_SIZE)+1:BlocL*BlOC_SIZE,1:P.point.nb);
                MeanRank=mean(CurrRanks,2);
                StdRank=std(CurrRanks,1,2);
                P.chip.cv(((BlocL-1)*BlOC_SIZE)+1:BlocL*BlOC_SIZE)=MeanRank./StdRank;
            end
            CurrRanks=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le',((BlocNb-1)*BlOC_SIZE)+1:P.chip.currProbeSetNb,1:P.point.nb);
            MeanRank=mean(CurrRanks,2);
            StdRank=std(CurrRanks,1,2);
            P.chip.cv(((BlocNb-1)*BlOC_SIZE)+1:end)=MeanRank./StdRank;
        end
        cd(P.dir.project)
        eval(sprintf('save %s P',P.project.name))

        h=figure;
        set(h,'name',sprintf('VARIATION COEFFICIENTS OF RANKS - %s',P.project.name));
        hist(P.chip.cv,100);
        xlabel('variation coefficient of ranks')
        ylabel('frequency')
        cd(P.dir.resProp)
        saveas(h,sprintf('cv_ranks_%s_%s',P.project.name,date),'png')
        
        
end


