% TRS_LOADPROJECT - load project of transcriptomic analysis
% verify that the directory has a correct structure (can 
% manage with files or directories that have be renamed)

% c) Michel Bellis
% arraymatic@gmail.com

% NO INPUT

% OUTPUT

% VERSIONS
%
% V01 - 2010 05 17- Refactoring existing version (2000-2010)

function trs_loadproject
global F K P S DataRanks

[ProjectFile,ProjectDir]=uigetfile('*.mat','Select the project file ');
if isequal(class(ProjectFile),'double')
    h=warndlg('Process canceled');
    waitfor(h)
else
    MatPos=findstr(upper(ProjectFile),'.MAT');
    ProjectName=ProjectFile(1:MatPos(end)-1);
    cd(ProjectDir);
    load(ProjectFile)
    if ~isequal(ProjectName,P.project.name)
        h=warndlg(sprintf('Project name has been changed from %s to %s',P.project.name,ProjectName));
        waitfor(h)        
        DirNames={'data','results'};
        for DirL=1:2
            if ~exist(sprintf('%s_%s',ProjectName,DirNames{DirL}),'dir')
                if exist(sprintf('%s_%s',P.project.name,DirNames{DirL}),'dir')
                    [Status,Message,MessageId]=movefile(sprintf('%s_%s',P.project.name,DirNames{DirL}),sprintf('%s_%s',ProjectName,DirNames{DirL}));
                else
                    h=warndlg(sprintf('Should exist %s_%s. If it exist at another location and you want to use it copy it there before continuing',P.project.name,DirNames{DirL}));
                    waitfor(h)
                    if ~exist(sprintf('%s_%s',ProjectName,DirNames{DirL}),'dir')
                        [Status,Message,MessageId]=mkdir(ProjectDir,sprintf('%s_%s',ProjectName,DirNames{DirL}));
                        if isequal(DirNames{DirL},'results')
                            [Status,Message,MessageId]=mkdir(fullfile(ProjectDir,sprintf('%s_%s',ProjectName,DirNames{DirL})),'calib');
                            [Status,Message,MessageId]=mkdir(fullfile(ProjectDir,sprintf('%s_%s',ProjectName,DirNames{DirL})),'tree');
                            [Status,Message,MessageId]=mkdir(fullfile(ProjectDir,sprintf('%s_%s',ProjectName,DirNames{DirL})),'comp');
                        end
                    end
                end
            end
        end

    end
     %update directories
     SaveIt=0;
    if ~isequal(P.dir.project,ProjectDir)
        P.dir.project=ProjectDir;
        SaveIt=1;
    end
    if ~isequal(P.project.name,ProjectName)               
        P.dir.data=fullfile(P.dir.project,sprintf('%s_data',ProjectName));
        P.dir.results=fullfile(P.dir.project,sprintf('%s_results',ProjectName));
        P.dir.resProp=fullfile(P.dir.project,sprintf('%s_results',ProjectName),'prop');
        P.dir.resCalib=fullfile(P.dir.project,sprintf('%s_results',ProjectName),'calib');
        P.dir.resTree=fullfile(P.dir.project,sprintf('%s_results',ProjectName),'tree');
        P.dir.resComp=fullfile(P.dir.project,sprintf('%s_results',ProjectName),'comp');        
        P.project.name=ProjectName;
        SaveIt=1;
    end
    if SaveIt
        cd(P.dir.project)
        eval(sprintf('save %s P',ProjectName));
    end
    %load data for displaying signals and keep them in memory if necessary
    %of possible
    DataRanks=[];
    try
        cd(P.dir.resProp)
        %do not load ifalready exist figure
        DisplayFig=1;
        DirInfo=dir;
        for DirL=1:length(DirInfo)
            if ~isempty(findstr(sprintf('ranks_heathmap_%s',P.project.name),DirInfo(DirL).name))
                DisplayFig=0;
                break
            end
        end
        if P.flag.loadData==0
            DataRanks=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le');
        end
        if DisplayFig
            if P.point.nb>500&P.flag.loadData
                DoIt=questdlg(sprintf('There exist %u points. Do you want to load data to display the first 2000 rank values ?',P.point.nb),'','no','yes','no')
            end
            if isequal(DoIt,'yes')
                if P.flag.loadData
                    DataRanks=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le');
                end
                h=figure;
                set(h,'name','upper left corner of Data')
                if size(P.point.nb)<=2000
                    image(DataRanks(1:2000,:)');
                    xlabel('probe set rank')
                    ylabel('point rank')
                else
                    image(DataRanks(1:2000,:));
                    ylabel('probe set rank')
                    xlabel('point rank')
                end
                if P.flag.testAlgo
                    AlgoSize=length(strmatch(P.point.algo{1},P.point.algo,'exact'));
                    set(gca,'ytick',[round(AlgoSize/2):AlgoSize:P.point.nb])
                    Algos={};
                    for PointL=1:P.point.nb
                        if isempty(strmatch(P.point.algo{PointL},Algos,'exact'))
                            Algos{end+1,1}=P.point.algo{PointL};
                        end
                    end
                    set(gca,'yticklabel',Algos)
                end
                set(gcf,'color',[1,1,1])
                set(gca,'clim',[0,100])                
                %P.flag.loadData=1
                if P.flag.loadData
                    DataRanks=[];                
                end
                h=msgbox(sprintf('eventually resize and save as \nranks_heatmap_%s_%s.png',P.project.name,date));
                waitfor(h)
            end
        end
    catch
        if isempty(DataRanks)&P.flag.loadData==0
            %loading process failed
            h=warndlg(sprintf('Data had not been loaded (may be not enough memory for %u data points',P.point.nb));
            waitfor(h)
            h=warndlg('Data are kept disk during this session');
            waitfor(h)
            trs_menucom('keep data on disk')
        end
    end
    %enable menu items
    set(get(F.gh.amcmenu.TRANSCRIPTOME,'children'),'enable','on')
    TrsChildren=get(F.gh.amcmenu.TRANSCRIPTOME,'children');
    for ChildL=1:length(TrsChildren)
        set(get(TrsChildren(ChildL),'children'),'enable','on')
        set(get(TrsChildren(ChildL),'children'),'checked','off')        
    end
    if P.flag.loadData
        set(F.gh.amcmenu.keep_data_on_disk,'checked','on')
        set(F.gh.amcmenu.keep_data_on_disk,'enable','off')
        set(F.gh.amcmenu.keep_data_in_memory,'checked','off')
        set(F.gh.amcmenu.keep_data_in_memory,'enable','on')
    else
        set(F.gh.amcmenu.keep_data_in_memory,'checked','on')
        set(F.gh.amcmenu.keep_data_in_memory,'enable','off')
        set(F.gh.amcmenu.keep_data_on_disk,'checked','off')
        set(F.gh.amcmenu.keep_data_on_disk,'enable','on')
    end 
    if P.flag.testAlgo==0
        set(F.gh.amcmenu.same_algorithm,'checked','on')
        set(F.gh.amcmenu.same_algorithm,'enable','off')
        set(F.gh.amcmenu.several_algorithms,'checked','off')
        set(F.gh.amcmenu.several_algorithms,'enable','on')
    else
        set(F.gh.amcmenu.several_algorithms,'checked','on')
        set(F.gh.amcmenu.several_algorithms,'enable','off')
        set(F.gh.amcmenu.same_algorithm,'checked','off')
        set(F.gh.amcmenu.same_algorithm,'enable','on')
    end 
    if isequal(P.par.analType,'network')
        set(F.gh.amcmenu.plot_signal_vs_rank,'enable','off')
        set(F.gh.amcmenu.plot_rank_vs_rank,'enable','off')
        set(F.gh.amcmenu.plot_algo_vs_algo,'enable','off')
    end
    set(F.gh.amcmenu.LOAD_A_PROJECT,'checked','on')
    set(F.gh.amcmenu.LOAD_A_PROJECT,'enable','off')
    set(F.gh.amcmenu.NEW_PROJECT,'enable','off')
    switch P.par.analType
        case {'transcription','network'}
            Children1=get(F.gh.amcmenu.TRANSCRIPTOME,'children');
    end
    P.tmp.menu={};
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
    switch P.par.analType
        case 'transcription'
            set(F.gh.amcmenu.ANALYSIS_OF_USER_GROUPS,'enable','on')
            set(F.gh.amcmenu.NETWORK_CONSTRUCTION,'enable','off')
        case 'network'
            set(F.gh.amcmenu.ANALYSIS_OF_USER_GROUPS,'enable','off')
            set(F.gh.amcmenu.NETWORK_CONSTRUCTION,'enable','on')
    end

end

