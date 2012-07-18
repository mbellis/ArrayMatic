%AMC_MENUCOM
% menu programme
% function amc_menucom triggered by uicontrols of amc_menu.fig

function varargout=amc_menucom(action,varargin)

global K C D F M N O P R S T X Data

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
    case 'keep data on disk'
        if P.flag.loadData
            clear Data
            P.flag.loadData=0;
            cd(P.dir.project)
            eval(sprintf('save %s P',P.project.name))
        end
    case 'keep data in memory'
        if P.flag.loadData==0
            Data=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le');
            P.flag.loadData=1;
            cd(P.dir.project)
            eval(sprintf('save %s P',P.project.name))
        end
    case 'make median comparison'
        UpdateS=questdlg('Do you want to recalculate calibration set if they exist ?','','yes','no','cancel','no');
        if ~isequal(UpdateS,'cancel')
            if isequal(UpdateS,'yes')
                UpdateSFlag=1;
            else
                UpdateSFlag=0;
            end
            trs_mediancomp(UpdateSFlag);
        end
end
