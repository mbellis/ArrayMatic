%PROG: functions used by programer
function prog(action)
global F K P
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
    case 'reload amcmenu'
        F.h.amcmenu=findobj('tag','amc_menu');
        if ~isempty(F.h.amcmenu)
            delete(F.h.amcmenu(1:end))
        end
        F.h.amcmenu=openfig('amc_menu');
        set(F.h.amcmenu,'visible','on');
        set(gcf,'pointer','watch')

        F.gh.amcmenu=guihandles(F.h.amcmenu);
        guidata(F.h.amcmenu,F.gh.amcmenu)
        set(0,'tag','')
        if isfield(P,'tmp')
            if isfield(P.tmp,'menu')
                for MenuL=1:length(P.tmp.menu)
                    eval(P.tmp.menu{MenuL})
                end
            end
        end
end