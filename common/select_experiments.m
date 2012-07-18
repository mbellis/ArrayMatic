function [ExpRanks,Ok]=select_experiments(SelectionMode,Title)
global P
Ok=1;
if nargin~=2
    h=warndlg('select_experiments needs two input parameters');
    waitfor(h)
    Ok=0;
end
if Ok
    if ~isequal(SelectionMode,'single')&&~isequal(SelectionMode,'multiple')
        h=errodlg('SelectionMode must be equal to single or multiple');
        waitfor(h)
        error('process canceled')
    end
    List={};
    for ExpL=1:length(P.exp.name)
            List{end+1,1}=sprintf('%03u - %s',ExpL,P.exp.name{ExpL});
    end
    [ExpRanks,Ok]=listdlg('liststring',List,'selectionmode',SelectionMode,'listsize',[400 300],'promptstring','select','name',Title);
    if Ok
    else
        if nargout==1
            h=errodlg('Canceled by user');
            waitfor(h)
            error('process canceled')
        else
            Ok=0;
        end

    end
end
