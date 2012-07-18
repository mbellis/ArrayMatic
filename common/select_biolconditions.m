function [BiolRanks,Ok]=select_biolconditions(SelectionMode,Title)
global P
Ok=1;
if nargin~=2
    h=warndlg('select_biolconditions needs two input parameters');
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
    for BiolL=1:P.biol.nb
            List{end+1,1}=sprintf('%03u - %s',BiolL,P.biol.name{BiolL});
    end
    [BiolRanks,Ok]=listdlg('liststring',List,'selectionmode',SelectionMode,'listsize',[400 300],'promptstring','select','name',Title);
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
