function [PointRanks,Ok]=select_points(SelectionMode,Title)
global P
Ok=1;
if nargin~=2
    if nargout==1
        h=errodlg('select_points needs two parameters');
        waitfor(h)
        error('process canceled')
    else
        h=warndlg('select_points needs two parameters');
        waitfor(h)
        Ok=0;
    end
end
if Ok
    if ~isequal(SelectionMode,'single')&&~isequal(SelectionMode,'multiple')
        h=errodlg('SelectionMode must be equal to single or multiple');
        waitfor(h)
        error('process canceled')
    end
    List={};
    Header='PointName';
    if isfield(P.point,'factorTypes')
        for FactorL=1:length(P.point.factorTypes)
            Header=sprintf('%s_%s',Header,P.point.factorNames{FactorL});
        end
    end
    PointRanks=[];
    if P.flag.testAlgo
        PointIndex=strmatch(P.point.algo{1},P.point.algo,'exact');
    else
        PointIndex=1:P.point.nb;
    end
    for PointL=1:length(PointIndex)
        if P.point.used(PointIndex(PointL))
            CurrPoint=P.point.name{PointIndex(PointL)};
            if isfield(P.point,'factorTypes')
                for FactorL=1:length(P.point.factorTypes)
                    switch P.point.factorTypes{FactorL}
                        case 'char'
                            CurrPoint=sprintf('%s_%c',CurrPoint,P.point.factorValues{FactorL}(PointIndex(PointL)));
                        case 'str'
                            CurrPoint=sprintf('%s_%s',CurrPoint,P.point.factorValues{FactorL}{PointIndex(PointL)});
                        case 'int'
                            CurrPoint=sprintf('%s_%u',CurrPoint,P.point.factorValues{FactorL}(PointIndex(PointL)));
                        case 'float'
                            CurrPoint=sprintf('%s_%f',CurrPoint,P.point.factorValues{FactorL}(PointIndex(PointL)));
                    end
                end
            end
            List{end+1,1}=CurrPoint;
            PointRanks=[PointRanks;PointIndex(PointL)];
        end
    end

    [PointIndex,Ok]=listdlg('liststring',List,'selectionmode',SelectionMode,'listsize',[400 300],'promptstring',Header,'name',Title);
    if Ok
        PointRanks=PointRanks(PointIndex);
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
