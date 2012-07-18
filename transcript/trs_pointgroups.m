function trs_pointgroups(Action)
global P
switch Action
    case 'make groups'
        Ok=1;
        MakeGrp=0;
        while Ok
            [PointIndex,Ok]=select_points('multiple','SELECT POINTS TO MAKE A GROUP');
            if Ok
                NameFlag=1;
                while NameFlag
                    GrpName=inputdlg('give the name of the current group','');
                    GrpName=GrpName{1};
                    if P.flag.testAlgo
                       UsedGrpName=sprintf('%s_%s',GrpName,P.point.algo{1});
                    else
                        UsedGrpName=GrpName;
                    end
                    if ~isempty(UsedGrpName)
                        if isfield(P,'grp')
                            if ~isempty(strmatch(UsedGrpName,P.grp.name,'exact'))
                                h=warndlg(sprintf('%s already used. Give another name',UsedGrpName));
                                waitfor(h)
                            else
                                NameFlag=0;
                                P.grp.name{end+1,1}=UsedGrpName;
                                P.grp.pointIndex{end+1,1}=PointIndex';
                            end
                        else
                            NameFlag=0;
                            P.grp.name{1}=UsedGrpName;
                            P.grp.pointIndex{1}=PointIndex';
                        end
                    end
                end
                MakeGrp=1;
                if P.flag.testAlgo
                    Algos=unique(P.point.algo);
                    Algos(strmatch(P.point.algo{1},Algos,'exact'))=[];
                    if ~isempty(Algos)
                        for AlgoL=1:length(Algos)
                            CurrAlgo=Algos{AlgoL};                                                        
                            AlgoPos=strmatch(CurrAlgo,P.point.algo,'exact');
                            CurrPointIndex=[];                            
                            for PointL=1:length(PointIndex)
                                NamePos=strmatch(P.point.name{PointIndex(PointL)},P.point.name,'exact');
                                CurrPointIndex=[CurrPointIndex,intersect(AlgoPos,NamePos)];
                            end
                            P.grp.name{end+1,1}=sprintf('%s_%s',GrpName,CurrAlgo);;
                            P.grp.pointIndex{end+1,1}=CurrPointIndex;                            
                        end
                    end
                end
            end
        end
        if MakeGrp
            cd(P.dir.project)
            eval(sprintf('save %s P;',P.project.name))
        end
    case 'delete groups'
        if isfield(P,'grp')
            [GrpIndex,Ok]=listdlg('liststring',P.grp.name,'selectionmode','multiple','listsize',[400 300],'promptstring','Select group(s) to be deleted','name','trs_groups');
            if Ok
                P.grp.name(GrpIndex)=[];
                P.grp.pointIndex(grpIndex)=[];
                if isempty(P.grp.name)
                    P=rmfield(P,'grp');
                end
                cd(P.dir.project)
                eval(sprintf('save %s P;',P.project.name))
            end
        end
end