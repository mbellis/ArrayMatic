function trs_prepcomp(Action)
global F P M
Action
switch Action
    case 'groups'
        %load M
        cd(P.dir.data)
        if exist('Comp.mat','file')
            load Comp
        else
            M{1}=[];
        end
        if isfield(P,'grp')
            %select comparisons to be made
            h=figure;
            F.h.prepcomp=h;
            set(h,'name','SELECT COMPARISONS TO BE MADE')
            if P.flag.testAlgo
                Labels={};
                for GrpL=1:length(P.grp.name)
                    EndPos=regexp(P.grp.name{GrpL},['_',P.point.algo{1},'$'],'start');
                    if ~isempty(EndPos)
                        Labels{end+1}=P.grp.name{GrpL}(1:EndPos-1);
                    end
                end
            else
                Labels=P.grp.name;
            end
            P.tmp.labels=Labels;
            GrpNb=length(Labels);
            a=axes;
            set(a,'xlim',[0,GrpNb+1])
            set(a,'ylim',[0,GrpNb+1])
            set(a,'xtick',1:GrpNb)
            set(a,'ytick',1:GrpNb)
            %    set(a,'xticklabel',Labels)
            set(a,'yticklabel',Labels)
            rotateXtickLabel(1:GrpNb,90,Labels)
            %trace lines
            for GrpL=1:GrpNb
                line([GrpL,GrpL],[0,GrpNb+1])
                line([0,GrpNb+1],[GrpL,GrpL])
            end
            set(gcf,'color',[1,1,1])
            %place comparison buttons
            hold on
            for GrpL1=1:GrpNb
                for GrpL2=1:GrpNb
                    if GrpL1~=GrpL2
                        p=plot(GrpL1,GrpL2,'ko','markersize',10);
                        set(p,'ButtonDownFcn','trs_prepcomp(''select comp'')');
                        set(p,'markerfacecolor',[0,0,0]);
                    end
                end
            end
            set(a,'box','on')
            %place OK button
            uicontrol('style','pushbutton','string','Ok','position',[20 20 60 20 ],'CallBack','trs_prepcomp(''ok'')')
            uicontrol('style','pushbutton','string','Cancel','position',[100 20 60 20 ],'CallBack','trs_prepcomp(''cancel'')')


        else
            h=warndlg('first do *** TRANSCRIPTOME/TRANSCRIPTOME ANALYSIS/make groups of points ***');
            waitfor(h)
        end
   
    case 'select comp'
        % coordinates of the clicked pixel
        FirstGrpRank=get(gco,'YData');
        SndGrpRank=get(gco,'XData');
        if P.flag.testAlgo
            CompName=sprintf('%s_vs_%s_%s',strrep(P.tmp.labels{FirstGrpRank},' ',':'),strrep(P.tmp.labels{SndGrpRank},' ',':'),P.point.algo{1});
        else
            CompName=sprintf('%s_vs_%s',strrep(P.tmp.labels{FirstGrpRank},' ',':'),strrep(P.tmp.labels{SndGrpRank},' ',':'));
        end
        if isequal(get(gco,'color'),[0,0,0])
            set(gco,'color',[1,0,0])
            set(gco,'markerfacecolor',[1,0,0]);
            if P.flag.testAlgo
                FirstGrpPos=strmatch(sprintf('%s_%s',P.tmp.labels{FirstGrpRank},P.point.algo{1}),P.grp.name,'exact');
                SndGrpPos=strmatch(sprintf('%s_%s',P.tmp.labels{SndGrpRank},P.point.algo{1}),P.grp.name,'exact');
            else
                FirstGrpPos=FirstGrpRank;
                SndGrpPos=SndGrpRank;
            end
            AddedComp=1;
            if isfield(M{1},'compName')
                Pos=strmatch(CompName,M{1}.compName,'exact');
                if ~isempty(Pos)
                    %assume that 
                    AddedComp=0;
                    if M{1}.made(Pos,1)
                        DoIt=questdlg(sprintf('comparison %s already exists and has been made',CompName),'','make comparison','don''t make comparison','don''t make comparison');
                    else
                        DoIt=questdlg(sprintf('comparison %s already exists and has not  been made',CompName),'','make comparison','don''t make comparison','make comparison');
                    end
                    if isequal(DoIt,'make comparison')
                        M{1}.make(Pos,1)=1;
                    else
                        M{1}.make(Pos,1)=0;
                    end
                    %update other algorithms
                    if P.flag.testAlgo
                        Algos=unique(P.point.algo);
                        %remove already processed algorithm
                        Algos(strmatch(P.point.algo{1},Algos,'exact'))=[];
                        if length(Algos)>1
                            for AlgoL=1:length(Algos)                              
                                CompName=sprintf('%s_vs_%s_%s',strrep(P.tmp.labels{FirstGrpRank},' ',':'),strrep(P.tmp.labels{SndGrpRank},' ',':'),Algos{AlgoL});
                                Pos=strmatch(CompName,M{1}.compName,'exact');
                                if isequal(DoIt,'make comparison')
                                    M{1}.make(Pos,1)=1;
                                else
                                    M{1}.make(Pos,1)=0;
                                end
                            end
                        end
                    end
                else
                    M{1}.compName{end+1,1}=CompName;
                    try
                    M{1}.firstGrpRank(end+1,1)=FirstGrpPos;
                    catch
                        error('process canceled')
                    end
                    M{1}.sndGrpRank(end+1,1)=SndGrpPos;
                    M{1}.firstPointRanks{end+1,1}=P.grp.pointIndex{FirstGrpPos};
                    M{1}.sndPointRanks{end+1,1}=P.grp.pointIndex{SndGrpPos};
                    M{1}.make(end+1,1)=1;
                    M{1}.made(end+1,1)=0;                  
                end
            else
                M{1}.compName{1}=CompName;
                M{1}.firstGrpRank=FirstGrpPos;
                M{1}.sndGrpRank=SndGrpPos;
                M{1}.firstPointRanks{1}=P.grp.pointIndex{FirstGrpPos};
                M{1}.sndPointRanks{1}=P.grp.pointIndex{SndGrpPos};
                M{1}.make=1;
                M{1}.made=0;
            end
            if P.flag.testAlgo & AddedComp
                Algos=unique(P.point.algo);
                %remove already processed algorithm
                Algos(strmatch(P.point.algo{1},Algos,'exact'))=[];
                if length(Algos)>1
                    for AlgoL=1:length(Algos)
                        FirstGrpPos=strmatch(sprintf('%s_%s',P.tmp.labels{FirstGrpRank},Algos{AlgoL}),P.grp.name,'exact');
                        SndGrpPos=strmatch(sprintf('%s_%s',P.tmp.labels{SndGrpRank},Algos{AlgoL}),P.grp.name,'exact');
                        M{1}.compName{end+1,1}=sprintf('%s_vs_%s_%s',strrep(P.tmp.labels{FirstGrpRank},' ',':'),strrep(P.tmp.labels{SndGrpRank},' ',':'),Algos{AlgoL});
                        M{1}.firstGrpRank(end+1,1)=FirstGrpPos;
                        M{1}.sndGrpRank(end+1,1)=SndGrpPos;
                        M{1}.firstPointRanks{end+1,1}=P.grp.pointIndex{FirstGrpPos};
                        M{1}.sndPointRanks{end+1,1}=P.grp.pointIndex{SndGrpPos};
                        M{1}.make(end+1,1)=1;
                        M{1}.made(end+1,1)=0;
                    end
                end
            end
        elseif isequal(get(gco,'color'),[1,0,0])
            set(gco,'color',[0,0,0])
            set(gco,'markerfacecolor',[0,0,0]);
            Pos=strmatch(CompName,M{1}.compName,'exact');
            if ~isempty(Pos)
                M{1}.make(Pos)=0;
                if P.flag.testAlgo
                    Algos=unique(P.point.algo);
                    %remove already processed algorithm
                    Algos(strmatch(P.point.algo{1},Algos,'exact'))=[];
                    if length(Algos)>1
                        for AlgoL=1:length(Algos)
                            CompName=sprintf('%s_vs_%s_%s',strrep(P.tmp.labels{FirstGrpRank},' ',':'),strrep(P.tmp.labels{SndGrpRank},' ',':'),Algos{AlgoL});
                            Pos=strmatch(CompName,M{1}.compName,'exact');                           
                            M{1}.make(Pos,1)=0;
                        end
                    end
                end
            end
        end

    case 'ok'
        %eliminate comparisons with make=made=0 (selected once or more and unselected latter)
        FieldNames=fieldnames(M{1});
        ClearComp=[];
        for CompL=length(M{1}.compName):-1:1
            if M{1}.make(CompL)==0&M{1}.made(CompL)==0
                ClearComp=[ClearComp;CompL];
            end
        end
        if ~isempty(ClearComp)
            M{1}=del_structitems(M{1},ClearComp);
        end
        %add other algorithm comparison if necessary
      
        cd(P.dir.data)
        save Comp M
        M=[];
        delete(F.h.prepcomp)
        F.h=rmfield(F.h,'prepcomp');

    case 'cancel'
        M=[];
        delete(F.h.prepcomp)
        F.h=rmfield(F.h,'prepcomp');
end