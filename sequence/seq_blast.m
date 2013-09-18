cd(fullfile('E:','sosma','collab','2012_delsert','séquences'))
%fid=fopen('blastn_SRR334212_02.txt','r');
%fid=fopen('SRR334212_plus_CDNA_01_gigas.txt','r');
fid=fopen('exons_SRR334212_plus_CDNA_02_gigas.txt','r');
CurrLine='';
NextQuery=1;
Queries={};
Scaffolds=[];
Scaffold={};
LineRank=0;
while NextQuery
    %find the next query
    while isempty(findstr(CurrLine,'Query='))
        CurrLine=fgetl(fid);
        if ~ischar(CurrLine)
            NextQuery=0;
            break
        end
    end
    %process the current query
    if NextQuery        
        CurrQuery=CurrLine(8:end);
        if isequal(CurrQuery,'SRR334212.3995242')
            'stop'
        end        
        QueryRank=strmatch(CurrQuery,Queries,'exact');
        if isempty(QueryRank)
            Queries{end+1,1}=CurrQuery;
            QueryRank=length(Queries);
        end
        NextScaffold=1;
        while NextScaffold
            %find the next scaffold
            while isempty(findstr(CurrLine,'> scaffold'))
                CurrLine=fgetl(fid);
                if ~ischar(CurrLine)
                    NextScaffold=0;
                    break
                end
            end
            %process the next Scaffold
            if NextScaffold             
                CurrScaffold=str2num(CurrLine(11:end));
                ScaffoldRank=find(Scaffolds==CurrScaffold);
                if isempty(ScaffoldRank)
                    Scaffolds(end+1,1)=CurrScaffold;
                    ScaffoldRank=length(Scaffolds);
                    Scaffold{ScaffoldRank}.queryRank=[];
                    Scaffold{ScaffoldRank}.queryPos=[];
                    Scaffold{ScaffoldRank}.scaffoldPos=[];
                    Scaffold{ScaffoldRank}.scaffoldStrand=[];                    
                end
               
                %find the next alignment
                NextScore=1;
                while NextScore
                    while isempty(findstr(CurrLine,'Score'))
                        CurrLine=fgetl(fid);
                        if ~ischar(CurrLine)
                            NextScore=0;
                            break
                        end
                    end
                    if NextScore
                        NextScaffoldQuery=1;
                        while NextScaffoldQuery
                            %find the next query
                            while isempty(findstr(CurrLine,'Query '))
                                CurrLine=fgetl(fid);
                                if ~ischar(CurrLine)
                                    NextScaffoldQuery=0;
                                    break
                                end
                            end
                            if NextScaffoldQuery
                                %process the current alignment
                                QueryPos=regexp(CurrLine,'\d+','match');
                                Scaffold{ScaffoldRank}.queryRank(end+1,1)=QueryRank;
                                Scaffold{ScaffoldRank}.queryPos(end+1,1)=str2num(QueryPos{1});
                                Scaffold{ScaffoldRank}.queryPos(end,2)=str2num(QueryPos{2});
                                %find the targeted scaffold position
                                CurrLine=fgetl(fid);
                                CurrLine=fgetl(fid);
                                if ~isequal('Sbjct',CurrLine(1:5))
                                    h=errordlg('wrong structure');
                                    waitfor(h)
                                else
                                    ScaffoldPos=regexp(CurrLine,'\d+','match');
                                    Scaffold{ScaffoldRank}.scaffoldPos(end+1,1)=str2num(ScaffoldPos{1});
                                    if str2num(ScaffoldPos{2})>Scaffold{ScaffoldRank}.scaffoldPos(end,1)
                                        Scaffold{ScaffoldRank}.scaffoldPos(end,2)=str2num(ScaffoldPos{2});
                                        Scaffold{ScaffoldRank}.scaffoldStrand(end+1,1)=1;
                                    else
                                        Scaffold{ScaffoldRank}.scaffoldPos(end,2)=Scaffold{ScaffoldRank}.scaffoldPos(end,1);
                                        Scaffold{ScaffoldRank}.scaffoldPos(end,1)=str2num(ScaffoldPos{2});
                                        Scaffold{ScaffoldRank}.scaffoldStrand(end+1,1)=-1;
                                    end


                                    %test the following lines
                                    while 1
                                        CurrLine=fgetl(fid);
                                        if ~ischar(CurrLine)
                                            NextQuery=0;
                                            NextScaffold=0;
                                            NextScore=0;
                                            NextScaffoldQuery=0;
                                            break
                                        else
                                            if isempty(findstr(CurrLine,'Score'))&isempty(findstr(CurrLine,'> scaffold'))&isempty(findstr(CurrLine,'Query='))
                                                if ~isempty(findstr(CurrLine,'Query '))
                                                    break
                                                end
                                            else
                                                NextScaffoldQuery=0;
                                                if ~isempty(findstr(CurrLine,'> scaffold'))
                                                    NextScore=0;
                                                elseif ~isempty(findstr(CurrLine,'Query='))
                                                    NextScore=0;
                                                    NextScaffold=0;
                                                end
                                                break
                                            end                                            
                                        end
                                    end
                                end
                            end % if NextScaffoldQuery
                        end %while NextScaffoldQuery
                    end %if NextScore
                end %while NextScore
            end %if NextScaffold
            %search eventual Scaffold
            if isempty(findstr(CurrLine,'Query=')) & isempty(findstr(CurrLine,'> scaffold'))
                while 1
                    CurrLine=fgetl(fid);
                    if ~ischar(CurrLine)
                        NextScaffold=0;
                        NextQuery=0;
                        break
                    else
                        if isempty(findstr(CurrLine,'> scaffold'))
                            if isempty(findstr(CurrLine,'Query='))
                                NextScaffold=0;
                                break
                            end
                        else
                            break
                        end
                    end
                end
            else
                if ~ischar(CurrLine)
                    NextScaffold=0;
                    NextQuery=0;
                    break
                else
                    if ~isempty(findstr(CurrLine,'Query='))
                        NextScaffold=0;
                    end
                end
            end
        end %while NextScaffold
    end %if NextQuery
end %while NextQuery
'stop'
%plot results
for ScaffoldL=1:length(Scaffolds)
    CurrScaffold=Scaffold{ScaffoldL};
    h=figure;
    set(gcf,'color',[1,1,1])
    set(h,'name',sprintf('scaffold%u',Scaffolds(ScaffoldL)))
    %find min and max
    MinPos=min(CurrScaffold.scaffoldPos(:,1));
    MaxPos=max(CurrScaffold.scaffoldPos(:,1));
    %plot queries
    QueryRank=[];
    for QueryL=1:size(CurrScaffold.queryPos,1)
        LinePos=CurrScaffold.queryRank(QueryL);
        if isempty(find(QueryRank==LinePos))
            QueryRank(end+1)=LinePos;
        end
        if CurrScaffold.scaffoldStrand(QueryL)==1
            line([CurrScaffold.scaffoldPos(QueryL,1)-MinPos+1,CurrScaffold.scaffoldPos(QueryL,2)-MinPos+1],[LinePos,LinePos],'color','b','linewidth',3)
        else
            line([CurrScaffold.scaffoldPos(QueryL,1)-MinPos+1,CurrScaffold.scaffoldPos(QueryL,2)-MinPos+1],[LinePos,LinePos],'color','r','linewidth',3)
        end
    end
    YLabel=Queries(QueryRank);
    set(gca,'ytick',QueryRank)
    set(gca,'yticklabel',YLabel)
    set(gca,'box','on')
    XLabel=str2num(get(gca,'xticklabel'))+MinPos;
    set(gca,'xticklabel',XLabel)
    YLim=get(gca,'ylim');
    YLim(2)=YLim(2)+1;
    set(gca,'ylim',YLim);
    title(sprintf('scaffold%u',Scaffolds(ScaffoldL)))


    h=figure;
    set(gcf,'color',[1,1,1])
    set(h,'name',sprintf('scaffold%u',Scaffolds(ScaffoldL)))
    [temp,SortIndex]=sort(CurrScaffold.scaffoldPos(:,1));
    ScaffoldPos=CurrScaffold.scaffoldPos(SortIndex,:);
    ScaffoldStrand=CurrScaffold.scaffoldStrand(SortIndex);
    for QueryL=1:size(CurrScaffold.queryPos,1)
        if ScaffoldStrand(QueryL)==1
            line([ScaffoldPos(QueryL,1)-MinPos+1,ScaffoldPos(QueryL,2)-MinPos+1],[QueryL,QueryL],'color','b','linewidth',3)
        else
            line([ScaffoldPos(QueryL,1)-MinPos+1,ScaffoldPos(QueryL,2)-MinPos+1],[QueryL,QueryL],'color','r','linewidth',3)
        end
    end
    set(gca,'box','on')
    set(gca,'xticklabel',XLabel)
    YLim=get(gca,'ylim');
    YLim(2)=YLim(2)+1;
    set(gca,'ylim',YLim);
    title(sprintf('scaffold%u',Scaffolds(ScaffoldL)))
end
