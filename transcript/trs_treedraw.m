function trs_treedraw(DataType,varargin)
global K P

cd(P.dir.data)
if exist('Tree.mat','file')
    load Tree
else
    errordlg('run prepare data for tree first')
end

BiolCondFlag=0;
WriteCoupleFlag=0;
if isequal(DataType,'points')
    if nargin==1
        TRank=1;
    else
        TRank=varargin{1};
    end

    %do not allow to select biol cond if Test Algo
    if P.flag.testAlgo==0
        Answer=questdlg('Do you want to select a subset of biological condition for drawing dendrogram ?','trs_treedraw','Yes','No','No');
        if isequal(Answer,'Yes')
            BiolCondFlag=1;
        end
    end

    Answer=questdlg('Do you want to display replicate information ?','trs_treedraw','Yes','No','Yes');
    if isequal(Answer,'Yes')
        WriteCoupleFlag=1;
    end
elseif  isequal(DataType,'biological conditions')
    if nargin==1
        TRank=2;
    else
        TRank=varargin{1};
    end
end

SavePngFlag=0;
Answer=questdlg('Do you want to save image of tree structure ?','trs_treedraw','Yes','No','No');
if isequal(Answer,'Yes')
    SavePngFlag=1;
end


Continue=0;
while Continue==0
    Answer=inputdlg('Enter the nb of leaf in tree figures','Input for display parameters',1,{'50'});
    TreeSize=str2num(Answer{1});
    if ~isempty(TreeSize)
        if TreeSize>1
            Continue=1;
        else
            warndlg('Nb of leef >=2 ','trs_treedraw')
        end
    else
        warndlg('Nb of leef >=2 ','trs_treedraw')
    end
end


Answer=questdlg('Do you want identical range for all tree figures (Yes) or adapted range for each (No)?','trs_treedraw','Yes','No','Yes');
if isequal(Answer,'Yes')
    Answer=inputdlg('give superior range limit (e.g ''200'' or ''>'' for full)','');
    if isequal(Answer{1},'>')
        XFullFlag=1;
    else
        XFullFlag=str2num(Answer{1});
    end
else
    XFullFlag=0;
end

if P.flag.testAlgo
    AlgoList=unique(P.point.algo);
    TreeNb=length(AlgoList)+1;
else
    TreeNb=1;
end
for TreeL=1:TreeNb
    %recover tree parameters
    if isequal(DataType,'points')
        ItemIndex=T{TRank}.pointIndex;
    else
        ItemIndex=T{TRank}.biolIndex;        
    end
    DistType=T{TRank}.distType;
    SelType=T{TRank}.selType;
    if isfield(T{TRank},'treeType')
        TreeType=T{TRank}.treeType;
    else
        TreeType='';
    end
    if isequal(SelType,'fdr')
    Pv=T{TRank}.parPv;
    Fdr=T{TRank}.parFdr;
    Sensitivity=T{TRank}.parS;
    end
    if isfield(T{TRank},'parTopSize')
        TopSize=T{TRank}.parTopSize;
    end
    MemDistances=squareform(T{TRank}.distances);
    %add empty rows for not used points
    Distances=zeros(P.point.nb);
    Distances(T{TRank}.pointIndex,T{TRank}.pointIndex)=MemDistances;
    clear MemDistances;
    LeafNames=T{TRank}.leafNames;
    LeafNb=length(LeafNames);

    % search couples of samples
    
    LeafName1=LeafNames;
    %write information on the pair of points that are the more similar in each
    %biological condition
    if WriteCoupleFlag==1
        BiolCondNb=P.biol.nb;
        for BiolCondL=1:BiolCondNb     
            if P.biol.used(BiolCondL)
            %keep only used points
            CurrItemIndex=P.biol.pointIndex{BiolCondL};
            UsedPoint=P.point.used(P.biol.pointIndex{BiolCondL});
            CurrItemIndex=CurrItemIndex(UsedPoint==1);
            PointNb=length(CurrItemIndex);
            if PointNb>1
                %recover the matrix of distances corresponding to the current biol condition               
                CurrDistances=Distances(CurrItemIndex,CurrItemIndex);               
                % replace 0 by Inf to find the real minimum (first diagonal == 0)
                if length(CurrItemIndex)>1
                    ZeroIndex=find(CurrDistances==0);
                    CurrDistances(ZeroIndex)=Inf;
                    [MinC CIndex]=min(CurrDistances);
                    [MinL LIndex]=min(MinC);
                    MinC=CurrItemIndex(CIndex(LIndex));
                    MinL=CurrItemIndex(LIndex);

                    %indicate each member of a couple by prefixing with
                    %'CondNb=> '
                    NameC=sprintf('%s_E%02u_C%03u_P%03u',P.point.name{MinC},P.point.expRank(MinC),P.point.biolRank(MinC),MinC);
                    Pos=strmatch(NameC,LeafNames,'exact');
                    if ~isempty(Pos)
                        if MinC~=MinL
                            if isfield(P.biol,'pairs')
                                if isempty(P.biol.pairs{BiolCondL})
                                    %some bolcond have been eliminated after a first round of automatic assignation
                                    LeafName1{Pos,1}=['C',sprintf('%02.f',BiolCondL),'>-< ',LeafNames{Pos}];
                                else
                                    if ~isempty(find(P.biol.pairs{BiolCondL}==MinC))
                                        LeafName1{Pos,1}=['C',sprintf('%02.f',BiolCondL),'>>>> ',LeafNames{Pos}];
                                    else
                                        %something wrong : currently detected couple does not correspond to registered one
                                        LeafName1{Pos,1}=['C',sprintf('%02.f',BiolCondL),'>?< ',LeafNames{Pos}];
                                    end
                                end
                            else
                                LeafName1{Pos,1}=['C',sprintf('%02.f',BiolCondL),'>>>> ',LeafNames{Pos}];
                            end
                        else
                            % identical samples
                            LeafName1{Pos,1}=['C',sprintf('%02.f',BiolCondL),'>=< ',LeafNames{Pos}];
                        end
                    end

                    NameL=[P.point.name{MinL},'_E',sprintf('%02.f',P.point.expRank(MinL)),'_C',sprintf('%03.f',P.point.biolRank(MinL)),'_P',sprintf('%03.f',MinL)];
                    Pos=strmatch(NameL,LeafNames,'exact');
                    if ~isempty(Pos)
                        if MinC~=MinL
                            if isfield(P.biol,'pairs')
                                if isempty(P.biol.pairs{BiolCondL})
                                    LeafName1{Pos,1}=['C',sprintf('%02.f',BiolCondL),'>-< ',LeafNames{Pos}];
                                else
                                    if ~isempty(find(P.biol.pairs{BiolCondL}==MinC))
                                        LeafName1{Pos,1}=['C',sprintf('%02.f',BiolCondL),'>>> ',LeafNames{Pos}];
                                    else
                                        LeafName1{Pos,1}=['C',sprintf('%02.f',BiolCondL),'>?< ',LeafNames{Pos}];
                                    end
                                end
                            else
                                LeafName1{Pos,1}=['C',sprintf('%02.f',BiolCondL),'>>> ',LeafNames{Pos}];
                            end
                        else
                            LeafName1{Pos,1}=['C',sprintf('%02.f',BiolCondL),'>=< ',LeafNames{Pos}];
                        end
                    end
                end
            elseif PointNb==1
                Name=[P.point.name{CurrItemIndex},'_E',sprintf('%02.f',P.point.expRank(CurrItemIndex)),'_C',sprintf('%03.f',P.point.biolRank(CurrItemIndex)),'_P',sprintf('%03.f',CurrItemIndex)];
                Pos=strmatch(Name,LeafNames,'exact');
                if ~isempty(Pos)
                    LeafName1{Pos,1}=['C',sprintf('%02.f',BiolCondL),' @ ',LeafNames{Pos}];
                end
            end
            end
        end
    end
    
    %add algorithm information if necessary
    if length(unique(P.point.algo))>1
        for LeafL=1:LeafNb
            LeafName1{LeafL}=sprintf('%s_%s',LeafName1{LeafL},P.point.algo{ItemIndex(LeafL)});
        end
    end
    %add usedPoint information if necessary
    if ~isempty(find(P.point.used==0))
        for LeafL=1:LeafNb
            LeafName1{LeafL}=sprintf('%s_U%u',LeafName1{LeafL},P.point.used(ItemIndex(LeafL)));
        end
    end



    if BiolCondFlag
        %construct the list of existing biological condition for the current chip
        BiolIndex=[];
        for BiolL=1:P.biol.nb
            if ~isempty(P.biol.pointIndex{BiolL})
                BiolIndex=[BiolIndex;BiolL];
            end
        end
        %select the biological conditions to be displayed
        BiolName={};
        for BiolL=1:length(BiolIndex)
            BiolName{BiolL,1}=sprintf('%04u - %s',BiolIndex(BiolL),P.biol.name{BiolIndex(BiolL)});
        end
        BiolSelType=questdlg('How do you want select biological conditions ?)','','on names','on ranks','on names');
        if isequal(BiolSelType,'on names')
            [BiolSelect OK]=listdlg('liststring',BiolName,'selectionmode','multiple','listsize',[300,900]);
            BiolSelect=BiolIndex(BiolSelect);
        else
            BiolSelect=inputdlg({'enter the biol cond ranks [1,3,...]'},'',1,{'[1,]'});
            BiolSelect=eval(BiolSelect{1});
            BiolSelect=BiolSelect';
        end
        BiolNb=length(BiolSelect);
        ClearIndex=[];
        for BiolL=1:P.biol.nb
            if isempty(find(BiolSelect==BiolL))
                if ~isempty(P.biol.pointIndex{BiolL})
                    for PointL=1:length(P.biol.pointIndex{BiolL})
                        PointPos=find(ItemIndex==P.biol.pointIndex{BiolL}(PointL));
                        if ~isempty(PointPos)
                            ClearIndex=[ClearIndex;PointPos];
                        end
                    end
                end
            end
        end
        if ~isempty(ClearIndex)
            LeafName1(ClearIndex)=[];
            LeafNb=length(LeafName1);
            if LeafNb<3
                h=errordlg('need at least 3 points !');
                waitfor(h)
                error('process canceled')
            else
                Distances(ClearIndex,:)=[];
                Distances(:,ClearIndex)=[];
                Distances(Distances==0)=eps;
                Distances=triu(Distances,1);
                Distances=Distances';
                Distances=Distances(:);
                Distances(Distances==0)=[];
                Linkage=linkage(Distances');
            end
        end
    end

    if TreeL>1
        CurrItemIndex=[];
        CurrPoints=strmatch(AlgoList{TreeL-1},P.point.algo,'exact');

        for PointL=1:length(CurrPoints)
            PointPos=find(ItemIndex==CurrPoints(PointL));
            if ~isempty(PointPos)
                CurrItemIndex=[CurrItemIndex;PointPos];
            end
        end
        ClearIndex=setdiff(ItemIndex,CurrItemIndex);
        if ~isempty(ClearIndex)
            LeafName1(ClearIndex)=[];
            LeafNb=length(LeafName1);
            if LeafNb<3
                h=errordlg('need at least 3 points !');
                waitfor(h)
                error('process canceled')
            else
                CurrDistances=Distances(CurrItemIndex,CurrItemIndex);
                CurrDistances(CurrDistances==0)=eps;
                CurrDistances=triu(CurrDistances,1);
                CurrDistances=CurrDistances';
                CurrDistances=CurrDistances(:);
                CurrDistances(CurrDistances==0)=[];
                Linkage=linkage(CurrDistances');
            end
        end
    end

    hfig=figure;
    set(hfig,'color',[1,1,1])

    if BiolCondFlag&&~isempty(ClearIndex)
        [h,Tree,v]=dendrogram(Linkage,0);
    elseif TreeL>1
        [h,Tree,v]=dendrogram(Linkage,0);
    else
        [h,Tree,v]=dendrogram(T{TRank}.linkage,0);
    end

    set(h,'linewidth',2)
    set(h,'color','k')
    XLimit=get(gca,'xlim');
    if LeafNb>=TreeSize
        set(gca,'ytick',[1:TreeSize]);
        set(gca,'xlim',XLimit);
        set(gca,'ylim',[0,TreeSize+1]);
        set(gca,'yticklabel',LeafName1(v((1:TreeSize))));
        set(gca,'position',[0.50 0.11 0.450 0.815])
    else
        set(gca,'ytick',[1:LeafNb]);
        set(gca,'ylim',XLimit);
        set(gca,'ylim',[0,LeafNb+1]);
        set(gca,'yticklabel',LeafName1(v));
        set(gca,'position',[0.50 0.11 0.450 0.815])
    end


    if isequal(SelType,'fdr')
        if isempty(TreeType)
             Title=sprintf('%s - %s - %s (%s)\nSelection at fdr<=%.2f,pv<=%.2f,s<=%.2f - Distances between %s',T{TRank}.distType,T{TRank}.selType,strrep(K.chip.name{P.chip.chipSetRank},'_','\_'),date,Fdr,Pv,Sensitivity,DataType);
        else
            Title=sprintf('%s - %s - %s (%s)\nSelection at fdr<=%.2f,pv<=%.2f,s<=%.2f - %u probe sets - Distances between %s (%s)',T{TRank}.distType,T{TRank}.selType,strrep(K.chip.name{P.chip.chipSetRank},'_','\_'),date,Fdr,Pv,Sensitivity,length(T{TRank}.psIndex),DataType,TreeType);
        end
    else
        if ~isempty(TreeType)
             Title=sprintf('%s - %s - %s (%s)\n%u probe sets - Distances between %s (%s)',T{TRank}.distType,T{TRank}.selType,strrep(K.chip.name{P.chip.chipSetRank},'_','\_'),date,length(T{TRank}.psIndex),DataType,TreeType);
        else
            if isequal(SelType,'point distances')
                Title=sprintf('%s - %s - %s (%s)\nfdr<=%.2f,pv<=%.2f,s<=%.2f - %u probe sets - Distances between %s (%s)',T{TRank}.distType,T{TRank}.selType,strrep(K.chip.name{P.chip.chipSetRank},'_','\_'),date,Fdr,Pv,Sensitivity,length(T{1}.psIndex),DataType,T{1}.treeType);
            else
                Title=sprintf('%s - %s - %s (%s)\nfirst %u of top lists - Distances between %s',T{TRank}.distType,T{TRank}.selType,strrep(K.chip.name{P.chip.chipSetRank},'_','\_'),date,DataType);
            end
        end
    end
    if TreeL>1
        Title=sprintf('%s_%s',AlgoList{TreeL-1},Title);
    end
    title(Title)

    Units=get(0,'units');
    set(0,'units','normalized')
    scrsz = get(0,'ScreenSize');
    eval(['set(0,''units'',''',Units,''')'])
    h1=get(gca,'parent');
    set(h1,'units','normalized')
    set(h1,'Position',[0 0 scrsz(3)  scrsz(4)*0.93])
    cd(P.dir.resTree)
    if XFullFlag==1
        set(gca,'xlim',XLimit);
    elseif XFullFlag==0
        set(gca,'xlimmode','auto');
    else
        XLimit=[XLimit(1),XFullFlag];
        set(gca,'xlim',XLimit);
    end
    if SavePngFlag==1
        if TreeL==1
            if isempty(TreeType)
                saveas(h1,sprintf('TS01_%s_%s_%s',DistType,SelType,date),'png')
            else
                saveas(h1,sprintf('TS01_%s_%s_%s_%s',DistType,SelType,TreeType,date),'png')
            end
        else
            saveas(h1,sprintf('TS01_top%uPs_%s_%sb',TopSize*2,AlgoList{TreeL-1},date),'png')
        end
    end
    DisplayNb=ceil(LeafNb/TreeSize);
    if DisplayNb>=2
        %     if DisplayNb>2
        %for DisplayLoop=[2:DisplayNb-1]
        Continue=1;
        Rank=1;
        while Continue==1
            Rank=inputdlg(sprintf('give the rank of sub tree to be displayed (>0 & <%u) or cancel (<=0 or >=%u)',DisplayNb+1,DisplayNb+1),'',1,{sprintf('%u',Rank+1)});
            if ~isempty(Rank)
                Rank=str2num(Rank{1});
            else
                Rank=0
            end
            if Rank<1|Rank>DisplayNb
                Answer=questdlg('do you want to cancel','','Yes','No','No');
                if isequal(Answer,'Yes')
                    Continue=0;
                end
            else
                if Rank<DisplayNb
                    set(gca,'ytick',[(Rank-1)*TreeSize+1:(Rank*TreeSize)]);
                else
                    set(gca,'ytick',[(DisplayNb-1)*TreeSize+1:LeafNb]);
                end
                if XFullFlag==0
                    set(gca,'xlimmode','auto');
                else
                    set(gca,'xlim',XLimit);
                end
                if Rank<DisplayNb
                    set(gca,'ylim',[(Rank-1)*TreeSize,(Rank*TreeSize)+1]);
                    set(gca,'yticklabel',LeafName1(v((Rank-1)*TreeSize+1:Rank*TreeSize)));
                else
                    set(gca,'ylim',[(DisplayNb-1)*TreeSize,(DisplayNb*TreeSize)+1]);
                    set(gca,'yticklabel',LeafName1(v((DisplayNb-1)*TreeSize+1:LeafNb)));
                end
                set(gca,'position',[0.50 0.11 0.450 0.815])

                if isequal(SelType,'fdr')
                    if isempty(TreeType)
                        Title=sprintf('%s - %s - %s (%s) - tree %u\nSelection at fdr<=%.2f,pv<=%.2f,s<=%.2f - Distances between %s',T{TRank}.distType,T{TRank}.selType,strrep(K.chip.name{P.chip.chipSetRank},'_','\_'),date,Rank,Fdr,Pv,Sensitivity,DataType);
                    else
                        Title=sprintf('%s - %s - %s (%s) - tree %u\nSelection at fdr<=%.2f,pv<=%.2f,s<=%.2f - %u probe sets - Distances between %s (%s)',T{TRank}.distType,T{TRank}.selType,strrep(K.chip.name{P.chip.chipSetRank},'_','\_'),date,Rank,Fdr,Pv,Sensitivity,length(T{TRank}.psIndex),DataType,TreeType);
                    end
                else
                    if ~isempty(TreeType)
                        Title=sprintf('%s - %s - %s (%s) - tree %u\n%u probe sets - Distances between %s (%s)',T{TRank}.distType,T{TRank}.selType,strrep(K.chip.name{P.chip.chipSetRank},'_','\_'),date,Rank,length(T{TRank}.psIndex),DataType,TreeType);
                    else
                        if isequal(SelType,'point distances')
                            Title=sprintf('%s - %s - %s (%s) - tree %u\nfdr<=%.2f,pv<=%.2f,s<=%.2f - %u probe sets - Distances between %s (%s)',T{TRank}.distType,T{TRank}.selType,strrep(K.chip.name{P.chip.chipSetRank},'_','\_'),date,Rank,Fdr,Pv,Sensitivity,length(T{1}.psIndex),DataType,T{1}.treeType);
                        else
                            Title=sprintf('%s - %s - %s (%s) - tree %u\nfirst %u of top lists - Distances between %s',T{TRank}.distType,T{TRank}.selType,strrep(K.chip.name{P.chip.chipSetRank},'_','\_'),date,Rank,DataType);
                        end
                    end
                end
                if TreeL>1
                    Title=sprintf('%s_%s_tree %u',AlgoList{TreeL-1},Title,Rank);
                end
                title(Title)
                %title([strrep(P.chip.name,'_',' '),' (',date,'): tree nb ',sprintf('%.f',Rank),' - Selection at p<=',sprintf('%01.3f',Pv),' - ',sprintf('%.f',length(SelIndex)),' genes - Distances between points (1 - Similarity) - Chip',CurrChip])
                h1=get(gca,'parent');
                set(h1,'units','normalized')
                set(h1,'Position',[0 0 scrsz(3)  scrsz(4)*0.93])
                if SavePngFlag==1
                    if TreeL==1
                        if isempty(TreeType)
                            saveas(h1,sprintf('TS%02u_%s_%s_%s',Rank,DistType,SelType,date),'png')
                        else
                            saveas(h1,sprintf('TS%02u_%s_%s_%s_%s',Rank,DistType,SelType,TreeType,date),'png')
                        end
                    else
                        saveas(h1,sprintf('TS%02u_top%uPs_%s_%sb',Rank,TopSize*2,AlgoList{TreeL-1},date),'png')
                    end
                end
            end
        end
    end
end


