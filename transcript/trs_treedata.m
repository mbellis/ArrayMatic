%calculate distances betweeen points or between biological conditions
%distances are calculated either by measuring euclidian distances between
%points or biological conditions in a multi-dimentional space
%or by the overlap between lists of probe set selected by
%fdr/pv/sensitivity or counted from the top (list ordered on ZVar)

function trs_treedata(DataType,varargin)
global K P DataRanks

cd(P.dir.data)
if exist('Tree.mat','file')
    load Tree
else
    T=[];
end
% selection parameters
if nargin==0
    h=errordlg('prepare_treedata needs DataType parameter');
    waitfor(h)
    error('process canceled')
end
Continue=0;
Fdr=1;
Sensitivity=1;
Pv=1;
TopSize=0;
CvLimit=0;


while Continue==0

    %INDICATE THE DISTANCE TYPE AND HOW ARE SELECTED THE POINTS OR THE TOP LIST
    if isequal(DataType,'points')|isequal(DataType,'biological conditions')
        if isequal(DataType,'points')
            if nargin==1
                TRank=1;
            else
                TRank=varargin{1};
            end
        elseif isequal(DataType,'biological conditions')
            if nargin==1
                TRank=2;
            else
                TRank=varargin{1};
            end
        end
        DistType=questdlg('Which type of distances do you want to use','','Euclidian distances','Overlap of top lists','Correlations','Euclidian distances');
        switch DistType
            case 'Euclidian distances'
                DistType='euclidian';
            case 'Correlations'
                DistType='correlation';
            case 'Overlap of top lists'
                DistType='overlap';
        end
        if isequal(DistType,'euclidian')|isequal(DistType,'correlation')
            if isequal(DataType,'points')
                SelType=questdlg('How do you want to select significative probe sets','','fdr,sensitivity,pv','rank coefficient of variation','fdr,sensitivity,pv');
            elseif isequal(DataType,'biological conditions')
                Continue=1;
                if ~isempty(T)
                    if ~isempty(T{1})
                        SelType='point';
                    else
                        h=errordlg('calculate tree data for points before');
                        waitfor(h)
                        error('process canceled')
                    end
                else
                    h=errordlg('calculate tree data for points before');
                    waitfor(h)
                    error('process canceled')
                end
            end
            if isequal(SelType,'fdr,sensitivity,pv')
                SelType='fdr';
                Prompt={'Enter the Fdr selection value:','Enter the Sensitivity selection value:','Enter the selection p-value:'};
                Def={'1','1','0.01'};
                DlgTitle='Input for selection parameters (<=Fdr&<=Sensitivity&<=Pv)';
            elseif isequal(SelType,'rank coefficient of variation')
                SelType='cv';
                Prompt={'Enter k in limit> mean(cv)+k*std(cv)'};
                Def={'1'};
                DlgTitle='Input for selection parameters';
            end
        elseif isequal(DistType,'overlap')
            SelType=questdlg('How do you want to select significative probe sets','','fdr,sensitivity,pv','ZVar top list size','fdr,sensitivity,pv');
            if isequal(SelType,'fdr,sensitivity,pv')
                SelType='fdr';
                Prompt={'Enter the Fdr selection value:','Enter the Sensitivity selection value:','Enter the selection p-value:'};
                Def={'1','1','0.01'};
                DlgTitle='Input for selection parameters (<=Fdr&<=Sensitivity&<=Pv)';
            elseif isequal(SelType,'ZVar top list size')
                SelType='zvar';
                Prompt={'Enter the size of the top list'};
                Def={'500'};
                DlgTitle='Input for selection parameter (TopList(ZVar))';
            end
        end

    else
        errordlg(sprintf('Unknown DataType %s.',DataType))
    end
    if ~isequal(SelType,'point')
        LineNb=1;
        Answer=inputdlg(Prompt,DlgTitle,LineNb,Def);
        if isequal(SelType,'fdr')
            Fdr=str2double(Answer{1});
            Sensitivity=str2double(Answer{2});
            Pv=str2double(Answer{3});
            if ~isempty(Fdr)&&~isempty(Sensitivity)&&~isempty(Pv)
                if Fdr>=0&&Fdr<=1&&Sensitivity>=0&&Sensitivity<=1&&Pv>=0&&Pv<=1
                    if Fdr==1&&Sensitivity==1&&Pv==1
                        warndlg('at least one of Fdr, Sensitivity or P-value must be <1','trs_treedata')
                    else
                        Continue=1;
                    end
                else
                    warndlg('Fdr, Sensitivity and P-value must be >=0 and <=1','trs_treedata')
                end
            else
                warndlg('Fdr, Sensitivity and P-value must filled with numbers','trs_treedata')
            end
        elseif  isequal(SelType,'cv')
            CvLimit=str2num(Answer{1});
            if CvLimit>=0&CvLimit<=5
                Continue=1;
                CvLimit=mean(P.chip.cv)+str2num(Answer{1})*std(P.chip.cv);
            else
                warndlg('select k between 0 and 5')
            end
        elseif isequal(SelType,'zvar')
            TopSize=str2double(Answer{1});
            if ~isempty(TopSize)
                if TopSize>1
                    Continue=1;
                else
                    warndlg('Top list size must be > 1','trs_treedata')
                end
            else
                warndlg('Top list size must be a number','trs_treedata')
            end
        end
    end
end %while Continue=0

if ~isequal(SelType,'point')
    TreeType='';
    if isequal(DistType,'euclidian')|isequal(DistType,'correlation')
        if isequal(DataType,'points')
            Answer=questdlg('Do you want to calculate distance between comparisons (pv or zvar) or points (signal or rank) ?','trs_treedata','points','comparisons','points');
            if isequal(Answer,'comparisons')
                Answer=questdlg('Do you want to calculate distances by using zvar or pv?','trs_treedata','zvar','pv','zvar');
                if isequal(Answer,'zvar')
                    TreeType='zvar';
                else
                    TreeType='pv';
                end
            else
                Answer=questdlg('Do you want to calculate distances by using log2(signal) or rank?','trs_treedata','signal','rank','signal');
                if isequal(Answer,'signal')
                    TreeType='log2(signal)';
                else
                    TreeType='rank';
                end
            end
        else
            Answer=questdlg('Do you want to calculate distances by using zvar or pv?','trs_treedata','zvar','pv','zvar');
            if isequal(Answer,'zvar')
                TreeType='zvar';
            else
                TreeType='pv';
            end
        end
        T{TRank}.treeType=TreeType;
    end
end

T{TRank}.distType=DistType;
T{TRank}.selType=SelType;
if isequal(SelType,'fdr')
    T{TRank}.parFdr=Fdr;
    T{TRank}.parS=Sensitivity;
    T{TRank}.parPv=Pv;
elseif isequal(SelType,'zvar')
    T{TRank}.parTopSize=TopSize;
elseif isequal(SelType,'cv')
    T{TRank}.parCv=CvLimit;
end



PsIndex=[];
% Answer=questdlg('Do you want to restrict trees to a particular probe set list ?','trs_treedata','No','Yes','No');
% if isequal(Answer,'Yes')
%     PsIndex=SELECTGENE(sprintf('Select probe sets used to draw %s dendrogram',DataType));
%     if isempty(PsIndex)
%         errordlg('Gene selection must not be empty')
%     end
%     T{TRank}.psIndex=PsIndex;
% end

if isequal(DataType,'points')
    Answer=questdlg('Do you want to restrict trees to a particular point list ?','trs_treedata','No','Yes','No');
    if isequal(Answer,'Yes')
        PointIndex=select_points('multiple','Select experimental points used to draw dendrogram');
        if isempty(PointIndex)
            errordlg('Point selection must not be empty')
        end
    else
        PointIndex=1:P.point.nb;
    end
    T{TRank}.pointIndex=PointIndex;
elseif isequal(DataType,'biological conditions')
    if isequal(P.par.analType,'network')
        cd(P.dir.data)
        if exist('UsedBiolIndex.mat','file')
            load UsedBiolIndex
        else
            h=errordlg('Make biol median comparison first');
            waitfor(h)
            error('process canceled')
        end
    else
        Answer=questdlg('Do you want to restrict trees to a particular biological condition list ?','trs_treedata','No','Yes','No');
        if isequal(Answer,'Yes')
            BiolIndex=select_biolconditions('Select biological conditions used to draw dendrogram');
            if isempty(BiolIndexIndex)
                errordlg('Biological selection  selection must not be empty')
            end
        else
            BiolIndex=1:P.biol.nb;
        end
    end
    T{TRank}.biolIndex=BiolIndex;
end

if ~isequal(SelType,'point')
    % search genes wich satisfy the selection condition (fdr,pv,s) and construct SelIndex
    % or constructs T{TRank.}incIndex and T{TRank}.decIndex used for dendrograms constructed upon 1 - similarity

    % initialisation of variables
    SelBindex=zeros(P.chip.currProbeSetNb,1);
    KeepBindex=zeros(P.chip.currProbeSetNb,1);
    if ~isempty(PsIndex)
        KeepBindex(PsIndex)=1;
    end

    if isequal(SelType,'fdr')
        FdrVal=ones(P.chip.currProbeSetNb,1);
        SensitivityVal=ones(P.chip.currProbeSetNb,1);
        PvVal=ones(P.chip.currProbeSetNb,1);
    end
    if isequal(SelType,'zvar')
        if isequal(DataType,'points')
            Tree.incIndex=cell(P.point.nb,1);
            Tree.decIndex=cell(P.point.nb,1);
        elseif isequal(DataType,'biological conditions')
            Tree.incIndex=cell(length(BiolIndex),1);
            Tree.decIndex=cell(length(BiolIndex),1);
        end
    end
    if isequal(DataType,'points')
        TotalItemNb=P.point.nb;
        ItemIndex=PointIndex;
        ResRank=1;
    elseif isequal(DataType,'biological conditions')
        TotalItemNb=length(BiolIndex);
        ItemIndex=BiolIndex;
        if isequal(P.par.analType,'network')
            ResRank=2;
        elseif isequal(P.par.analType,'transcriptome')
            h=errordlg('not yet implemented');
            waitfor(h)
            error('process canceled')
        end
    end

    ItemNb=length(ItemIndex);
    for ItemL=1:ItemNb
        ItemL
        if isequal(DataType,'points')
            CurrItem=ItemIndex(ItemL);
        elseif isequal(DataType,'biological conditions')
            if isequal(P.par.analType,'network')
                CurrItem=ItemL;
            else
                CurrItem=ItemIndex(ItemL);
            end
        end
        %load data
        if isequal(SelType,'fdr')
            SelStr='';
            if Fdr<1
                FdrVal=load_data(sprintf('Fdr_%02u.float32le',ResRank),P.dir.data,P.chip.currProbeSetNb,TotalItemNb,'single','ieee-le',1:P.chip.currProbeSetNb,CurrItem);
                SelStr='FdrVal';
            end
            if Sensitivity<1
                SensitivityVal=load_data(sprintf('Sensitivity_%02u.float32le',ResRank),P.dir.data,P.chip.currProbeSetNb,TotalItemNb,'single','ieee-le',1:P.chip.currProbeSetNb,CurrItem);
                if isempty(SelStr)
                    SelStr='SensitivityVal';
                end
            end
            if Pv<1
                PvVal=load_data(sprintf('Pv_%02u.float32le',ResRank),P.dir.data,P.chip.currProbeSetNb,TotalItemNb,'single','ieee-le',1:P.chip.currProbeSetNb,CurrItem);
                if isempty(SelStr)
                    SelStr='PvVal';
                end
            end
            CurrIndex=find(abs(FdrVal)<=Fdr&abs(SensitivityVal)<=Sensitivity&abs(PvVal)<=Pv);

            if isequal(DistType,'overlap')
                if ~isempty(PsIndex)
                    SelIndex=intersect(CurrIndex,PsIndex);
                else
                    SelIndex=CurrIndex;
                end
                eval(sprintf('Tree.incIndex{ItemL}=SelIndex(find(%s(SelIndex>=0)));',SelStr));
                eval(sprintf('Tree.decIndex{ItemL}=SelIndex(find(%s(SelIndex<0)));',SelStr));
            else
                SelBindex(CurrIndex)=1;
            end

        elseif isequal(SelType,'cv')
            SelBindex=P.chip.cv>=CvLimit;
        elseif isequal(SelType,'zvar')
            ZVar=load_data(sprintf('ZVar_%02u.float32le',ResRank),P.dir.data,P.chip.currProbeSetNb,TotalItemNb,'single','ieee-le',1:P.chip.currProbeSetNb,CurrItem);
            [IncIndex,DecIndex]=topindexes(ZVar,TopSize,'max',1);
            if ~isempty(PsIndex)
                IncIndex=intersect(IncIndex,PsIndex);
                DecIndex=intersect(DecIndex,PsIndex);
            end
            Tree.incIndex{ItemL}=IncIndex;
            Tree.decIndex{ItemL}=DecIndex;
        end

        % SelIndex used for recover zvar, pv, signal or rank
        if ~isequal(DistType,'overlap')
            if ~isempty(PsIndex)
                SelBindex=SelBindex&KeepBindex;
            end
            SelIndex=find(SelBindex);
            if length(SelIndex)<2
                errordlg('SelIndex must have at least two elements')
            end
            T{TRank}.psIndex=SelIndex;
        end



        % Recover Ranks
        if ~isequal(DistType,'overlap')
            if isequal(DataType,'points')
                Rank=[];
                for ItemL=1:ItemNb
                    CurrItem=ItemIndex(ItemL);
                    if P.flag.loadData
                        CurrRank=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le',1:P.chip.currProbeSetNb,CurrItem);
                        CurrRank=CurrRank(SelIndex);
                    else
                        CurrRank=DataRanks(SelIndex,CurrItem);
                    end
                    Rank=[Rank,CurrRank];
                end

                % eliminates probe sets which are not measured in at least one point (based on Rank==-1,or -2)
                ClearPos=[];
                for LineL=1:size(Rank,1);
                    if ~isempty(find(Rank(LineL,:)<0))
                        ClearPos=[ClearPos;LineL];
                    end
                end
                if ~isempty(ClearPos)
                    %h=warndlg(sprintf('%u genes to be deleted (absent on some point) out of a total of %u',length(ClearPos),size(Rank,2)));
                    %waitfor(h);
                    Rank(ClearPos,:)=[];
                end

                if isequal(TreeType,'rank')
                    Values=Rank;
                    Rank=[];
                end


                if isequal(TreeType,'log2(signal)')
                    try
                        Values=log2(interp1(P.chip.refRank,P.chip.refSignal,Rank));
                    catch
                        Values=log2(interp1(K.chip.ref.rank,K.chip.ref.signal,Rank));
                    end
                    if ~isempty(find(isnan(Values)))
                        h=warndlg(sprintf('%u NaN values in Signal of chip %s!',length(find(isnan(Signal))),F.list.array));
                        waitfor(h)
                        h=errordlg('Process canceled');
                        waitfor(h)
                    end
                end
            end


            if isequal(TreeType,'zvar')|isequal(TreeType,'pv')
                ItemNb=length(ItemIndex);
                for ItemL=1:ItemNb
                    if isequal(DataType,'points')
                        CurrItem=ItemIndex(ItemL);
                    elseif isequal(DataType,'biological conditions')
                        if isequal(P.par.analType,'network')
                            CurrItem=ItemL;
                        else
                            CurrItem=ItemIndex(ItemL);
                        end
                    end
                    if isequal(TreeType,'zvar')
                        CurrVal=load_data(sprintf('ZVar_%02u.float32le',ResRank),P.dir.data,P.chip.currProbeSetNb,ItemNb,'single','ieee-le',1:P.chip.currProbeSetNb,CurrItem);
                    else
                        CurrVal=load_data(sprintf('Pv_%02u.float32le',ResRank),P.dir.data,P.chip.currProbeSetNb,ItemNb,'single','ieee-le',1:P.chip.currProbeSetNb,CurrItem);
                    end
                    CurrVal=CurrVal(SelIndex);
                    Values=[Values,CurrVal];
                end
                Values(ClearPos,:)=[];
            end
        end
    end
end


if isequal(DataType,'points')
    LeafNames=P.point.name(PointIndex);
    LeafNb=length(LeafNames);
    for LeafL=1:LeafNb
        LeafNames{LeafL}=[LeafNames{LeafL},'_E',sprintf('%02.f',P.point.expRank(PointIndex(LeafL))),'_C',sprintf('%03.f',P.point.biolRank(PointIndex(LeafL))),'_P',sprintf('%03.f',PointIndex(LeafL))];
    end
elseif isequal(DataType,'biological conditions')
    LeafNames=P.biol.name(BiolIndex);
    LeafNb=length(LeafNames);
    for LeafL=1:LeafNb
        LeafNames{LeafL}=[LeafNames{LeafL},'_E',sprintf('%02.f',P.point.expRank(P.biol.pointIndex{BiolIndex(LeafL)}(1))),'_C',sprintf('%03.f',BiolIndex(LeafL))];
    end
end
T{TRank}.leafNames=LeafNames;

if isequal(DistType,'overlap')
    Y=single(zeros(1,ItemNb*(ItemNb-1)/2));
    tic
    t=0;
    YPos=0;
    for ItemL1=1:ItemNb-1
        for ItemL2=ItemL1+1:ItemNb
            YPos=YPos+1;
            if isequal(SelType,'zvar')
                Y(YPos)=1-(length(intersect(Tree.incIndex{ItemL1},Tree.incIndex{ItemL2}))+length(intersect(Tree.decIndex{ItemL1},Tree.decIndex{ItemL2})))/(2*T{TRank}.parTopSize);
            else
                Y(YPos)=1-(length(intersect(Tree.incIndex{ItemL1},Tree.incIndex{ItemL2}))+length(intersect(Tree.decIndex{ItemL1},Tree.decIndex{ItemL2})))/(sqrt(length(Tree.incIndex{ItemL1})*length(Tree.incIndex{ItemL2}))+min(length(Tree.decIndex{ItemL1}),length(Tree.decIndex{ItemL2})));
                %Y(YPos)=1-(length(intersect(Tree.incIndex{ItemL1},Tree.incIndex{ItemL2}))+length(intersect(Tree.decIndex{ItemL1},Tree.decIndex{ItemL2})))/(sqrt(length(Tree.incIndex{ItemL1})*length(Tree.incIndex{ItemL2}))+sqrt(length(Tree.decIndex{ItemL1})*length(Tree.decIndex{ItemL2})));
            end
        end
        t=t+toc;
        sprintf('%u items processed in %u min: estimate that process will end in %u min',ItemL1,round(t/60),round((ItemNb-1-ItemL1)*t/(ItemL1*60)))
        tic
    end
elseif isequal(DistType,'Correlation')
    Y=single(zeros(1,ItemNb*(ItemNb-1)/2));
    tic
    t=0;
    YPos=0;
    for ItemL1=1:ItemNb-1
        for ItemL2=ItemL1+1:ItemNb
            CurrCorr=corrcoef(Values(:,ItemL1),Values(:,ItemL2));
            Y(YPos)=1-CurrCorr(1,2);
        end
    end
    t=t+toc;
    sprintf('%u items processed in %u min: estimate that process will end in %u min',ItemL1,round(t/60),round((ItemNb-1-ItemL1)*t/(ItemL1*60)))
    tic
elseif isequal(SelType,'point')
    ItemNb=length(BiolIndex);
    Y=single(zeros(1,ItemNb*(ItemNb-1)/2));
    YPos=0;
    PointDistances=squareform(T{1}.distances);
    for ItemL1=1:ItemNb-1
        PointRanks1=P.biol.pairs{BiolIndex(ItemL1)};
        for ItemL2=ItemL1+1:ItemNb
            PointRanks2=P.biol.pairs{BiolIndex(ItemL2)};
            YPos=YPos+1;
            Y(YPos)=mean([PointDistances(PointRanks1(1),PointRanks2(1)),...
                PointDistances(PointRanks1(1),PointRanks2(2)),...
                PointDistances(PointRanks1(2),PointRanks2(1)),...
                PointDistances(PointRanks1(2),PointRanks2(2))]);
        end
    end
else
    Y=sqrt(dist2(Values',Values'));
end

T{TRank}.distances=Y;
T{TRank}.linkage=linkage(Y);

h=figure;
set(gcf,'color',[1,1,1])
if ~isfield(T{TRank},'treeType')
    Title=sprintf('Distance distribution between %s - %s - %s',DataType,T{TRank}.distType,T{TRank}.selType);
    FileName=sprintf('%s_distances_hist_%s_%s_%s',strrep(DataType,' ','_'),T{TRank}.selType,P.project.name,date);
else
    Title=sprintf('Distance distribution between %s - %s - %s - %s ',DataType,T{TRank}.distType,T{TRank}.selType,T{TRank}.treeType);
    FileName=sprintf('%s_distances_hist_%s_%s_%s_%s',strrep(DataType,' ','_'),T{TRank}.selType,T{TRank}.treeType,P.project.name,date);
end
set(h,'name',Title)
hist(T{TRank}.distances,100);
title(Title)
xlabel(sprintf('distances between %s',DataType))
ylabel('frequency')
set_figsize('780px')
cd(P.dir.resTree)
saveas(h,FileName,'png')

h=figure;
set(gcf,'color',[1,1,1])
if ~isfield(T{TRank},'treeType')
    Title=sprintf('Heathmap of distances between %s - %s - %s',DataType,T{TRank}.distType,T{TRank}.selType);
    FileName=sprintf('%s_distances_heathmap_%s_%s_%s',strrep(DataType,' ','_'),T{TRank}.selType,P.project.name,date);
else
    Title=sprintf('Heathmap of distribution between %s - %s - %s - %s ',DataType,T{TRank}.distType,T{TRank}.selType,T{TRank}.treeType);
    FileName=sprintf('%s_distances_heathmap_%s_%s_%s_%s',strrep(DataType,' ','_'),T{TRank}.selType,T{TRank}.treeType,P.project.name,date);
end
set(h,'name',Title)
if max(T{TRank}.distances)<1
    image(squareform(T{TRank}.distances)*70)
else
    image(squareform(T{TRank}.distances)/5)
end
title(Title)
xlabel(sprintf('distances between %s',DataType))
set_figsize('960px')
cd(P.dir.resTree)
saveas(h,FileName,'png')




cd(P.dir.data)
save Tree T
T=[];

% h=figure;
% set(gcf,'color',[1,1,1])
% set(h,'name','BIOL COND DISTANCES COMPARISON')
% Labels={'points','fdr sqrt','fdr min','top list'};
% Pos=0;
% TList=[2,3,5];
% for T1=[1,2]
%     for T2=T1+1:3
%         Pos=Pos+1;
%         subplot(1,3,Pos)
%         plot(T{TList(T1)}.distances,T{TList(T2)}.distances,'b.','markersize',3)
%         xlabel(Labels{TList(T1)-1})
%         ylabel(Labels{TList(T2)-1})
%     end
% end