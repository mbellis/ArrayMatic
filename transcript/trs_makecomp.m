function trs_makecomp(CompType,varargin)
global P S DataRanks
switch CompType
    case 'user groups'
        CompRank=1;
        ResRank=2;
        cd(P.dir.data)
        if exist('Comp.mat','file')
            load Comp
        else
            h=errodlg('make *** TRANSCRIPTION/USER GROUPS ANALYSIS/prepare group comparisons ***');
            waitfor(h)
            error('process canceled')
        end
        MakeIndex=find(M{CompRank}.make==1);
    case 'network'
        ResRank=2;
end
PsNb=P.chip.currProbeSetNb;

%% MAKE CALIBRATION SETS

%load calib set
cd(P.dir.data)
if exist(sprintf('CalibSet_%02u.mat',ResRank),'file')
    eval(sprintf('load CalibSet_%02u',ResRank))
end
%noise_distribution parameters
RankThreshold=[0,0];
ClearIndex=[];
NormType='quantile';
AnalyseType='transcriptome';
CalibType='idem';
SizerFittingDegree=7;
DisplayFlag=0;
ClearIndexFlag=0;
%set to zero not to save after each calcul in
%noise_distrihution
SaveFlag=0;

switch CompType

    case 'user groups'
        if length(MakeIndex)>0
            
            SNb=0;
            t=0;
            tic
            for MakeL=1:length(MakeIndex)
                SNb=SNb+1;
                %         waitbar(0,w,sprintf('Calculating Calibration Sets for the %uth comparison among %u',MakeL,length(MakeIndex)))
                %fill default comparison scheme
                Grp={};
                Grp{1}=M{CompRank}.firstPointRanks{MakeIndex(MakeL)};
                Grp{2}=M{CompRank}.sndPointRanks{MakeIndex(MakeL)};
                GrpNb(1)=length(Grp{1});
                GrpNb(2)=length(Grp{2});
                [MinNb,MinGrp]=min(GrpNb);
                MaxNb=max(GrpNb);
                if MinGrp==1
                    MaxGrp=2;
                else
                    MaxGrp=1;
                end
                %calculate point order and CompScheme only once
                MakeScheme=1;
                if isfield(M{CompRank},'firstOrder')
                    if length(M{CompRank}.firstOrder)>=MakeIndex(MakeL)
                        if ~isempty(M{CompRank}.firstOrder{MakeIndex(MakeL)})
                            MakeScheme=0;
                            if MinGrp==1
                                MinOrder=M{CompRank}.firstOrder{MakeIndex(MakeL),1};
                                MaxOrder=M{CompRank}.sndOrder{MakeIndex(MakeL),1};
                            else
                                MinOrder=M{CompRank}.sndOrder{MakeIndex(MakeL),1};
                                MaxOrder=M{CompRank}.firstOrder{MakeIndex(MakeL),1};
                            end
                            CompScheme=M{CompRank}.compScheme{MakeIndex(MakeL)};
                        end
                    end
                end
                if MakeScheme
                    if P.flag.testAlgo
                        %keep the same order for all identical comparisons but with
                        %different algorithms
                        CompName=regexp(M{CompRank}.compName{MakeIndex(MakeL)},'^\w+(?=_\w+$)','match');
                        CompPos=strmatch(CompName{1},M{CompRank}.compName);
                    end

                    Combinations=zeros(MinNb,MinNb);
                    MinOrder=randperm(MinNb);
                    Circle=[MinOrder,MinOrder];
                    MaxOrder=randperm(MaxNb);
                    if MinGrp==1
                        M{CompRank}.firstOrder{MakeIndex(MakeL),1}=MinOrder;
                        M{CompRank}.sndOrder{MakeIndex(MakeL),1}=MaxOrder;
                        if P.flag.testAlgo&length(CompPos)>1
                            for CompL=2:length(CompPos)
                                M{CompRank}.firstOrder{CompPos(CompL),1}=MinOrder;
                                M{CompRank}.sndOrder{CompPos(CompL),1}=MaxOrder;
                            end
                        end
                    else
                        M{CompRank}.firstOrder{MakeIndex(MakeL),1}=MaxOrder;
                        M{CompRank}.sndOrder{MakeIndex(MakeL),1}=MinOrder;
                        if P.flag.testAlgo&length(CompPos)>1
                            for CompL=2:length(CompPos)
                                M{CompRank}.firstOrder{CompPos(CompL),1}=MaxOrder;
                                M{CompRank}.sndOrder{CompPos(CompL),1}=MinOrder;
                            end
                        end
                    end
                    for CombL=1:MinNb
                        Combinations(CombL,:)=Circle(CombL:CombL+MinNb-1);
                    end
                    CompScheme={};
                    CurrPos=0;
                    while CurrPos<=GrpNb(MaxGrp)-GrpNb(MinGrp)
                        CurrPos=CurrPos+1;
                        for CombL=1:MinNb
                            if MinGrp==1
                                CompScheme{1,end+1}=[Combinations(CombL,:);CurrPos:CurrPos+MinNb-1];
                            else
                                CompScheme{1,end+1}=[CurrPos:CurrPos+MinNb-1;Combinations(CombL,:)];
                            end
                        end
                    end
                    M{CompRank}.compScheme{MakeIndex(MakeL),1}=CompScheme;
                    if P.flag.testAlgo&length(CompPos)>1
                        for CompL=2:length(CompPos)
                            M{CompRank}.compScheme{CompPos(CompL),1}=CompScheme;
                        end
                    end
                end

                %make calibration sets


                %maximum number of calib set to be done
                %         WaitSteps=sum(GrpNb);
                %         WaitStep=0;
                for GrpL=1:2
                    for RoundL=1:MaxNb-MinNb+1
                        for SetL=1:MinNb
                            if MinGrp==GrpL
                                %all MinNb in MinGrp are used, with random order
                                if SetL<MinNb
                                    CalibRanks=sort([Grp{GrpL}(MinOrder(SetL)),Grp{GrpL}(MinOrder(SetL+1))]);
                                else
                                    CalibRanks=sort([Grp{GrpL}(MinOrder(SetL)),Grp{GrpL}(MinOrder(1))]);
                                end
                            else
                                %MinNb in MaxGrp in random order starting at RoundL are used
                                if SetL+RoundL-1<MaxNb
                                    CalibRanks=sort([Grp{GrpL}(MaxOrder(SetL+RoundL-1)),Grp{GrpL}(MaxOrder(SetL+RoundL))]);
                                else
                                    CalibRanks=sort([Grp{GrpL}(MaxOrder(SetL+RoundL-1)),Grp{GrpL}(MaxOrder(RoundL))]);
                                end
                            end
                            %do calibration if necessary
                            DoS=1;

                            
                                if length(S)>=3
                                    if ~isempty(S{3})
                                        if isfield(S{3},'testSurf')
                                            Pos=find(S{3}.position(:,1)==CalibRanks(1)&S{3}.position(:,2)==CalibRanks(2));
                                            if ~isempty(Pos)
                                                if ~isempty(S{3}.testSurf(Pos))
                                                    DoS=0;
                                                end                                                
                                            end
                                        end
                                    end
                                end
                            

                            if DoS
                                if P.flag.loadData==0
                                    HL=DataRanks(:,CalibRanks(1));
                                else
                                    HL=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le',1:P.chip.currProbeSetNb,CalibRanks(1));
                                end
                                if P.flag.loadData==0
                                    BL=DataRanks(:,CalibRanks(2));
                                else
                                    BL=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le',1:P.chip.currProbeSetNb,CalibRanks(2));
                                end

                                SaveCalibSet=1;
                                [RankGrid,Grid,ZVarGrid,ZVar]=noise_distribution(HL,BL,RankThreshold,CalibRanks(1),CalibRanks(2),ClearIndex,NormType,AnalyseType,CalibType,SizerFittingDegree,DisplayFlag);
                                fill_s(NormType,CalibType,ClearIndexFlag,CalibRanks,Grid,ZVar,ZVarGrid,ResRank,SaveFlag)
                            end
                        end
                    end
                end
                if SNb==100
                    t=t+toc;
                    cd(P.dir.data)
                    eval(sprintf('save CalibSet_%02u S',ResRank))
                    SNb=0;
                    sprintf('I evaluate that the process will end in around %u min',round(t/60*(length(MakeIndex)-MakeL)/MakeL))
                    tic
                end
            end
            if SNb>1
                cd(P.dir.data)                
                eval(sprintf('save CalibSet_%02u S',ResRank))
            end
            cd(P.dir.data)
            save Comp M
        end


    case 'network'
        SNb=0;
        %maximum number of calib set to be done
        t=0;
        tic
        SaveIt=0;
        for MakeL=1:length(P.biol.pairs)
            SNb=SNb+1;
            %do as if all biological conditions are used even if
            %P.biol.used is set to 0 or if P.biol.pairs=0, whatever is
            %the reason (e.g. no replicate).
            %If in a future update of the data, a not used biological
            %condition becames usable, its place in S then already
            %reserved.
            if ~isempty(P.biol.pairs{MakeL})
                CalibRanks=sort(P.biol.pairs{MakeL});
                %do calibration if necessary
                DoS=1;
                if length(S)>=3
                    if ~isempty(S{3})
                        if isfield(S{3},'testSurf')
                            Pos=find(S{3}.position(:,1)==CalibRanks(1)&S{3}.position(:,2)==CalibRanks(2));
                            if ~isempty(Pos)
                                if ~isempty(S{3}.testSurf(Pos))
                                    DoS=0;
                                end
                            end
                        end
                    end
                end


                if DoS
                    SaveIt=1;
                    HL=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le',1:P.chip.currProbeSetNb,CalibRanks(1));
                    BL=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le',1:P.chip.currProbeSetNb,CalibRanks(2));
                    [RankGrid,Grid,ZVarGrid,ZVar]=noise_distribution(HL,BL,RankThreshold,CalibRanks(1),CalibRanks(2),ClearIndex,NormType,AnalyseType,CalibType,SizerFittingDegree,DisplayFlag);
                    %prevent analysis of identical data
                    fill_s(NormType,CalibType,ClearIndexFlag,CalibRanks,Grid,ZVar,ZVarGrid,ResRank,SaveFlag)
                end
            end
            if SNb==100
                t=t+toc;
                cd(P.dir.data)
                if SaveIt
                    eval(sprintf('save CalibSet_%02u S',ResRank))
                    SaveIt=0;
                end
                SNb=0;
                sprintf('%u biol cond processed, %u min elapsed: I evaluate that the process will end in around %u min',MakeL,round(t/60),round(t*(length(P.biol.pairs)-MakeL)/(MakeL*60)))
                tic
                
            end
        end
        if SNb>1&SaveIt
            cd(P.dir.data)
            eval(sprintf('save CalibSet_%02u S',ResRank))
        end
end % of switch


%% MAKE COMPARISONS
switch CompType
    case 'user groups'
        if length(MakeIndex)>0
            %load existing comparison
            ResNb=0;
            if isfield(M{CompRank},'column')
                for ResL=1:length(M{CompRank})
                    if M{CompRank}.made(ResL)
                        ResNb=max(ResNb,M{CompRank}.column(ResL));
                    end
                end
            end
            if P.flag.loadData==0
                if ResNb==0
                    ZVar=[];
                    Fdr=[];
                    Sensitivity=[];
                    Pv=[];
                    FC=[];
                else
                    [ZVar,Fdr,Sensitivity,Pv,FC,TotalVar]=loadsave_comp('load',ResRank,0,P.dir.data,1,1,1,1,1,P.chip.currProbeSetNb,ResNb,0);
                end
            end
            %rdam parameters
            CalibType='idem';
            NormType='quantile';
            AnalyseType='transcriptome';
            %DataRanks are filled here if P.flag.loadData=1, so do not load data again in rdam;
            LoadDataFlag=0;
            SingleCalibPointFlag=0;
            SingleCalibCurveFlag=0;
            CalibUpdateSFlag=0;
            CalibSaveFlag=1;
            DisplayFlag=0;
            SizerFittingDegree=7;
            ComparisonFlag=1;
            CalibSchemeFlag=0;
            cd(P.dir.data)
            for MakeL=1:length(MakeIndex)
                Grp={};
                Grp{1}=M{CompRank}.firstPointRanks{MakeIndex(MakeL)}(M{CompRank}.firstOrder{MakeIndex(MakeL)});
                Grp{2}=M{CompRank}.sndPointRanks{MakeIndex(MakeL)}(M{CompRank}.sndOrder{MakeIndex(MakeL)});
                GrpNb(1)=length(Grp{1});
                GrpNb(2)=length(Grp{2});

                %calculate mean signal
                Signal{1}=[];
                Signal{2}=[];
                for GrpL=1:2
                    for PointL=1:GrpNb(GrpL)
                        Signal{GrpL}=[Signal{GrpL},load_data('DataSignals.float32le',P.dir.data,PsNb,P.point.nb,'single','ieee-le',1:PsNb,Grp{GrpL}(PointL))];
                    end
                    Signal{GrpL}=mean(Signal{GrpL},2);
                    Signal{GrpL}(Signal{GrpL}<=0)=0.01;
                end
                CurrFC=Signal{1}./Signal{2};

                CompScheme=M{CompRank}.compScheme{MakeIndex(MakeL)};
                if P.flag.loadData==0
                    TGRankList=Grp{1};
                    CGRankList=Grp{2};
                    eval(sprintf('[CurrZVar,CurrPv,CurrPpv,CurrFdr,CurrSensitivity,CurrTotalVar]=rdam(CompScheme,TGRankList,CGRankList,%u,[0,0],''%s'',[],''%s'',''%s'',%u,%u,%u,%u,%u,%u,%u,%u,%u);',...
                        LoadDataFlag,CalibType,NormType,AnalyseType,SizerFittingDegree,SingleCalibPointFlag,...
                        SingleCalibCurveFlag,CalibUpdateSFlag,CalibSaveFlag,DisplayFlag,ComparisonFlag,...
                        ResRank,CalibSchemeFlag));

                else
                    %DataRanks is a global variable and is filled here only with
                    %used data ...
                    DataRanks=[];
                    for PointL=1:GrpNb(1)
                        DataRanks(:,PointL)=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le',1:P.chip.currProbeSetNb,Grp{1}(PointL));
                    end
                    for PointL=1:GrpNb(2)
                        DataRanks(:,PointL+GrpNb(1))=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le',1:P.chip.currProbeSetNb,Grp{2}(PointL));
                    end
                    %... so TG and CG RankLists are continuous and start at 1
                    TGRankList=1:GrpNb(1);
                    CGRankList=1+GrpNb(1):GrpNb(2)+GrpNb(1);
                    TGRankList4S=Grp{1};
                    CGRankList4S=Grp{2};
                    CalibSchemeFlag=2;
                    MakeL
                    eval(sprintf('[CurrZVar,CurrPv,CurrPpv,CurrFdr,CurrSensitivity,CurrTotalVar]=rdam(CompScheme,TGRankList,CGRankList,%u,[0,0],''%s'',[],''%s'',''%s'',%u,%u,%u,%u,%u,%u,%u,%u,%u,TGRankList4S,CGRankList4S);',...
                        LoadDataFlag,CalibType,NormType,AnalyseType,SizerFittingDegree,SingleCalibPointFlag,...
                        SingleCalibCurveFlag,CalibUpdateSFlag,CalibSaveFlag,DisplayFlag,ComparisonFlag,...
                        ResRank,CalibSchemeFlag));

                end
                if ~isempty(CurrZVar)
                    if isfield(M{CompRank},'column')
                        try
                            ResNb=M{CompRank}.column(MakeIndex(MakeL),1);
                        catch
                            ResNb=ResNb+1;
                        end
                    else
                        ResNb=ResNb+1;
                    end
                    %change FC
                    DecPos=find(CurrZVar<0);
                    CurrFC(DecPos)=-1./CurrFC(DecPos);
                    M{CompRank}.column(MakeIndex(MakeL),1)=ResNb;
                    if P.flag.loadData
                        loadsave_comp('save',ResRank,ResNb,P.dir.data,1,1,1,1,1,P.chip.currProbeSetNb,CurrZVar,CurrFdr,CurrSensitivity,CurrPv,CurrFC,CurrTotalVar);
                    else
                        ZVar(:,ResNb)=CurrZVar;
                        Sensitivity(:,ResNb)=CurrSensitivity;
                        Fdr(:,ResNb)=CurrFdr;
                        Pv(:,ResNb)=CurrPv;
                        FC(:,ResNb)=CurrFC;
                        TotalVar.inc(ResNb,1)=CurrTotalVar(1);
                        TotalVar.dec(ResNb,1)=CurrTotalVar(2);
                    end
                    M{1}.make(MakeIndex(MakeL))=0;
                    M{1}.made(MakeIndex(MakeL))=1;
                else
                    cd(P.dir.data)
                    save Comp M
                    h=errordlg('CurrZVar is empty');
                    waitfor(h)
                    error('process canceled')
                end
            end
            if P.flag.loadData==0
                loadsave_comp('save',ResRank,0,P.dir.data,1,1,1,1,1,P.chip.currProbeSetNb,ZVar,Fdr,Sensitivity,Pv,FC,TotalVar);
            end
            cd(P.dir.data)
            save Comp M
            M=[];
        end
        %do not clear S before since it is not loaded in rdam
        S=[];
    case 'network'
        BlocList1=varargin{1};
        BlocList2=varargin{2};
        IndexRank=varargin{3};
        BlocIndex=[];
        BlocRank=0;
        Continue=1;
        %verify that current index corresponds
        if sum(P.net.biolIndex)>0
            if ~isequal(P.net.biolIndex,P.net.biolIndexes{IndexRank})
                %search currently used BiolIndex
                CurrIndexRank=0;                
                for BiolL=1:length(P.net.biolIndexes)
                    if isequal(P.net.biolIndex,P.net.biolIndexes{BiolL})
                        CurrIndexRank=BiolL;
                        break
                    end
                end
                if CurrIndexRank==0
                    h=warndlg(sprintf('Index selected previously does not correspond to the one in position %u',IndexRank));
                    waitfor(h)
                    SelIndexRank=questdlg('Which index do you want to use','',num2str(IndexRank),'none','none');
                else
                    h=warndlg(sprintf('Index selected previously correspond to the one in position %u (not %u)',CurrIndexRank,IndexRank));
                    waitfor(h)
                    SelIndexRank=questdlg('Which index do you want to use','',num2str(CurrIndexRank),num2str(IndexRank),'none','none');
                end
                if isequal(SelIndexRank,'none')|isempty(SelIndexRank)%questdlg has been closed by user
                    Continue=0;
                else
                    IndexRank=str2num(SelIndexRank);
                end
            end
        end

        if Continue
            %DataRanks are filled here if P.flag.loadData=1, so do not load data again in rdam;
            LoadDataFlag=0;
            SingleCalibPointFlag=0;
            SingleCalibCurveFlag=0;
            CalibUpdateSFlag=0;
            CalibSaveFlag=1;
            DisplayFlag=0;
            SizerFittingDegree=7;
            ComparisonFlag=1;
            CalibSchemeFlag=2;
            CompScheme={[1,2;1,2],[1,2;2,1]};
            TGRankList=1:2;
            CGRankList=3:4;



            cd(P.dir.data)
            fid=fopen(sprintf('log_%s_%s.txt',P.project.name,date),'w');
            BlocNb=length(BlocList1)*length(BlocList2);
            MadeNb=0;
            t=0;
            BiolList=find(P.net.biolIndexes{IndexRank});

            for BlocL1=1:length(BlocList1)
                Bloc1=BlocList1(BlocL1);
                for BlocL2=1:length(BlocList2)
                    Bloc2=BlocList2(BlocL2);
                    if Bloc2>=Bloc1
                        if BlocRank>0
                            MadeNb=MadeNb+1;
                            t=t+toc;
                            cd(P.dir.data)
                            %RawNb=size(BlocIndex,1);
                            %ColNb=size(BlocIndex,2);
                            eval(sprintf('save bloc%u BlocIndex',BlocRank))
                            %save_data(BlocIndex,sprintf('bloc%u.uint16le',BlocRank),P.dir.data,'r+','uint16','ieee-le',RawNb,[],1:ColNb)
                            sprintf('%u bloc processed (currrent=%u), %uH elapsed: I evaluate that the process will end in around %uH',MadeNb,BlocRank,round(t/3600),round(t*(BlocNb-MadeNb)/(MadeNb*3600)))
                            tic
                        else
                            tic
                        end
                        BlocRank=P.net.blocGrid(Bloc1,Bloc2);
                        cd(P.dir.data)
                        %if exist(sprintf('bloc%u.uint16le',BlocRank),'file')
                        if exist(sprintf('bloc%u.mat',BlocRank),'file')
                            %BlocIndex=load_data(sprintf('bloc%u.uint16le',BlocRank),P.dir.data,0,0,'uint16','ieee-le');
                            eval(sprintf('load bloc%u',BlocRank))
                        else
                            BlocIndex=[];
                        end
                        %comparison between selected biol cond
                        BiolIndex1=find(BiolList>(Bloc1-1)*100&BiolList<=Bloc1*100);
                        BiolList1=BiolList(BiolIndex1);
                        BiolIndex2=find(BiolList>(Bloc2-1)*100&BiolList<=Bloc2*100);
                        BiolList2=BiolList(BiolIndex2);
                        if length(BiolList1)>0
                            for BiolL1=1:length(BiolList1)
                                Biol1=BiolList1(BiolL1);
                                if length(BiolList2)>0
                                    for BiolL2=1:length(BiolList2)
                                        Biol2=BiolList2(BiolL2);
                                        if Biol1~=Biol2
                                            CompPos=[];
                                            if ~isempty(BlocIndex)
                                                CompPos=find(BlocIndex(:,1)==Biol1&BlocIndex(:,2)==Biol2);
                                            end
                                            if isempty(CompPos)
                                                DataRanks=[];
                                                DataRanks(:,1)=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le',1:P.chip.currProbeSetNb,P.biol.pairs{Biol1}(1));
                                                DataRanks(:,2)=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le',1:P.chip.currProbeSetNb,P.biol.pairs{Biol1}(2));
                                                DataRanks(:,3)=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le',1:P.chip.currProbeSetNb,P.biol.pairs{Biol2}(1));
                                                DataRanks(:,4)=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le',1:P.chip.currProbeSetNb,P.biol.pairs{Biol2}(2));
                                                if length(find(DataRanks(:,1)==DataRanks(:,3)))~=length(DataRanks(:,1))&...
                                                        length(find(DataRanks(:,1)==DataRanks(:,4)))~=length(DataRanks(:,1))&...
                                                        length(find(DataRanks(:,2)==DataRanks(:,3)))~=length(DataRanks(:,1))&...
                                                        length(find(DataRanks(:,2)==DataRanks(:,4)))~=length(DataRanks(:,1))
                                                    %... so TG and CG RankLists are continuous and start at 1
                                                    TGRankList4S=P.biol.pairs{Biol1};
                                                    CGRankList4S=P.biol.pairs{Biol2};
                                                    eval(sprintf('[CurrZVar,CurrPv,CurrPpv,CurrFdr,CurrSensitivity,CurrTotalVar]=rdam(CompScheme,TGRankList,CGRankList,%u,[0,0],''%s'',[],''%s'',''%s'',%u,%u,%u,%u,%u,%u,%u,%u,%u,TGRankList4S,CGRankList4S);',...
                                                        LoadDataFlag,CalibType,NormType,AnalyseType,SizerFittingDegree,SingleCalibPointFlag,...
                                                        SingleCalibCurveFlag,CalibUpdateSFlag,CalibSaveFlag,DisplayFlag,ComparisonFlag,...
                                                        ResRank,CalibSchemeFlag));
                                                    if ~isempty(CurrZVar)
                                                        BlocIndex=[BlocIndex;[Biol1,Biol2,length(find(CurrFdr>0&CurrFdr<=0.001)),length(find(CurrFdr<0&CurrFdr>=-0.001)),...
                                                            length(find(CurrFdr>0&CurrFdr<=0.01)),length(find(CurrFdr<0&CurrFdr>=-0.01)),...
                                                            length(find(CurrFdr>0&CurrFdr<=0.1)),length(find(CurrFdr<0&CurrFdr>=-0.1))]];
                                                        save_data(CurrFdr,sprintf('bloc%u.float32le',BlocRank),P.dir.data,'r+','single','ieee-le',P.chip.currProbeSetNb,[],size(BlocIndex,1))
                                                    end
                                                else
                                                    fprintf(fid,'Net-identical points in %u and %u\n',Biol1,Biol2)
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        %systematic intra experiment comparisons
                        BiolList1=((Bloc1-1)*100)+1:Bloc1*100;
                        BiolList2=((Bloc2-1)*100)+1:Bloc2*100;
                        for BiolL1=1:length(BiolList1)
                            Biol1=BiolList1(BiolL1);
                            if Biol1<=P.biol.nb
                                if length(P.biol.pairs{Biol1})==2
                                    Exp1=P.point.expRank(P.biol.pointIndex{Biol1}(1));
                                    for BiolL2=1:length(BiolList2)
                                        Biol2=BiolList2(BiolL2);
                                        if Biol2<=P.biol.nb
                                            if length(P.biol.pairs{Biol2})==2
                                                Exp2=P.point.expRank(P.biol.pointIndex{Biol2}(1));
                                                if Biol1~=Biol2&Exp1==Exp2
                                                    CompPos=[];
                                                    if ~isempty(BlocIndex)
                                                        CompPos=find(BlocIndex(:,1)==Biol1&BlocIndex(:,2)==Biol2);
                                                    end
                                                    if isempty(CompPos)
                                                        DataRanks=[];
                                                        DataRanks(:,1)=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le',1:P.chip.currProbeSetNb,P.biol.pairs{Biol1}(1));
                                                        DataRanks(:,2)=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le',1:P.chip.currProbeSetNb,P.biol.pairs{Biol1}(2));
                                                        DataRanks(:,3)=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le',1:P.chip.currProbeSetNb,P.biol.pairs{Biol2}(1));
                                                        DataRanks(:,4)=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le',1:P.chip.currProbeSetNb,P.biol.pairs{Biol2}(2));
                                                        if length(find(DataRanks(:,1)==DataRanks(:,3)))~=length(DataRanks(:,1))&...
                                                                length(find(DataRanks(:,1)==DataRanks(:,4)))~=length(DataRanks(:,1))&...
                                                                length(find(DataRanks(:,2)==DataRanks(:,3)))~=length(DataRanks(:,1))&...
                                                                length(find(DataRanks(:,2)==DataRanks(:,4)))~=length(DataRanks(:,1))

                                                            %... so TG and CG RankLists are continuous and start at 1
                                                            TGRankList4S=P.biol.pairs{Biol1};
                                                            CGRankList4S=P.biol.pairs{Biol2};
                                                            eval(sprintf('[CurrZVar,CurrPv,CurrPpv,CurrFdr,CurrSensitivity,CurrTotalVar]=rdam(CompScheme,TGRankList,CGRankList,%u,[0,0],''%s'',[],''%s'',''%s'',%u,%u,%u,%u,%u,%u,%u,%u,%u,TGRankList4S,CGRankList4S);',...
                                                                LoadDataFlag,CalibType,NormType,AnalyseType,SizerFittingDegree,SingleCalibPointFlag,...
                                                                SingleCalibCurveFlag,CalibUpdateSFlag,CalibSaveFlag,DisplayFlag,ComparisonFlag,...
                                                                ResRank,CalibSchemeFlag));
                                                            if ~isempty(CurrZVar)
                                                                BlocIndex=[BlocIndex;[Biol1,Biol2,length(find(CurrFdr>0&CurrFdr<=0.001)),length(find(CurrFdr<0&CurrFdr>=-0.001)),...
                                                                    length(find(CurrFdr>0&CurrFdr<=0.01)),length(find(CurrFdr<0&CurrFdr>=-0.01)),...
                                                                    length(find(CurrFdr>0&CurrFdr<=0.1)),length(find(CurrFdr<0&CurrFdr>=-0.1))]];
                                                                save_data(CurrFdr,sprintf('bloc%u.float32le',BlocRank),P.dir.data,'r+','single','ieee-le',P.chip.currProbeSetNb,[],size(BlocIndex,1))
                                                            end
                                                        else
                                                            fprintf(fid,'Exp-identical points in %u and %u\n',Biol1,Biol2)
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            if BlocRank>0
                MadeNb=MadeNb+1;
                t=t+toc;
                cd(P.dir.data)
                eval(sprintf('save bloc%u BlocIndex',BlocRank))
                %save_data(BlocIndex,sprintf('bloc%u.uint16le',BlocRank),P.dir.data,'r+','uint16','ieee-le')
                sprintf('%u bloc processed (currrent=%u), %uH elapsed: I evaluate that the process will end in around %uH',MadeNb,BlocRank,round(t/3600),round(t*(BlocNb-MadeNb)/(MadeNb*3600)))
            end
            fclose(fid)
        end
end