%TRS_MEDIANCOMP
function trs_mediancomp(UpdateSFlag)
global P S DataRanks
%make MedianPoint if does not exist create it
cd(P.dir.data)
if exist('MedianRanks.mat','file')
    load MedianRanks
else
    trs_makemedian
end

%load CalibSet
cd(P.dir.data)
ResRank=1;
if exist('CalibSet_01.mat','file')
    load CalibSet_01
else
    S=[];
end

%load all Results if P.flag.loadData==0 (otherwise write results directly
%in files)
if P.flag.loadData==0
    [ZVar,Fdr,Sensitivity,Pv,Fc]=loadsave_comp('load',1,1:P.point.nb,P.dir.data,1,1,1,1,0,P.chip.currProbeSetNb,P.point.nb,1);
end

%if P.flag.loadData==0, add temporarily MedianRanks
if P.flag.loadData==0
    if P.point.currNb==P.point.nb
        DataRanks=[DataRanks,MedianRanks];
        P.point.currNb=P.point.nb+1;
    end
else
    DataRanks=zeros(length(MedianRanks),2);
    DataRanks(:,2)=MedianRanks;   
end

%CONSTRUCT CALIB SET

%noise_distribution parameters
RankThreshold=[0,0];
ClearIndex=[];
NormType='quantile';
AnalyseType='transcriptome';
CalibType='idem';
SizerFittingDegree=7;
DisplayFlag=0;
%do not save S at the end of each calibration set calculus
SaveFlag=0;
ClearIndexFlag=0;
UpdateSFlag=0;

BL=MedianRanks;
SNb=0;
tic
t=0;
SaveSFlag=0;
for PointL=1:P.point.nb
    if P.point.used(PointL)
        SNb=SNb+1;
        CalibRanks=[PointL,P.point.nb+1];
        %do calibration if necessary
        DoS=1;
        if ~UpdateSFlag
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
        end
        if DoS
            if P.flag.loadData==0
                HL=DataRanks(:,PointL);
            else
                HL=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le',1:P.chip.currProbeSetNb,PointL);
            end
            [RankGrid,Grid,ZVarGrid,ZVar]=noise_distribution(HL,BL,RankThreshold,CalibRanks(1),CalibRanks(2),ClearIndex,NormType,AnalyseType,CalibType,SizerFittingDegree,DisplayFlag);            
            fill_s(NormType,CalibType,ClearIndexFlag,CalibRanks,Grid,ZVar,ZVarGrid,ResRank,SaveFlag)
            SaveSFlag=1;
        end
        if SNb==100
            t=t+toc;
            SNb=0;
            sprintf('%u comparisons processed, %u min elapsed: I evaluate that the process will end in around %u min',PointL,round(t/60),round(t*(P.point.nb-PointL)/(PointL*60)))
            tic
        end
    end
end
if SaveSFlag
    cd(P.dir.data)
    save CalibSet_01 S
end


%MAKE COMPARISONS
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
CalibSchemeFlag=1;
cd(P.dir.data)
if P.flag.loadData==0
    ZVar=single(zeros(P.chip.currProbeSetNb,P.point.nb));
    Sensitivity=ZVar;
    Fdr=ZVar;
    Pv=ZVar;    
end

for PointL=1:P.point.nb
    if P.point.used(PointL)
        PointL
        CompScheme='{[1;1]}';
        CalibScheme=sprintf('{[%u;%u;3]}',PointL,P.point.nb+1);
        if P.flag.loadData==0
            TGRankList=PointL;
            CGRankList=P.point.nb+1;
        else
            DataRanks(:,1)=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le',1:P.chip.currProbeSetNb,PointL);
            TGRankList=1;
            CGRankList=2;
        end
        CurrZVar=[];
        eval(sprintf('[CurrZVar,CurrPv,CurrPpv,CurrFdr,CurrSensitivity,CurrTotalVar]=rdam(%s,[%u],[%u],%u,[0,0],''%s'',[],''%s'',''%s'',%u,%u,%u,%u,%u,%u,%u,%u,%u,%s);',...
            CompScheme,TGRankList,CGRankList,LoadDataFlag,CalibType,NormType,AnalyseType,SizerFittingDegree,SingleCalibPointFlag,...
            SingleCalibCurveFlag,CalibUpdateSFlag,CalibSaveFlag,DisplayFlag,ComparisonFlag,...
            ResRank,CalibSchemeFlag,CalibScheme));
        if ~isempty(CurrZVar)
            if P.flag.loadData
                loadsave_comp('save',ResRank,PointL,P.dir.data,1,1,1,1,0,P.chip.currProbeSetNb,CurrZVar,CurrFdr,CurrSensitivity,CurrPv,CurrTotalVar);
            else
                ZVar(:,PointL)=CurrZVar;
                Sensitivity(:,PointL)=CurrSensitivity;
                Fdr(:,PointL)=CurrFdr;
                Pv(:,PointL)=CurrPv;
                TotalVar.inc(PointL,1)=Values(1);
                TotalVar.dec(PointL,1)=Values(2);
            end
        else
            h=errordlg('CurrZVar is empty');
            waitfor(h)
            error('process canceled')
        end
    else
        if P.flag.loadData
            CurrZVar=single(zeros(P.chip.currProbeSetNb,1));
            CurrSensitivity=CurrZVar;
            CurrFdr=single(ones(P.chip.currProbeSetNb,1));
            CurrPv=CurrFdr;
            CurrTotalVar=[0,0];
            loadsave_comp('save',ResRank,PointL,P.dir.data,1,1,1,1,0,P.chip.currProbeSetNb,CurrZVar,CurrFdr,CurrSensitivity,CurrPv,CurrTotalVar);
        end
    end
end

if P.flag.loadData==0
    if P.point.currNb==P.point.nb+1
        DataRanks(:,end)=[];
        P.point.currNb=P.point.nb;
    end
    loadsave_comp('save',ResRank,0,P.dir.data,1,1,1,1,P.chip.currProbeSetNb,ZVar,Fdr,Sensitivity,Pv);
else
    DataRanks=[];
end
S=[];




