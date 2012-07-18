%TRS_BIOLMEDIANCOMP makes comparison between biological conditions and median biological conditions which
%                   is the biological condition containing the median point

%INPUT PARAMETERS
%none

%EXTERNAL FILES
%Dataranks

%OUTPUT PARAMETERS
%none


%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%
%                          c) Michel Bellis                                                %
%                          michel.bellis@crbm.cnrs.fr                                      %
%            Affiliation:  CNRS (Centre National de la Recherche Scientifique - France)    %                               
%  Bioinformatic Project:  ARRAYMATIC => http://code.google.com/p/arraymatic               %
%        Code Repository:  GITHUB => http://github.com/mbellis                             %
%          Personal Page:  http://bns.crbm.cnrs.fr                                         %
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
%  THIS CODE IS DISTRIBUTED UNDER THE CeCILL LICENSE, WHICH IS COMPATIBLE WITH       %
%  THE GNU GENERAL PUBLIC LICENCE AND IN ACCORDANCE WITH THE EUROPEAN LEGISLATION.   %
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%


function trs_biolmediancomp()
global DataRanks P S

if isequal(P.par.analType,'network')
    cd(P.dir.data)
    load MedianRanks
    ResRank=2;
    DataRanks=[];
    DataRanks(:,3)=MedianRanks;
    CompScheme='{[1,2;1,1]}';
    TGRankList='1:2';
    CGRankList=3;
elseif isequal(P.par.analType,'transcriptome')
    h=errordlg('not yet implemented for simple transcriptome studies');
    waitfor(h)
    errror('process canceled')
else
    h=errordlg(sprintf('not yet implemented for %s studies',P.par.analType));
    waitfor(h)
    errror('process canceled')
end

if exist('CalibSet_02.mat','file')
    cd(P.dir.data)
    load CalibSet_02
else
    S=[];
end


fid=fopen(sprintf('UsedBiolCond4Comp_%s_%s.txt',P.project.name,date),'w');
fprintf(fid,'ExpRank\tExp\tBiolRank\tBioCond\n');
BiolIndex=[];
for BiolL=1:P.biol.nb
    if length(P.biol.pairs{BiolL})==2&P.biol.used(BiolL)&P.biol.corrCoeff{BiolL}>=P.biol.corrLimits(1)&P.biol.corrCoeff{BiolL}<=P.biol.corrLimits(2)
        BiolIndex=[BiolIndex;BiolL];
        ExpRank=P.point.expRank(P.biol.pointIndex{BiolL}(1));
        fprintf(fid,'%u\t%s\t%u\t%s\n',ExpRank,P.exp.name{ExpRank},BiolL,P.biol.name{BiolL});
    end
end
fclose(fid)
save UsedBiolIndex BiolIndex

RankThreshold=[0,0];
ClearIndex=[];
NormType='quantile';
AnalyseType='transcriptome';
CalibType='idem';
LoadDataFlag=0;
SingleCalibPointFlag=0;
SingleCalibCurveFlag=0;
CalibUpdateSFlag=0;
CalibSaveFlag=1;
DisplayFlag=0;
SizerFittingDegree=7;
ComparisonFlag=1;
CalibSchemeFlag=1;




for BiolL=1:length(BiolIndex)
    BiolL
    CurrBiol=BiolIndex(BiolL);
    DataRanks(:,1)=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le',1:P.chip.currProbeSetNb,P.biol.pairs{CurrBiol}(1));
    DataRanks(:,2)=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le',1:P.chip.currProbeSetNb,P.biol.pairs{CurrBiol}(2));
    CalibScheme=sprintf('{[%u,%u;%u,%u;3,3]}',P.biol.pairs{CurrBiol}(1),P.biol.pairs{CurrBiol}(1),P.biol.pairs{CurrBiol}(2),P.biol.pairs{CurrBiol}(2));
    CurrZVar=[];
    eval(sprintf('[CurrZVar,CurrPv,CurrPpv,CurrFdr,CurrSensitivity,CurrTotalVar]=rdam([%s],%s,[%u],%u,[0,0],''%s'',[],''%s'',''%s'',%u,%u,%u,%u,%u,%u,%u,%u,%u,%s);',...
        CompScheme,TGRankList,CGRankList,LoadDataFlag,CalibType,NormType,AnalyseType,SizerFittingDegree,SingleCalibPointFlag,...
        SingleCalibCurveFlag,CalibUpdateSFlag,CalibSaveFlag,DisplayFlag,ComparisonFlag,...
        ResRank,CalibSchemeFlag,CalibScheme));
    if ~isempty(CurrZVar)
        loadsave_comp('save',ResRank,BiolL,P.dir.data,1,1,1,1,0,P.chip.currProbeSetNb,CurrZVar,CurrFdr,CurrSensitivity,CurrPv,CurrTotalVar);
    end
end

cd(P.dir.data)
eval(sprintf('load TotalVar_%02u',ResRank));
h=figure;
set(gcf,'color',[1,1,1])
set(h,'name','TOTALVAR 10%FDR BIOL COND VS MEDIAN POINT')
plot(TotalVar.inc,TotalVar.dec,'b+','markersize',3)
xlabel('nb of increase')
ylabel('nb of idecrease')
title('variation at 10% fdr for biol cond vs median point')
cd(P.dir.resComp)
saveas(h,sprintf('biolcond_vs_medianpoint_var10pcfdr_%s_%s',P.project.name,date),'png')



DataRanks=[];
S=[];