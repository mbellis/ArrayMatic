% TRS_POINT_PROPERTIES
% calculate some properties on the signal distribution of a particular
% point
%
% INPUT
% Signal: signals of one point
% PointRank: point rank
% OUTPUT
% fill fields of the global variable P.point
% minSignal
% maxSignal

% VERSION
% V01   03/2010 Refactoring existing version



function trs_point_properties(Signal,PointRank)
global P



P.point.minSignal(PointRank,1)=min(Signal);
P.point.maxSignal(PointRank,1)=max(Signal);


NaNIndex=find(isnan(Signal));


P.point.nanNb(PointRank,1)=length(NaNIndex);
P.point.nullNb(PointRank,1)=length(find(Signal==0));
P.point.negNb(PointRank,1)=length(find(Signal<0));


%number of signals wich have shared values
Signal=sort(Signal);
DiffSignal=diff(Signal);
P.point.nullDiffNb(PointRank,1)=length(find(DiffSignal==0));

%number of unique signal values
UniqueVal=unique(Signal);
P.point.diversity(PointRank,1)=length(UniqueVal);

%search for the longest run of identical values
if length(UniqueVal)==length(Signal)
    RunNb=0;
else
    ValRep=histc(Signal,UniqueVal);
    [RunNb,ValPos]=max(ValRep);
    RunSignal=UniqueVal(ValPos);
    P.point.runVal(PointRank,1)=RunSignal;
    P.point.runNb(PointRank,1)=RunNb;
    P.point.negRunNb(PointRank,1)=length(find(Signal<RunSignal));
end

%consider that more than 500 identical signals indicates a threshold
if RunNb>=500
    P.point.threshVal(PointRank,1)=RunSignal;
    P.point.threshNb(PointRank,1)=RunNb;
    P.point.negThreshNb(PointRank,1)=P.point.negRunNb(PointRank);
else
    P.point.threshVal(PointRank,1)=min(Signal)-1;
    P.point.threshNb(PointRank,1)=0;
    P.point.negThreshNb(PointRank,1)=0;    
end
