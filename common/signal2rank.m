%SIGNAL2RANK converts a series of numeric values into a series of relative ranks on
%a range [0:100]

%INPUT
% Signal: numeric values to be transformed into ranks
% SortedFlag: indicates if the data are already sorted
% SingleFlag: indicates if the format of rank must be single (otherwise double)
% varargin{1} = Threshold: the value under or equal which all signal will have the
% same rank = 0
% varargin{2} = Saturated:   the value over or equal which all signal will have the
% same rank = 100
% Rank -2 is reserved to genes absent on array (in case of comparison between two different array type
% Rank - 1 is reserved to NaN values
% Rank 0 is reserved for absent genes (lower or equal to Threshold)
% Rank 100 is reserved for saturated values
% Rank -2 is given outside this programm


function [Rank,Success]=signal2rank(Signal,SortedFlag,SingleFlag,varargin)

Success=1;
if min(size(Signal))>1|ndims(Signal)>2
    h=warndlg('Signal must be a vector','signal2rank');
    waitfor(h)
    Success=0;
end
if Success
    %Rank must have the same size than Signal
    FlipFlag=0;
    if size(Signal,2)>1
        Signal=Signal';
        FlipFlag=1;
    end
    
    ThresholdNb=0;
    SaturatedNb=0;

    if nargin>=4
        Threshold=varargin{1};
        ThresholdNb=length(find(Signal<=Threshold));
    end
    if nargin==5
        Saturated=varargin{2};
        SaturatedNb=length(find(Signal>=Saturated));
    end

    if SortedFlag==0
        [Temp DirSortIndex]=sort(Signal);
        [Temp RevSortIndex]=sort(DirSortIndex);
    end

    NanIndex=find(isnan(Signal));
    NanNb=length(NanIndex);
    
    TotalNb=length(Signal);
    % all the values equal or less than the Threshold values are counted as one
    % and only one value (0)
    UpThresholdNb=TotalNb-ThresholdNb-SaturatedNb-NanNb;
    NumericNb=TotalNb-NanNb;

    if nargin==3
        % construction of the sorted relative position vector for actual and numeric values
        Rank=[[1:NumericNb]'*100/(NumericNb+1) ; ones(NanNb,1)*1];
    else
        % all the values equal or less than the Thresholdold values are counted as one
        % and only one value (0)
        Rank=[zeros(ThresholdNb,1); [1:UpThresholdNb]'*100/(UpThresholdNb+1); ones(SaturatedNb,1)*100 ; ones(NanNb,1)*-1];
    end
    % reorder them to have relative position vector in original order
    if SingleFlag
        Rank=single(Rank);
    end
    if SortedFlag==0
        Rank=Rank(RevSortIndex);
    end
    if FlipFlag
        Rank=Rank';
    end   
end