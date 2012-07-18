%TOPINDEXES : find the ordered index of toplist
%INPUT
%Data: a vector of numerical data
%ListSize: the size of the top list
%TopType: if TopType='min', Top list is searched from the smallest value
%         if TopType='max', Top list is searched from the greatest value
%AbsFlag: if AbsFlag=1, absolute values of negative and positive values are processed
%                       independantly
%         if AbsFlag=0, values are processed as such
%OUTPUT
%either one or two indexes
%if AbsFlag=0, only one index (Index{2}=[]) which is min or max toplist
%if AbsFlag=1, two indexes
%       either min(abs(positive values)) in Index{1} and min(abs(negative
%       values)) in Index{2} if TopType='min'
%       or max(abs(positive values)) in Index{1} and max(abs(negative
%       values)) in Index{2}if TopType='min'
function [FirstIndex,SndIndex]=topindexes(Data,ListSize,TopType,AbsFlag)
if size(Data,1)>1&&size(Data,2)>1
    errordlg('Data must be a vector in toplist.m')
else
    if size(Data,2)>1
        Data=Data';
    end
end

FirstIndex=[];
SndIndex=[];
if AbsFlag
    PosIndex=find(Data>=0);
    PosData=abs(Data(PosIndex));
    [Temp SortIndex]=sort(PosData);
    PosIndex=PosIndex(SortIndex);

    NegIndex=find(Data<0);
    NegData=abs(Data(NegIndex));
    [Temp SortIndex]=sort(NegData);
    NegIndex=NegIndex(SortIndex);

    if isequal(TopType,'min')
        FirstIndex=PosIndex(1:ListSize);
        SndIndex=NegIndex(1:ListSize);        
    else
        FirstIndex=flipud(PosIndex(end-ListSize+1:end));
        SndIndex=flipud(NegIndex(end-ListSize+1:end));
    end
else
    [Temp SortIndex]=sort(Data);
    if isequal(TopType,'min')
        FirstIndex=SortIndex(1:ListSize);
    else
        FirstIndex=flipud(SortIndex(end-ListSize+1:end));
    end
end
