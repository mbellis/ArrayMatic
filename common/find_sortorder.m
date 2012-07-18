% FIND_SORTORDER - Find the index of test list items in a reference list
% lists must have unique values

% c) Michel Bellis
% arraymatic@gmail.com

% INPUT
% RefList: the reference list
% TestList: the testlist
% DataType: {'string','numeric')
% FindAllFlag: {0/1} indicates if all the items of the test list must be found
%               in the reference list
%
% OUTPUT
% SortIndex: the index of test items as they are ordered in the reference list
% FoundIndex: the index of reference items have been found in the test list
% SortIndex and FoundIndex have the same size

function [SortIndex,FoundIndex]=find_sortorder(RefList,TestList,DataType,FindAllFlag)
%test unicity
if length(RefList)~=length(unique(RefList))
    errordlg('reference list has repeated values')
    error('process canceled')
end
if length(TestList)~=length(unique(TestList))
    errordlg('test list has repeated values')
    error('process canceled')
end
if FindAllFlag
    if length(RefList)==length(TestList)
        [SortedRefList RefSortIndex]=sort(RefList);
        [SortedTestList TestSortIndex]=sort(TestList);
        if sum(cellfun(@isequal,SortedTestList,SortedRefList))~=length(SortedRefList)
            errordlg('reference and test list does not have the same values')
            error('process canceled')
        else
            [Temp,ReverseIndex]=sort(RefSortIndex);
            SortIndex=TestSortIndex(ReverseIndex);
        end
        FoundIndex=1:length(TestList);
    end    
end
if FindAllFlag==0 || length(RefList)~=length(TestList)
    SortIndex=[];
    %index of found reference items
    FoundIndex=[];
    if isequal(DataType,'string')
        for RefL=1:length(RefList)
            CurrPos=strmatch(RefList{RefL},TestList,'exact');
            if ~isempty(CurrPos)
                FoundIndex=[FoundIndex;RefL];
                SortIndex=[SortIndex;CurrPos];
            end
        end
    else
        for RefL=1:length(RefList)
            CurrPos=find(TestList==RefList{RefL});
            if ~isempty(CurrPos)
                FoundIndex=[FoundIndex;RefL];
                SortIndex=[SortIndex;CurrPos];
            end
        end
    end
    if FindAllFlag
        %control that all test items have been found
        if length(SortIndex)~=length(TestList)
            errordlg(sprintf('%u test items had not been found in reference list',length(TestList)-length(SortIndex)))
            error('process canceled')
        end
    end
end

