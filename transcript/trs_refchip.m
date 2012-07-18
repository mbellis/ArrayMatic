function [Rank,Signal]=RefChip(Signal)

%do not allow negative signal (pb with log)
if min(Signal)<=0
    Signal=Signal+min(Signal)+1;
end
Rank=signal2rank(Signal,0,1);
NanIndex=find(Signal==-1);
Signal(NanIndex)=[];
Rank(NanIndex)=[];
NegIndex=find(Rank<=0);
Signal(NegIndex)=[];
Rank(NegIndex)=[];
%order the set and sample the ref
[Rank SortIndex]=sort(Rank);
Signal=Signal(SortIndex);
GeneNb=length(Signal);
Step=round(GeneNb/100);
Signal=[Signal(1:Step:GeneNb-1);Signal(end)];
Rank=[Rank(1:Step:GeneNb-1);Rank(end)];
%eliminate Rank doublons
if length(unique(Rank))<length(Rank)
    RankDiff=diff(Rank);
    NullPos=find(RankDiff==0);
    Rank[NullPos)=[];
    Signal[NullPos)=[];
end
%enforce a point at 0 and 100 to be sure that no error occur with interp1
if Rank(1)>0
    Rank(1)=0;
end
if Rank(end)<100
    Rank(end)=100;
end