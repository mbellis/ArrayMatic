%=============
% MCL_SORTCLU
%=============

% MCL_SORTCLU - sorts several columns of a clustering result by keeping correspondance between
% probe sets (each column may correspond to a clustering result made with different
% paramters)
% The first column must be ordered according to the cluster ranks (if it is not the case, set
% SortFirstFlag to 1)

% PARAMETERS
%           Clu : the clustering result to be sorted (size = Nb of probe sets x Nb of clustering results)
%         Index : probe set order
%        ColorNb: number of colors usd to display clusters (used if NestedFlag==1)
%  MoveZeroFlag : indicates if the non clusterized probe sets (clusters 0) must be moved to
%                 the end
% SortFirstFlag : indicates if the first column must be sorted
% NestedFlag : indicates that a cluster at a given level of corr limit is a subset of
%                another cluster defined at a lower corr limit 

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

function [SClu,Index]=mcl_sortclu(Clu,Index,ColorNb,MoveZeroFlag,SortFirstFlag,NestedFlag)

MyColors=[0.70,0.70,0.70;
0.65,0.81,0.89;
0.12,0.47,0.71;
0.70,0.87,0.54;
0.20,0.63,0.17;
0.98,0.60,0.60;
0.89,0.10,0.11;
0.99,0.75,0.44;
1.00,0.50,0.00;
0.79,0.70,0.84;
0.42,0.24,0.60;
1.00,1.00,0.60];



if size(Index,1)==1
    ColumnFlag=0;
    Index=Index';
elseif size(Index,2)==1
    ColumnFlag=1;
else
    h=errordlg('Index must be a one dimentional vector');
    waitfor(h)
    error('process canceled')
end
PsNb=size(Clu,1);
ColNb=size(Clu,2);
if SortFirstFlag
    [tmp,SortIndex]=sort(Clu(:,1));
    Clu=Clu(SortIndex,:);
    Index=Index(SortIndex);
end
%assumes that first column is sorted
if MoveZeroFlag
    ZeroPos=find(Clu(:,1)==0);
    if ~isempty(ZeroPos)
        ZeroNb=length(ZeroPos);
        if ZeroPos(end)~=PsNb
            if ZeroPos(1)==1
                Clu=[Clu(ZeroNb+1:PsNb,:);Clu(1:ZeroNb,:)];
                Index=[Index(ZeroNb+1:PsNb,:);Index(1:ZeroNb,:)];
            else
                Clu=[Clu(1:ZeroPos-1,:);Clu(ZeroPos,:);Clu(ZeroPos(end)+1:end,:)];
                Index=[Index(1:ZeroPos-1,:);Index(ZeroPos,:);Index(ZeroPos(end)+1:end,:)];
            end
        end
    end
end
SClu=zeros(size(Clu));
SClu(:,1)=Clu(:,1);

%process second column
%process each cluster of the first column
for CluL=min(SClu(:,1)):max(SClu(:,1))
    PsPos=find(SClu(:,1)==CluL);
    if ~isempty(PsPos)
        %reconstruct the second column and keep in register all probe sets in other columns
        CluVal2=unique(Clu(PsPos,2));
        CluNb2=histc(Clu(PsPos,2),CluVal2);
        [CluNb2,SortIndex]=sort(CluNb2,'descend');
        CluVal2=CluVal2(SortIndex);
        StartPos=PsPos(1);
        MemIndex=Index;
        for CluL2=1:length(CluVal2)
            SClu(StartPos:StartPos+CluNb2(CluL2)-1,2)=CluVal2(CluL2);
            %reorder other columns
            PsPos1=find(Clu(PsPos,2)==CluVal2(CluL2));
            SClu(StartPos:StartPos+CluNb2(CluL2)-1,3:ColNb)=Clu(PsPos(PsPos1),3:ColNb);
            %update Index
            Index(StartPos:StartPos+CluNb2(CluL2)-1)=MemIndex(PsPos(PsPos1));
            StartPos=StartPos+CluNb2(CluL2);
        end        
        if ~isempty(find(CluVal2==0)) & MoveZeroFlag
            ZeroPos=find(SClu(PsPos,2)==0);
            if ~isempty(ZeroPos)
                ZeroNb=length(ZeroPos);
                if ZeroPos(end)~=length(PsPos)
                    if ZeroPos(1)==1
                        SClu(PsPos,:)=[SClu(PsPos(ZeroNb+1):PsPos(end),:);SClu(PsPos(1:ZeroNb),:)];
                        Index(PsPos)=[Index(PsPos(ZeroNb+1):PsPos(end));Index(PsPos(1:ZeroNb))];
                    else
                        SClu(PsPos,:)=[SClu(PsPos(1):PsPos(ZeroPos(1)-1),:);SClu(PsPos(ZeroPos(end)+1:end),:);SClu(PsPos(ZeroPos),:)];
                        Index(PsPos,:)=[Index(PsPos(1):PsPos(ZeroPos(1)-1));Index(PsPos(ZeroPos(end)+1:end));Index(PsPos(ZeroPos))];

                    end
                end
            end
        end
        %change color of the first most populated cluster
        if CluVal2(1)~=CluL & CluVal2(1)~=0
            ModifPos=find(SClu(:,2)==CluL);
            SClu(PsPos(1):PsPos(1)+CluNb2(1)-1,2)=CluL;
            SClu(ModifPos,2)=CluVal2(1);
        end
    end
end



    



%process other columns
%process each cluster of the first column
for CluL=1:max(SClu(:,1))
    %for CluL=1
    %position of cluster CluL in the first column
    PsPos=find(SClu(:,1)==CluL);
    %process each other column
    %for ColL1=2
    for ColL1=2:ColNb-1                
        %process each subdivision of the current column
        % (a cluster may be, at this step splitted in several division)
        % use the preceeding divisions in case of merging of several division with the same
        CluPos1=diff(SClu(PsPos,ColL1));
        CluPos1=find(CluPos1);
        if ColL1>2
            % cluster rank (specifically with rank = 0)
            CluPos0=diff(SClu(PsPos,ColL1-1));
            CluPos0=find(CluPos0);
            CommonLimit=intersect(CluPos0,CluPos1);
            if ~isempty(CommonLimit)
                UnionPsPos=zeros(length(PsPos),1);
                UnionPsPos(CluPos0)=1;
                UnionPsPos(CluPos1)=1;
                CluPos1=find(UnionPsPos);
            end
        end
        StartPos1=PsPos(1);
        if isempty(CluPos1)
            %only one cluster
            EndCluL1=1;
            CluPos1=length(PsPos);
        else
            EndCluL1=length(CluPos1)+1;
        end
        for CluL1=1:EndCluL1
            %position of the cluster CluL1 in the current column
            if CluL1==1
                PsPos1=[1:CluPos1(1)];
            elseif CluL1==EndCluL1
                if CluPos1(end)==length(PsPos)-1
                    PsPos1=length(PsPos);
                else
                    PsPos1=[CluPos1(CluL1-1)+1:length(PsPos)];
                end
            else
                PsPos1=[CluPos1(CluL1-1)+1:CluPos1(CluL1)];
            end
            %find cluster in the next column at position of the cluster CluL1
            CluVal2=unique(SClu(PsPos(PsPos1),ColL1+1));
            %reorder
            CluNb2=histc(SClu(PsPos(PsPos1),ColL1+1),CluVal2);
            [CluNb2,SortIndex]=sort(CluNb2,'descend');
            CluVal2=CluVal2(SortIndex);
            %reorder clusters in the next column
            %recover positions of each cluster in the next column
            StartPos2=StartPos1(1);
            MemSClu=SClu;
            MemIndex=Index;
            for CluL2=1:length(CluVal2)
                PsPos2=find(MemSClu(PsPos(PsPos1),ColL1+1)==CluVal2(CluL2));                
                SClu(StartPos2:StartPos2+CluNb2(CluL2)-1,ColL1+1)=CluVal2(CluL2);
                %update Index
                CurrPsPos=PsPos(PsPos1);
                CurrPsPos=CurrPsPos(PsPos2);
                Index(StartPos2:StartPos2+CluNb2(CluL2)-1)=MemIndex(CurrPsPos);
                %reorder other columns
                if ColL1<ColNb-1                    
                    SClu(StartPos2:StartPos2+CluNb2(CluL2)-1,ColL1+2:ColNb)=MemSClu(CurrPsPos,ColL1+2:ColNb);
                end
                StartPos2=StartPos2+CluNb2(CluL2);
            end
            if ~isempty(find(CluVal2==0)) & MoveZeroFlag
                ZeroPos=find(SClu(PsPos(PsPos1),ColL1+1)==0);
                if ~isempty(ZeroPos)
                    ZeroNb=length(ZeroPos);
                    if ZeroPos(end)~=length(PsPos1)
                        if ZeroPos(1)==1
                            SClu(PsPos(PsPos1),:)=[SClu(PsPos(PsPos1(ZeroNb+1)):PsPos(PsPos1(end)),:);SClu(PsPos(PsPos1(1:ZeroNb)),:)];
                            Index(PsPos(PsPos1))=[Index(PsPos(PsPos1(ZeroNb+1)):PsPos(PsPos1(end)));Index(PsPos(PsPos1(1:ZeroNb)))];
                        else
                            SClu(PsPos(PsPos1),:)=[SClu(PsPos(PsPos1(1)):PsPos(PsPos1(ZeroPos(1)-1)),:);SClu(PsPos(PsPos1(ZeroPos(end)+1)):PsPos(PsPos1(end)),:);SClu(PsPos(PsPos1(ZeroPos)),:)];
                            Index(PsPos(PsPos1))=[Index(PsPos(PsPos1(1)):PsPos(PsPos1(ZeroPos(1)-1)));Index(PsPos(PsPos1(ZeroPos(end)+1)):PsPos(PsPos1(end)));Index(PsPos(PsPos1(ZeroPos)))];
                        end
                    end
                end
            end
            StartPos1=StartPos1+length(PsPos1);
        end
        %change color of the first most populated cluster, but don't if the current most
        %populated cluster is less than CluL (to prevent form modifying previous change)
        if SClu(PsPos(1),ColL1+1)>CluL & SClu(PsPos(1),ColL1+1)~=0
            ModifPos1=find(SClu(:,ColL1+1)==SClu(PsPos(1),ColL1+1));
            ModifPos2=find(SClu(:,ColL1+1)==CluL);
            SClu(ModifPos2,ColL1+1)=SClu(PsPos(1),ColL1+1);
            SClu(ModifPos1,ColL1+1)=SClu(PsPos(1),1);
        end       
    end
end


%modify the second cluster in order to have a different color than the first one

for CluL=1:max(SClu(:,1))
    PsPos=find(SClu(:,1)==CluL);
    %control color of second cluster
    for ColL1=2:ColNb-1
        CluPos1=diff(SClu(PsPos,ColL1));
        CluPos1=find(CluPos1);

        if isempty(CluPos1)
            %only one cluster
            EndCluL1=1;
            CluPos1=length(PsPos);
        else
            EndCluL1=length(CluPos1)+1;
        end


        CluVal1=SClu(PsPos(1),ColL1);
        if length(CluPos1)>1
            CluVal2=SClu(PsPos(CluPos1(2)),ColL1);
        else
            CluVal2=0;
        end
        if length(CluPos1)>2
            CluVal3=SClu(PsPos(CluPos1(3)),ColL1);
        else
            CluVal3=0;
        end

        if mod(CluVal1,ColorNb)==mod(CluVal2,ColorNb) & CluVal2~=0
            ModifPos1=find(SClu(PsPos,ColL1)==CluVal2);
            if CluVal3==0
                CluVal4=CluVal2+1;
            else
                CluVal4=CluVal2+1;
                while mod(CluVal4,ColorNb)==mod(CluVal3,ColorNb) | mod(CluVal4,ColorNb)==mod(CluVal1,ColorNb)
                    CluVal4=CluVal4+1;
                end
            end
            SClu(PsPos(ModifPos1),ColL1)=CluVal4;
        end

        if NestedFlag
            %modify cluster rank (cluster rank are not longer comparable in a single column (a single
            %cluster can be splitted in two different clusters and cluster with the same rank ùay be
            %different
            if EndCluL1>2
                for CluL1=2:EndCluL1
                    %position of the cluster CluL1 in the current column
                    if CluL1==1
                        PsPos1=[1:CluPos1(1)];
                    elseif CluL1==EndCluL1
                        if CluPos1(end)==length(PsPos)-1
                            PsPos1=length(PsPos);
                        else
                            PsPos1=[CluPos1(CluL1-1)+1:length(PsPos)];
                        end
                    else
                        PsPos1=[CluPos1(CluL1-1)+1:CluPos1(CluL1)];
                    end
                    CluVal=(SClu(PsPos(PsPos1(1)),ColL1));
                    SearchPos=PsPos(PsPos1);
                    if SClu(SearchPos(1),ColL1+1)~=CluVal & SClu(SearchPos(1),ColL1+1)~=0
                        ModifPos1=find(SClu(SearchPos,ColL1+1)==SClu(SearchPos(1),ColL1+1));
                        SClu(SearchPos(ModifPos1),ColL1+1)=CluVal;
                    end
                end
            end
        end
    end
end

if ColumnFlag==0
    Index=Index';
end