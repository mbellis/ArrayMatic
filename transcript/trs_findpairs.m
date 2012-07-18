
function trs_findpairs()
global P
cd(P.dir.data)
TRANK=1;
if exist('Tree.mat','file')
    load Tree
else
    errordlg('run *** TRANSCRIPTOME/DENDROGRAM/prepare data for tree *** first')
end
Distances=squareform(T{TRANK}.distances);
for PointL=1:P.point.nb
    Distances(PointL,PointL)=Inf;
end
% find couple of points

cd(P.dir.resTree)
fid=fopen(sprintf('modif_scouple_pv_%s.txt',date),'w');
if fid==-1
    errordlg('modif_scouple.txt not open. Process canceled');
end
fprintf(fid,'date : %s\n',date)


%find abnormal distances
NullPos=find(Distances<=5);
if ~isempty(NullPos)
    fprintf(fid,'point separated by a distance <=5\n\n');
    fprintf(fid,'already not used\tfirst point\tsnd point\tfirst biol cond\tsnd biol cond\n');
    NullIndex=[];
    for PosL=1:length(NullPos)
        PointRank1=ceil(NullPos(PosL)/P.point.nb);
        if mod(NullPos(PosL),P.point.nb)==0
            PointRank2=P.point.nb;
        else
            PointRank2=NullPos(PosL)-P.point.nb*(PointRank1-1);
        end
        if PointRank1<PointRank2
            NullIndex=[NullIndex;[PointRank1,PointRank2,P.point.biolRank(PointRank1),P.point.biolRank(PointRank2)]];
            if P.point.used(PointRank1)
                fprintf(fid,'0\t%u\t%u\t%u\t%u\n',PointRank1,PointRank2,P.point.biolRank(PointRank1),P.point.biolRank(PointRank2));
            else
                fprintf(fid,'1\t%u\t%u\t%u\t%u\n',PointRank1,PointRank2,P.point.biolRank(PointRank1),P.point.biolRank(PointRank2));
            end
        end
    end
    fprintf(fid,'\n\n');
end

            
        


%update information
%mark as not used points which are identical to others
ClearedBiol=0;
for PosL=1:size(NullIndex)
    if NullIndex(PosL,1)~=0&NullIndex(PosL,2)~=0
        PointRank1=NullIndex(PosL,1);
        PointRank2=NullIndex(PosL,2);
        BiolRank1=P.point.biolRank(PointRank1);
        BiolRank2=P.point.biolRank(PointRank2);
        %delete the point belonging to the biol condition having the greatest
        %nomber of used points
        %try to keep pairs of point in at least one biol conditions
        PointNb=length(find((NullIndex(:,3)==NullIndex(PosL,3)&NullIndex(:,4)==NullIndex(PosL,4))|(NullIndex(:,4)==NullIndex(PosL,3)&NullIndex(:,3)==NullIndex(PosL,4))))
        if PointNb>1
            Pos=find(NullIndex(:,3)==NullIndex(PosL,3)&NullIndex(:,4)==NullIndex(PosL,4));
            PointRanks1=NullIndex(Pos,1);
            PointRanks2=NullIndex(Pos,2);
            Pos=find(NullIndex(:,4)==NullIndex(PosL,3)&NullIndex(:,3)==NullIndex(PosL,4));
            PointRanks1=[PointRanks1;NullIndex(Pos,2)];
            PointRanks2=[PointRanks2;NullIndex(Pos,1)];
            if length(find(P.point.used(P.biol.pointIndex{BiolRank1})))>length(find(P.point.used(P.biol.pointIndex{BiolRank2})))
                PointRanks=PointRanks1;
                BiolRank=BiolRank1;
            else
                PointRanks=PointRanks2;
                BiolRank=BiolRank2;
            end
        else
            if length(find(P.point.used(P.biol.pointIndex{BiolRank1})))>length(find(P.point.used(P.biol.pointIndex{BiolRank2})))
                PointRanks=PointRank1;
                BiolRank=BiolRank1;
            else
                PointRanks=PointRank2;
                BiolRank=BiolRank2;
            end
        end
        for PointL=1:length(PointRanks)
            Pos=find(NullIndex(:,1)==PointRanks(PointL));
            NullIndex(Pos,1)=0;
            Pos=find(NullIndex(:,2)==PointRanks(PointL));
            NullIndex(Pos,2)=0;
            P.point.used(PointRanks(PointL))=0;
            Pos=find(P.biol.pointIndex{BiolRank}==PointRanks(PointL));
            if  ~isempty(Pos)
                P.biol.pointIndex{BiolRank}(Pos)=[];
                if isempty(P.biol.pointIndex{BiolRank})
                    P.biol.used(BiolRank)=0;
                    if ClearedBiol==0
                        fprintf(fid,'\n\nnot used empty biol cond\n\n');
                        ClearedBiol=1;
                    end
                    fprintf(fid,'%u\n',BiolRank)
                end
            end
        end
    end
end
if ClearedBiol
    fprintf(fid,'\n\n');
end

fprintf(fid,'pairs of points for each biol cond\n\n')

BiolCondNb=P.biol.nb;
for BiolCondL=1:BiolCondNb
    %keep only used points
    CurrPointIndex=P.biol.pointIndex{BiolCondL};
    UsedPoint=P.point.used(CurrPointIndex);
    CurrPointIndex=CurrPointIndex(UsedPoint==1);
    PointNb=length(CurrPointIndex);
    if PointNb>1
        %recover the matrix of distances corresponding to the current biol condition
        CurrDistances=Distances(CurrPointIndex,CurrPointIndex);
        % replace 0 by Inf to find the real minimum (first diagonal == 0)

        [MinC CIndex]=min(CurrDistances);
        [MinL LIndex]=min(MinC);
        MinC=CurrPointIndex(CIndex(LIndex));
        MinL=CurrPointIndex(LIndex);

        NewPairFlag=0;
        if isfield(P.biol,'pairs')
            if length(P.biol.pairs)>=BiolCondL
                if isempty(P.biol.pairs{BiolCondL})
                    NewPairFlag=1;
                elseif P.biol.pairs{BiolCondL}(1)~=min(MinC,MinL)|P.biol.pairs{BiolCondL}(2)~=max(MinC,MinL)
                    NewPairFlag=2;
                end
            else
                NewPairFlag=1;
            end
        else
            NewPairFlag=1;
        end
        if MinC~=MinL
            if NewPairFlag==1
                %fill with current pair
                fprintf(fid,'     C%04u %s : [] => [%u %u]\n',BiolCondL,P.biol.name{BiolCondL},MinC,MinL);
            elseif NewPairFlag==2
                %replace by current pair
                fprintf(fid,'!!!  C%04u %s : [%u,%u] => [%u,%u]\n',BiolCondL,P.biol.name{BiolCondL},P.biol.pairs{BiolCondL}(1),P.biol.pairs{BiolCondL}(2),min(MinC,MinL),max(MinC,MinL));
            end
            P.biol.pairs{BiolCondL,1}=sort([MinC,MinL]);
        else
            if NewPairFlag==1
                %do not fill pair
                fprintf(fid,'***  C%04u %s : same points [%u,%u] => []\n',BiolCondL,P.biol.name{BiolCondL},min(MinC,MinL),max(MinC,MinL));
            elseif NewPairFlag==2
                %keep existing pair
                fprintf(fid,'***  C%04u %s : same points [%u,%u] => [%u,%u]\n',BiolCondL,P.biol.name{BiolCondL},min(MinC,MinL),max(MinC,MinL),P.biol.pairs{BiolCondL}(1),P.biol.pairs{BiolCondL}(2));
            end
        end
    else
        if length(P.biol.pairs{BiolCondL})==2
            fprintf(fid,'!!!  C%04u %s : [%u,%u] => []\n',BiolCondL,P.biol.name{BiolCondL},P.biol.pairs{BiolCondL}(1),P.biol.pairs{BiolCondL}(2));
            P.biol.pairs{BiolCondL}=[];
        end
    end
end
fclose(fid)
cd(P.dir.project)
eval(sprintf('save %s P',P.project.name))
