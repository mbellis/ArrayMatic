%==========================%
% FUNCTION GSEA_MULTIPLEPS %
%==========================%

% GSEA_MULTIPLE test the effect of multiplicity of probe sets on GSEA score
%
%INPUT PARAMETERS
%
% 1        Species: species
% 2       ChipRank: chip model rank
% 3        Suffix: suffix appended to chip model name (e.g. '_2to22_1p' in 'm8_2to22_1p')
% 4  NetRanks: list of networks
% 5   ProbeNbLimit: minimal number of probes targeting a gene
% 6     PvCorrRank: pv(overlap) is calculated for corr limit >[0,40,50,60]. PvCorrRank
% 7  FileName: GSEA lists
% 8  NrBiolRank: rank of set of non redundant biological conditions (P.biol.nr.scoupleIndex{NrBiolRank})

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

%gsea_multipleps('mouse',8,'_1p',[7:21],1,1,'c3.tft.v3.1.symbols_gene2psrank')
%gsea_multipleps('mouse',27,'_2to22_1p',[2:22],1,1,'c3.tft.v3.1.symbols_gene2psrank',9)

function gsea_multipleps(Species,ChipRank,Suffix,NetRanks,ProbeNbLimit,PvCorrRanks,FileName,NrBiolRank)
global K

%complete K.dir
if ~isfield(K.dir,'mldata')
    K.dir.mldata='/home/mbellis/sosma/data/psawn/mldata';
    K.dir.pydata='/home/mbellis/sosma/data/psawn/pydata';
    K.dir.rawdata='/home/mbellis/sosma/data/psawn/rawdata';
    K.dir.point='/home/mbellis/sosma/raydata/mlf/point';
end

%load point info
cd(K.dir.point)
eval(sprintf('load m%u_point',ChipRank))
%find Point to be used
Point=[];
ScoupleIndex=P.biol.nr{1}.scoupleindex{NrBiolRank};
for BiolL=1:length(ScoupleIndex)
    Point=[Point,P.biol.scouple{1}{ScoupleIndex(BiolL)}(1)];
end
PointNb=length(Point);

NetNb=length(NetRanks);
%load PsPair
DirName=fullfile(K.dir.mldata,Species,sprintf('m%u%s',ChipRank,Suffix));
cd(DirName)
[PairGeneId,PairGeneName,PsName1,PsName2,PsPair(:,1),PsPair(:,2),PsPair(:,3),PsPair(:,4),PsPair(:,5),...
    PsPair(:,6),PsPair(:,7),PsPair(:,8),PsPair(:,9),PsPair(:,10),PsPair(:,11),PsPair(:,12),PsPair(:,13),PsPair(:,14),...
    PsPair(:,15),PsPair(:,16),Mark,PsPair(:,17),PsPair(:,18),PsPair(:,19),PsPair(:,20),...
    PsPair(:,21),PsPair(:,22),PsPair(:,23),PsPair(:,24),PsPair(:,25),PsPair(:,26),PsPair(:,27),...
    PsPair(:,28),PsPair(:,29),PsPair(:,30),PsPair(:,31),PsPair(:,32),PsPair(:,33),PsPair(:,34),...
    PsPair(:,35),PsPair(:,36),PsPair(:,37),PsPair(:,38),PsPair(:,39),PsPair(:,40),PsPair(:,41),...
    PsPair(:,42),PsPair(:,43),PsPair(:,44),PsPair(:,45),PsPair(:,46),PsPair(:,47),PsPair(:,48),...
    PsPair(:,49),PsPair(:,50),PsPair(:,51),PsPair(:,52),PsPair(:,53),PsPair(:,54),PsPair(:,55),...
    PsPair(:,56),PsPair(:,57),PsPair(:,58),PsPair(:,59),PsPair(:,60),PsPair(:,61),PsPair(:,62)]...
    =textread(fullfile('.',sprintf('m%u_n%u_netnb%u_probenb%u_pvcorr%u_pspair.txt',ChipRank,NetRanks(1),NetNb,ProbeNbLimit,PvCorrRanks)),...
    [repmat('%s',1,4),repmat('%u',1,16),'%c',repmat('%u',1,12),repmat('%d',1,3),repmat('%u',1,20),repmat('%d',1,3),repmat('%u',1,8)],'delimiter','\t');

PairNb=size(PsPair,1);

%load ps assignation
 cd(K.dir.chip)
 [PsName,PsGeneId,PsGeneName,PsClass,PsProbeNb]=textread(sprintf('m%u_gene.txt',ChipRank),'%s%s%s%u%u','delimiter','\t');

%load GSEA list (load GenePos containing for each GSEA list and for several chip model
% the targeting probe sets
cd(K.dir.list)
eval(sprintf('load %s.mat',FileName))

%load ps ranks
cd (fullfile('/usr/data/net',sprintf('m%u',ChipRank),sprintf('m%u_data',ChipRank)))
load DataRanks


%recover chip information
ChipPos=find(K.chip.rank==ChipRank);
PsNb=K.chip.probesetNb(ChipPos);

MinProbeNb=7;
%find pairs of probesets with different similarity starting with the less similar ones
PairSim=zeros(PairNb,1);
Pos=find(PsPair(:,3)==0&PsPair(:,9)>=MinProbeNb&PsPair(:,12>=MinProbeNb));
PairSim(Pos)=1;
Pos=find((PsPair(:,3)==1|PsPair(:,3)==2)&PsPair(:,9)>=MinProbeNb&PsPair(:,12>=MinProbeNb));
PairSim(Pos)=2;
for SimL=2:5
    Pos=find((PsPair(:,1+SimL)==1|PsPair(:,1+SimL)==2)&(PsPair(:,2+SimL)==1|PsPair(:,2+SimL)==2)&PsPair(:,9)>=MinProbeNb&PsPair(:,12>=MinProbeNb));
    PairSim(Pos)=SimL+1;
end

%construct 6 list of genes with different scheme of probe set merging
%

%detect genes that are targeted by mutliple probe sets with all the same
%similarity when taken two by two (clique of multiple ps with the same similarity)

[MultiGeneId,MultiIdx]=unique(PairGeneId);
MultiGeneName=PairGeneName(MultiIdx);
[temp,UniGeneIdx]=setdiff(PsGeneId,union(MultiGeneId,'-'));
UniGeneName=unique(PsGeneName(UniGeneIdx));
Uni=zeros(length(UniGeneName),1);
for GeneL=1:length(UniGeneName)
    PsPos=strmatch(UniGeneName{GeneL},PsGeneName,'exact');
    if length(PsPos)==1
        Uni(GeneL)=PsPos;
    else
        GeneL
    end
end
UniGeneName=upper(UniGeneName);
ClearPos=find(Uni==0);
Uni(ClearPos)=[];
UniGeneName(ClearPos)=[];




MultiPs=cell(1,6);
for SimL=1:6
    MultiPs{SimL}=cell(length(MultiGeneId),1);
end

for GeneL=1:length(MultiGeneId)
    CurrPs=strmatch(MultiGeneId{GeneL},PsGeneId,'exact');
    if isempty(CurrPs)
        sprintf('gene %s (%u) is not multiple',MultiGeneId{GeneL},GeneL)
    else
    CurrPsNb=length(CurrPs);
    CurrPairNb=CurrPsNb*(CurrPsNb-1)/2;
    %fill ps matrix with similarity values corresponding to each pair of ps
    SimVal=zeros(CurrPsNb);
    for PsL1=1:CurrPsNb-1
        for PsL2=PsL1+1:CurrPsNb
            PairPos=find(PsPair(:,1)==CurrPs(PsL1)&PsPair(:,2)==CurrPs(PsL2));
            SimVal(PsL1,PsL2)=PairSim(PairPos);
            SimVal(PsL2,PsL1)=PairSim(PairPos);
        end
    end
    %find existing similarity values
    USimVal=setdiff(unique(SimVal),0);
    if length(USimVal)==1
        MultiPs{USimVal}{GeneL}=CurrPs;
    else
        %find the multiple ps with the largest SimVal
        UsedPos=[];
        for SimL=length(USimVal):-1:1
            %clear ps already used
            if ~isempty(UsedPos)
                SimVal(UsedPos,:)=0;
                SimVal(:,UsedPos)=0;
            end
            SimValPos=find(SimVal==USimVal(SimL));
            if ~isempty(SimValPos)
                %tag the pair of probe sets with the current similarity
                Bin=zeros(CurrPsNb);
                Bin(SimValPos)=1;
                RowNb=sum(Bin);
                ColNb=sum(Bin,2);
                if USimVal(SimL)==1
                    %similarity = 0 => all probe set kept
                    ComPos=union(find(RowNb),find(ColNb));
                    MultiPs{USimVal(SimL)}{GeneL}=CurrPs(ComPos);
                    UsedPos=[UsedPos,ComPos];
                else
                    Nb=setdiff(unique(union(RowNb,ColNb)),0);
                    for NbL=length(Nb):-1:1
                        RowPos=find(RowNb==Nb(NbL));
                        ColPos=find(ColNb==Nb(NbL));
                        ComPos=setdiff(intersect(RowPos,ColPos),UsedPos);
                        if length(ComPos)>0
                            MultiPs{USimVal(SimL)}{GeneL}=CurrPs(ComPos);
                            UsedPos=[UsedPos,ComPos];
                        end
                    end
                end
            end
        end
    end
    end
end

%construct definitive gene lists with their associate list of merged probe set
AllGene=cell(length(MultiGeneId),1);
for SimL=1:6    
    for GeneL=1:length(MultiGeneId)
        AllGene{GeneL}=[AllGene{GeneL};MultiPs{SimL}{GeneL}];
    end
end
for SimL1=6:-1:1
    Multi{SimL1}=cell(length(MultiGeneId),1);
    Rand{SimL1}=cell(length(MultiGeneId),1);
    for GeneL=1:length(MultiGeneId)
        if ~isempty(AllGene{GeneL})
            Multi{SimL1}{GeneL}=MultiPs{SimL1}{GeneL};
            if SimL1<6
                for SimL2=SimL1:6
                    Multi{SimL1}{GeneL}=[Multi{SimL1}{GeneL};Multi{SimL2}{GeneL}];
                end
                Multi{SimL1}{GeneL}=unique(Multi{SimL1}{GeneL});
            end
            Rand{SimL1}{GeneL}=setdiff(AllGene{GeneL},Multi{SimL1}{GeneL});
        end
    end
end
%clear empty slots
ClearPos=[];
for GeneL=1:length(MultiGeneId)
    if isempty(AllGene{GeneL})
        ClearPos(end+1)=GeneL;
    end
end
AllGene(ClearPos)=[];
AllGeneName=MultiGeneName;
AllGeneName(ClearPos)=[];
for SimL=1:6
    ClearPos=[];
    for GeneL=1:length(MultiGeneId)
        if isempty(Multi{SimL}{GeneL})
            ClearPos(end+1)=GeneL;
        end
    end
    Multi{SimL}(ClearPos)=[];
    MGeneName{SimL}=MultiGeneName;
    MGeneName{SimL}(ClearPos)=[];    
end
for SimL=1:6
    ClearPos=[];
    for GeneL=1:length(MultiGeneId)
        if isempty(Rand{SimL}{GeneL})
            ClearPos(end+1)=GeneL;
        end
    end
    Rand{SimL}(ClearPos)=[];
    RGeneName{SimL}=MultiGeneName;
    RGeneName{SimL}(ClearPos)=[];    
end

AllGeneName=upper(AllGeneName);
for SimL=1:6
    MGeneName{SimL}=upper(MGeneName{SimL});
    RGeneName{SimL}=upper(RGeneName{SimL});    
end




GsNb=length(GenePos);
DataNb=PointNb;
UMultiPsNb=[];
MaxEs={};
MinEs={};
MaxEsPos={};
MinEsPos={};
MeanEs={};
MeanEsPos={};
MeanPerc1Nb={};
MeanPerc5Nb={};
MaxPerc1Nb={};
MaxPerc5Nb={};
MinPerc1Nb={};
MinPerc5Nb={};

for GsL=1:GsNb
    GsL
    %indicate position of the genes of the curren gene set in the ordered list of genes
    %existing in each similarity level list of genes
    CurrGene=Gene{GsL};
    UMultiPsNb(GsL,1)=length(intersect(UniGeneName,CurrGene));
    SumGene=0;
    for SimL=1:6
        Gs{SimL}=zeros(length(Uni)+length(Multi{SimL})+length(Rand{SimL}),1);
        SumGene=SumGene+length(Gs{SimL});
        [temp,UniPos,temp]=intersect(UniGeneName,CurrGene);
        Gs{SimL}(UniPos)=1;
        [temp,MultiPos,temp]=intersect(MGeneName{SimL},CurrGene);
        UMultiPsNb(GsL,SimL+1)=length(MultiPos);
        Gs{SimL}(MultiPos+length(Uni))=1;
        [temp,RandPos,temp]=intersect(RGeneName{SimL},CurrGene);
        Gs{SimL}(RandPos+length(Uni)+length(Multi{SimL}))=1;
    end
    if GsL==1
        %construct p-value table
        GeneNb=round(SumGene/6);
        Values=zeros(GeneNb,1);
        Perc1=[];
        Perc5=[];
        Range=[[10:10:200],[250:50:1000],[1250:250:5000]];
        for GsNbL=Range
            Es=zeros(1000,1);
            EsPos=zeros(1000,1);
            for RandL=1:1000
                RandPos=ceil(rand(GsNbL,1)*GeneNb);
                CurrValues=Values;
                CurrValues(RandPos)=1;
                [Es(RandL),EsPos(RandL)]=ES_RUNNING_SUM(CurrValues);
            end
            Es=sort(Es);
            Perc1(end+1,1)=Es(999);
            Perc5(end+1,1)=Es(950);
        end
        Perc5=interp1(Range,Perc5,[10:5000]);
        Perc1=interp1(Range,Perc1,[10:5000]);
    end




    MaxEs{GsL}=[];
    MinEs{GsL}=[];
    MaxEsPos{GsL}=[];
    MinEsPos{GsL}=[];
    MeanEs{GsL}=[];
    MeanEsPos{GsL}=[];
    MinPerc1Nb{GsL}=zeros(DataNb,6);
    MinPerc5Nb{GsL}=zeros(DataNb,6);
    MaxPerc1Nb{GsL}=zeros(DataNb,6);
    MaxPerc5Nb{GsL}=zeros(DataNb,6);
    MeanPerc1Nb{GsL}=zeros(DataNb,6);
    MeanPerc5Nb{GsL}=zeros(DataNb,6);
    for DataL=1:DataNb
    %for DataL=1:50
        %recalculate ranks according to the merging schemes indicated in each similarity
        %level
        CurrData=DataRanks(:,Point(DataL));
        for SimL=1:6
            MinSignal=zeros(length(Uni)+length(Multi{SimL})+length(Rand{SimL}),1);
            MaxSignal=MinSignal;
            MeanSignal=MinSignal;
            MinSignal(1:length(Uni))=CurrData(Uni);
            MeanSignal(1:length(Uni))=CurrData(Uni);
            MaxSignal(1:length(Uni))=CurrData(Uni);
            Offset=length(Uni);
            for GeneL=1:length(Multi{SimL})
                CurrSignal=CurrData(Multi{SimL}{GeneL});
                MinSignal(GeneL+Offset)=min(CurrSignal);
                MaxSignal(GeneL+Offset)=max(CurrSignal);
                MeanSignal(GeneL+Offset)=mean(CurrSignal);
            end
            if SimL>1
                Offset=length(Uni)+length(Multi{SimL});
                for GeneL=1:length(Rand{SimL})
                    CurrSignal=CurrData(Rand{SimL}{GeneL});
                    PsPos=randperm(length(CurrSignal));
                    CurrSignal=CurrSignal(PsPos(1));
                    MinSignal(GeneL+Offset)=CurrSignal;
                    MaxSignal(GeneL+Offset)=CurrSignal;
                    MeanSignal(GeneL+Offset)=CurrSignal;
                end
            end
            %sort Data in descending order to calculate ES
            [temp SortIndex]=sort(MinSignal,'descend');
            [temp DataIndex]=sort(SortIndex);
            CurrGs=Gs{SimL}(DataIndex);
            [MinEs{GsL}(DataL,SimL),MinEsPos{GsL}(DataL,SimL)]=ES_RUNNING_SUM(CurrGs);
            if MinEs{GsL}(DataL,SimL)>=Perc1(length(find(CurrGs)))
                MinPerc1Nb{GsL}(DataL,SimL)=2;
            end
             if MinEs{GsL}(DataL,SimL)>=Perc5(length(find(CurrGs)))
                MinPerc5Nb{GsL}(DataL,SimL)=2;
            end
            [temp SortIndex]=sort(MeanSignal,'descend');
            [temp DataIndex]=sort(SortIndex);
            CurrGs=Gs{SimL}(DataIndex);
            [MaxEs{GsL}(DataL,SimL),MaxEsPos{GsL}(DataL,SimL)]=ES_RUNNING_SUM(CurrGs);
             if MaxEs{GsL}(DataL,SimL)>=Perc1(length(find(CurrGs)))
                MaxPerc1Nb{GsL}(DataL,SimL)=4;
             end
             if MaxEs{GsL}(DataL,SimL)>=Perc5(length(find(CurrGs)))
                MaxPerc5Nb{GsL}(DataL,SimL)=4;
            end

            [temp SortIndex]=sort(MaxSignal,'descend');
            [temp DataIndex]=sort(SortIndex);
            CurrGs=Gs{SimL}(DataIndex);
            [MeanEs{GsL}(DataL,SimL),MeanEsPos{GsL}(DataL,SimL)]=ES_RUNNING_SUM(CurrGs);
             if MeanEs{GsL}(DataL,SimL)>=Perc1(length(find(CurrGs)))
                MeanPerc1Nb{GsL}(DataL,SimL)=1;
             end
             if MeanEs{GsL}(DataL,SimL)>=Perc5(length(find(CurrGs)))
                MeanPerc5Nb{GsL}(DataL,SimL)=1;
            end
        end
    end

%      h=figure;
%      set(gcf,'color',[1,1,1])
%      set(h,'name',sprintf('GsL%u',GsL))
%      subplot(1,2,1)
%      image((MaxPerc1Nb{GsL}+MinPerc1Nb{GsL}+MeanPerc1Nb{GsL})*11)
%      subplot(1,2,2)
%      image((MaxPerc5Nb{GsL}+MinPerc5Nb{GsL}+MeanPerc5Nb{GsL})*11)


end




'stop'
save gsea_multipleps

%correction MinGseaPos and MinGsea







%% ES RUNNING SUM
    function [Es,EsPos]=ES_RUNNING_SUM(Values)
        %one running sum from start to end
        Es=zeros(size(Values));
        HitPos=find(Values);
        HitNb=length(HitPos);
        MissPos=find(Values==0);
        MissNb=length(MissPos);
        if MissNb==0
            Es(:)=-1/HitNb;
        elseif HitNb==0
            Es(:)=-1/MissNb;
        else
            Es(MissPos)=-1/MissNb;
            Es(HitPos)=1/HitNb;
        end
        Es=cumsum(Es);
        [Es EsPos]=max(Es);
        EsPos=EsPos(1);




      

        
        
%         for DataL=1:DataNb
%         %sort Data in descending order to calculate ES
%         [temp DataSortIndex]=sort(DataRanks(:,DataL),'descend');
%         for SimL=1:6
%             % use probe set with the smallest signal
%             MinPos=zeros(PsNb,1);
%             MinPos(DataSortIndex(UniPs))=1;
%             for GeneL=1:length(GsPsRanks)
%                 [CurrPs CurrPos,temp]=intersect(GsPsRanks{GeneL},MultiPs{SimL}{GeneL});
%                 if ~isempty(CurrPs)
%                     [temp PsMinPos]=min(DataSortIndex(CurrPs));
%                     MinPos(GsPsRanks{GeneL}(CurrPos(PsMinPos)))=1;
%                 end
%             end
%             [MinEs(DataL,SimL)),MinEsPos(DataL,SimL))]=ES_RUNNING_SUM(MinPos);
% 
%             % use probe set with the smallest signal
%             MaxPos=zeros(PsNb,1);
%             MaxPos(DataSortIndex(UniPs))=1;
%             for GeneL=1:length(GsPsRanks)
%                 [CurrPs CurrPos,temp]=intersect(GsPsRanks{GeneL},MultiPs{SimL}{GeneL});
%                 if ~isempty(CurrPs)                
%                     [temp PsMaxPos]=max(DataSortIndex(CurrPs));
%                     MaxPos(GsPsRanks{GeneL}(CurrPos(PsMaxPos)))=1;
%                 end
%             end
%             [MaxEs(DataL,SimL)),MaxEsPos(DataL,SimL))]=ES_RUNNING_SUM(MaxPos);
%             %random choice
%             for RandL=1:100
%                 RandPos=zeros(PsNb,1);
%                 RandPos(DataSortIndex(UniPs))=1;
%                 for GeneL=1:length(GsPsRanks)
%                     [CurrPs CurrPos,temp]=intersect(GsPsRanks{GeneL},MultiPs{SimL}{GeneL});
%                     if ~isempty(CurrPs)
%                         PsPos=randperm(length(CurrPs));
%                         RandPos(GsPsRanks{GeneL}(CurrPos(PsPos(1))))=1;
%                     end
%                 end
%                 [RandEs{DataL}(SimL,RandL),RandEsPos{DataL}(SimL,RandL)]=ES_RUNNING_SUM(RandPos);
%             end
%         end
% 
%         h=figure;
%         set(gca,'color',[1,1,1])
%         set(h,'name',sprintf('ES GS %u Data %u',GsL,DataL))
%         for SimL=1:6
%             subplot(sort(2,3,SimL))
%             hold on
%             plot(sort(sort(RandEs{DataL}(SimL,:)))
%             h=line([0,100],[MinEs(DataL,SimL)),MinEs(DataL,SimL))]);
%             set(h,'color','b')
%             h=line([0,100],[MaxEs(DataL,SimL)),MaxEs(DataL,SimL))]);
%             set(h,'color','r')
%         end
% 
%         h=figure;
%         set(gca,'color',[1,1,1])
%         set(h,'name',sprintf('ESPos GS %u Data %u',GsL,DataL))
%         for SimL=1:6
%             subplot(sort(2,3,SimL))
%             hold on
%             plot(sort(sort(RandEsPos{DataL}(SimL,:)))
%             h=line([0,100],[MinEsPos(DataL,SimL)),MinEsPos(DataL,SimL))]);
%             set(h,'color','b')
%             h=line([0,100],[MaxEsPos(DataL,SimL)),MaxEsPos(DataL,SimL))]);
%             set(h,'color','r')
%         end
% 
%     end