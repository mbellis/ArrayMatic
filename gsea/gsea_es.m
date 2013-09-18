%==========================%
% FUNCTION GSEA_MULTIPLEPS %
%==========================%

% GSEA_ES calculate Enrichment Score andp-value (multiple probe sets are not merged)
%
%INPUT PARAMETERS
%
% 1        Species: species
% 2       ChipRank: chip model rank
% 3  FileName: GSEA lists
% 4  NrBiolRank: rank of set of non redundant biological conditions (P.biol.nr.scoupleIndex{NrBiolRank})
% 5 SingleFlag: uses only probe sets with a one to one relationship with genes (class SS)

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

%gsea_es('mouse',8,'_1p',[7:21],1,1,'c3.tft.v3.1.symbols_gene2psrank')
%gsea_es('mouse',27,'c3.tft.v3.1.symbols_gene2psrank_clique',9,0)
%gsea_es('mouse',27,'c3.tft.v3.1.symbols_gene2psrank_clique',9,1)

function [EsVal,EsPos,Perc1Gs,Perc5Gs,BiolName]=gsea_es(ChipRank,FileName,NrBiolRank,SingleFlag)
global K

%load GSEA list (load GenePos containing for each GSEA list and for several chip model
% the targeting probe sets
if ischar(FileName)
    cd(K.dir.list)
    eval(sprintf('load %s.mat',FileName))
    GsPs=SelPs;
else
    GsPs=FileName;
end

%load eventually PsMatrix to find probe sets of class SS
if SingleFlag
    cd(K.dir.chip)
    [PsId,GeneId,GeneName,Class,PNb]=textread(sprintf('m%u_gene.txt',ChipRank),'%s%s%s%u%u','delimiter','\t');
    SIndex=find(Class==1);    
end
    
    
%load point info
K.dir.point='/home/mbellis/sosma/raydata/mlf/point';
cd(K.dir.point)
eval(sprintf('load m%u_point',ChipRank))
%find Point to be used
Point=[];
BiolName={};
if NrBiolRank>0
    ScoupleIndex=P.biol.nr{1}.scoupleindex{NrBiolRank};
    for BiolL=1:length(ScoupleIndex)
        Point=[Point,P.biol.scouple{1}{ScoupleIndex(BiolL)}(1)];
        BiolName{end+1,1}=P.biol.name{ScoupleIndex(BiolL)};
    end
else
    %use all SCouples
    for BiolL=1:length(P.biol.scoupleindex{1})
        Point=[Point,P.biol.scouple{1}{P.biol.scoupleindex{1}(BiolL)}(1)];
        BiolName{end+1,1}=P.biol.name{P.biol.scoupleindex{1}(BiolL)};
    end
end
PointNb=length(Point);


%load ps ranks
cd (fullfile('/usr/data/net',sprintf('m%u',ChipRank),sprintf('m%u_data',ChipRank)))
load DataRanks


%recover chip information
ChipPos=find(K.chip.rank==ChipRank);
PsNb=K.chip.probesetNb(ChipPos);


GsNb=length(GsPs);
DataNb=PointNb;
EsVal{1}=zeros(GsNb,DataNb);
EsPos{1}=zeros(GsNb,DataNb);
Perc1Gs{1}=zeros(GsNb,DataNb);
Perc5Gs{1}=zeros(GsNb,DataNb);
EsVal{2}=zeros(GsNb,DataNb);
EsPos{2}=zeros(GsNb,DataNb);
Perc1Gs{2}=zeros(GsNb,DataNb);
Perc5Gs{2}=zeros(GsNb,DataNb);
GsPsNb=zeros(GsNb,2);
for GsL=1:GsNb
    GsL
    %indicate position of the genes of the current gene set in the ordered list of genes
    %existing in each similarity level list of genes
    CurrPs=GsPs{GsL};
    GsPsNb(GsL,1)=length(CurrPs);
    if SingleFlag
        CurrPs=intersect(CurrPs,SIndex);
        GsPsNb(GsL,2)=length(CurrPs);
    end
    if GsL==1
        %construct p-value table
        Values=zeros(PsNb,1);
        Perc1=[];
        Perc5=[];
        Range=[[10:10:200],[250:50:1000],[1250:250:5000]];
        for GsNbL=Range
            Es=zeros(1000,1);            
            for RandL=1:1000
                RandPos=ceil(rand(GsNbL,1)*PsNb);
                CurrValues=Values;
                CurrValues(RandPos)=1;
                [CurrEs,tmp]=ES_RUNNING_SUM(CurrValues,0);
                Es(RandL)=CurrEs(1);
            end
            Es=sort(Es);
            Perc1(end+1,1)=Es(999);
            Perc5(end+1,1)=Es(950);
        end
        Perc5=interp1(Range,Perc5,[10:5000]);
        Perc1=interp1(Range,Perc1,[10:5000]);
    end
    if SingleFlag
        CurrPsNb=GsPsNb(GsL,2);
    else
        CurrPsNb=GsPsNb(GsL,1);
    end    
    PercPos=find(Range<=CurrPsNb);
    if ~isempty(PercPos)
        PercPos=PercPos(end);
    else
        PercPos=1;
    end
    if length(CurrPs)>0
        for DataL=1:DataNb
            [temp SortIndex]=sort(DataRanks(:,Point(DataL)));
            %search from the highest signal
            SortIndex=flipud(SortIndex);
            GsBindex=zeros(PsNb,1);
            GsBindex(CurrPs)=1;
            [CurrEsVal,CurrEsPos]=ES_RUNNING_SUM(GsBindex(SortIndex),1);
            for LoopL=1:2
                if abs(CurrEsVal(LoopL))>=Perc1(PercPos)
                    Perc1Gs{LoopL}(GsL,DataL)=1;
                end
                if abs(CurrEsVal(LoopL))>=Perc5(PercPos)
                    Perc5Gs{LoopL}(GsL,DataL)=1;
                end
                EsVal{LoopL}(GsL,DataL)=CurrEsVal(LoopL);
                EsPos{LoopL}(GsL,DataL)=CurrEsPos(LoopL);
            end
        end
    end
end

% cd('/home/mbellis/sosma/sosresult/gsea')
% if SingleFlag
%     eval(sprintf('save %s_singleres EsVal EsPos Perc1Gs Perc5Gs Perc1 Perc5 GsPsNb',FileName))
% else
%     eval(sprintf('save %s_res EsVal EsPos Perc1Gs Perc5Gs Perc1 Perc5 GsPsNb',FileName))
% end


%% ES RUNNING SUM
function [CurrEs,CurrEsPos]=ES_RUNNING_SUM(Values,FlipFlag)
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

%ES from top
[MaxEs EsPos]=max(Es);
EsPos=EsPos(1);
CurrEs(1)=MaxEs;
CurrEsPos(1)=EsPos;
if FlipFlag
    %find Es and EsPos if the values are flipped (ES from bottom)
    [MinEs EsPos]=min(Es);
    EsPos=EsPos(1);
    CurrEs(2)=MinEs;
    CurrEsPos(2)=EsPos;
end




