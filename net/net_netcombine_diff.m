%==============================%
% FUNCTION NET_NETCOMBINE_DIFF %
%==============================%


% NET_NETCOMBINE make a decision to construct correspondance between several set of chips
% using the results of NETCOMP_DIFFSP_DO

%INPUT PARAMETERS

% 1 ChipRank
% 2 Postfix
% 3 StartNet
% 4 NetNb
% 5 ProbeNbLimit
% 6 PvCorr
% 7 CompCorr
% 8 TestNet

%EXTERNAL FILES

%OUTPUT PARAMETERS

%%% !!!! ADD IdemFlag for m2 !!!!

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



% net_netcombine_diff([91,92],[3,4],60,[0,0])
% net_netcombine_diff([93,6],[3,63],60,[0,1],{'','_1p'},[0,6],[0,15],[0,1],[0,1])
% net_netcombine_diff([92,6],[4,63],60,[0,1],{'','_1p'},[0,6],[0,15],[0,1],[0,1])

% net_netcombine_diff([2,5],[80,123],60,[1,1],{'_1p','_1p'},[12,17],[21,35],[1,1],[1,1])
% net_netcombine_diff([2,6],[80,63],60,[1,1],{'_1p','_1p'},[12,6],[21,15],[1,1],[1,1])
% net_netcombine_diff([2,8],[80,228],60,[1,1],{'_1p','_1p'},[12,7],[21,15],[1,1],[1,1])
% net_netcombine_diff([2,27],[80,164],60,[1,1],{'_1p','_12to22_1p'},[12,12],[21,11],[1,1],[1,1])

% net_netcombine_diff([3,5],[86,123],60,[1,1],{'_1p','_1p'},[6,17],[18,35],[1,1],[1,1])
% net_netcombine_diff([3,6],[86,63],60,[1,1],{'_1p','_1p'},[6,6],[18,15],[1,1],[1,1])
% net_netcombine_diff([3,8],[86,228],60,[1,1],{'_1p','_1p'},[6,7],[18,15],[1,1],[1,1])
% net_netcombine_diff([3,27],[86,164],60,[1,1],{'_1p','_12to22_1p'},[6,12],[18,11],[1,1],[1,1])

% net_netcombine_diff([5,6],[123,63],60,[1,1],{'_1p','_1p'},[17,6],[35,15],[1,1],[1,1])
% net_netcombine_diff([6,8],[63,228],60,[1,1],{'_1p','_1p'},[6,7],[15,15],[1,1],[1,1])
% net_netcombine_diff([6,27],[63,164],60,[1,1],{'_1p','_2to22_1p'},[6,2],[15,21],[1,1],[1,1])




function net_netcombine_diff(ChipRank,NetRank,CompCorr,PsPairFlag,varargin)
global K

if nargin==9
    Postfix=varargin{1};
    StartNet=varargin{2};
    NetNb=varargin{3};
    ProbeNbLimit=varargin{4};
    PvCorr=varargin{5};
end

if ~isfield(K.dir,'mldata')
    K.dir.mldata='/home/mbellis/sosma/data/psawn/mldata';
    K.dir.pydata='/home/mbellis/sosma/data/psawn/pydata';
    K.dir.rawdata='/home/mbellis/sosma/data/psawn/rawdata';
    K.dir.point='/home/mbellis/sosma/raydata/mlf/point';
end

%recover chip information
ChipNb=length(ChipRank);
for ChipL=1:ChipNb
    ChipPos(ChipL)=strmatch(sprintf('m%u',ChipRank(ChipL)),K.chip.myName,'exact');
    Species{ChipL}=K.chip.species{ChipPos(ChipL)};
    cd(K.dir.chip)
    PsNb(ChipL)=K.chip.probesetNb(ChipPos(ChipL));
end

%load correspondance between probe set list order and network order
NetPs={};
ModifFlag=[0,0];
cd(K.dir.chip)
for ChipL=1:2
    if exist(sprintf('m%u_ps2net.mat',ChipRank(ChipL)))
        ModifFlag(ChipL)=1;
        load(sprintf('m%u_ps2net.mat',ChipRank(ChipL)))
        NetPs{ChipL}=Ps2Net;
    end
end

%load correspondance between probe sets found in a previous step
cd(K.dir.chip)
eval(sprintf('load m%u_m%u_combineps_corr%u.mat',ChipRank(1),ChipRank(2),CompCorr))


PsPair={};
PairSim={};
for ChipL=1:2
    if PsPairFlag(ChipL)
        %load PsPair
        PairGeneId={};
        PairGeneName={};
        PsName1={};
        PsName2={};
        
        DirName=fullfile(K.dir.mldata,Species{ChipL},sprintf('m%u%s',ChipRank(ChipL),Postfix{ChipL}));
        cd(DirName)
        [PairGeneId{ChipL},PairGeneName{ChipL},PsName1{ChipL},PsName2{ChipL},PsPair{ChipL}(:,1),PsPair{ChipL}(:,2),PsPair{ChipL}(:,3),PsPair{ChipL}(:,4),PsPair{ChipL}(:,5),...
            PsPair{ChipL}(:,6),PsPair{ChipL}(:,7),PsPair{ChipL}(:,8),PsPair{ChipL}(:,9),PsPair{ChipL}(:,10),PsPair{ChipL}(:,11),PsPair{ChipL}(:,12),PsPair{ChipL}(:,13),PsPair{ChipL}(:,14),...
            PsPair{ChipL}(:,15),PsPair{ChipL}(:,16),Mark,PsPair{ChipL}(:,17),PsPair{ChipL}(:,18),PsPair{ChipL}(:,19),PsPair{ChipL}(:,20),...
            PsPair{ChipL}(:,21),PsPair{ChipL}(:,22),PsPair{ChipL}(:,23),PsPair{ChipL}(:,24),PsPair{ChipL}(:,25),PsPair{ChipL}(:,26),PsPair{ChipL}(:,27),...
            PsPair{ChipL}(:,28),PsPair{ChipL}(:,29),PsPair{ChipL}(:,30),PsPair{ChipL}(:,31),PsPair{ChipL}(:,32),PsPair{ChipL}(:,33),PsPair{ChipL}(:,34),...
            PsPair{ChipL}(:,35),PsPair{ChipL}(:,36),PsPair{ChipL}(:,37),PsPair{ChipL}(:,38),PsPair{ChipL}(:,39),PsPair{ChipL}(:,40),PsPair{ChipL}(:,41),...
            PsPair{ChipL}(:,42),PsPair{ChipL}(:,43),PsPair{ChipL}(:,44),PsPair{ChipL}(:,45),PsPair{ChipL}(:,46),PsPair{ChipL}(:,47),PsPair{ChipL}(:,48),...
            PsPair{ChipL}(:,49),PsPair{ChipL}(:,50),PsPair{ChipL}(:,51),PsPair{ChipL}(:,52),PsPair{ChipL}(:,53),PsPair{ChipL}(:,54),PsPair{ChipL}(:,55),...
            PsPair{ChipL}(:,56),PsPair{ChipL}(:,57),PsPair{ChipL}(:,58),PsPair{ChipL}(:,59),PsPair{ChipL}(:,60),PsPair{ChipL}(:,61),PsPair{ChipL}(:,62)]...
            =textread(fullfile('.',sprintf('m%u_n%u_netnb%u_probenb%u_pvcorr%u_pspair.txt',ChipRank(ChipL),StartNet(ChipL),NetNb(ChipL),ProbeNbLimit(ChipL),PvCorr(ChipL))),...
            [repmat('%s',1,4),repmat('%u',1,16),'%c',repmat('%u',1,12),repmat('%d',1,3),repmat('%u',1,20),repmat('%d',1,3),repmat('%u',1,8)],'delimiter','\t');
           
        PairSim{ChipL}=zeros(length(PsPair{ChipL}),1);
        for PairL=1:length(PsPair{ChipL})
            CurrSim=PsPair{ChipL}(PairL,3:7);
            CurrSim=find(CurrSim);
            if isempty(CurrSim)
                %similarity = 0% (1)
                PairSim{ChipL}(PairL)=1;
            else
                %similarity = 1% (2), 25%(3), 50%(4), 75%(5), 100%(6)
                PairSim{ChipL}(PairL)=CurrSim(end)+1;
            end
        end
    else
        PsPair{ChipL}=Sim{ChipL}(:,[1:2]);
        PairSim{ChipL}=Sim{ChipL}(:,3);
    end
end




GeneName={};
GeneId={};
PsRanks=[];
StartPos=0;
Pv=[];
UGeneId=unique(SelGeneId{1});
%recover similarity matrix on each chip for all probe sets targeting the same gene
%construct one-to-one probe set correspondance
Sim{1}=cell(length(UGeneId),1);
Sim{2}=cell(length(UGeneId),1);
SimPbNb=ones(length(UGeneId),1)*-1;
for GeneL=1:length(UGeneId)
    GenePos=strmatch(UGeneId{GeneL},SelGeneId{1},'exact');
    %only one targeting probe set
    if length(GenePos)==1
        GeneName{end+1,1}=SelGeneName{1}{GenePos};
        GeneId{end+1,1}=SelGeneId{1}{GenePos};
        PsRanks(end+1,:)=SelRanks(GenePos,:);
        EndPos=length(GeneName);
        Pv(end+1,1)=SelVal(GenePos,1);
        
    else
        %several targeting probe sets => select only one corresponding probe set in the
        %second chip
        %probe set ranks for the first chip
        CurrPsRank=unique(SelRanks(GenePos,1));
        %p-values
        CurrVal=SelVal(GenePos,:);
        for PsL=1:length(CurrPsRank)
            CurrPos=find(SelRanks(GenePos,1)==CurrPsRank(PsL)&SelRanks(GenePos,1)>0);
            if ~isempty(CurrPos)
                MinPos=find(CurrVal(CurrPos,2)==1);
                if isempty(MinPos)
                    MinPos=1;
                else
                    if length(MinPos)>1
                        MinPos=MinPos(1);
                    end
                end
                GeneName{end+1,1}=SelGeneName{1}{GenePos(CurrPos(MinPos))};
                GeneId{end+1,1}=SelGeneId{1}{GenePos(1)};
                PsRanks(end+1,:)=[CurrPsRank(PsL),SelRanks(GenePos(CurrPos(MinPos)),2)];
                Pv(end+1,1)=CurrVal(CurrPos(MinPos),1);
                %do not reuse current probe set of the second chip
                ClearPos=find(SelRanks(GenePos,2)==PsRanks(end,2));
                SelRanks(GenePos(ClearPos),1)=-SelRanks(GenePos(ClearPos),1);
                SelRanks(GenePos(ClearPos),2)=-SelRanks(GenePos(ClearPos),2);
            end
        end
        %control similarity for each possible pair of aternative probe sets targeting the
        %current gene
        EndPos=length(GeneName);
        CurrPairNb=EndPos-StartPos;
        if CurrPairNb>1
            for ChipL=1:2
                Sim{ChipL}{GeneL}=ones(CurrPairNb)*-1;
                for PsL1=1:CurrPairNb-1
                    for PsL2=PsL1+1:CurrPairNb
                        try
                            PairPos=find(PsPair{ChipL}(:,1)==min(PsRanks(StartPos+PsL1,ChipL),PsRanks(StartPos+PsL2,ChipL))&PsPair{ChipL}(:,2)==max(PsRanks(StartPos+PsL1,ChipL),PsRanks(StartPos+PsL2,ChipL)));
                            if ~isempty(PairPos)
                                Sim{ChipL}{GeneL}(PsL1,PsL2)=PairSim{ChipL}(PairPos);
                            end
                        catch
                        end
                    end
                end
            end
            %discrepancy between similarities of the same pair of probe set between the two chips
            PbPos=find(((Sim{1}{GeneL}==0|Sim{1}{GeneL}==1)&(Sim{2}{GeneL}==5|Sim{2}{GeneL}==6))|((Sim{2}{GeneL}==0|Sim{2}{GeneL}==1)&(Sim{1}{GeneL}==5|Sim{1}{GeneL}==6)));
            SimPbNb(GeneL)=length(PbPos);
        end
    end
    StartPos=EndPos;
end

PbNbVal=unique(SimPbNb);
fprintf('nb of pairs with discrepancy\n')
[PbNbVal,histc(SimPbNb,PbNbVal)]

[Pv,SortIndex]=sort(Pv);
PsRanks=PsRanks(SortIndex,:);
GeneName=GeneName(SortIndex);
GeneId=GeneId(SortIndex);

%add probe sets with similarity equal to 75 or 100% to be merged

UGeneId=unique(GeneId);
for ChipL=1:2
    PsRank=cell(length(PsRanks),1);
    MultiNb=zeros(length(PsRanks),1);
    for PsL1=1:length(PsRanks)    
        Pos=find(PsPair{ChipL}(:,1)==PsRanks(PsL1,ChipL)|PsPair{ChipL}(:,2)==PsRanks(PsL1,ChipL));
        if isempty(Pos)
            PsRank{PsL1}=PsRanks(PsL1,ChipL);
        else
            %find pairs that contains the current probe set and that have a similarity >=75%
            SelPos=find(PairSim{ChipL}(Pos)>=5);            
            if ~isempty(SelPos)
                SelPsRank=unique([PsPair{ChipL}(Pos(SelPos),1);PsPair{ChipL}(Pos(SelPos),2)]);
                UsedPsRank=[];
                for PsL2=1:length(SelPsRank)
                    if ~isempty(find(PsRanks(:,ChipL)==SelPsRank(PsL2)))
                        UsedPsRank(end+1,1)=SelPsRank(PsL2);
                    end
                end
                SelPsRank=setdiff(SelPsRank,UsedPsRank);
                if isempty(SelPsRank)
                    PsRank{PsL1}=PsRanks(PsL1,ChipL);
                else
                    PsRank{PsL1}=[PsRanks(PsL1,ChipL);SelPsRank];
                    MultiNb(PsL1)=length(PsRank{PsL1});
                end
            else
                PsRank{PsL1}=PsRanks(PsL1,ChipL);
            end
        end
    end
    
    %construct similarity
    Sim=[];
    for GeneL=1:length(UGeneId)
        GenePos=strmatch(UGeneId{GeneL},GeneId,'exact');
        if length(GenePos)>1
            CurrPsRank=unique(PsRanks(GenePos,ChipL));
            for PsL1=1:length(CurrPsRank)-1
                for PsL2=PsL1+1:length(CurrPsRank)
                    PairPos=find(PsPair{ChipL}(:,1)==CurrPsRank(PsL1)&PsPair{ChipL}(:,2)==CurrPsRank(PsL2));
                    if ~isempty(PairPos)
                        Sim(end+1,:)=[GenePos(PsL1),GenePos(PsL2),PairSim{ChipL}(PairPos)];
                    end
                end
            end
        end
    end
    cd(K.dir.chip)
    if ChipL==1
        eval(sprintf('save m%un%u_m%un%u_combinedps_corr%u.mat Pv PsRank PsRanks Sim GeneName GeneId',ChipRank(1),NetRank(1),ChipRank(2),NetRank(2),CompCorr))
    else
        eval(sprintf('save m%un%u_m%un%u_combinedps_corr%u.mat Pv PsRank PsRanks Sim GeneName GeneId',ChipRank(2),NetRank(2),ChipRank(1),NetRank(1),CompCorr))
    end
end

cd(K.dir.chip)
eval(sprintf('load m%utom%u_normcurve OldVal',ChipRank(2),ChipRank(1)));
%PlotSize=max(500,floor(length(PsRanks)/10));
PlotSize=300;
for LoopL=1:3
    h=figure;
    set(gcf,'color',[1,1,1])
    if LoopL==1
        PsPos=1:PlotSize;
    elseif LoopL==2
        PsPos=round(length(PsRanks)/2-PlotSize/2):round(length(PsRanks)/2+PlotSize/2);
    else
        PsPos=round(length(PsRanks)-PlotSize+1):length(PsRanks);
    end
    for ChipL=1:2
        CurrPsRank=PsRanks(PsPos,ChipL);
        if ModifFlag(ChipL)
            CurrPsRank=NetPs{ChipL}(CurrPsRank);
        end
        cd(fullfile(K.dir.net,sprintf('m%u',ChipRank(ChipL)),sprintf('n%u',NetRank(ChipL))))   
        Corr{ChipL}=load_data(sprintf('c_m%u_n%u.4mat',ChipRank(ChipL),NetRank(ChipL)),'./',PsNb(ChipL),...
            PsNb(ChipL),'uint8','ieee-le',CurrPsRank,CurrPsRank);   
        if ChipL==2
            CurrC=Corr{ChipL};
            MaxC=max(max(Corr{ChipL}));
            MaxPos=find(OldVal<=MaxC);
            for ValL=1:min(MaxPos(end),99)
                Pos=find(CurrC>=OldVal(ValL)&CurrC<OldVal(ValL+1));
                Corr{ChipL}(Pos)=ValL;
            end
            if MaxPos(end)<99
                Pos=find(CurrC>=OldVal(MaxPos(end)+1)&CurrC<100);
                Corr{ChipL}(Pos)=MaxPos(end)+1;
            end
        end
        subplot(1,2,ChipL)
        image(Corr{ChipL})
        %h=pcolor(double(Corr{ChipL}));
        %set(h,'linestyle','none')
        title(sprintf('m%un%u',ChipRank(ChipL),NetRank(ChipL)))
    end
    set(gcf,'position',[5        1379        1010         399])
    Corr{1}=Corr{1}(:);
    Corr{2}=Corr{2}(:);
    [Corr{1},SortIndex]=sort(Corr{1});
    figure
    set(gcf,'color',[1,1,1])    
    hold on
    plot(double(Corr{2}(SortIndex))+rand(length(Corr{2}),1)-0.5,'b.','markersize',3)
    plot(Corr{1},'k.')
    set(gca,'box','on')
    set(gca,'ylim',[0,60])
end


