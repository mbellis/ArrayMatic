%==========================%
% FUNCTION NET_COPYCOMBINE %
%==========================%


% NET_NETCOMBINE make a decision to construct correspondance between several set of chips
% using the results of NET_NETCOMP

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


%net_copycombine([27,8],2,'_1p',7,15,1,1,[60,60],[164,228],'m5n123_m27n164_combinedps_corr60')




function net_copycombine(ChipRank,FirstChipPos,Postfix,StartNet,NetNb,ProbeNbLimit,PvCorr,CompCorr,TestNet,CombinedFile)
global K


ChipNb=length(ChipRank);
if ~isfield(K.dir,'mldata')
    K.dir.mldata='/home/mbellis/sosma/data/psawn/mldata';
    K.dir.pydata='/home/mbellis/sosma/data/psawn/pydata';
    K.dir.rawdata='/home/mbellis/sosma/data/psawn/rawdata';
    K.dir.point='/home/mbellis/sosma/raydata/mlf/point';
end

for ChipL=1:ChipNb
    CurrChipPos=strmatch(sprintf('m%u',ChipRank(ChipL)),K.chip.myName,'exact');
    Species{ChipL}=K.chip.species{CurrChipPos};
    PsNb(ChipL)=K.chip.probesetNb(CurrChipPos);
end

if FirstChipPos==1
    ChipPos(1)=1;
    ChipPos(2)=2;
else
    ChipPos(1)=2;
    ChipPos(2)=1;
end

cd(K.dir.chip)
FileName=sprintf('m%u_m%u_commonps.mat',min(ChipRank),max(ChipRank));
if exist(FileName,'file')
    load(FileName)
    if ChipRank(1)>ChipRank(2)
        Temp=ComPsRank;
        ComPsRank(:,1)=ComPsRank(:,2);
        ComPsRank(:,2)=Temp(:,1);
        clear Temp
    end
    TransPsRank=zeros(max(ComPsRank(:,1)),2);
    TransPsRank(ComPsRank(:,1),1)=ComPsRank(:,1);
    TransPsRank(ComPsRank(:,1),2)=ComPsRank(:,2);
end

cd(K.dir.chip)
eval(sprintf('load %s',CombinedFile))
Sim=[];

%load PsPair
DirName=fullfile(K.dir.mldata,Species{2},sprintf('m%u%s',ChipRank(2),Postfix));
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
    =textread(fullfile('.',sprintf('m%u_n%u_netnb%u_probenb%u_pvcorr%u_pspair.txt',ChipRank(2),StartNet,NetNb,ProbeNbLimit,PvCorr)),...
    [repmat('%s',1,4),repmat('%u',1,16),'%c',repmat('%u',1,12),repmat('%d',1,3),repmat('%u',1,20),repmat('%d',1,3),repmat('%u',1,8)],'delimiter','\t');
cd(K.dir.chip)

ModifFlag=[0,0];
NetPs={};
for ChipL=1:2
    if exist(sprintf('m%u_ps2net.mat',ChipRank(2)))
        ModifFlag(ChipL)=1;
        load(sprintf('m%u_ps2net.mat',ChipRank(2)))
        NetPs{ChipL}=Ps2Net;
    end
end

%search for similarity between pairs of alternative probe sets
PairSim=zeros(length(PsPair),1);
for PairL=1:length(PsPair)
    CurrSim=PsPair(PairL,3:7);
    CurrSim=find(CurrSim);
    if isempty(CurrSim)
        PairSim(PairL)=1;
    else
        PairSim(PairL)=CurrSim(end)+1;
    end
end

%update PsRanks
for PsL=1:length(PsRanks)
    PsRanks(PsL,ChipPos(2))=TransPsRank(PsRanks(PsL,ChipPos(1)),2);
end

%eval(sprintf('save m%u_combinedps_corr%u.mat Pv PsRank Sim GeneName GeneId',ChipRank(ChipL),CompCorr))

%add probe sets with similarity equal to 75 or 100% to be merged
UGeneName=unique(GeneName);

PsRank=cell(length(PsRanks),1);
MultiNb=zeros(length(PsRanks),1);
for PsL1=1:length(PsRanks)
    Pos=find(PsPair(:,1)==PsRanks(PsL1,ChipPos(2))|PsPair(:,2)==PsRanks(PsL1,ChipPos(2)));
    if isempty(Pos)
        PsRank{PsL1}=PsRanks(PsL1,ChipPos(2));
    else
        SelPos=find(PairSim(Pos)>=5);
        if ~isempty(SelPos)
            SelPsRank=unique([PsPair(Pos(SelPos),1);PsPair(Pos(SelPos),2)]);
            UsedPsRank=[];
            for PsL2=1:length(SelPsRank)
                if ~isempty(find(PsRanks(:,ChipPos(2))==SelPsRank(PsL2)))
                    UsedPsRank(end+1,1)=SelPsRank(PsL2);
                end
            end
            SelPsRank=setdiff(SelPsRank,UsedPsRank);
            if isempty(SelPsRank)
                PsRank{PsL1}=PsRanks(PsL1,ChipPos(2));
            else
                PsRank{PsL1}=[PsRanks(PsL1,ChipPos(2));SelPsRank];
                MultiNb(PsL1)=length(PsRank{PsL1});
            end
        else
            PsRank{PsL1}=PsRanks(PsL1,ChipPos(2));
        end
    end
end

%construct similarity
Sim=[];
for GeneL=1:length(UGeneName)
    GenePos=strmatch(UGeneName{GeneL},GeneName,'exact');
    if length(GenePos)>1
        CurrPsRank=unique(PsRanks(GenePos,ChipPos(2)));
        for PsL1=1:length(CurrPsRank)-1
            for PsL2=PsL1+1:length(CurrPsRank)
                PairPos=find(PsPair(:,1)==CurrPsRank(PsL1)&PsPair(:,2)==CurrPsRank(PsL2));
                if ~isempty(PairPos)
                    Sim(end+1,:)=[CurrPsRank(PsL1),CurrPsRank(PsL2),PairSim(PairPos)];
                end
            end
        end
    end
end
cd(K.dir.chip)
eval(sprintf('save m%u_%s_combinedps_corr%u.mat Pv PsRank PsRanks Sim GeneName GeneId',ChipRank(2),regexp(CombinedFile,'.+(?=combined)','match');CompCorr(2)))



PlotSize=floor(length(PsRanks)/10);
Corr={};
for LoopL=1:3
    h=figure;
    set(gcf,'color',[1,1,1])
    if LoopL==1
        PsPos=1:PlotSize;
    elseif LoopL==2
        PsPos=(PlotSize*5)+1:PlotSize*6;
    else
        PsPos=length(PsRanks)-PlotSize+1:length(PsRanks);
    end
    for ChipL=1:2
        cd(fullfile(K.dir.net,sprintf('m%u',ChipRank(ChipL)),sprintf('n%u',TestNet(ChipL))))
        if ModifFlag(ChipL)==0
            Corr{ChipL}=load_data(sprintf('c_m%u_n%u.4mat',ChipRank(ChipL),TestNet(ChipL)),'./',PsNb(ChipL),...
                PsNb(ChipL),'uint8','ieee-le',PsRanks(PsPos,ChipPos(ChipL)),PsRanks(PsPos,ChipPos(ChipL)));
        else
            Corr{ChipL}=load_data(sprintf('c_m%u_n%u.4mat',ChipRank(ChipL),TestNet(ChipL)),'./',PsNb(ChipL),...
                PsNb(ChipL),'uint8','ieee-le',NetPs{ChipL}(PsRanks(PsPos,ChipPos(ChipL))),NetPs{ChipL}(PsRanks(PsPos,ChipPos(ChipL))));
        end
        subplot(1,2,ChipL)
        image(Corr{ChipL})
        %h=pcolor(double(Corr));
        %set(h,'linestyle','none')
        title(sprintf('m%un%u',ChipRank(ChipL),TestNet(ChipL)))
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
