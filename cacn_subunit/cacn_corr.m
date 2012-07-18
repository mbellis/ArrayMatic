% Ensembl and Gene IDs of calcium channel subunits
ModelRank=27;
NetRanks=[2:22];
NetNb=length(NetRanks);
ModelPos=strmatch(sprintf('m%u',ModelRank),K.chip.myName,'exact');
PsNb=K.chip.probesetNb(ModelPos);
ProjectDir='/home/mbellis/sosma/raydataraw/2011_CACN';

cd(ProjectDir)
if exist(sprintf('cacn_m%u_net%u_netnb%u.m',ModelRank,NetRanks(1),NetNb),'file')
    eval(sprintf('load cacn_m%u_net%u_netnb%u',ModelRank,NetRanks(1),NetNb))
else

    MmEnsId={
        'ENSMUSG00000034656';...
        'ENSMUSG00000004113';...
        'ENSMUSG00000051331';...
        'ENSMUSG00000015968';...
        'ENSMUSG00000004110';...
        'ENSMUSG00000031142';...
        'ENSMUSG00000020866';...
        'ENSMUSG00000024112';...
        'ENSMUSG00000022416';...
        'ENSMUSG00000026407';...
        'ENSMUSG00000040118';...
        'ENSMUSG00000010066';...
        'ENSMUSG00000021991';...
        'ENSMUSG00000041460';...
        'ENSMUSG00000020882';...
        'ENSMUSG00000057914';...
        'ENSMUSG00000003352';...
        'ENSMUSG00000017412';...
        'ENSMUSG00000020722';...
        'ENSMUSG00000019146';...
        'ENSMUSG00000066189';...
        'ENSMUSG00000020723';...
        'ENSMUSG00000040373';...
        'ENSMUSG00000078815';...
        'ENSMUSG00000069806';...
        'ENSMUSG00000053395'};


    HsEnsId={
        'ENSG00000141837';...
        'ENSG00000148408';...
        'ENSG00000151067';...
        'ENSG00000157388';...
        'ENSG00000198216';...
        'ENSG00000102001';...
        'ENSG00000006283';...
        'ENSG00000196557';...
        'ENSG00000100346';...
        'ENSG00000081248';...
        'ENSG00000153956';...
        'ENSG00000007402';...
        'ENSG00000157445';...
        'ENSG00000151062';...
        'ENSG00000067191';...
        'ENSG00000165995';...
        'ENSG00000167535';...
        'ENSG00000182389';...
        'ENSG00000108878';...
        'ENSG00000166862';...
        'ENSG00000006116';...
        'ENSG00000075461';...
        'ENSG00000075429';...
        'ENSG00000130433';...
        'ENSG00000105605';...
        'ENSG00000142408'};

    GeneId={
        'CACNA1A';...
        'CACNA1B';...
        'CACNA1C';...
        'CACNA1D';...
        'CACNA1E';...
        'CACNA1F';...
        'CACNA1G';...
        'CACNA1H';...
        'CACNA1I';...
        'CACNA1S';...
        'CACNA2D1';...
        'CACNA2D2';...
        'CACNA2D3';...
        'CACNA2D4';...
        'CACNB1';...
        'CACNB2';...
        'CACNB3';...
        'CACNB4';...
        'CACNG1';...
        'CACNG2';...
        'CACNG3';...
        'CACNG4';...
        'CACNG5';...
        'CACNG6';...
        'CACNG7';...
        'CACNG8'};

    %search index in Gene List (Genes.name)
    cd('/home/mbellis/sosma/data/psawn/mldata/mouse/m27_2to22')
    eval(sprintf('load m%u_n%u_netnb%u_probenb1_newps_stat_netprc100_pvcorr1',ModelRank,NetRanks(1),NetNb))


    GenePos=zeros(length(MmEnsId),1);
    for GeneL=1:length(MmEnsId)
        CurrGenePos=strmatch(MmEnsId{GeneL},Genes.name,'exact');
        if ~isempty(CurrGenePos)
            GenePos(GeneL)=CurrGenePos;
        end
    end

    %find probe set that target genes
    PsPos=cell(length(MmEnsId),1);
    for GeneL=1:length(MmEnsId)
        if GenePos(GeneL)>0
            PsPos{GeneL}=find(PsMatrix(:,1)==GenePos(GeneL));
        else
            PsPos{GeneL}=[];
        end
    end
    %find targetted GOT
    TargGot=cell(length(MmEnsId),1);
    for GeneL=1:length(MmEnsId)
        TargGot{GeneL}=[];
        if ~isempty(PsPos{GeneL})
            if length(PsPos{GeneL})>1
                for PsL=1:length(PsPos{GeneL})
                    TargGot{GeneL}=[TargGot{GeneL},PsMatrix(PsPos{GeneL}(PsL),10)];
                end
            end
        end
    end
    %multiple ps
    MultiPs=[];
    for GeneL=1:length(MmEnsId)
        if length(PsPos{GeneL})>1
            MultiPs(end+1)=GeneL;
        end
    end

    %exist pivot ?
    for MultiL=1:length(MultiPs)
        MultiPs(MultiL)
        GeneL=MultiPs(MultiL);
        for PsL=1:length(PsPos{GeneL})
            PsMatrix(PsPos{GeneL}(PsL),12)
        end
    end

    GeneRanks=[];
    PsRanks=[];
    for GeneL=1:length(MmEnsId)
        if ~isempty(PsPos{GeneL})
            for PsL=1:length(PsPos{GeneL})
                GeneRanks(end+1,1)=GeneL;
            end
            PsRanks=[PsRanks;PsPos{GeneL}];
        end
    end


    for NetL=1:length(NetRanks)
        DataDir=sprintf('/home/mbellis/array1/sosma/net/m%03u/n%05u',ModelRank,NetRanks(NetL));
        Corr{NetL}=load_data(sprintf('c_m%u_n%u.4mat',ModelRank,NetRanks(NetL)),DataDir,PsNb,PsNb,'uint8','ieee-le',PsRanks);
        Anti{NetL}=load_data(sprintf('a_m%u_n%u.4mat',ModelRank,NetRanks(NetL)),DataDir,PsNb,PsNb,'uint8','ieee-le',PsRanks);
    end


    cd(ProjectDir)
    eval(sprintf('save cacn_m%u_net%u_netnb%u Corr Anti PsRanks GeneRanks MmEnsId GeneId GenePos PsPos PsMatrix',ModelRank,NetRanks(1),NetNb))
end

%CONSTRUCT CACNMAT
CorrMat=uint8(zeros(PsNb,length(PsRanks)));
AntiMat=uint8(zeros(PsNb,length(PsRanks)));
IdemMat=uint8(zeros(PsNb,length(PsRanks)));
for NetL=1:length(NetRanks)
    CurrCorr=Corr{NetL};
    CurrAnti=Anti{NetL};
    %filter Corr and Anti
    CurrCorr(find(CurrCorr<10))=0;
    CurrAnti(find(CurrAnti<10))=0;
    for PsL=1:length(PsRanks)
        %Pos=find(CurrCorr(PsL,:)>=CurrAnti(PsL,:)&CurrCorr(PsL,:)>0);
        Pos=find(CurrCorr(PsL,:)>CurrAnti(PsL,:));
        CorrMat(Pos,PsL)=CorrMat(Pos,PsL)+1;
        Pos=find(CurrCorr(PsL,:)<CurrAnti(PsL,:));
        AntiMat(Pos,PsL)=AntiMat(Pos,PsL)+1;
        Pos=find(CurrCorr(PsL,:)==CurrAnti(PsL,:)&CurrCorr(PsL,:)>0);
        IdemMat(Pos,PsL)=IdemMat(Pos,PsL)+1;
    end
end


%REPRODUCIBILITY IN x NETWORKS
NetNbLimits=[1,5,10,15,16,17,18,19,20,21];

Label={};
for PsL=1:length(PsRanks)
    Label{end+1}=sprintf('%s-%u',GeneId{GeneRanks(PsL)},PsRanks(PsL));
end

for LimitL=1:length(NetNbLimits)
    CurrLimit=NetNbLimits(LimitL);
    for TypeL=1:2
        if TypeL==1
            CurrMat=CorrMat;
        else
            CurrMat=AntiMat;
        end
        CurrMat(find(CurrMat<CurrLimit))=0;
        CurrMat(find(CurrMat))=1;
        ClearPos=find(sum(CurrMat,2)==0);
        CurrMat(ClearPos,:)=[];

        for i=1:size(CurrMat,1)
            CurrMat(i,find(CurrMat(i,:)))=sum(CurrMat(i,:));
        end
        %order
        for PsL=size(CurrMat,2):-1:1
            [Temp SortOrder]=sort(CurrMat(:,PsL));
            SortOrder=flipud(SortOrder);
            CurrMat=CurrMat(SortOrder,:);
        end

        h=figure;
        if TypeL==1
            set(h,'name',sprintf('INTERACTANTS OF CALCIUM CHANEL SUBUNITS - CORR - >=%u NETWORKS',CurrLimit))
        else
            set(h,'name',sprintf('INTERACTANTS OF CALCIUM CHANEL SUBUNITS - CORR - >=%u NETWORKS',CurrLimit))
        end
        set(gcf,'color',[1,1,1])
        MatPlot=CurrMat*floor(255/size(CurrMat,2));
        image(MatPlot');       
        set(gca,'ytick',[1:length(Label)])
        set(gca,'yticklabel',Label)
        %set(gca,'ytick',[])
        if TypeL==1
            title(sprintf('INTERACTANTS OF CALCIUM CHANEL SUBUNITS - CORR - >=%u NETWORKS',CurrLimit))
        else
            title(sprintf('INTERACTANTS OF CALCIUM CHANEL SUBUNITS - ANTI - >=%u NETWORKS',CurrLimit))
        end
        map=colormap;   
        map(1,:)=[1,1,1];
        colormap(map)        
        set(gcf,'units','pixel')
        Position=get(gcf,'position');
        Position(3)=1100;
        Position(4)=500;
        set(gcf,'position',Position);
        set(gca,'tickdir','out')
        
    end
end


%STATISTICS ON DIFFERENT CATEGORIES
%CONSTRUCT SUCCESSIVE LISTS
PsLists=cell(21,1);
PsNbs=zeros(21,1);
for LimitL=1:21
    CurrMat=CorrMat;
    CurrMat(find(CurrMat<LimitL))=0;
    CurrMat(find(CurrMat))=1;
    PsLists{LimitL}=find(sum(CurrMat,2)>0);
    PsNbs(LimitL)=length(PsLists{LimitL});
end

h=figure;
set(gcf,'color',[1,1,1])
plot([1:21],PsNbs,'b-')
hold on
plot([1:21],PsNbs,'r+')
xlabel('Nb of networks where exist a correlation with one of the CACN subunits')
ylabel('Frequency of probe sets correlated with one of the CACN subunits')

% PS CLASSES
PsClass{1}=zeros(21,7);
PsClass{2}=zeros(21,7);
SupPsPos=[];
for LimitL=21:-1:1    
    PsPos=setdiff(PsLists{1},SupPsPos);
    for ClassL=1:7
        PsClass{1}(LimitL,ClassL)=length(find(PsMatrix(PsLists{LimitL},17)==ClassL))*100/length(PsLists{LimitL});
        PsClass{2}(LimitL,ClassL)=length(find(PsMatrix(PsPos,17)==ClassL))*100/length(PsPos);
    end
    SupPsPos=union(SupPsPos,PsLists{LimitL});
end

h=figure;
set(gcf,'color',[1,1,1])
Colors=colors(colormap,7);
hold on
for ClassL=1:7
    plot([1:21],PsClass{1}(:,ClassL),'color',Colors(ClassL,:))
end
for ClassL=1:7
    plot([1:21],PsClass{2}(:,ClassL),'color',Colors(ClassL,:),'linestyle','-.')
end
set(gca,'box','on')
legend('su','sm','mu','mm','cx','hx','0')
xlabel('Nb of networks where exist a correlation with one of the CACN subunits')
ylabel('ps class frequency')

% PROBE NB
PsProbeNb{1}=zeros(21,11);
PsProbeNb{2}=zeros(21,11);
SupPsPos=[];
for LimitL=21:-1:1
    PsPos=setdiff(PsLists{1},SupPsPos);
    for ProbeNbL=1:10
        PsProbeNb{1}(LimitL,ProbeNbL)=length(find(PsMatrix(PsLists{LimitL},5)==ProbeNbL))*100/length(PsLists{LimitL});
        PsProbeNb{2}(LimitL,ProbeNbL)=length(find(PsMatrix(PsPos,5)==ProbeNbL))*100/length(PsPos);
    end
    PsProbeNb{1}(LimitL,11)=length(find(PsMatrix(PsLists{LimitL},5)>=11))*100/length(PsLists{LimitL});
    PsProbeNb{2}(LimitL,11)=length(find(PsMatrix(PsPos,5)>=11))*100/length(PsPos);
    SupPsPos=union(SupPsPos,PsLists{LimitL});
end

h=figure;
set(gcf,'color',[1,1,1])
Colors=colors(colormap,11);
subplot(1,2,1)
hold on
for ProbeNbL=1:10
    plot([1:21],PsProbeNb{1}(:,ProbeNbL),'color',Colors(ProbeNbL,:))
end
for ProbeNbL=1:10
    plot([1:21],PsProbeNb{2}(:,ProbeNbL),'color',Colors(ProbeNbL,:),'linestyle','-.')
end
set(gca,'box','on')
legend('1','2','3','4','5','6','7','8','9','10')
xlabel('Nb of networks where exist a correlation with one of the CACN subunits')
ylabel('probe nb frequency')

subplot(1,2,2)
hold on
for ProbeNbL=11
    plot([1:21],PsProbeNb{1}(:,ProbeNbL),'color',Colors(ProbeNbL,:))
end
for ProbeNbL=11
    plot([1:21],PsProbeNb{2}(:,ProbeNbL),'color',Colors(ProbeNbL,:),'linestyle','-.')
end
set(gca,'box','on')
legend('11')
xlabel('Nb of networks where exist a correlation with one of the CACN subunits')
ylabel('probe nb frequency')

% TRANSCRIPT FREQUENCY
PsTrsNb{1}=zeros(21,11);
PsTrsNb{2}=zeros(21,11);
SupPsPos=[];
for LimitL=21:-1:1
    PsPos=setdiff(PsLists{1},SupPsPos);
    for TrsNbL=0:9
        PsTrsNb{1}(LimitL,TrsNbL+1)=length(find(PsMatrix(PsLists{LimitL},8)==TrsNbL))*100/length(PsLists{LimitL});
        PsTrsNb{2}(LimitL,TrsNbL+1)=length(find(PsMatrix(PsPos,8)==TrsNbL))*100/length(PsPos);
    end
    PsTrsNb{1}(LimitL,11)=length(find(PsMatrix(PsLists{LimitL},8)>=10))*100/length(PsLists{LimitL});
    PsTrsNb{2}(LimitL,11)=length(find(PsMatrix(PsPos,8)>=10))*100/length(PsPos);    
    SupPsPos=union(SupPsPos,PsLists{LimitL});
end

h=figure;
set(gcf,'color',[1,1,1])
Colors=colors(colormap,11);
hold on
for TrsNbL=1:11
    plot([1:21],PsTrsNb{1}(:,TrsNbL),'color',Colors(TrsNbL,:))
end
for TrsNbL=1:11
    plot([1:21],PsTrsNb{2}(:,TrsNbL),'color',Colors(TrsNbL,:),'linestyle','-.')
end

set(gca,'box','on')
legend('0','1','2','3','4','5','6','7','8','9','10')
xlabel('Nb of networks where exist a correlation with one of the CACN subunits')
ylabel('transcript nb frequency')


%FILTERED LIST

GoodPsLists=cell(21,1);
GoodPsNbs=zeros(21,1);
for LimitL=1:21
    CurrMat=CorrMat;
    CurrMat(find(CurrMat<LimitL))=0;
    CurrMat(find(CurrMat))=1;
    GoodPsLists{LimitL}=find(sum(CurrMat,2)>0&(PsMatrix(:,17)==1|PsMatrix(:,17)==3|PsMatrix(:,17)==4|PsMatrix(:,17)==6)&PsMatrix(:,5)>=8);
    GoodPsNbs(LimitL)=length(GoodPsLists{LimitL});
end

h=figure;
set(gcf,'color',[1,1,1])
plot([1:21],PsNbs,'b-')
hold on
plot([1:21],PsNbs,'r+')
plot([1:21],GoodPsNbs,'b-')
plot([1:21],GoodPsNbs,'ro')
xlabel('Nb of networks where exist a correlation with one of the CACN subunits')
ylabel('Frequency of probe sets correlated with one of the CACN subunits')
label('all','filtered')


%GENE LIST
%eload Genes
cd('/home/mbellis/sosma/data/psawn/mldata/mouse/m27_2to22')
eval(sprintf('load m%u_n%u_netnb%u_probenb1_newps_stat_netprc100_pvcorr1',ModelRank,NetRanks(1),NetNb))

GoodGenePosLists=cell(21,1);
for LimitL=1:21
    GoodGenePosLists{LimitL}=PsMatrix(GoodPsLists{LimitL},1);
end

%ONE NETWORK
GenePos=GoodGenePosLists{1};
for i=2:21
GenePos=setdiff(GenePos,GoodGenePosLists{i});
end

GeneNames{1}=Genes.name(GenePos);
KeepPos=strmatch('ENS',GeneNames{1});
GeneNames{1}=GeneNames{1}(KeepPos);

%TWO TO FIVE NETWORKS
GenePos=GoodGenePosLists{2};
for i=[6:21]
GenePos=setdiff(GenePos,GoodGenePosLists{i});
end
GeneNames{2}=Genes.name(GenePos);
KeepPos=strmatch('ENS',GeneNames{2});
GeneNames{2}=GeneNames{2}(KeepPos);

%SIX TO TEN NETWORKS
GenePos=GoodGenePosLists{6};
for i=[11:21]
GenePos=setdiff(GenePos,GoodGenePosLists{i});
end
GeneNames{3}=Genes.name(GenePos);
KeepPos=strmatch('ENS',GeneNames{3});
GeneNames{3}=GeneNames{3}(KeepPos);

%11 TO 15 NETWORKS
GenePos=GoodGenePosLists{11};
for i=[16:21]
GenePos=setdiff(GenePos,GoodGenePosLists{i});
end
GeneNames{4}=Genes.name(GenePos);
KeepPos=strmatch('ENS',GeneNames{4});
GeneNames{4}=GeneNames{4}(KeepPos);

%16 TO 18 NETWORKS
GenePos=GoodGenePosLists{16};
for i=[19:21]
GenePos=setdiff(GenePos,GoodGenePosLists{i});
end
GeneNames{5}=Genes.name(GenePos);
KeepPos=strmatch('ENS',GeneNames{5});
GeneNames{5}=GeneNames{5}(KeepPos);

for i=1:4
    for j=i+1:5
        length(unique(intersect(GeneNames{i},GeneNames{j})))
    end
end

% NO VISIBLE DIFFERENCE AFTER FILTRATION

% MultiPs=zeros(length(PsRanks),1);
% for i=1:length(MultiPs)
%     if length(find(GeneRanks==GeneRanks(i)))>1
%         MultiPs(i)=1;
%     end
% end
% 

% Label={};
% for PsL=1:length(PsRanks)
%     if MultiPs(PsL)==0
%         Label{end+1}=sprintf('%s-%u',GeneId{GeneRanks(PsL)},PsRanks(PsL));        
%     else
%         Label{end+1}=sprintf('!-%s-%u',GeneId{GeneRanks(PsL)},PsRanks(PsL));
%     end
% end
% 
% NetNbLimits=[1,5,10,15,16,17,18,19,20,21];
% 
% for LimitL=1:length(NetNbLimits)
%     CurrLimit=NetNbLimits(LimitL);
%     CurrMat=CorrMat(GoodPsLists{CurrLimit},:);
%     CurrMat(find(CurrMat<CurrLimit))=0;
%     CurrMat(find(CurrMat))=1;
%     for i=1:size(CurrMat,1)
%         CurrMat(i,find(CurrMat(i,:)))=sum(CurrMat(i,:));
%     end
%     %order
%     for PsL=size(CurrMat,2):-1:1
%         [Temp SortOrder]=sort(CurrMat(:,PsL));
%         SortOrder=flipud(SortOrder);
%         CurrMat=CurrMat(SortOrder,:);
%     end
% 
%     h=figure;
%     if TypeL==1
%         set(h,'name',sprintf('INTERACTANTS OF CALCIUM CHANEL SUBUNITS - CORR - >=%u NETWORKS',CurrLimit))
%     else
%         set(h,'name',sprintf('INTERACTANTS OF CALCIUM CHANEL SUBUNITS - CORR - >=%u NETWORKS',CurrLimit))
%     end
%     set(gcf,'color',[1,1,1])
%     MatPlot=CurrMat*floor(255/size(CurrMat,2));
%     image(MatPlot');
%     set(gca,'ytick',[1:length(Label)])
%     set(gca,'yticklabel',Label)
%     %set(gca,'ytick',[])
%     title(sprintf('INTERACTANTS OF CALCIUM CHANEL SUBUNITS - CORR - >=%u NETWORKS',CurrLimit))    
%     map=colormap;
%     map(1,:)=[1,1,1];
%     colormap(map)
%     set(gcf,'units','pixel')
%     Position=get(gcf,'position');
%     Position(3)=1100;
%     Position(4)=500;
%     set(gcf,'position',Position);
%     set(gca,'tickdir','out')
% end