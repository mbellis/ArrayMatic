%NET_PROPERTIES - Allows to find the component of covariation matrix.
% The whole structure is constructed by ascendant combination of cliques
% giving this schema : cliques => quasi-cliques (clusters of cliques) => regions (clusters of
% quasi-cliques)

%SUB FONCTIONS:

%   CONNECTIVITY
% load Connect_Climit_Alimit_mxxx_nxxxxx.mat with CLimit=ALimit=[0:10:60)
% contains five columns giving for each probe set
% 1 - the connectivity of CORR calculaed on CORR>CLimit
% 2 - the connectivity of ANTI calculaed on ANTI>ALimit
% 3 - the connectivity of CORR calculaed on CORR>CLimit and ANTI>ALimit
% 4 - the connectivity of CORR calculaed on CORR>CLimit and ANTI==0
% 4 - the connectivity of ANTI calculaed on CORR==0 and ANTI>ALimit
% display frequency of several class of connectivity ([0:100:2000]) for all CLimit values
%display ANTI connectivity vs CORR connectivity for the each CLimit
% posibility of removing probe set that are present in a second network
% (e.g. the 'intersect' network)

%SEARCH_CLIQUES
%

%INPUT PARAMETERS


%EXTERNAL SOURCES

%EXTERNAL APPLICATION
%need cliquer installed for finding cliques
%need to code of cliquer in order to mak it stop after 60s (very large
%network can be loaded and no algorithm can find the greatest clique)
%http://www.tkk.fi/~pat/cliquer.html


%OUTPUT PARAMETERS

% c) Michel Bellis
% CNRS
% michel.bellis@crbm.cnrs.fr
% arraymatic@gmail.com
% http://code.google.com/p/arraymatic
% distributed under the CeCILL license, which is compatible with the GNU
% General Public Licence and adapted to the European legislation.

function net_properties(Action,ModelRank,NetRank,CorrLimit,varargin)
global K
% %STUDY NETWORK STRUCTURE
BLOCSIZE=5000;
ChipPos=find(K.chip.rank==ModelRank);
PsNb=K.chip.probesetNb(ChipPos);
BlocNb=ceil(PsNb./BLOCSIZE);

switch Action
    case 'connectivity'
        NetRank1=NetRank;
        if nargin==4
            CONNECTIVITY(ModelRank,NetRank1)
        elseif nargin==5
            NetRank2=varargin{1};
            CONNECTIVITY(ModelRank,NetRank1,NetRank2)
        end
        
    case 'MCL properties'        
        CLUSTER_PROPERTIES(ModelRank,NetRank,BLOCSIZE,BlocNb,'Mcl',PsNb,CorrLimit)        
    case 'search cliques'
        DisplayFlag=0;
        MyAlgoFlag=0;
        if nargin>=5
            DisplayFlag=varargin{1};
        end
        if nargin>=6
            MyAlgoFlag=varargin{2};
        end
        SEARCH_CLIQUES(ModelRank,NetRank,PsNb,CorrLimit,DisplayFlag,MyAlgoFlag)
        
    case 'search adjacent cliques'
        InterNetRank=NetRank;
        UnionNetRank=varargin{1};               
        UnionCorrLimit=varargin{2};       
        SEARCH_ADJACENT_CLIQUES(ModelRank,InterNetRank,CorrLimit,UnionNetRank,PsNb,UnionCorrLimit)
        
    case 'cliques properties'
        CLUSTER_PROPERTIES(ModelRank,NetRank,BLOCSIZE,BlocNb,'cliques',PsNb,CorrLimit)
        
    case 'calculate cliques overlap'        
        CliNetRank=varargin{1};
        CALCULATE_OVERLAP(ModelRank,NetRank,CliNetRank,CorrLimit,BLOCSIZE,BlocNb,PsNb)
        
    case 'calculate MCL overlap'
        CALCULATE_OVERLAP(ModelRank,NetRank,varargin{1},varargin{2},'MCL',BLOCSIZE,BlocNb,PsNb)
        %net_properties('calculate overlap',8,1,6,8,54)
        
    case 'calculate MCL anti corr'
        if nargin==4
            CALCULATE_ANTI_CORR(ModelRank,NetRank,BLOCSIZE,BlocNb,'MCL',PsNb)
        else
            CALCULATE_ANTI_CORR(ModelRank,NetRank,BLOCSIZE,BlocNb,'MCL',PsNb,varargin{1},varargin{2})
        end
        
    case 'calculate cliques anti corr'
        if nargin==4
            CALCULATE_ANTI_CORR(ModelRank,NetRank,'cliques',PsNb,CorrLimit)
        else
            CALCULATE_ANTI_CORR(ModelRank,NetRank,'cliques',PsNb,CorrLimit,varargin{1},varargin{2})
        end
        
    case 'search quasi-cliques'
        %old version
        %FIND_CLIQUES_V1(DIRECTORY1,ModelRank,NetRank,PsNb)
        if nargin==4
            SEARCH_QUASI_CLIQUES(ModelRank,NetRank,PsNb,CorrLimit)
        elseif nargin==7
            RefNetRank=varargin{1};
            RefCliqueSize=varargin{2};
            RefCorrLimit=varargin{3};
            SEARCH_QUASI_CLIQUES(ModelRank,NetRank,PsNb,CorrLimit,RefNetRank,RefCliqueSize,RefCorrLimit);
        end
        
    case 'clusterize quasi-cliques'
        if nargin==4
        CLUSTERIZE_QUASI_CLIQUES(ModelRank,NetRank,PsNb,CorrLimit)
        elseif nargin==8
            RefNetRank=varargin{1};
            RefCliqueSize=varargin{2};
            RefCorrLimit=varargin{3};
            KeepRefFlag=varargin{4};
            CLUSTERIZE_QUASI_CLIQUES(ModelRank,NetRank,PsNb,CorrLimit,RefNetRank,RefCliqueSize,RefCorrLimit,KeepRefFlag)
        end
        
    case 'compare networks'        
        NetRank2=varargin{1};
        CliNetRank=varargin{2};
        COMPARE_NETWORKS(ModelRank,NetRank,NetRank2,CliNetRank)
end


%% CONNECTIVITY
function CONNECTIVITY(ModelRank,NetRank1,varargin)
global K
NetRank2=0;
Continue=1;
if nargin==3
    NetRank2=varargin{1};
    Directory2=fullfile(K.dir.net,sprintf('m%u',ModelRank),sprintf('n%u',NetRank2));
    cd(Directory2)
    if ~exist('stat','dir')
        Continue=0;
        h=warndlg(sprintf('RUN NETSTAT for m%un%u before',ModelRank,NetRank2));
        waitfor(h)
    else
        cd('stat')
        if ~exist(sprintf('Connect_C00_A00_m%u_n%u.mat',ModelRank,NetRank2),'file')
            Continue=0;
            h=warndlg(sprintf('RUN NETSTAT for m%un%u before',ModelRank,NetRank2));
            waitfor(h)
        end
    end
end
Directory1=fullfile(K.dir.net,sprintf('m%u',ModelRank),sprintf('n%u',NetRank1));
cd(Directory1)
if ~exist('stat','dir')
    Continue=0;
    h=warndlg(sprintf('RUN NETSTAT for m%un%u before',ModelRank,NetRank1));
    waitfor(h)
else
    cd('stat')
    if ~exist(sprintf('Connect_C00_A00_m%u_n%u.mat',ModelRank,NetRank1),'file')
        Continue=0;
        h=warndlg(sprintf('RUN NETSTAT for m%un%u before',ModelRank,NetRank1));
        waitfor(h)
    end
end
if Continue
    Colors=colors(colormap,7);
    Pos=[];
    for Limit=0:10:60
        cd(Directory1)
        cd('stat')
        if exist(sprintf('Connect_C%02u_A%02u_m%u_n%u.mat',Limit,Limit,ModelRank,NetRank1),'file')
            load(sprintf('Connect_C%02u_A%02u_m%u_n%u.mat',Limit,Limit,ModelRank,NetRank1))
            CConn1=Connect(:,1);
            AConn1=Connect(:,2);
            if NetRank2>0
                cd(Directory2)
                cd('stat')
                if exist(sprintf('Connect_C%02u_A%02u_m%u_n%u.mat',Limit,Limit,ModelRank,NetRank2),'file')
                    load(sprintf('Connect_C%02u_A%02u_m%u_n%u.mat',Limit,Limit,ModelRank,NetRank2))
                    CConn2=Connect(:,1);
                    Pos=find(CConn2);
                else
                    Pos=[];
                end
            end
            if Limit==0
                h1=figure;
                set(gcf,'color',[1,1,1])
                set(h1,'name','CONNECTIVITY PLOT')
                subplot(1,2,1)
                hold on
                subplot(1,2,2)
                hold on
                h2=figure;
                set(gcf,'color',[1,1,1])
                set(h2,'name','ANTI VS CORR CONNECTIVITY')
                hold on
            end
            figure(h1)
            N=histc(CConn1,0:100:2000);
            subplot(1,2,1)
            plot(0:100:2000,N,'O','color',Colors(Limit/10+1,:))
            plot(0:100:2000,N,'color',Colors(Limit/10+1,:))
            subplot(1,2,2)
            N=histc(AConn1,0:100:2000);
            plot(0:100:2000,N,'+','color',Colors(Limit/10+1,:))
            plot(0:100:2000,N,'color',Colors(Limit/10+1,:))
            if Limit==0
                subplot(1,2,1)
                xlabel('CORR(O) connectivity')
                ylabel('frequency')
                set(gca,'box','on')
                set(gca,'yscale','log')
                title(sprintf('CORR connectivity of m%un%u',ModelRank,NetRank1));
                subplot(1,2,2)
                xlabel('ANTI(+) connectivity')
                ylabel('frequency')
                set(gca,'box','on')
                set(gca,'yscale','log')
                title(sprintf('ANTI connectivity of m%un%u',ModelRank,NetRank1));
            end
            figure(h2)
            subplot(2,4,Limit/10+1)
            if ~isempty(Pos)
                %remove data corresponding to second network
                CConn1(Pos)=[];
                AConn1(Pos)=[];
            end
            plot(CConn1,AConn1,'.','color',Colors(Limit/10+1,:),'markersize',3)
            xlabel('CORR connectivity')
            ylabel('ANTI connectivity')
            title(sprintf('CORR>%u - ANTI>%u',Limit,Limit))
        end
    end
end
'stop'

Colors=colors(colormap,7);
h=figure;
set(h,'color',[1,1,1])
set(h,'name',sprintf('COMPARISON CONNECTIVITY m%un%u VS m%un%u',ModelRank,NetRank1,ModelRank,NetRank2))
PlotPos=0;
for Limit=0:10:60
    cd(Directory1)
    cd('stat')
    load(sprintf('Connect_C%02u_A%02u_m%u_n%u.mat',Limit,Limit,ModelRank,NetRank1))
    CConn1=Connect(:,1);
    AConn1=Connect(:,2);

    cd(Directory2)
    cd('stat')
    load(sprintf('Connect_C%02u_A%02u_m%u_n%u.mat',Limit,Limit,ModelRank,NetRank2))
    CConn2=Connect(:,1);
    AConn2=Connect(:,2);

    PlotPos=PlotPos+1;
    subplot(4,2,PlotPos)
    plot(CConn1,CConn2,'b.','markersize',3,'color',Colors(PlotPos,:))
    set(gca,'xlim',[0,1500])
    set(gca,'ylim',[0,1500])
    xlabel(sprintf('m%un%u connectivity',ModelRank,NetRank1))
    ylabel(sprintf('m%un%u connectivity',ModelRank,NetRank2))
    title(sprintf('CORR>%u',Limit))
end

cd(Directory1)
cd('stat')
load(sprintf('MeanCorr_m%u_n%u.mat',ModelRank,NetRank1))
MeanCorr1=MeanCorr;
cd(Directory2)
cd('stat')
load(sprintf('MeanCorr_m%u_n%u.mat',ModelRank,NetRank2))
MeanCorr2=MeanCorr;
h=figure;
set(h,'color',[1,1,1])
set(h,'name',sprintf('COMPARISON MEAN(CORR) m%un%u VS m%un%u',ModelRank,NetRank1,ModelRank,NetRank2))
subplot(1,3,1)
plot(MeanCorr1,MeanCorr2,'b.','markersize',3)
set(gca,'xlim',[0,100])
set(gca,'ylim',[0,100])
xlabel(sprintf('m%un%u mean(CORR)',ModelRank,NetRank1))
ylabel(sprintf('m%un%u mean(CORR)',ModelRank,NetRank2))
subplot(1,3,2)
hold on
Hist1=histc(MeanCorr1,[0:100]);
plot([0:100],Hist1/sum(Hist1),'r')
Hist2=histc(MeanCorr2,[0:100]);
plot([0:100],Hist2/sum(Hist2),'m')
xlabel('mean(CORR')
ylabel('dist')
subplot(1,3,3)
hold on
plot([0:100],cumsum(Hist1)/sum(Hist1),'r')
plot([0:100],cumsum(Hist2)/sum(Hist2),'m')
xlabel('mean(CORR')
ylabel('dist')





%% CLUSTER_PROPERTIES
function CLUSTER_PROPERTIES(ModelRank,NetRank,BlocSize,BlocNb,CluType,PsNb,CorrLimit)
global K
CLIQUE_SIZE=5;

NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',NetRank));
cd(NetDir)
if ~exist('stat','dir')
    Continue=0;
    h=warndlg(sprintf('RUN NETSTAT for m%un%u before',ModelRank,NetRank2));
    waitfor(h)
else
    cd('stat')
    if ~exist(sprintf('Connect_C%02u_A%02u_m%u_n%u.mat',CorrLimit,CorrLimit,ModelRank,NetRank),'file')
        Continue=0;
        h=warndlg(sprintf('RUN NETSTAT for m%un%u before',ModelRank,NetRank));
        waitfor(h)
    else
        eval(sprintf('load Connect_C%02u_A%02u_m%u_n%u.mat',CorrLimit,CorrLimit,ModelRank,NetRank))
        CConn=Connect(:,1);
    end
end


DoIt=1;
cd(NetDir)
switch CluType
    case 'MCL'
        DispType{1}='clusters';
        DispType{2}='cluster';

        if exist(sprintf('m%un%u_mcl_properties.mat',ModelRank,NetRank),'file')
            load(sprintf('m%un%u_mcl_properties.mat',ModelRank,NetRank))
            eval(sprintf('load m%un%u_%s_overlap.mat',ModelRank,NetRank,CluType));
            DoIt=0;
        else
            cd(Directory)
            if exist(sprintf('m%u_n%u_c0_i2_.mat',ModelRank,NetRank),'file')
                load(sprintf('m%u_n%u_c0_i2_.mat',ModelRank,NetRank));
            else
                h=errordlg('Do MCL FIRST');
                waitfor(h)
                error('process canceled')
            end
        end
    case 'cliques'
        DispType{1}='cliques';
        DispType{2}='clique';
%         if exist(sprintf('m%un%u_cliques_properties_%02u.mat',ModelRank,NetRank,CorrLimit),'file')
%             load(sprintf('m%un%u_cliques_properties_%02u.mat',ModelRank,NetRank,CorrLimit))
%             DoIt=0;
%         else
            if exist(sprintf('m%un%u_cliques_%02u.mat',ModelRank,NetRank,CorrLimit),'file')
                load(sprintf('m%un%u_cliques_%02u.mat',ModelRank,NetRank,CorrLimit));
            else
                h=errordlg('Do CLIQUES FIRST');
                waitfor(h)
                error('process canceled')
            end
%         end
end
CluNb=max(Clu);


if isequal(CluType,'cliques')

    h=figure;
    set(h,'name',sprintf('CONNECTIVITY vs  %s',CluType))
    plot(Clu,CConn,'k.','markersize',3)
    hold on
    plot(1:CluNb,CluSize,'r.')
    set(gca,'yscale','log')
    set(gcf,'color',[1,1,1])
    switch CluType
        case 'MCL'
            title(sprintf('connectivity of probe sets belonging to MCL clusters of m%un%u',ModelRank,NetRank));
            xlabel(' MCL clusters ordered by their rank')
            ylabel('connectivity of each probe set (black), size of MCL clusters (red)')
        case 'cliques'
            title(sprintf('connectivity of probe sets belonging to cliques of m%un%u',ModelRank,NetRank));
            xlabel(' cliques ordered by their rank')
            ylabel('connectivity of each probe set (black), size of cliques (red)')
    end
end
ConnTotal=zeros(CluNb,1);
for CluL=1:CluNb
    ConnTotal(CluL)=sum(CConn(Clu==CluL));
end


h=figure;
set(h,'name',sprintf('CONNECTIVITY vs  %s',CluType))
plot(Clu,CConn,'k.','markersize',3)
hold on
plot(1:CluNb,CluSize,'r.')
set(gca,'yscale','log')
set(gcf,'color',[1,1,1])
switch CluType
    case 'MCL'
        title(sprintf('connectivity of probe sets belonging to MCL clusters of m%un%u',ModelRank,NetRank));
        xlabel(' MCL clusters ordered by their rank')
        ylabel('connectivity of each probe set (black), size of MCL clusters (red)')
    case 'cliques'
        title(sprintf('connectivity of probe sets belonging to cliques of m%un%u',ModelRank,NetRank));
        xlabel(' cliques ordered by their rank')
        ylabel('connectivity of each probe set (black), size of cliques (red)')
end


h=figure;
set(h,'name',sprintf('EDGE NB vs  %s',CluType))
plot(Clu,CConn,'k.','markersize',3)
hold on
plot(1:CluNb,CluSize,'r')
plot(1:CluNb,ConnTotal,'g')
plot(1:CluNb,(ConnTotal/2)./((CluSize-1).*CluSize/2),'m')
set(gca,'yscale','log')
set(gcf,'color',[1,1,1])
switch CluType
    case 'MCL'
        title(sprintf('edge MCL clusters of m%un%u',ModelRank,NetRank));
        xlabel('MCL clusters ordered by size')
        ylabel(sprintf('cluster size in red, total nb of edges in green,/n ratio with edge nb in a cluster in magenta, ps connectivity in black)'))
    case 'cliques'
        title(sprintf('edge nb in cliques of m%un%u',ModelRank,NetRank));
        xlabel('cliques ordered by size')
        ylabel(sprintf('clique size in red, total nb of edges in green,\n ratio with edge nb in a clique in magenta, ps connectivity in black)'))
end

title(sprintf('edge nb vs  %s of m%un%u',CluType,ModelRank,NetRank));




%MORE INFORMATION ON  CLUSTERING (OUT and IN edges)
if DoIt

    %load(sprintf('m%un%u_connectivity.mat',ModelRank,NetRank))

    InClu=uint16(zeros(PsNb,CluNb));
    for BlocL=1:BlocNb
        if BlocL<BlocNb
            C=load_data(sprintf('c_m%u_n%u.4mat',ModelRank,NetRank),'./',PsNb,PsNb,'uint8','ieee-le',1:PsNb,((BlocL-1)*BlocSize)+1:BlocL*BlocSize);
        else
            C=load_data(sprintf('c_m%u_n%u.4mat',ModelRank,NetRank),'./',PsNb,PsNb,'uint8','ieee-le',1:PsNb,((BlocL-1)*BlocSize)+1:PsNb);
        end

        ColNb=size(C,2);
        C=C>0;
        %eliminate sel reference
        for ColL=1:ColNb
            C(((BlocL-1)*BlocSize)+ColL,ColL)=0;
        end
        for ColL=1:ColNb
            %the clusters to which belong the probe sets correlated to the current
            %probe set
            CurrClu=Clu(C(:,ColL)>0);
            CurrClu(CurrClu==0)=[];
            if ~isempty(CurrClu)
                %the number of probe set in each cluster
                CurrPsNb=histc(CurrClu,unique(CurrClu));
                %write for the current probe set the number of probe sets
                %linked to it in each cluster
                InClu(ColL+((BlocL-1)*BlocSize),unique(CurrClu))=CurrPsNb;
            end
        end
    end
    clear C

    %Number of external clusters
    OutCluNb=zeros(CluNb,1);
    %calculate In and Out edge nb
    OutEdgeNb=zeros(CluNb,1);
    InEdgeNb=zeros(CluNb,1);
    %The external cluster wich reveive the greatest number of link
    FirstOutClu=zeros(CluNb,1);
    %The percentage of out link that target the first external cluster
    FirstPrcOverlaping=zeros(CluNb,1);
    %The percentage of nodes that target the first external cluster
    FirstPrcOutputNodes=zeros(CluNb,1);
    %The percentage of nodes of the the first out external that are targeted
    FirstPrcInputNodes=zeros(CluNb,1);
    for CluL=1:CluNb
        PsPos=find(Clu==CluL);
        EdgeNb=sum(InClu(PsPos,:),1);        
        OutCluNb(CluL)=length(find(EdgeNb))-1;
        InEdgeNb(CluL)=max(1,EdgeNb(CluL)/2);
        EdgeNb(CluL)=0;
        OutEdgeNb(CluL)=sum(EdgeNb);
        if OutEdgeNb(CluL)>0
            [MaxEdgeNb MaxPos]=max(EdgeNb);
            OutPos=find(Clu==MaxPos);
            FirstOutClu(CluL)=MaxPos;
            FirstPrcOverlaping(CluL)=100*MaxEdgeNb/(length(PsPos)*length(OutPos));
            %percentage of nodes in the current cluster that target the first out cluster
            TargetingNodeNb=length(find(InClu(PsPos,MaxPos)));
            FirstPrcOutputNodes(CluL)=TargetingNodeNb*100/CluSize(CluL);
            %percentage of nodes in the first out cluster that are targeted by
            %the current cluster

            TargetedNodeNb=length(find(InClu(OutPos,CluL)));
        end
    end



    %INTER CLUSTER OVERLAPING
    %VERY LONG !!!
    CluOverlap=single(zeros(max(Clu)));
    for CluL=1:max(Clu)
        CluOverlap(CluL,CluL)=100;
    end
    %     for CluL1=1:max(Clu)-1
    %         CluL1
    %         for CluL2=CluL1+1:max(Clu)
    %             PsPos1=find(Clu==CluL1);
    %             EdgeNb=sum(InClu(PsPos1,CluL2));
    %             if EdgeNb>0
    %                 CluOverlap(CluL1,CluL2)=100*EdgeNb/(length(PsPos1)*length(find(Clu==CluL2)));
    %                 CluOverlap(CluL2,CluL1)=CluOverlap(CluL1,CluL2);
    %             end
    %         end
    %     end


    for CluL1=1:max(Clu)-1
        PsPos1=find(Clu==CluL1);
        EdgeNbs=sum(InClu(PsPos1,CluL1+1:max(Clu)));
        CluOverlap(CluL1,CluL1+1:max(Clu))=single(100*EdgeNbs./(length(PsPos1)*CluSize(CluL1+1:max(Clu))'));
        CluOverlap(CluL1+1:max(Clu),CluL1)=CluOverlap(CluL1,CluL1+1:max(Clu))';
    end


    cd(NetDir)
    eval(sprintf('save m%un%u_%s_overlap_%02u.mat CluOverlap',ModelRank,NetRank,CluType,CorrLimit));
    clear CluOverlap


    cd(NetDir)
    eval(sprintf('save m%un%u_%s_properties_%02u.mat Clu* In* Out* First*',ModelRank,NetRank,CluType,CorrLimit));

    eval(sprintf('load m%un%u_%s_overlap_%02u.mat',ModelRank,NetRank,CluType,CorrLimit));

    %DISTIBUTION OF CONNECTIVITY
    ConnVal=unique(CConn);
    BinNb=histc(CConn,ConnVal);
    h=figure;
    set(h,'name',sprintf('Connectivity distribution (m%un%u)',ModelRank,NetRank))
    plot(ConnVal(2:end),BinNb(2:end),'r.')
    hold on
    plot(ConnVal(2:end),BinNb(2:end),'r.')
    set(gca,'yscale','log')
    xlabel('connectivity')
    ylabel('frequency')
    title('dist connectivity (>0)')
    set(gcf,'color',[1,1,1])


    %DISTRIBUTION OF MAX(OVERLAP) AND <MAX(OVERLAP)
    for CluL=1:CluNb
        CluOverlap(CluL,CluL)=0;
    end
    CluSizeByPs=zeros(PsNb,1);
    for PsL=1:PsNb
        if Clu(PsL)>0
            CluSizeByPs(PsL)=CluSize(Clu(PsL));
        end
    end

    %selection on probe set connectivity
    CONN1=51;
    CONN2=max(CConn);
    PsPos=find(CConn>=CONN1&CConn<=CONN2);
    MaxOverlap=[];
    OverlapVal=zeros(1,101);
    for PsL=1:length(PsPos)
        CurrClu=Clu(PsPos(PsL));
        if CurrClu>0
            Val=CluOverlap(CurrClu,:);
            [MaxVal,MaxPos]=max(Val);
            MaxOverlap=[MaxOverlap;MaxVal];
            Val(MaxPos)=0;
            OverlapVal=OverlapVal+histc(Val,[0:100]);
        end
    end

    h=figure;
    set(h,'name',sprintf('overlap for cluster containing ps with Conn1=%u Conn2=%u (m%un%u)',CONN1,CONN2,ModelRank,NetRank))
    subplot(1,2,1)
    BinNb=histc(MaxOverlap,[0:100]);
    bar([1:100],BinNb(2:end)*100/sum(BinNb(2:end)),'histc')
    set(gca,'box','on')
    set(gca,'xlim',[0,100])
    xlabel('overlap')
    ylabel('frequency')
    title('dist max overlap (>0)')
    subplot(1,2,2)
    bar([1:100],OverlapVal(2:end)*100/sum(OverlapVal(2:end)),'histc')
    set(gca,'box','on')
    set(gca,'xlim',[0,100])
    xlabel('overlap')
    ylabel('frequency')
    title('dist other overlaps (>0)')
    set(gcf,'color',[1,1,1])

    %selection on Clustr size
    CLUSIZE1=51;
    CLUSIZE2=max(CluSize);
   

    CluPos=find(CluSize>=CLUSIZE1&CluSize<=CLUSIZE2);
    MaxOverlap=[];
    OverlapVal=zeros(1,101);
    PsNbClu=0;
    for CluL=1:length(CluPos)
        CurrClu=CluPos(CluL);
        PsNbClu=PsNbClu+length(find(Clu==CurrClu));
        if CurrClu>0
            Val=CluOverlap(CurrClu,:);
            Val(CurrClu)=0;
            [MaxVal,MaxPos]=max(Val);
            MaxOverlap=[MaxOverlap;MaxVal];
            Val(MaxPos)=0;
            OverlapVal=OverlapVal+histc(Val,[0:100]);
        end
    end
    MaxPos=find(MaxOverlap>0);
    length(find(MaxOverlap(MaxPos)>=30))*100/length(MaxPos)
    length(find(MaxOverlap(MaxPos)>=40))*100/length(MaxPos)
    length(find(MaxOverlap(MaxPos)>=50))*100/length(MaxPos)

    h=figure;
    set(gcf,'color',[1,1,1])
    set(h,'name',sprintf('overlap for cluster of sizes >==%u and <=%u (m%un%u)',CLUSIZE1,CLUSIZE2,ModelRank,NetRank))
    hold on
    BinNb=histc(MaxOverlap,[0:100]);
    plot([0:100],cumsum(BinNb)/sum(BinNb),'g.')
    plot([0:100],cumsum(BinNb)/sum(BinNb),'g')
    plot([0:100],OverlapVal/sum(OverlapVal),'r.')
    plot([0:100],OverlapVal/sum(OverlapVal),'r')
    set(gca,'box','on')
    xlabel('overlap value')
    ylabel('cumulative frequency (max(green) and others(red))')
    set(gca,'ylim',[0,1])
end

h=figure;
Pos2=find(CluSize==2);
Pos2=Pos2(1)-1;
set(h,'name',sprintf('IN AND OUT CONNECTIVITY vs  %s',CluType))
hold on
plot(1:CluNb,OutEdgeNb,'g');
plot(1:CluNb,InEdgeNb,'c');
plot(1:CluNb,CluSize,'r');
% plot(1:CluNb,OutCluNb,'k+');
% plot(1:CluNb,FirstOutClu,'ko');
plot(1:CluNb,(InEdgeNb)./((CluSize-1).*CluSize/2),'m');
line([0,Pos2],[1,1],'color','b','linestyle',':')
set(gca,'yscale','log')
set(gca,'box','on')

switch CluType
    case 'MCL'
        title(sprintf('In and Out connectivity in MCL clusters of m%un%u',ModelRank,NetRank));
        xlabel(' clusters ordered by size')
        ylabel(sprintf('cliques size in red, nb of in and out edges in cyan and green,\n ratio of nb of in edges with edge nb in a cluster in magenta'))
    case 'cliques'
        title(sprintf('In and Out connectivity in cliques of m%un%u',ModelRank,NetRank));
        xlabel(' cliques ordered by size')
        ylabel(sprintf('cliques size in red, nb of in and out edges in cyan and green,\n ratio of nb of in edges with edge nb in a cluster in magenta'))
end

set(gcf,'color',[1,1,1])

h=figure;
set(gcf,'color',[1,1,1])
Pos2=find(CluSize==2);
Pos2=Pos2(1)-1;
set(h,'name','RELATION WITH OTHER CLUSTERS')
hold on
plot(1:CluNb,CluSize,'r');
plot(1:CluNb,(InEdgeNb)./((CluSize-1).*CluSize/2),'m');
plot(1:CluNb,OutCluNb,'b+','markersize',3);
plot(1:CluNb,FirstOutClu,'ro','markersize',3);
line([0,Pos2],[1,1],'color','b','linestyle',':')
set(gca,'yscale','log')
set(gca,'box','on')
switch CluType
    case 'MCL'
        title(sprintf('In and Out connectivity in MCL clusters of m%un%u',ModelRank,NetRank));
        xlabel(' clusters ordered by size')
        ylabel('cluster size (red line),nb of out cluster (b+), first out cluster (ro)')
    case 'cliques'
        title(sprintf('In and Out connectivity in cliques of m%un%u',ModelRank,NetRank));
        xlabel(' cliques ordered by size')
        ylabel('clique size (red line),nb of out cliques (b+), first out clique (ro)')
end


for TypeL=1:3
    h=figure;
    Pos2=find(CluSize==2);
    Pos2=Pos2(1)-1;
    switch TypeL
        case 1
            set(h,'name',sprintf('PERCENTAGE OF EDGES POINTING IN THE FIRST OUT %s',upper(DispType{1})))
        case 2
            set(h,'name',sprintf('PERCENTAGE OF NODES TARGETING THE FIRST OUT %s',upper(DispType{1})))
        case 3
            set(h,'name',sprintf('PERCENTAGE OF NODES OF THE FIRST OUT %s',upper(DispType{1})))
    end
    hold on
    plot(1:CluNb,CluSize,'r');
    plot(1:CluNb,OutCluNb,'b+','markersize',3);
    plot(1:CluNb,FirstOutClu,'ro','markersize',3);
    if TypeL==1
        plot(1:CluNb,FirstPrcOverlaping/100,'g');
        title(sprintf('In and Out connectivity in %s of m%un%u',DispType{1},ModelRank,NetRank));
        xlabel(sprintf('%s ordered by size',DispType{1}))
        ylabel(sprintf('fraction of overlaping with the first out %s (green),\n size in red,nb of out %s (b+), first out %s (ro)',DispType{2},DispType{2},DispType{2}))
    elseif TypeL==2
        plot(1:CluNb,FirstPrcOutputNodes/100,'g');
        title(sprintf('Output nodes in %s of m%un%u',DispType{1},ModelRank,NetRank));
        xlabel(sprintf('%s ordered by size',DispType{1}))
        ylabel(sprintf('fraction of overlaping with the first out %s (green),\n size in red,nb of out %s (b+), first out %s (ro)',DispType{2},DispType{2},DispType{2}))
    else
        plot(1:CluNb,FirstPrcInputNodes/100,'g');
        title(sprintf('Input nodes in %s of m%un%u',DispType{1},ModelRank,NetRank));
        xlabel(sprintf('%s ordered by size',DispType{1}))
        ylabel(sprintf('fraction of nodes of the first out %s that are targeted (green),\n size in red,nb of out %s (b+), first out %s (ro)',DispType{2},DispType{2},DispType{2}))
    end
    line([0,Pos2],[1,1],'color','b','linestyle',':')
    set(gca,'yscale','log')
    set(gca,'box','on')

    set(gcf,'color',[1,1,1])
end

%IMAGES OF CLUSTERS
Pos=find(CluSize>=20);
PsRanks=[];
for CluL=1:length(Pos)
    PsRanks=[PsRanks;find(Clu==Pos(CluL))];
end
NetRank=41;
NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',NetRank));
cd(NetDir)
C=load_data(sprintf('c_m%u_n%u.4mat',ModelRank,NetRank),'./',PsNb,PsNb,'uint8','ieee-le',PsRanks,PsRanks);
C=triu(C);
A=load_data(sprintf('a_m%u_n%u.4mat',ModelRank,NetRank),'./',PsNb,PsNb,'uint8','ieee-le',PsRanks,PsRanks);
C=C+tril(A,-1);
clear A
figure
image(C);
title(sprintf('anti.corr cliques m%un%u',ModelRank,NetRank))
set(gcf,'color',[1,1,1])
%% CALCULATE ANTI CORR
function CALCULATE_ANTI_CORR(ModelRank,NetRank,CluType,PsNb,CorrLimit,varargin)
global K
if nargin==5
    NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',NetRank));
    CluType=lower(CluType);
    cd(NetDir)
    DoIt=1;
%     if exist(sprintf('m%un%u_%s_anticorr_%02u.mat',ModelRank,NetRank,CluType,CorrLimit),'file')
%         load(sprintf('m%un%u_%s_anticorr_%02u.mat',ModelRank,NetRank,CluType,CorrLimit))
%         DoIt=0;
%     else
        if exist(sprintf('m%un%u_%s_properties_%02u.mat',ModelRank,NetRank,CluType,CorrLimit),'file')
            load(sprintf('m%un%u_%s_properties_%02u.mat',ModelRank,NetRank,CluType,CorrLimit))
            clear CluSize CluSizeByPs First* In* Out*
        else
            h=errordlg(sprintf('Do %s PROPERTIES',CluType));
            waitfor(h)
            error('process canceled')
        end
        if exist(sprintf('m%un%u_%s_overlap_%02u.mat',ModelRank,NetRank,CluType,CorrLimit),'file')
            load(sprintf('m%un%u_%s_overlap_%02u.mat',ModelRank,NetRank,CluType,CorrLimit))
        else
            h=errordlg(sprintf('Do %s OVERLAP',CluType));
            waitfor(h)
            error('process canceled')
        end

%     end
elseif nargin==7
    CliModelRank=varargin{1};
    CliNetRank=varargin{2};

    NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',NetRank));
    CliNetDir=fullfile(K.dir.net,sprintf('m%03u',CliModelRank),sprintf('n%05u',CliNetRank));
    CluType=lower(CluType);
    DoIt=1;
    cd(NetDir)
    if exist(sprintf('m%un%u_%s_anticorr_cli_m%un%u.mat',ModelRank,NetRank,CluType,CliModelRank,CliNetRank),'file')
        load(sprintf('m%un%u_%s_anticorr_m%un%u.mat',ModelRank,NetRank,CluType,CliModelRank,CliNetRank))
        DoIt=0;
    else
        cd(CliNetDir)
        if exist(sprintf('m%un%u_%s_properties.mat',CliModelRank,CliNetRank,CluType),'file')
            load(sprintf('m%un%u_%s_properties.mat',CliModelRank,CliNetRank,CluType))
            clear CluSize CluSizeByPs First* In* Out*
        else
            h=errordlg(sprintf('Do %s PROPERTIES',CluType));
            waitfor(h)
            error('process canceled')
        end
    end
end
if DoIt

    %INTER CLUSTER MEAN ANTI/CORR
    %VERY LONG !!!
    CluAntiCorr=single(zeros(max(Clu)));
    cd(NetDir)
    for CluL1=1:max(Clu)-1
        CluL1
        PsPos=find(Clu==CluL1);
        C=load_data(sprintf('c_m%u_n%u.4mat',ModelRank,NetRank),'./',PsNb,PsNb,'uint8','ieee-le',1:PsNb,PsPos);
        A=load_data(sprintf('a_m%u_n%u.4mat',ModelRank,NetRank),'./',PsNb,PsNb,'uint8','ieee-le',1:PsNb,PsPos);
        for CluL2=CluL1+1:max(Clu)
            if CluOverlap(CluL1,CluL2)>0
                CluAntiCorr(CluL1,CluL2)=mean(mean(C(find(Clu==CluL2),:)));
                CluAntiCorr(CluL2,CluL1)=mean(mean(A(find(Clu==CluL2),:)));
            end
        end
    end
    cd(NetDir)
    if nargin==5
        eval(sprintf('save m%un%u_%s_anticorr_%02u CluAntiCorr',ModelRank,NetRank,CluType,CorrLimit));
    else
        eval(sprintf('save m%un%u_%s_anticorr_%02u_cli_m%un%u CluAntiCorr',ModelRank,NetRank,CluType,CorrLimit,CliModelRank,CliNetRank))
    end
end
h=figure;
set(h,'name',sprintf('M%uN%u -  ANTI CORR BETWEEN CLIQUES',ModelRank,NetRank))
image(CluAntiCorr)
title(sprintf('M%uN%u -  CORR BETWEEN CLIQUES',ModelRank,NetRank))
xlabel(sprintf('M%uN%u -  ANTI BETWEEN CLIQUES',ModelRank,NetRank))
set(gcf,'color',[1,1,1])
set(gca,'xlim',[0,500])
set(gca,'ylim',[0,500])

AntiCorr=CluAntiCorr;
RAntiCorr=AntiCorr';
M=zeros(size(RAntiCorr));
M(find(AntiCorr>RAntiCorr))=AntiCorr(find(AntiCorr>RAntiCorr))-RAntiCorr(find(AntiCorr>RAntiCorr));
%M(find(AntiCorr>RAntiCorr))=C(find(AntiCorr>RAntiCorr))-AntiCorr(find(AntiCorr>RAntiCorr));
h=figure;
set(h,'name',sprintf('M%uN%u -  MAX DIFF ANTI CORR BETWEEN CLIQUES',ModelRank,NetRank))
image(M)
title(sprintf('M%uN%u -  CORR-ANTI BETWEEN CLIQUES',ModelRank,NetRank))
xlabel(sprintf('M%uN%u -  ANTI-CORR BETWEEN CLIQUES',ModelRank,NetRank))
set(gcf,'color',[1,1,1])

U=triu(AntiCorr);
L=tril(AntiCorr);
L=L';
h=figure;
set(h,'name',sprintf('M%uN%u - BETWEEN CLIQUES ANTI vs CORR ',ModelRank,NetRank))
plot(U(U>L),L(U>L),'r.')
hold on
plot(U(U<=L),L(U<=L),'b.')
set(gcf,'color',[1,1,1])
xlabel('CORR BETWEEN CLIQUES')
ylabel('ANTI BETWEEN CLIQUES')

DispTYpe=questdlg('do you want to display anti-corr for clusters of cliques','','yes','no','yes');
if isequal(DispTYpe,'yes')
    [FileName,FileDir]=uigetfile(sprintf('m%un%u*.mat',ModelRank,NetRank),'select a file');
    cd(FileDir)
    load(FileName)
    LastRound=length(Cliques);
    CurrCliques=[];
    for CliL=1:length(Cliques{LastRound});
        CurrCliques=[CurrCliques;Cliques{LastRound}{CliL}];
    end
    CliNb=length(CurrCliques);
    if CliNb>length(CluAntiCorr)
        h=figure;
        set(h,'name',sprintf('M%uN%u - BETWEEN CLIQUES ANTI vs CORR FOR CLUSTERS',ModelRank,NetRank))
        image(CluAntiCorr(CurrCliques,CurrCliques))
        title(sprintf('M%uN%u - BETWEEN CLIQUES CORR FOR CLUSTERS',ModelRank,NetRank))
        xlabel(sprintf('M%uN%u - BETWEEN CLIQUES ANTI FOR CLUSTERS',ModelRank,NetRank))
        set(gcf,'color',[1,1,1])
    else
        h=warndlg('%s can be processed with %u (different nb of cliques)',FileName,sprintf('save m%un%u_%s_anticorr.mat CluAntiCorr',ModelRank,NetRank,CluType));
        waitfor(h)
    end

end



%% CALCULATE OVERLAP
function CALCULATE_OVERLAP(ModelRank,NetRank,CliNetRank,CorrLimit,BlocSize,BlocNb,PsNb)
global K
tic
%load definition of cliques constructed on model network
NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',CliNetRank));
cd(NetDir)
%load definition of cliques
cd(NetDir)
if exist(sprintf('m%un%u_cliques_%02u.mat',ModelRank,CliNetRank,CorrLimit),'file')
    load(sprintf('m%un%u_cliques_%02u.mat',ModelRank,CliNetRank,CorrLimit))
else
    h=errordlg('Do search cliques');
    waitfor(h)
    error('process canceled')
end



if exist(sprintf('m%un%u_cliques_properties_%02u.mat',ModelRank,CliNetRank,CorrLimit),'file')
    load(sprintf('m%un%u_cliques_properties_%02u.mat',ModelRank,CliNetRank,CorrLimit))
    clear CluSize CluSizeByPs First* In* Out*
else
    h=errordlg(sprintf('Do %s PROPERTIES',CluType));
    waitfor(h)
    error('process canceled')
end
% if exist(sprintf('m%un%u_cliques_overlap_%02u.mat',ModelRank,CliNetRank,CorrLimit),'file')
%     load(sprintf('m%un%u_cliques_overlap_%02u.mat',ModelRank,CliNetRank,CorrLimit))
%     clear CluSize CluSizeByPs First* In* Out*
% else
%     h=errordlg(sprintf('Do %s PROPERTIES',CluType));
%     waitfor(h)
%     error('process canceled')
% end
% ModelCluOverlap=CluOverlap;
% clear CluOverlap

CluNb=max(Clu);

%load CORR values for the processed network in order to calculate overlap between cliques defined on the model network
NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',NetRank));
cd(NetDir)
%calulate for each probe set the number of links that comes from each clique
InClu=uint16(zeros(PsNb,CluNb));
for BlocL=1:BlocNb
    if BlocL<BlocNb
        C=load_data(sprintf('c_m%u_n%u.4mat',ModelRank,NetRank),'./',PsNb,PsNb,'uint8','ieee-le',1:PsNb,((BlocL-1)*BlocSize)+1:BlocL*BlocSize);
    else
        C=load_data(sprintf('c_m%u_n%u.4mat',ModelRank,NetRank),'./',PsNb,PsNb,'uint8','ieee-le',1:PsNb,((BlocL-1)*BlocSize)+1:PsNb);
    end
    ColNb=size(C,2);
    C=C>0;
    %eliminate self reference
    for ColL=1:ColNb
        C(((BlocL-1)*BlocSize)+ColL,ColL)=0;
    end
    for ColL=1:ColNb
        %the clusters to which belong the probe sets correlated to the current
        %probe set
        CurrClu=Clu(C(:,ColL)>0);
        CurrClu(CurrClu==0)=[];
        if ~isempty(CurrClu)
            %the number of probe set in each cluster
            CurrPsNb=histc(CurrClu,unique(CurrClu));
            %write for the current probe set the number of probe sets
            %linked to it in each cluster
            InClu(ColL+((BlocL-1)*BlocSize),unique(CurrClu))=CurrPsNb;
        end
    end
end
clear C


CluSize=zeros(max(Clu),2);
for CluL=1:max(Clu)
    CluSize(CluL,:)=[CluL,length(find(Clu==CluL))];
end

%CluOVerlap
CluOverlap=single(zeros(max(Clu)));
for CluL=1:max(Clu)
    CluOverlap(CluL,CluL)=100;
end
for CluL1=1:max(Clu)-1
    PsPos1=find(Clu==CluL1);
    EdgeNbs=sum(InClu(PsPos1,CluL1+1:max(Clu)));
    CluOverlap(CluL1,CluL1+1:max(Clu))=single(100*EdgeNbs./(length(PsPos1)*CluSize(CluL1+1:max(Clu),2)'));
    CluOverlap(CluL1+1:max(Clu),CluL1)=CluOverlap(CluL1,CluL1+1:max(Clu))';
end

%     CluOverlap=single(zeros(max(Clu)));
%     for CluL=1:max(Clu)
%         CluOverlap(CluL,CluL)=100;
%     end
%     for CluL1=1:max(Clu)-1
%         CluL1
%         for CluL2=CluL1+1:max(Clu)
%             PsPos1=find(Clu==CluL1);
%             EdgeNb=sum(InClu(PsPos1,CluL2));
%             if EdgeNb>0
%                 CluOverlap(CluL1,CluL2)=100*EdgeNb/(length(PsPos1)*length(find(Clu==CluL2)));
%                 CluOverlap(CluL2,CluL1)=CluOverlap(CluL1,CluL2);
%             end
%         end
%     end
toc
cd(NetDir)
eval(sprintf('save m%un%u_cliques_overlap_with_m%un%u_%02u CluOverlap InClu;',ModelRank,NetRank,ModelRank,CliNetRank,CorrLimit))

%% COMPARE NETWORKS
function COMPARE_NETWORKS(ModelRank,NetR1,NetR2,CliNetRank)
%net_properties('compare networks',8,1,6,54,54)
global K


DISP_SIZE=1000;
CLIQUE_SIZE=5;

NetRank(1)=NetR1;
NetRank(2)=NetR2;
Overlap1=0;
Overlap2=0;

NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',CliNetRank));
cd(NetDir)
eval(sprintf('load m%un%u_quasi-cliques_clusters_size%u',ModelRank,CliNetRank,CLIQUE_SIZE))
clear CluOverlap Clusters
CliqueIndex=[];
for CliL=1:length(Cliques{end})
    CliqueIndex=[CliqueIndex;Cliques{end}{CliL}];
end


for NetL=1:2
    NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',NetRank(NetL)));
    cd(NetDir)
    if NetRank(NetL)==CliNetRank
        FileName=sprintf('m%un%u_cliques_overlap.mat',ModelRank,NetRank(NetL));
    else
        FileName=sprintf('m%un%u_cliques_overlap_cli_m%un%u.mat',ModelRank,NetRank(NetL),ModelRank,CliNetRank);
    end
    if exist(FileName,'file')
        if NetL==1
            load(FileName)
            clear InClu
            Overlap1=1;
            Overlap{1}=CluOverlap(CliqueIndex,CliqueIndex);
            clear CluOverlap
        else
            if Overlap1==1
                load(FileName)
                clear InClu
                Overlap2=1;
                Overlap{2}=CluOverlap(CliqueIndex,CliqueIndex);
                clear CluOverlap
            end
        end
    end
    if NetRank(NetL)==CliNetRank
        FileName=sprintf('m%un%u_cliques_anticorr.mat',ModelRank,NetRank(NetL));
    else
        FileName=sprintf('m%un%u_cliques_anticorr_cli_m%un%u.mat',ModelRank,NetRank(NetL),ModelRank,CliNetRank);
    end

    if exist(FileName,'file')
        load(FileName)
        AntiCorr{NetL}=CluAntiCorr(CliqueIndex,CliqueIndex);
        clear CluAntiCorr
    else
        h=errordlg(sprintf('%s does not exist',FileName));
        waitfor(h)
        error('process canceled')
    end
end
CliNb=length(AntiCorr{1});
if Overlap2
    h=figure;
    subplot(2,2,1)
    image(Overlap{1})
    title(sprintf('OVERLAP M%uN%u',ModelRank,NetRank(1)))
    subplot(2,2,2)
    h=pcolor([double(flipud(Overlap{1})-flipud(Overlap{2})),zeros(CliNb,1);zeros(1,CliNb+1)]);
    title(sprintf('OVERLAP M%uN%u - OVERLAP M%uN%u',ModelRank,NetRank(1),ModelRank,NetRank(2)))
    set(h,'linestyle','none')
end
if Overlap2==0
    h=figure;
end
if Overlap2
    subplot(2,2,3)
else
    subplot(1,2,1)
end
image(AntiCorr{1})
title(sprintf('ANTI/CORR M%uN%u',ModelRank,NetRank(1)))
if Overlap2
    subplot(2,2,4)
else
    subplot(1,2,2)
end
h=pcolor([double(flipud(AntiCorr{1})-flipud(AntiCorr{2})),zeros(CliNb,1);zeros(1,CliNb+1)]);
title(sprintf('ANTI/CORR M%uN%u - ANTI/CORR M%uN%u',ModelRank,NetRank(1),ModelRank,NetRank(2)))
set(h,'linestyle','none')
set(gcf,'color',[1,1,1])





ACDiff=AntiCorr{1}-AntiCorr{2};
DiffAnti=triu(ACDiff',1);
DiffCorr=triu(ACDiff,1);
APos=find(abs(DiffCorr)<abs(DiffAnti));
CPos=find(abs(DiffCorr)>abs(DiffAnti));

ODiff=Overlap{1}-Overlap{2};
DiffOverlap=triu(ODiff,1);
DiffOverlap=DiffOverlap(:);
DiffOverlap(DiffOverlap==0)=[];

figure
set(gcf,'color',[1,1,1])

CurrDiffCorr=DiffCorr;
CurrDiffCorr(APos)=0;
CurrDiffCorr=CurrDiffCorr(:);
CurrDiffCorr(CurrDiffCorr==0)=[];
subplot(2,3,1)
hist(CurrDiffCorr,200);
title('CORR>ANTI DIFFERENCE')
xlabel('corr difference')
ylabel('frequency of non null difference')

CurrDiffAnti=DiffAnti;
CurrDiffAnti(CPos)=0;
CurrDiffAnti=CurrDiffAnti(:);
CurrDiffAnti(CurrDiffAnti==0)=[];
subplot(2,3,2)
hist(CurrDiffAnti,200);
title('ANTI>CORR DIFFERENCE')
xlabel('anti difference')
ylabel('frequency of non null difference')

subplot(2,3,3)
hist(DiffOverlap,200);
title('OVERLAP DIFFERENCE')
xlabel('overlap difference')
ylabel('frequency of non null difference')



subplot(2,3,4)
hold on
plot(DiffCorr(CPos),DiffAnti(CPos),'r.','markersize',3)
plot(DiffCorr(APos),DiffAnti(APos),'b.','markersize',3)
title('ANTI vs CORR DIFFERENCES')
xlabel('corr difference')
ylabel('anti difference')
set(gca,'box','on')
set(gca,'xlim',[-100,100])
set(gca,'ylim',[-100,100])

subplot(2,3,5)
hold on
plot(ODiff(CPos),DiffCorr(CPos),'r.','markersize',3)
plot(ODiff(APos),DiffCorr(APos),'b.','markersize',3)
title('CORR vs OVERLAP DIFFERENCES')
xlabel('overlap difference')
ylabel('corr difference')
set(gca,'box','on')
set(gca,'xlim',[-100,100])
set(gca,'ylim',[-100,100])


subplot(2,3,6)
hold on
plot(ODiff(CPos),DiffAnti(CPos),'r.','markersize',3)
plot(ODiff(APos),DiffAnti(APos),'b.','markersize',3)
title('ANTI vs OVERLAP DIFFERENCES')
xlabel('overlap difference')
ylabel('anti difference')
set(gca,'box','on')
set(gca,'xlim',[-100,100])
set(gca,'ylim',[-100,100])

figure
set(gcf,'color',[1,1,1])
h=pcolor([double(flipud(Overlap{1})-flipud(Overlap{2})),zeros(CliNb,1);zeros(1,CliNb+1)]);
title(sprintf('OVERLAP M%uN%u - OVERLAP M%uN%u',ModelRank,NetRank(1),ModelRank,NetRank(2)))
set(h,'linestyle','none')

figure
set(gcf,'color',[1,1,1])
h=pcolor([double(flipud(AntiCorr{1})-flipud(AntiCorr{2})),zeros(CliNb,1);zeros(1,CliNb+1)]);
title(sprintf('ANTI/CORR M%uN%u - ANTI/CORR M%uN%u',ModelRank,NetRank(1),ModelRank,NetRank(2)))
set(h,'linestyle','none')
set(gcf,'color',[1,1,1])










%% SEARCH CLIQUES
function SEARCH_CLIQUES(ModelRank,NetRank,PsNb,CorrLimit,DisplayFlag,MyAlgoFlag)
global K
%SEARCH CLIQUES
NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',NetRank));
cd(NetDir)
Continue=1;
if ~exist('stat','dir')
    Continue=0;
    h=warndlg(sprintf('RUN NETSTAT for m%un%u before',ModelRank,NetRank));
    waitfor(h)
else
    cd('stat')
    if ~exist(sprintf('Connect_C00_A00_m%u_n%u.mat',ModelRank,NetRank),'file')
        Continue=0;
        h=warndlg(sprintf('RUN NETSTAT for m%un%u before',ModelRank,NetRank));
        waitfor(h)
    else
        load(sprintf('Connect_C%02u_A%02u_m%u_n%u.mat',CorrLimit,CorrLimit,ModelRank,NetRank))
        LoadedPsRanks=[1:PsNb]';
        CurrConn=Connect(:,1);
    end
end
if Continue==1
    CliqueRank=0;
    CluSize=[];
    %Clique membership. Generic name Clu is used to allow using other
    %functions
    Clu=zeros(PsNb,1);
    %indicates if the current probe set when used as hub for agregating a
    %clique belongs or not to that clique
    % 0 = the probebe set does not belong to a clique or belongs to a clique agregated around another probe
    % set
    % 1 = the probe set belongs to  a clique that has bee agregated around
    % itself
    IsSeed=zeros(PsNb,1);
    CliqueSize=zeros(PsNb,1);
    [Clu,CluSize,CliqueSize,IsSeed]=CLIQUES(ModelRank,NetRank,LoadedPsRanks,CurrConn,PsNb,CorrLimit,Clu,CluSize,CliqueSize,IsSeed,CliqueRank,MyAlgoFlag,DisplayFlag);
    CliqueSize=zeros(PsNb,1);
    for i=1:PsNb
        if Clu(i)>0
            CliqueSize(i)=CluSize(Clu(i));
        end
    end
    %reorder on clique size
    [Temp,SortOrder]=sort(CluSize,'descend');
    [Temp NewRank]=sort(SortOrder);
    for PsL=1:length(Clu)
        if Clu(PsL)>0
            Clu(PsL)=NewRank(Clu(PsL));
        end
    end
    CluSize=zeros(max(Clu),1);
    for CluL=1:max(Clu)
        CluSize(CluL)=length(find(Clu==CluL));
    end
    cd(NetDir)
    eval(sprintf('save m%un%u_cliques_%02u Clu CluSize CliqueSize IsSeed',ModelRank,NetRank,CorrLimit))
end
%% COMPARE CLIQUES
function COMPARE_CLIQUES()
cd(K.dir.net)
[RefFileName,FileDir]=uigetfile('*.mat','Load reference cliques');
cd(FileDir)
load(RefFileName)
RefClu=Clu;
RefCluSize=CluSize;
RefCliqueSize=CliqueSize;
[FileName,FileDir]=uigetfile('*.mat','Load cliques to be compared');
cd(FileDir)
load(FileName)
[RefClu,SortIndex]=sort(RefClu);
RefCliqueSize=RefCliqueSize(SortIndex);
Clu=Clu(SortIndex);
CliqueSize=CliqueSize(SortIndex);
PsNb=length(Clu);

RefFileName=strrep(RefFileName(1:end-4),'_',' ');
FileName=strrep(FileName(1:end-4),'_',' ');
h=figure;
set(gcf,'color',[1,1,1])
set(h,'name','CLIQUE COMPARISON')
subplot(2,2,1)
hold on
plot(Clu,'g.','markersize',3)
plot(RefClu,'k.')
xlabel(sprintf('%s (black) %s (green)',RefFileName,FileName))
ylabel('Clu Rank')
title('Clu Ranks')
set(gca,'box','on')
subplot(2,2,2)
hold on
plot(RefClu,Clu,'m.','markersize',3)
title('CluRanks')
set(gca,'box','on')
xlabel(RefFileName)
ylabel(FileName)
set(gca,'box','on')
subplot(2,2,3)
hold on
plot(CliqueSize,'g.','markersize',3)
plot(RefCliqueSize,'k.')
title('Clu Size')
xlabel(sprintf('%s (black) %s (green)',RefFileName,FileName))
ylabel('CluSize')
set(gca,'box','on')
subplot(2,2,4)
plot(RefCliqueSize,CliqueSize,'m.','markersize',3)
title('Clu Size')
set(gca,'box','on')
xlabel(RefFileName)
ylabel(FileName)

%% SEARCH ADJACENT CLIQUES
function SEARCH_ADJACENT_CLIQUES(ModelRank,InterNetRank,InterCorrLimit,UnionNetRank,PsNb,UnionCorrLimit)
%net_properties('search adjacent cliques',8,54,0,55)
global K
%Seach cliques in union network that are adjacent to cliques found in intersect network:
%For the current intersect clique, all the probe set that are linked to all the probe sets of the intersect cliques are selected
%and cliques are searched among these probe sets. So each of the resulting cliques and the current intersect clique form a clique.
%When all the intersect cliques (above a given size, in general equal to 4), are processed,
%cliques are searched in the residual probe sets.

CLIQUE_SIZE=5;
NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',InterNetRank));
cd(NetDir)
eval(sprintf('load m%un%u_cliques_%02u',ModelRank,InterNetRank,InterCorrLimit))
% %patch (to do update search cliques results)
% if size(CluSize,2)==1
%     CluSize=[[1:length(CluSize)]',CluSize];
% end
% 
% %patch
% CliqueSize=zeros(PsNb,1);
% for i=1:PsNb
%     if Clu(i)>0
%         CliqueSize(i)=CluSize(Clu(i),2);
%     end
% end
%keep information on intersect cliques
IntersectClu=Clu;
IntersectCluSize=CluSize;
IntersectCliqueSize=CliqueSize;
%Cliques ranks are colinear (CluSize(:,1)=1:length(CluSize)
IntersectCluRanks=find(CluSize>=CLIQUE_SIZE);
IntersectIsSeed=IsSeed;
UnionIsSeed=IsSeed;

NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',UnionNetRank));
cd(NetDir)
%eliminate small intersect cliques to allow their probe sets to be clusterized in union cliques
Clu(find(CliqueSize<CLIQUE_SIZE))=0;
CluSize(find(CluSize<CLIQUE_SIZE),:)=[];
CliqueSize(find(CliqueSize<CLIQUE_SIZE))=0;
MaxClu=max(Clu);

% %renumbering of cliques
% CurrCliques=setdiff(unique(Clu),0);
% for PsL=1:length(Clu)
%     if Clu(PsL)>0
%         Clu(PsL)=find(CurrCliques==Clu(PsL));
%     end
% end

UnionClu=Clu;
UnionCluSize=CluSize;
UnionCliqueSize=CliqueSize;


UnionParentClu=[];
UnionCluRank=max(Clu);


for CluL=1:length(IntersectCluRanks)
    sprintf('Clul: %u\n',CluL)
    if length(setdiff(unique(UnionClu),0))~=max(UnionClu)
        'stop'
    end
    %find the probe sets of the union network that are linked to all the probe
    %sets of the current clique of the intersection network
    CurrClu=find(IntersectClu==IntersectCluRanks(CluL));
    if isempty(CurrClu)
        'stop'
    end
    Corr=load_data(sprintf('c_m%u_n%u.4mat',ModelRank,UnionNetRank),NetDir,PsNb,PsNb,'uint8','ieee-le',1:PsNb,CurrClu);
    Anti=load_data(sprintf('a_m%u_n%u.4mat',ModelRank,UnionNetRank),NetDir,PsNb,PsNb,'uint8','ieee-le',1:PsNb,CurrClu);
    Corr=Corr>UnionCorrLimit&Corr>Anti;
    clear Anti
    %find probe set that are linked to all the probe set
    %of the current intersect cliques
    PsRanks=find(sum(Corr,2)==size(Corr,2));
    %remove all probe sets found in previous steps
    PsRanks=setdiff(PsRanks,find(UnionClu));
    if ~isempty(PsRanks)
        %find cliques inside thes PsRanks
        Continue=1;
    else
        Continue=0;
    end
    while Continue
        if length(PsRanks)<=1
            Continue=0;
        elseif length(PsRanks)==2
            CurrCluster=PsRanks;
            UnionCluRank=UnionCluRank+1;
            UnionClu(CurrCluster)=UnionCluRank;
            UnionParentClu(end+1,:)=[UnionCluRank,CluL];
            UnionCluSize(end+1,1)=length(CurrCluster);
            UnionCliqueSize(CurrCluster)=length(CurrCluster);
            Continue=0;
        else
            %load sub network
            Corr=load_data(sprintf('c_m%u_n%u.4mat',ModelRank,UnionNetRank),NetDir,PsNb,PsNb,'uint8','ieee-le',PsRanks,PsRanks);
            Anti=load_data(sprintf('a_m%u_n%u.4mat',ModelRank,UnionNetRank),NetDir,PsNb,PsNb,'uint8','ieee-le',PsRanks,PsRanks);
            Corr=Corr>UnionCorrLimit&Corr>Anti;
            clear Anti

            Conn=sum(Corr,2);
            [Temp,SortIndex]=sort(Conn);
            %write ASCII file for cliquer program
            fid=fopen(sprintf('m%un%u_union_cliquer.txt',ModelRank,UnionNetRank),'w');
            fprintf(fid,'c rank %u\n',CluL);
            fprintf(fid,'p edge %u %u\n',length(Corr),max(1,uint32((sum(sum(Corr))-length(Corr))/2)));
            for PsL=1:length(Corr)
                Pos=find(Corr(SortIndex(PsL),:));
                Pos=Pos(find(Pos>SortIndex(PsL)));
                if ~isempty(Pos)
                    for PosL=1:length(Pos)
                        fprintf(fid,'e %u %u\n',SortIndex(PsL),Pos(PosL));
                    end
                end
            end
            sprintf('length(Corr): %u\n',length(Corr))
            fclose(fid);
            PsPos=[];
            eval(sprintf('! /usr/local/cliquer/cl -su %s/m%03u/n%05u/m%un%u_union_cliquer.txt>%s/m%03u/n%05u/m%un%u_union_cliquer_res.txt;',K.dir.net,ModelRank,UnionNetRank,ModelRank,UnionNetRank,K.dir.net,ModelRank,UnionNetRank,ModelRank,UnionNetRank))
            eval(sprintf('PsPos=load(''m%un%u_union_cliquer_res.txt'');',ModelRank,UnionNetRank));
            %The maximal clique found (can be inferior to the real greatest clique because calculus is stoped after 60"
            if isempty(PsPos)
                Continue=0;
            else
                CurrCluster=PsRanks(PsPos);
                UnionCluRank=UnionCluRank+1;
                UnionClu(CurrCluster)=UnionCluRank;
                UnionParentClu(end+1,:)=[UnionCluRank,CluL];
                UnionCluSize(end+1,1)=length(CurrCluster);
                UnionCliqueSize(CurrCluster)=length(CurrCluster);
                PsRanks(PsPos)=[];
                if isempty(PsRanks)
                    Continue=0;
                end
            end
        end
    end
end %of Clul


%reorder cliques by size
Union=MaxClu+1:max(UnionClu);
SortOrder=[1:MaxClu]';
[Temp,USortOrder]=sort(UnionCluSize(Union),'descend');
SortOrder=[SortOrder;USortOrder+MaxClu];
[Temp NewRank]=sort(SortOrder);

MemClu=Clu;
for PsL=1:length(UnionClu)
    if UnionClu(PsL)>0
        UnionClu(PsL)=NewRank(UnionClu(PsL));
    end
end
MemUnionCluSize=UnionCluSize;
UnionCluSize=zeros(max(UnionClu),1);
for CluL=1:max(UnionClu)
    UnionCluSize(CluL)=length(find(UnionClu==CluL));
end
for PsL=1:length(UnionCliqueSize)
    if UnionClu(PsL)>0
        UnionCliqueSize(PsL)=length(find(UnionClu==UnionClu(PsL)));
    end
end
MemUnionParentClu=UnionParentClu;
for i=1:length(UnionParentClu)
    UnionParentClu(i,1)=NewRank(UnionParentClu(i,1));
    UnionParentClu(i,2)=NewRank(UnionParentClu(i,2));
end
UnionMaxClu=max(UnionClu);

%process residuals probe set
PsRanks=find(UnionCliqueSize<CLIQUE_SIZE);
UnionClu(PsRanks)=0;
Corr=load_data(sprintf('c_m%u_n%u.4mat',ModelRank,UnionNetRank),NetDir,PsNb,PsNb,'uint8','ieee-le',PsRanks,PsRanks);
Corr=uint8(Corr>UnionCorrLimit);
CurrConn=sum(Corr,2);
clear Corr


FoundPsRanks=setdiff(find(UnionClu),PsRanks);
LoadedPsRanks=PsRanks;
CliqueRank=max(UnionClu);
ClearedPsRanks=[CliqueRank:max(UnionParentClu(:,1))];
for PsL=1:length(ClearedPsRanks)
    Pos=find(UnionParentClu(:,1)==ClearedPsRanks(PsL));
    if ~isempty(Pos)
        UnionParentClu(Pos,:)=[];
    end
end


[UClu,UCluSize,UCliqueSize,UIsSeed]=CLIQUES(ModelRank,UnionNetRank,LoadedPsRanks,CurrConn,PsNb,UnionCorrLimit,UnionClu,UnionCluSize,UnionCliqueSize,UnionIsSeed,CliqueRank,0,0,FoundPsRanks);
%CLIQUES(ModelRank,NetRank,LoadedPsRanks,CurrConn,PsNb,CorrLimit,Clu,CluSize,CliqueSize,IsSeed,CliqueRank,MyAlgoFlag,DisplayFlag,varargin)
cd(NetDir)
Clu=UClu;

%reorder cliques by size
Union=CliqueRank+1:max(Clu);
SortOrder=[1:CliqueRank]';
[Temp,USortOrder]=sort(UCluSize(Union),'descend');
SortOrder=[SortOrder;USortOrder+CliqueRank];
[Temp NewRank]=sort(SortOrder);

MemClu=Clu;
for PsL=1:length(Clu)
    if Clu(PsL)>0
        Clu(PsL)=NewRank(Clu(PsL));
    end
end

CluSize=zeros(max(Clu),1);
for CluL=1:max(Clu)
    CluSize(CluL)=length(find(Clu==CluL));
end

CliqueSize=zeros(PsNb,1);
for PsL=1:length(CliqueSize)
    if Clu(PsL)>0
        CliqueSize(PsL)=length(find(Clu==Clu(PsL)));
    end
end

eval(sprintf('save m%un%u_cliques_%02u Clu CliqueSize CluSize UnionParentClu',ModelRank,UnionNetRank,UnionCorrLimit))



%% SEARCH QUASI-CLIQUES
function SEARCH_QUASI_CLIQUES(ModelRank,NetRank,PsNb,CorrLimit,varargin)
global K
%net_properties('search quasi-cliques',8,55,0,8,54,5,0)
%net_properties('search quasi-cliques',11,42,0,41,5,0)

RefFlag=0;
if nargin==7
    %load already calculated quasi-clique
    %cliques to be merged are merged together (to form quasi-clique)
    %and not with these loaded quasi-cliques
    RefNetRank=varargin{1};
    RefCliqueSize=varargin{2};
    RefCorrLimit=varargin{3};
    RefFlag=1;
end

DISPLAY_IT=1;
CLIQUE_SIZE=5;
FACTOR=2;

DoIt=1;
% if exist(sprintf('m%un%u_quasi-cliques_size%u_%02u.mat',ModelRank,NetRank,CLIQUE_SIZE,CorrLimit),'file')
%     load(sprintf('m%un%u_quasi-cliques_size%u_%02u.mat',ModelRank,NetRank,CLIQUE_SIZE,CorrLimit))
%     DoIt=0;
% else
if RefFlag
    RefNetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',RefNetRank));
    cd(RefNetDir)

    if exist(sprintf('m%un%u_cliques_%02u.mat',ModelRank,RefNetRank,RefCorrLimit),'file')
        load(sprintf('m%un%u_cliques_%02u.mat',ModelRank,RefNetRank,RefCorrLimit))
        Pos=find(CluSize>=RefCliqueSize);
        MaxRefClu=max(Pos);
        clear Cl* IsSeed
    else
        h=errordlg(sprintf('m%un%u_cliques_%02u.mat',ModelRank,RefNetRank,RefCorrLimit));
        waitfor(h)
        error('process canceled')
    end
    if exist(sprintf('m%un%u_quasi-cliques_size%u_%02u.mat',ModelRank,RefNetRank,RefCliqueSize,RefCorrLimit),'file')
        load(sprintf('m%un%u_quasi-cliques_size%u_%02u.mat',ModelRank,RefNetRank,RefCliqueSize,RefCorrLimit))
        clear QCliOverlap
    else
        h=errordlg(sprintf('m%un%u_quasi-cliques_size%u_%02u.mat',ModelRank,RefNetRank,RefCliqueSize,RefCliqueSize,RefCorrLimit));
        waitfor(h)
        error('process canceled')
    end
end

NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',NetRank));
cd(NetDir)
if exist(sprintf('m%un%u_cliques_%02u.mat',ModelRank,NetRank,CorrLimit),'file')
    load(sprintf('m%un%u_cliques_%02u.mat',ModelRank,NetRank,CorrLimit))
else
    h=errordlg(sprintf('m%un%u_cliques_%02u.mat does not exist. Run find cliques.',ModelRank,NetRank,CorrLimit));
    waitfor(h)
    error('process canceled')
end
if exist(sprintf('m%un%u_cliques_properties_%02u.mat',ModelRank,NetRank,CorrLimit),'file')
    load(sprintf('m%un%u_cliques_properties_%02u.mat',ModelRank,NetRank,CorrLimit))
else
    h=errordlg(sprintf('m%un%u_cliques_properties_%02u.mat does not exist. Run cliques properties.',ModelRank,NetRank,CorrLimit));
    waitfor(h)
    error('process canceled')
end
if exist(sprintf('m%un%u_cliques_anticorr_%02u.mat',ModelRank,NetRank,CorrLimit),'file')
    load(sprintf('m%un%u_cliques_anticorr_%02u.mat',ModelRank,NetRank,CorrLimit))
else
    h=errordlg(sprintf('m%un%u_cliques_anticorr_%02u.mat does not exist. Run ANTI CORR.',ModelRank,NetRank,CorrLimit));
    waitfor(h)
    error('process canceled')
end
if exist(sprintf('m%un%u_cliques_overlap_%02u.mat',ModelRank,NetRank,CorrLimit),'file')
    load(sprintf('m%un%u_cliques_overlap_%02u.mat',ModelRank,NetRank,CorrLimit))
else
    h=errordlg(sprintf('m%un%u_cliques_overlap_%02u.mat does not exist. Run ANTI CORR.',ModelRank,NetRank,CorrLimit));
    waitfor(h)
    error('process canceled')
end


%remove Clique<CLIQUE_SIZE
% CluList=1:length(CluSize);
% ClearPos=find(CliqueSize<CLIQUE_SIZE);
% if RefFlag
%     for CluL=1:MaxRefClu
%         ClearPos=[ClearPos;find(Clu==CluL)];
%     end
%     ClearPos=unique(ClearPos);
% end
% Clu(ClearPos)=0;
% ClearPos=find(CluSize<CLIQUE_SIZE);
% 
% if RefFlag
%     ClearPos=unique([ClearPos;[1:MaxRefClu]']);
% end
% CluSize(ClearPos)=[];
% CluList(ClearPos)=[];
% 
% CluAntiCorr=CluAntiCorr(CluList,CluList);
% CluOverlap=CluOverlap(CluList,CluList);


% end
%find LIMITS
if RefFlag
    %don't take in account ref clique that have another distributions
    CluSize(1:MaxRefClu)=0;
end
Limits=zeros(max(CluSize),2);

CluSizes=setdiff(unique(CluSize),0);
for CluL1=1:length(CluSizes)
    %sprintf('size nb %u',CluL1)
    %position of the cliques that have a size equal to the current size
    CluPos=find(CluSize==CluSizes(CluL1));
    %overlap values between the current cliques and all other cliques
    Val=CluOverlap(CluPos,:);
    %don't count the overlap inside a clique (=100)
    Val(:,CluPos)=0;
    %position of maximal overlap value in each clique
    [MaxOverlap,MaxPos]=max(Val,[],2);
    %eliminate the maximal value to calculate statistics on other values
    Val(:,MaxPos)=0;
    OverlapVal=Val(:);
    if length(MaxOverlap)==1
        Limits(CluSizes(CluL1),1)=MaxOverlap-1;
    else
        Limits(CluSizes(CluL1),1)=median(MaxOverlap);
    end
    OverlapVal=sort(OverlapVal);
    Limits(CluSizes(CluL1),2)=OverlapVal(round(length(OverlapVal)*0.8));
end

CluLimits=zeros(max(Clu),2);
for CluL=1:max(Clu)
    %     Pos=find(CluList==CluL);
    %     if ~isempty(Pos)
    %     CluLimits(CluL,1)=Limits(CluSize(Pos),1);
    %     CluLimits(CluL,2)=Limits(CluSize(Pos),2);
    %     end
    if CluSize(CluL)>0
    CluLimits(CluL,1)=Limits(CluSize(CluL),1);
    CluLimits(CluL,2)=Limits(CluSize(CluL),2);
    end
end


%ANTI cs CORR for cliques
if DISPLAY_IT
    h=figure;
    set(gcf,'color',[1,1,1])
    set(h,'name',sprintf('Limits for merging cliques into quasi-cliques (m%un%u)',ModelRank,NetRank))
    hold on
    Pos=find(Limits(:,1));
    plot(CluSizes,Limits(Pos,1),'g.')
    plot(CluSizes,Limits(Pos,1),'g')
    plot(CluSizes,Limits(Pos,2),'r.')
    plot(CluSizes,Limits(Pos,2),'r')
    title(sprintf('Limits for merging cliques into quasi-cliques (m%un%u)',ModelRank,NetRank))
    set(gca,'box','on')
    xlabel('smallest clique size')
    ylabel('inclusion (median(max)-green) and exclusion (80th perc - red) limits')



    h=figure;
    set(h,'name',sprintf('M%uN%u - ANTI VS CORR FOR CLIQUES',ModelRank,NetRank))
    PosCorr=find(triu(CluAntiCorr));
    PosAnti=find(tril(CluAntiCorr));
    subplot(1,3,1)
    hist(CluAntiCorr(PosAnti),100);
    set(gca,'xlim',[0,50])
    title('ANTI DISTRIBUTION')

    subplot(1,3,2)
    hist(CluAntiCorr(PosCorr),100);
    set(gca,'xlim',[0,50])
    title('CORR DISTRIBUTION')

    CluCorrAnti=CluAntiCorr';
    subplot(1,3,3)
    plot(CluAntiCorr(PosCorr),CluCorrAnti(PosCorr),'g.','markersize',3)
    set(gca,'xlim',[0,80])
    title('ANTI vs CORR')
    set(gcf,'color',[1,1,1])

    h=figure;
    set(h,'name',sprintf('M%uN%u - SELECTED CANDIDATES',ModelRank,NetRank))
    set(gcf,'color',[1,1,1])

end
%MERGE QUASI CLIQUES
if DoIt


    %cliques
    %FIRST ROUND AGREGATION IN ORDER OF DECREASING CLIQUE SIZE

    CliFid=fopen(sprintf('m%un%u_quasi-cliques_size%u_%s.txt',ModelRank,NetRank,CLIQUE_SIZE,date),'w');
    fprintf(CliFid,'\tbefore cleaning\t\tafter cleaning\t\t\n');
    fprintf(CliFid,'round\t#q\t#c\t#q\t#c\t#next\n');

    SizeFid=fopen(sprintf('m%un%u_size_quasi-cliques_size%u_%s.txt',ModelRank,NetRank,CLIQUE_SIZE,date),'w');
    fprintf(SizeFid,'round\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t>10\n');


    
    clear InClu;
    FullCluOverlap=CluOverlap;
    if RefFlag==0
        QClique={};
    end


    MergedClu=[];
    DisplayRank=0;

    if RefFlag
        %eliminate reference to cliques that are not used
        for QCliL=1:length(QClique);
            Pos=find(QClique{QCliL}>MaxRefClu);
            if ~isempty(Pos)
                QClique{QCliL}(Pos)=[];
            end
        end
        for QCliL=length(QClique):-1:1
            if isempty(QClique{QCliL})
                QClique(QCliL)=[];
            end
        end                        
        QCliqueRank=length(QClique);
        StartQCli=1+QCliqueRank;
        StartCli=1+MaxRefClu;
    else
        QCliqueRank=0;
        StartQCli=1;
        StartCli=1;
    end
    AllCliques=[];
    for i=1:length(QClique)
        AllCliques=[AllCliques;QClique{i}];
    end
    
    NotSearchedClu=ones(length(CluSize),1);  
    %mark as searched all cliques that are <CLIQUE_SIZE and which rank is <=MaxCliPos
    LessMinSize=find(CluSize<CLIQUE_SIZE);
    if RefFlag
        %mark as searched all cliques already merged in the reference network
        LessMinSize=unique([LessMinSize;[1:MaxRefClu]']);
    end
    NotSearchedClu(LessMinSize)=0;
    NotSearchedNb=length(find(NotSearchedClu));
    Continue=1;
    RoundL=0;
    NotUsedClu=ones(length(CluSize),1);


    while Continue

        RoundL=RoundL+1;
        MemCluOverlap=CluOverlap;

        MergedNb=[];
        CurrMergedClu=[];

        while ~isempty(find(NotSearchedClu))           
            UseIt=0;
            CurrCluPos=find(NotSearchedClu);
            %use the largest clique not yet searched as initiator of a
            %quasi-clique
            CurrClu=CurrCluPos(1);
            %set up the inclusion and exclusion limits which depend on the
            %clique size

            MergeIt=1;
            %don'use curr clique if it has significative overlap with previously
            %merged cliques
            if ~isempty(CurrMergedClu)
                for CluL=1:length(CurrMergedClu)
                    %if CluOverlap(CurrClu,CurrMergedClu(CluL))>OVERLAP_OUT
                    %if CluOverlap(CurrClu,CurrMergedClu(CluL))>FIND_LIMITS(min(length(find(Clu==CurrClu)),length(find(Clu==CurrMergedClu(CluL)))),OVERLAP_IN,OVERLAP_OUT,SIZE_TEST,'OUT')
                    if CluOverlap(CurrClu,CurrMergedClu(CluL))>Limits(min(length(find(Clu==CurrClu)),length(find(Clu==CurrMergedClu(CluL)))),2)
                        MergeIt=0;
                        break
                    end
                end
            end
            if MergeIt
                %find potential candidate for merging
                CurrVal=CluOverlap(CurrClu,:);
                %eliminate 100 value
                CurrVal(CurrClu)=0;
                %eliminate merged clu
                CurrVal(MergedClu)=0;
                %Candidates=find(CurrVal>=OVERLAP_IN);
                %Candidates=find(CurrVal>=min(OVERLAP_IN));
                Candidates=find(CurrVal'>=CluLimits(:,1));
                Candidates=intersect(Candidates,find(NotSearchedClu));
                CurrVal=CurrVal(Candidates);
                KeepIndex=[];
                for CandL=1:length(Candidates)
                    %if (CluAntiCorr(min(CurrClu,Candidates(CandL)),max(CurrClu,Candidates(CandL)))>=FACTOR*CluAntiCorr(max(CurrClu,Candidates(CandL)),min(CurrClu,Candidates(CandL))))&CurrVal(CandL)>=FIND_LIMITS(min(length(find(Clu==CurrClu)),length(find(Clu==Candidates(CandL)))),OVERLAP_IN,OVERLAP_OUT,SIZE_TEST,'IN')
                    if (CluAntiCorr(min(CurrClu,Candidates(CandL)),max(CurrClu,Candidates(CandL)))>=FACTOR*CluAntiCorr(max(CurrClu,Candidates(CandL)),min(CurrClu,Candidates(CandL))))&CurrVal(CandL)>=Limits(min(length(find(Clu==CurrClu)),length(find(Clu==Candidates(CandL)))),1)
                        KeepIndex=[KeepIndex,CandL];
                    end
                end
                Candidates= Candidates(KeepIndex);
                if ~isempty(Candidates)
                    if DISPLAY_IT
                        DisplayRank=DisplayRank+1;
                        if DisplayRank<=9
                            figure(h);
                            subplot(3,3,DisplayRank)
                            hold on
                            plot(CluOverlap(CurrClu,Candidates),'m')
                            plot(CluOverlap(CurrClu,Candidates),'b.')
                            set(gca,'box','on')
                            title(sprintf('Clique %u',CurrClu))
                            xlabel('candidate cliques (b)')
                            ylabel('overlaping (m), selected cliques(r)')
                            Position=get(gcf,'position');
                            Position(3)=370;
                            Position(4)=270;
                            set(gcf,'position',Position)
                            set(gcf,'color',[1,1,1])
                        end
                    end
                end
                Merged=[];
                if ~isempty(Candidates)

                    for CandL=1:length(Candidates)

                        MergeIt=0;
                        %the current candidate must have a significative overlap with the seeding clique
                        %if CluOverlap(Candidates(CandL),CurrClu)<OVERLAP_IN
                        %if (CluAntiCorr(min(CurrClu,Candidates(CandL)),max(CurrClu,Candidates(CandL)))>=FACTOR*CluAntiCorr(max(CurrClu,Candidates(CandL)),min(CurrClu,Candidates(CandL))))&CluOverlap(Candidates(CandL),CurrClu)>=FIND_LIMITS(min(length(find(Clu==CurrClu)),length(find(Clu==Candidates(CandL)))),OVERLAP_IN,OVERLAP_OUT,SIZE_TEST,'IN')
                        if (CluAntiCorr(min(CurrClu,Candidates(CandL)),max(CurrClu,Candidates(CandL)))>=FACTOR*CluAntiCorr(max(CurrClu,Candidates(CandL)),min(CurrClu,Candidates(CandL))))&CluOverlap(Candidates(CandL),CurrClu)>=Limits(min(length(find(Clu==CurrClu)),length(find(Clu==Candidates(CandL)))),1)
                            MergeIt=1;
                        end
                        if MergeIt
                            if ~isempty(Merged)
                                % with all cliques of the quasi clique being constructed
                                for MergeL=1:length(Merged)
                                    %if CluOverlap(Candidates(CandL),Merged(MergeL))<OVERLAP_IN
                                    %if (CluAntiCorr(min(CurrClu,Candidates(CandL)),max(CurrClu,Candidates(CandL)))<FACTOR*CluAntiCorr(max(CurrClu,Candidates(CandL)),min(CurrClu,Candidates(CandL))))|CluOverlap(Candidates(CandL),Merged(MergeL))<FIND_LIMITS(min(length(find(Clu==Merged(MergeL))),length(find(Clu==Candidates(CandL)))),OVERLAP_IN,OVERLAP_OUT,SIZE_TEST,'IN')
                                    if (CluAntiCorr(min(CurrClu,Candidates(CandL)),max(CurrClu,Candidates(CandL)))<FACTOR*CluAntiCorr(max(CurrClu,Candidates(CandL)),min(CurrClu,Candidates(CandL))))|CluOverlap(Candidates(CandL),Merged(MergeL))<Limits(min(length(find(Clu==Merged(MergeL))),length(find(Clu==Candidates(CandL)))),1)
                                        MergeIt=0;
                                        break
                                    end
                                end
                            end
                        end
                        if MergeIt
                            if ~isempty(CurrMergedClu)
                                %the current candidate must have a not
                                %significative overlap with all the cliques
                                %already merged in quasi cliques
                                for CluL=1:length(CurrMergedClu)
                                    %if CluOverlap(Candidates(CandL),CurrMergedClu(CluL))>OVERLAP_OUT
                                    %if CluOverlap(Candidates(CandL),CurrMergedClu(CluL))>FIND_LIMITS(min(length(find(Clu==CurrMergedClu(CluL))),length(find(Clu==Candidates(CandL)))),OVERLAP_IN,OVERLAP_OUT,SIZE_TEST,'OUT')
                                    if CluOverlap(Candidates(CandL),CurrMergedClu(CluL))>Limits(min(length(find(Clu==CurrMergedClu(CluL))),length(find(Clu==Candidates(CandL)))),2)
                                        MergeIt=0;
                                        break
                                    end
                                end
                            end
                            if MergeIt
                                UseIt=1;
                                Merged=[Merged;Candidates(CandL)];
                            end
                        end
                    end
                else
                    %no other clique agregated => quasiclique will contain only
                    %the intiating clique
                    %UseIt=1;
                    NotSearchedClu(CurrClu)=0;
                    NotUsedClu(CurrClu)=1;
                end
            end
            if UseIt
                QCliqueRank=QCliqueRank+1;
                MergedNb(end+1,1)=length(Merged);
                %order candidate by size
                if ~isempty(Merged)
                    if DISPLAY_IT &  DisplayRank<=9
                        figure(h)
                        figure(h);
                        subplot(3,3,DisplayRank)
                        for CliL=1:length(Merged)
                            Pos=find(Candidates==Merged(CliL));
                            plot(Pos,CluOverlap(CurrClu,Candidates(Pos)),'ro','markersize',10)
                        end
                    end
                    Merged=[CurrClu;Merged];
                    [CandSize SizeOrder]=sort(CluSize(Merged),'descend');
                    Merged=Merged(SizeOrder);
                    if ~isempty(intersect(AllCliques,Merged))
                        QCliqueRank
                        Merged
                        intersect(AllCliques,Merged)
                        'stop'
                    else
                        AllCliques=[AllCliques;Merged];
                        QClique{QCliqueRank,1}=Merged;
                    end
                    CluOverlap(QClique{QCliqueRank},QClique{QCliqueRank})=0;
                    CurrMergedClu=[CurrMergedClu;QClique{QCliqueRank}];
                    MergedClu=[MergedClu;QClique{QCliqueRank}];
                else
                    if ~isempty(intersect(AllCliques,CurrClu))
                        QCliqueRank
                        CurrClu
                        intersect(AllCliques,CurrClu)
                        'stop'
                    else
                        AllCliques=[AllCliques;CurrClu];
                        QClique{QCliqueRank,1}=CurrClu;
                    end
                    CluOverlap(CurrClu,CurrClu)=0;
                    Merged=CurrClu;
                    CurrMergedClu=[CurrMergedClu;CurrClu];
                    MergedClu=[MergedClu;CurrClu];
                end
                NotSearchedClu(QClique{QCliqueRank})=0;
                NotUsedClu(QClique{QCliqueRank})=0;
            else
                NotSearchedClu(CurrClu)=0;
                NotUsedClu(CurrClu)=1;
            end
        end

        %SECOND ROUND
        %CLEANING : remove cliques that have significative overlap between several
        %quasi cliques




        %INTER CLIQUE OVERLAPING
        CluOverlap=MemCluOverlap;
        % clear MemCluOverlap;
        CliNb=[];
        CliList=[];
        if RefFlag
            CurrStartCli=1;
        else
            CurrStartCli=1;
        end
            
        for QCliL=CurrStartCli:length(QClique)
            CliNb(QCliL-CurrStartCli+1,1)=length(QClique{QCliL});
            CliList=[CliList;QClique{QCliL}];            
        end
        QCliLimit=cumsum(CliNb);
        TotalCliNb=sum(CliNb);
        BeforeQNb=length(QClique);
        BeforeCNb=TotalCliNb;
        if UseIt==0
            AfterQNb=length(QClique);
            AfterCNb=TotalCliNb;
        end
        if length(CliList)~=length(unique(CliList));
            'stop'
        end
        if UseIt
            CliOverlap=FullCluOverlap(CliList,CliList);
            CliLimits=CluLimits(CliList,:);
            for CliL=1:TotalCliNb
                CliOverlap(CliL,CliL)=100;
            end

            if ~isempty(CliList)
                if DISPLAY_IT & RoundL==3
                    h1=figure;
                    set(h1,'name',sprintf('MERGED CLIQUES AT ROUND %u',RoundL));
                    image(CliOverlap)
                    set(gca,'xtick',cumsum(CliNb)+0.5)
                    set(gca,'xticklabel','')
                    set(gca,'ytick',cumsum(CliNb)+0.5)
                    set(gca,'yticklabel','')
                    set(gcf,'color',[1,1,1])
                    set(gca,'tickdir','out')
                    title(sprintf('merged cliques before cleaning (size>=%u, round %u)',CLIQUE_SIZE,RoundL))
                    set(gcf,'color',[1,1,1])
                end
            end


            %SEARCH CLIQUES TO BE CLEARED FROM QCLIQUES

            CliRank=StartCli-1;
            Pointer=CliRank;
            QCliDelete=[];
            QCliPosDelete=[];
            for QCliL=StartQCli:length(QClique)
                CurrQuasiCli=QClique{QCliL};
                for CliL=1:length(CurrQuasiCli)
                    RefCli=CurrQuasiCli(CliL);
                    CliRank=CliRank+1;
                    try
                        CurrVal=CliOverlap(CliRank,:);
                    catch
                        'stop'
                    end

                    %find clique outside the current quasi clique that
                    %have significative overlap with the current clique

                    %TestPos=find(CurrVal>OVERLAP_OUT);
                    %TestPos=find(CurrVal>min(OVERLAP_OUT));
                    TestPos=find(CurrVal'>CliLimits(:,2));
                    TestCli=CliList(TestPos);
                    %remove cliques that are in the current QClique
                    ClearPos=find(TestPos<Pointer+length(QClique{QCliL})+1);
                    TestPos(ClearPos)=[];
                    TestCli(ClearPos)=[];

                    if ~isempty(TestPos)
                        %process each clique
                        %mean overlap of the current clique inside the current quasi
                        %clique
                        Val=CliOverlap(CliRank,Pointer+1:Pointer+length(CurrQuasiCli));
                        %remove self overlap = 100
                        Val(CliL)=[];
                        %list of mean overlap values
                        MeanOverlap=mean(Val);
                        %list of concerned quasi cliques
                        QCli=QCliL;
                        %list of positions in quasi cliques
                        QCliPos=CliL;
                        for PosL=1:length(TestPos)
                            %verify that the position must be processed
                            %if CliOverlap(CliRank,TestPos)>FIND_LIMITS(min(length(find(Clu==RefCli)),length(find(Clu==TestCli(PosL)))),OVERLAP_IN,OVERLAP_OUT,SIZE_TEST,'OUT')
                            if CliOverlap(CliRank,TestPos)>Limits(min(length(find(Clu==RefCli)),length(find(Clu==TestCli(PosL)))),2)
                                %find the quasi clique of the currently processed clique that have
                                %significative overlap wwith the current clique
                                CurrQCli=find(QCliLimit>=TestPos(PosL));
                                CurrQCli=CurrQCli(1);
                                InQCliPos=find(QClique{CurrQCli}==TestCli(PosL));
                                Val=CliOverlap(TestPos(PosL),QCliLimit(CurrQCli-1)+1:QCliLimit(CurrQCli-1)+length(QClique{CurrQCli}));
                                Val(InQCliPos)=0;
                                MeanOverlap=[MeanOverlap;mean(Val)];
                                QCli=[QCli;CurrQCli];
                                QCliPos=[QCliPos;InQCliPos];
                            end
                        end
                        %find min value
                        [MaxOverlap,MaxPos]=max(MeanOverlap);
                        if MaxPos==1
                            %delete the max position
                            QCli(MaxPos)=[];
                            QCliPos(MaxPos)=[];
                            %and keep other cliques in list to be deleted
                            if ~isempty(QCli)
                                for i=1:length(QCli)
                                    if isempty(find(QCliDelete==QCli(i)&QCliPosDelete==QCliPos(i)))
                                        QCliDelete=[QCliDelete;QCli(i)];
                                        QCliPosDelete=[QCliPosDelete;QCliPos(i)];
                                    end
                                end
                            end
                        else
                            if isempty(find(QCliDelete==QCli(1)&QCliPosDelete==QCliPos(1)))
                                QCliDelete=[QCliDelete;QCli(1)];
                                QCliPosDelete=[QCliPosDelete;QCliPos(1)];
                            end
                        end
                    end
                end
                Pointer=Pointer+length(CurrQuasiCli);
            end

            %DELETE MULTI OVERLAPING CLIQUES
            %order the delete information
            [QCliPosDelete,SortOrder]=sort(QCliPosDelete);
            QCliDelete=QCliDelete(SortOrder);
            [QCliDelete,SortOrder]=sort(QCliDelete);
            QCliPosDelete=QCliPosDelete(SortOrder);

            if ~isempty(QCliDelete)
                for CliL=length(QCliDelete):-1:1
                    %record the deleted clique as not used
                    NotUsedClu(QClique{QCliDelete(CliL)}(QCliPosDelete(CliL)))=1;
                    QClique{QCliDelete(CliL)}(QCliPosDelete(CliL))=[];
                end


                %eliminate empty quasi cliques
                KeepIt=[];
                for QCliL=1:length(QClique)
                    if length(QClique{QCliL})>0
                        KeepIt=[KeepIt;QCliL];
                    end
                end
                Mem=QClique;
                QClique={};
                for QCliL=1:length(KeepIt)
                    QClique{end+1,1}=Mem{KeepIt(QCliL)};
                end

                %recount cliques
                CliNb=[];
                CliList=[];
                for QCliL=1:length(QClique)
                    CliNb(QCliL,1)=length(QClique{QCliL});
                    CliList=[CliList;QClique{QCliL}];
                end
                if length(CliList)~=length(unique(CliList));
                    'stop'
                end
                TotalCliNb=sum(CliNb);
            end
            AfterQNb=length(QClique);           
            AfterCNb=TotalCliNb;

            %         %ADD ALL THE NOT USED CLIQUES IN THE LAST POSITION (it is not a quasi
            %         %clique)
            %         NextPos=find(NotUsedClu);
            %
            %         %UPDATE NotSearched to search only cliques that where not merged
            %         NotSearchedClu=zeros(max(Clu),1);
            %         NotSearchedClu(NextPos)=1;
            %         NotSearchedClu(LessMinSize)=0;
        end

        % TEST IF MERGING MUST BE STOPPED

        if length(find(NotSearchedClu))==NotSearchedNb % stop merging
            Continue=0;
            AddClu=find(NotUsedClu);
            if RefFlag
                AddClu=setdiff(AddClu,[1:MaxRefClu]);
            end
            if ~isempty(AddClu)
                for CluL=1:length(AddClu)
                    QClique{end+1,1}=AddClu(CluL);
                end
            end
            fprintf(CliFid,'%u\t%u\t%u\t%u\t%u\t%u\n',RoundL,BeforeQNb,BeforeCNb,AfterQNb,AfterCNb,size(AddClu));

            %calculate ps nb in two classes of clusters with size >1 or =1
            CliList=[];
            CluSize=[];
            for QCliL=1:length(QClique)
                CliList=[CliList;QClique{QCliL}];
                if length(CliList)~=length(unique(CliList))
                    'stop'
                end
                CluSize(end+1,1)=length(QClique{QCliL});
            end

            Pos=find(CluSize>1);
            ProbesetNb=0;
            CliNb=0;
            for CluL=1:length(Pos)
                CliNb=CliNb+length(QClique{CluL});
                for CliL=1:length(QClique{CluL})
                    ProbesetNb=ProbesetNb+length(find(Clu==(QClique{CluL}(CliL))));
                end
            end
            fprintf(CliFid,'\n');
            fprintf(CliFid,'C size\t#q\t#c\t#ps\n');
            fprintf(CliFid,'>1\t%u\t%u\t%u\n',length(Pos),CliNb,ProbesetNb);

            Pos=find(CluSize==1);
            ProbesetNb=0;
            CliNb=0;
            for CluL=1:length(Pos)
                CliNb=CliNb+length(QClique{CluL});
                for CliL=1:length(QClique{CluL})
                    ProbesetNb=ProbesetNb+length(find(Clu==(QClique{CluL}(CliL))));
                end
            end
            fprintf(CliFid,'1\t%u\t%u\t%u\n',length(Pos),CliNb,ProbesetNb);


            fclose(CliFid);
            fclose(SizeFid);


            if length(CliList)~=length(unique(CliList))
                h=errordlg('CliList has doublons');
                waitfor(h)
                error('process canceled')
            end

            CliOverlap=FullCluOverlap(CliList,CliList);
            h=figure;
            set(h,'name',sprintf('OVERLAP OF MERGED CLIQUES OF SIZE >= %u IN M%uN%u',CLIQUE_SIZE,ModelRank,NetRank))
            image(CliOverlap);
            title(sprintf('Overlap of cliques of size >= %u in m%un%u',CLIQUE_SIZE,ModelRank,NetRank))
            set(gcf,'color',[1,1,1])


            %display relation between qcliques
            QCliOverlap=ones(length(QClique))*100;
            for QCliL1=1:length(QClique)-1
                for QCliL2=QCliL1+1:length(QClique)
                    %calculate mean overlap between two qcliques
                    Overlap=FullCluOverlap(QClique{QCliL1},QClique{QCliL2});
                    QCliOverlap(QCliL1,QCliL2)=mean(Overlap(:));
                    QCliOverlap(QCliL2,QCliL1)=QCliOverlap(QCliL1,QCliL2);
                end
            end
            h=figure;
            set(h,'name',sprintf('MEAN OVERLAP OF GROUPS OF CLIQUES OF SIZE >= %u IN M%uN%u',CLIQUE_SIZE,ModelRank,NetRank))
            image(QCliOverlap);
            title(sprintf('mean overlap of groups of cliques of size >= %u in m%un%u',CLIQUE_SIZE,ModelRank,NetRank))
            set(gcf,'color',[1,1,1])

            eval(sprintf('save m%un%u_quasi-cliques_size%u_%02u QClique CliOverlap QCliOverlap',ModelRank,NetRank,CLIQUE_SIZE,CorrLimit))


        else  % START ANOTHER ROUND
            NextPos=find(NotUsedClu);
            NotSearchedNb=length(find(NotSearchedClu));
            fprintf(CliFid,'%u\t%u\t%u\t%u\t%u\t%u\n',RoundL,BeforeQNb,BeforeCNb,AfterQNb,AfterCNb,length(NextPos));
            %statisitc on clu size
            CliLengths=zeros(1,11);
            for QCliL=StartQCli:length(QClique)
                for CliL=1:length(QClique{QCliL})
                    CliLength=length(find(Clu==QClique{QCliL}(CliL)));
                    if CliLength>10
                        CliLengths(11)=CliLengths(11)+1;
                    else
                        CliLengths(CliLength)=CliLengths(CliLength)+1;
                    end
                end
            end
            fprintf(SizeFid,'%u',RoundL);
            for i=1:11
                fprintf(SizeFid,'\t%u',CliLengths(i));
            end
            fprintf(SizeFid,'\n');



            %     %RECALCULATE INTER CLIQUE OVERLAPING FOR ALL CLIQUES
            %
            CliNb=[];
            CliList=[];
            CliSize=[];
            for QCliL=1:length(QClique)
                CliNb(QCliL,1)=length(QClique{QCliL});
                CliList=[CliList;QClique{QCliL}];
                for CliL=1:length(QClique{QCliL})
                    CliSize=[CliSize;length(find(Clu==QClique{QCliL}(CliL)))];
                end
            end
            if length(CliList)~=length(unique(CliList));
                'stop'
            end
            if length(CliList)~=length(unique(CliList))
                h=errordlg('CliList has doublons');
                waitfor(h)
                error('process canceled')
            end

            StartQCli=length(QClique)+1;
            StartCli=length(CliList)+1;

            %don't use information of already clusterized cliques to allow
            %clustering of resting cliques
            CluOverlap(CliList,CliList)=0;


            if DISPLAY_IT & RoundL==3
                if ~isempty(NextPos)
                    CliSize=[];
                    CliNb(end+1,1)=length(NextPos);
                    CliList=[CliList;NextPos];
                    for CliL=1:length(NextPos)
                        CliSize=[CliSize;length(find(Clu==NextPos(CliL)))];
                    end
                end
                TotalCliNb=sum(CliNb);
                CliOverlap=zeros(TotalCliNb);
                for CliL=1:TotalCliNb
                    CliOverlap(CliL,CliL)=100;
                end

                CliOverlap=FullCluOverlap(CliList,CliList);
                h=figure;
                set(h,'name',sprintf('OVERLAP BETWEEN CLIQUES IN QUASI-CLIQUES (SIZE>=%u,ROUND %u)',CLIQUE_SIZE,RoundL));
                image(CliOverlap)
                set(gca,'xtick',cumsum(CliNb)+0.5)
                set(gca,'xticklabel','')
                set(gca,'ytick',cumsum(CliNb)+0.5)
                set(gca,'yticklabel','')
                set(gcf,'color',[1,1,1])
                set(gca,'tickdir','out')
                title(sprintf('OVERLAP BETWEEN CLIQUES IN QUASI-CLIQUES (SIZE>=%u,ROUND %u)',CLIQUE_SIZE,RoundL))
                set(gcf,'color',[1,1,1])


                h=figure;
                set(h,'name',sprintf('ANTI CORR BETWEENE CLIQUES IN QUASI-CLIQUES,ROUND %u)',RoundL));
                image(CluAntiCorr(CliList,CliList))
                title(sprintf('ANTI CORR BETWEEN CLIQUES IN QUASI-CLIQUES (ROUND %u)',RoundL))
                set(gcf,'color',[1,1,1])

            end
        end
    end
end %if DoIt
if DISPLAY_IT
    

    First=271;
    End=400;

    PsRanks=[];
    CliNb=0;
    for QCliL=First:End
        CliNb=CliNb+length(QClique{QCliL});
        for CluL=1:length(QClique{QCliL})
            PsRanks=[PsRanks;find(Clu==(QClique{QCliL}(CluL)))];
        end
    end
    length(PsRanks)
    CliNb

    a=find(CluSize(:)>=5);
    PsRanks=[];
    CliNb=0;
    for CliL=1:length(a)
        PsRanks=[PsRanks;find(Clu==a(CliL))];
    end

    length(PsRanks)



    C=load_data(sprintf('c_m%u_n%u.4mat',ModelRank,NetRank),'./',PsNb,PsNb,'uint8','ieee-le',PsRanks,PsRanks);
    C=triu(C);
    A=load_data(sprintf('a_m%u_n%u.4mat',ModelRank,NetRank),'./',PsNb,PsNb,'uint8','ieee-le',PsRanks,PsRanks);
    C=C+tril(A,-1);
    clear A
    h=figure;
    set(h,'name',sprintf('QUASI CLIQUES >= %s',CLIQUE_SIZE));
    xlabel('ANTI of ps in quasi-cliques')
    title(sprintf('CORR of ps in quasi-cliques of m%un%u',ModelRank,NetRank))

    image(C);
    set(gcf,'color',[1,1,1])

end



%% CLUSTERIZE QUASI-CLIQUES
function CLUSTERIZE_QUASI_CLIQUES(ModelRank,NetRank,PsNb,CorrLimit,varargin)
%net_properties('clusterize quasi-cliques',8,54,0)
%net_properties('clusterize quasi-cliques',8,55,50,54,5,0,1)
%net_properties('clusterize quasi-cliques',11,22,50,21,5,0,1)
%net_properties('clusterize quasi-cliques',11,42,0,41,5,0,1)
global K
RefFlag=0;
if nargin==8
    %load already calculated quasi-clique
    %cliques to be merged are merged together (to form quasi-clique)
    %and not with these loaded quasi-cliques
    RefNetRank=varargin{1};
    RefCliqueSize=varargin{2};
    RefCorrLimit=varargin{3};
    KeepRefFlag=varargin{4};
    if KeepRefFlag
        RefFlag=1;
    end
end

DISPLAY_IT=1;
CLIQUE_SIZE=5;

DoIt=1;
% if exist(sprintf('m%un%u_quasi-cliques_size%u_%02u.mat',ModelRank,NetRank,CLIQUE_SIZE,CorrLimit),'file')
%     load(sprintf('m%un%u_quasi-cliques_size%u_%02u.mat',ModelRank,NetRank,CLIQUE_SIZE,CorrLimit))
%     DoIt=0;
% else
if RefFlag
    RefNetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',RefNetRank));
    cd(RefNetDir)    
    if exist(sprintf('m%un%u_quasi-cliques_clusters_size%u_%02u.mat',ModelRank,RefNetRank,RefCliqueSize,RefCorrLimit),'file')
        load(sprintf('m%un%u_quasi-cliques_clusters_size%u_%02u.mat',ModelRank,RefNetRank,RefCliqueSize,RefCorrLimit))
        RefCliques=Cliques;
        RefCluOverlap=CluOverlap;
        RefClusters=Clusters;
        clear Cliques CluOverlap Clusters
    else
        h=errordlg(sprintf('m%un%u_quasi-cliques_clusters_size%u_%02u.mat does not exist',ModelRank,RefNetRank,RefCliqueSize,RefCorrLimit));
        waitfor(h)
        error('process canceled')
    end
end


NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',NetRank));
cd(NetDir)
DoIt=1;
% if exist(sprintf('m%un%u_quasi-cliques_clusters_size%u_%02u.mat',ModelRank,NetRank,CLIQUE_SIZE,CorrLimit),'file')
%     load(sprintf('m%un%u_quasi-cliques_clusters_size%u_%02u.mat',ModelRank,NetRank,CLIQUE_SIZE,CorrLimit))
%     DoIt=0;
% else
if exist(sprintf('m%un%u_quasi-cliques_size%u_%02u.mat',ModelRank,NetRank,CLIQUE_SIZE,CorrLimit),'file')
    load(sprintf('m%un%u_quasi-cliques_size%u_%02u.mat',ModelRank,NetRank,CLIQUE_SIZE,CorrLimit))
else
    h=errordlg('Do QUASI CLIQUES');
    waitfor(h)
    error('process canceled')
end
%reorder CliOverlap on cliques order (presently ordered on quasi cliques)
Cli=[];
for i=1:length(QClique)
    Cli=[Cli;QClique{i}];
end
[Temp,SortOrder]=sort(Cli);
CliOverlap=CliOverlap(SortOrder,SortOrder);
if RefFlag
    Cli=[];
    for i=1:length(RefCliques{end})
        Cli=[Cli;RefCliques{end}{i}];
    end
    figure
    image(CliOverlap(Cli,Cli))
end
if exist(sprintf('m%un%u_cliques_anticorr_%02u.mat',ModelRank,NetRank,CorrLimit),'file')
    load(sprintf('m%un%u_cliques_anticorr_%02u.mat',ModelRank,NetRank,CorrLimit))
else
    h=errordlg('Do CALCULATE ANTI CORR');
    waitfor(h)
    error('process canceled')
end
%CluAntiCorr is ordered on cliques
CliAntiCorr=CluAntiCorr;

if exist(sprintf('m%un%u_cliques_%02u.mat',ModelRank,NetRank,CorrLimit),'file')
    load(sprintf('m%un%u_cliques_%02u.mat',ModelRank,NetRank,CorrLimit))
else
    h=errordlg('Search cliques');
    waitfor(h)
    error('process canceled')
end

% end
LastClu=max(Clu);
Pos=find(CliqueSize<CLIQUE_SIZE);
CliqueSize(Pos)=0;
Clu(Pos)=0;
MaxClu=max(Clu);
CliOverlap=CliOverlap(1:MaxClu,1:MaxClu);
CliAntiCorr=CliAntiCorr(1:MaxClu,1:MaxClu);
CluSize=CluSize(1:MaxClu);

if DoIt
    
    %update QClique & QCliqueOverlap
    ClearPos=[];
    for QCliL=1:length(QClique)        
        if isempty(setdiff(QClique{QCliL},[MaxClu+1:LastClu]))
            ClearPos=[ClearPos,QCliL];
        end
    end
    if ~isempty(ClearPos)
        QClique(ClearPos)=[];
        QCliOverlap(ClearPos,:)=[];
        QCliOverlap(:,ClearPos)=[];
    end


    %quasi-cliques constructed by FIND QUASI-CLIQUES are
    %considered as the result of first found of clustering
    if RefFlag
        %find max ref clique
        MaxRefCli=0;
        for CliL=1:length(RefCliques{end})
            MaxRefCli=max(MaxRefCli,max(RefCliques{end}{CliL}));
        end
        %find position of max ref clique in QClique
        for QCliL=1:length(QClique)
            if max(QClique{QCliL})==MaxRefCli
                RefQCliquePos=QCliL;
                break
            end
        end
        
        CluNb=length(QClique)-RefQCliquePos;
        Clusters{1}=cell(CluNb,1);
        for i=1:CluNb
            Clusters{1}{i}=i+RefQCliquePos;
        end
    else
        CluNb=length(QClique);
        Clusters{1}=cell(CluNb,1);
        for i=1:CluNb
            Clusters{1}{i}=i;
        end
    end

    % overlap between cluster of first round
    CluOverlap{1}=QCliOverlap;
    if RefFlag
        %clear overlap of reference cliques
        CluOverlap{1}(1:RefQCliquePos,:)=[];
        CluOverlap{1}(:,1:RefQCliquePos)=[];
        %clear reference cliques in Qclique
        QClique(1:RefQCliquePos)=[];
        CluNb=length(QClique);
    end
    for CluL=1:length(CluOverlap{1})
        CluOverlap{1}(CluL,CluL)=0;
    end

    
    if RefFlag
        CluPos=find(CluSize>=CLIQUE_SIZE&CluSize<=max(CluSize));
        MaxOverlap=[];
        OverlapVal=zeros(1,101);
        MaxAnti=[];
        AntiVal=zeros(1,101);
        MaxCorr=[];
        CorrVal=zeros(1,101);

        PsNbClu=0;
        for CluL=1:length(CluPos)
            CurrClu=CluPos(CluL);
            PsNbClu=PsNbClu+length(find(Clu==CurrClu));
            if CurrClu>0
                Val=CliOverlap(CurrClu,:);
                Val(CurrClu)=0;
                [MaxVal,MaxPos]=max(Val);
                MaxOverlap=[MaxOverlap;MaxVal];
                Val(MaxPos)=0;
                OverlapVal=OverlapVal+histc(Val,[0:100]);

                Val=CliAntiCorr(CurrClu,CurrClu:end);
                Val(CurrClu)=0;
                [MaxVal,MaxPos]=max(Val);
                MaxCorr=[MaxCorr;MaxVal];
                Val(MaxPos)=0;
                CorrVal=CorrVal+histc(Val,[0:100]);

                Val=CliAntiCorr(CurrClu:end,CurrClu);
                Val=Val';
                Val(CurrClu)=0;
                [MaxVal,MaxPos]=max(Val);
                MaxAnti=[MaxAnti;MaxVal];
                Val(MaxPos)=0;
                AntiVal=AntiVal+histc(Val,[0:100]);

            end
        end
        MaxPos=find(MaxOverlap>0);
        length(find(MaxOverlap(MaxPos)>=30))*100/length(MaxPos)
        length(find(MaxOverlap(MaxPos)>=40))*100/length(MaxPos)
        length(find(MaxOverlap(MaxPos)>=50))*100/length(MaxPos)

        h=figure;
        set(gcf,'color',[1,1,1])
        set(h,'name',sprintf('overlap for cluster of sizes >==%u and <=%u (m%un%u)',CLIQUE_SIZE,max(CluSize),ModelRank,NetRank))
        hold on
        BinNb=histc(MaxOverlap,[0:100]);
        plot([0:100],cumsum(BinNb)/sum(BinNb),'g.')
        plot([0:100],cumsum(BinNb)/sum(BinNb),'g')
        plot([0:100],OverlapVal/sum(OverlapVal),'r.')
        plot([0:100],OverlapVal/sum(OverlapVal),'r')
        set(gca,'box','on')
        xlabel('overlap value')
        ylabel('cumulative frequency (max(green) and others(red))')
        set(gca,'ylim',[0,1])

        h=figure;
        set(gcf,'color',[1,1,1])
        set(h,'name',sprintf('corr for cluster of sizes >==%u and <=%u (m%un%u)',CLIQUE_SIZE,max(CluSize),ModelRank,NetRank))
        hold on
        BinNb=histc(MaxCorr,[0:100]);
        plot([0:100],cumsum(BinNb)/sum(BinNb),'g.')
        plot([0:100],cumsum(BinNb)/sum(BinNb),'g')
        plot([0:100],CorrVal/sum(CorrVal),'r.')
        plot([0:100],CorrVal/sum(CorrVal),'r')
        set(gca,'box','on')
        xlabel('corr value')
        ylabel('cumulative frequency (max(green) and others(red))')
        set(gca,'ylim',[0,1])
        h=figure;
        set(gcf,'color',[1,1,1])
        set(h,'name',sprintf('anti for cluster of sizes >==%u and <=%u (m%un%u)',CLIQUE_SIZE,max(CluSize),ModelRank,NetRank))
        hold on
        BinNb=histc(MaxAnti,[0:100]);
        plot([0:100],cumsum(BinNb)/sum(BinNb),'g.')
        plot([0:100],cumsum(BinNb)/sum(BinNb),'g')
        plot([0:100],AntiVal/sum(AntiVal),'r.')
        plot([0:100],AntiVal/sum(AntiVal),'r')
        set(gca,'box','on')
        xlabel('anti value')
        ylabel('cumulative frequency (max(green) and others(red))')
        set(gca,'ylim',[0,1])
    end



    
    
    %transform absolue rank of cliques into index
    %recover list of cliques
    CliList=[];
    if RefFlag
        RefCliList=[];
        for CluL=1:CluNb
            %start numberinf of cliques at 1
            CliList=[CliList;QClique{CluL}-MaxRefCli];                        
        end

    else
        for CluL=1:CluNb
            CliList=[CliList;QClique{CluL}];
        end
    end
    CliList=sort(CliList);
    CliRanks=CliList;

    %reorder CliAntiCorr to match current quasi-cliques order
    MemCliOverlap=CliOverlap;
    MemCliAntiCorr=CliAntiCorr;
    if RefFlag
        CliAntiCorr=CliAntiCorr(CliList+MaxRefCli,CliList+MaxRefCli);        
        CliOverlap=CliOverlap(CliList+MaxRefCli,CliList+MaxRefCli);        
    else
        CliAntiCorr=CluAntiCorr(CliList,CliList);
        CliOverlap=CliOverlap(CliList,CliList);
    end

    %add index
    CliList=[CliList,[1:length(CliList)]'];
    CliIndex=zeros(max(CliList(:,1)),1);
    CliIndex(CliList(:,1))=CliList(:,2);

    %replace rank by indexes in Clusters and Cliques
    % composition of clusters indexed on clusters of the previous rank (apart
    % for the first round where it is indexed on quasi-cliques
    Cliques{1}=cell(CluNb,1);
    if RefFlag
        for CluL=1:CluNb
            Cliques{1}{CluL,1}=CliIndex(QClique{CluL}-MaxRefCli);
        end

    else
        for CluL=1:CluNb
            Cliques{1}{CluL,1}=CliIndex(QClique{CluL});
        end
    end

    % compositions of clusters indexed on cliques

    %calculate  a matrix of anti.corr for  quasi-cliques
    %split anti and corr in separate matrix
    CliAnti=tril(CliAntiCorr)+tril(CliAntiCorr)';
    CliCorr=triu(CliAntiCorr)+triu(CliAntiCorr)';
    CluAnti{1}=zeros(length(CluNb));
    CluCorr{1}=zeros(length(CluNb));
    for CluL1=1:CluNb-1
        for CluL2=CluL1+1:CluNb
            CluAnti{1}(CluL1,CluL2)=mean(mean(CliAnti(Cliques{1}{CluL1},Cliques{1}{CluL2})));
            CluAnti{1}(CluL2,CluL1)=CluAnti{1}(CluL1,CluL2);
            CluCorr{1}(CluL1,CluL2)=mean(mean(CliCorr(Cliques{1}{CluL1},Cliques{1}{CluL2})));
            CluCorr{1}(CluL2,CluL1)=CluCorr{1}(CluL1,CluL2);
        end
    end
  

    %construct reference clique and anticorr anti/corr matrix
    %used to display the effect of rearrangement at each round
    %CluAntiCorr is already sorted
    RefOverlap=CliOverlap;
    RefAntiCorr=CliAntiCorr;

%     CurrOverlap=CluOverlap{1};
%     CurrOverlap(find(CurrOverlap==0))=[];
%     CurrOverlap=sort(CurrOverlap);
%     OVERLAP_LIMIT=CurrOverlap(round(length(CurrOverlap)*0.95))
%     CurrCorr=CluCorr{1};
%     CurrCorr(find(CurrCorr==0))=[];
%     CurrCorr=sort(CurrCorr);
%     CORR_LIMIT=CurrCorr(round(length(CurrCorr)*0.95))
%     CurrAnti=CluAnti{1};
%     CurrAnti(find(CurrAnti==0))=[];
%     CurrAnti=sort(CurrAnti);
%     ANTI_LIMIT=CurrAnti(round(length(CurrAnti)*0.95))
%     FACTOR=CORR_LIMIT/ANTI_LIMIT
    
    MaxOverlap=[];
    MaxAnti=[];
    MaxCorr=[];
    Percentile=0.5;
    AntiPercentile=0.7;
    Delta=0.05;
    
    
     for CluL=1:length(CluOverlap{1})
            Val=CluOverlap{1}(CluL,:);
            if RefFlag
                Val(CurrClu)=0;
            end
            [MaxVal,MaxPos]=max(Val);
            MaxOverlap=[MaxOverlap;MaxVal];
            
            Val=CluCorr{1}(CluL,:);
            if RefFlag
                Val(CurrClu)=0;
            end
            [MaxVal,MaxPos]=max(Val);
            MaxCorr=[MaxCorr;MaxVal];            
            
            Val=CluAnti{1}(CluL,:);
            if RefFlag
                Val(CurrClu)=0;
            end
            %[MaxVal,MaxPos]=max(Val);
            %MaxAnti=[MaxAnti;MaxVal];
            MaxAnti=[MaxAnti;Val(MaxPos)];
     end
     BinNb=histc(MaxOverlap,[0:100]);
     BinNb=cumsum(BinNb)/sum(BinNb);
     OVERLAP_LIMIT=find(BinNb>=Percentile);
     OVERLAP_LIMIT=OVERLAP_LIMIT(1);
     
     BinNb=histc(MaxCorr,[0:100]);
     BinNb=cumsum(BinNb)/sum(BinNb);
     CORR_LIMIT=find(BinNb>=Percentile);
     CORR_LIMIT=CORR_LIMIT(1);
     
     BinNb=histc(MaxAnti,[0:100]);
     BinNb=cumsum(BinNb)/sum(BinNb);
     ANTI_LIMIT=find(BinNb>=AntiPercentile);
     ANTI_LIMIT=ANTI_LIMIT(1);
     
     FACTOR=CORR_LIMIT/ANTI_LIMIT;
     
     Limits=[1,length(CluOverlap{1}),Percentile,OVERLAP_LIMIT,CORR_LIMIT,FACTOR]
     


    Continue1=1;
    Round=2;
    LastStep1=0;
    LastStep2=0;
   
    while Continue1
        if Round==2            
%             if RefFlag
%                 OVERLAP_LIMIT=90;
%                 CORR_LIMIT=50;
%                 FACTOR=1;
%                 DELTA=5;
%             else
%                 OVERLAP_LIMIT=10;
%                 CORR_LIMIT=10;
%                 FACTOR=2;
%                 DELTA=1;
%             end
%         elseif OVERLAP_LIMIT>=5
%             OVERLAP_LIMIT=5;
%             CORR_LIMIT=5;
%             FACTOR=2;
        end
        CluNb=length(CluOverlap{Round-1});

        Clusters{Round}={};
        Cliques{Round}={};
        PrevCluIndex=[];
        AddedCliques=[];
        Continue2=1;
        CluBindex=ones(CluNb,1);
        sprintf('Round %u\n',Round)
        while Continue2
            if LastStep1
                Continue2=0;
                SelectFlag=1;
                while SelectFlag
                    CluBindex(PrevCluIndex)=0;
                    CluIndex=inputdlg('give the list of agregated cliques','');
                    if isempty(CluIndex)
                        SelectFlag=0;
                    else
                        CluIndex=str2num(CluIndex{1});
                        Clusters{Round}{end+1,1}=CluIndex(1);
                        PrevCluIndex=[PrevCluIndex;CluIndex(1)];
                        CluPos=length(Clusters{Round});
                        Cliques{Round}{CluPos,1}=Cliques{Round-1}{CluIndex(1)};
                        if length(CluIndex)>1
                            for PosL2=2:length(CluIndex)
                                Clusters{Round}{CluPos,1}=[Clusters{Round}{CluPos,1};CluIndex(PosL2)];
                                Cliques{Round}{CluPos,1}=[Cliques{Round}{CluPos,1};Cliques{Round-1}{CluIndex(PosL2)}];
                                PrevCluIndex=[PrevCluIndex;CluIndex(PosL2)];
                            end
                        end
                        AddedCliques=[AddedCliques;Cliques{Round}{CluPos}];
                        if length(AddedCliques)~=length(unique(AddedCliques))
                                                                'stop'
%                             h=errordlg('cliques doublon');
%                             waitfor(h)
%                             error('process canceled')
                        end
                    end

                end
            else
                %FIRST STEP FIND QUASI-CLIQUES THAT ARE USED TO SEED CLUSTERS
                %those that have no significative interaction with precedent
                %not processed quasi-cliques
                CluBindex(PrevCluIndex)=0;
                CluIndex=find(CluBindex);

                %The first quasi-cliques belongs to seeds by definition
                if length(CluIndex)==0
                    if isempty(Clusters{Round})
                        Continue1=0;
                    end
                    Continue2=0;
                else
                    SeedList=[CluIndex(1)];
                    if length(CluIndex)>1
                        for CluL=2:length(CluIndex)
                            try
                            CurrCorr=CluCorr{Round-1}(CluIndex(CluL),CluIndex(1:CluL-1));
                            catch
                                'stop'
                            end
                            CurrAnti=CluAnti{Round-1}(CluIndex(CluL),CluIndex(1:CluL-1));
                            CurrOverlap=CluOverlap{Round-1}(CluIndex(CluL),CluIndex(1:CluL-1));
                            if length(find(CurrOverlap>=OVERLAP_LIMIT&CurrCorr>FACTOR*CurrAnti&CurrCorr>=CORR_LIMIT))==0
                                SeedList=[SeedList;CluIndex(CluL)];
                            end
                        end
                        length(SeedList)
                    end
%%%%%%%                  
                    if length(SeedList)==0
                        Continue2=0;
                    else
                        CurrOverlap=CluOverlap{Round-1};
                        %do not take account of previously processed clusters
                        CurrOverlap(PrevCluIndex,:)=0;
                        CurrAnti=CluAnti{Round-1};
                        CurrAnti(PrevCluIndex,:)=0;
                        CurrAnti(:,PrevCluIndex)=0;
                        CurrCorr=CluCorr{Round-1};
                        CurrCorr(PrevCluIndex,:)=0;
                        CurrCorr(:,PrevCluIndex)=0;                        

                        %take upper triangular matrix not to have
                        %max(Overlap)< current position
                        Overlap=triu(CurrOverlap,1);

                        %calculate FACTOR,OVERLAP_LIMIT a CORR_LIMIT

                        %¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹
                        % GET RID OF NON SIGNIFICATIVE OVERLAP OR CORR
                        Overlap(Overlap<OVERLAP_LIMIT|CurrCorr<=FACTOR*CurrAnti|CurrCorr<CORR_LIMIT)=0;
                        %¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹

                        %recover clusters that have no overlap with others
                        % => the cluster which has the aximal overlap with them are themselve
                        [MaxVal,MaxPos]=max(Overlap);
                        NullPos=find(MaxVal==0);
                        MaxPos(NullPos)=NullPos;
                        for CluL=1:length(SeedList)
                            if isempty(find(PrevCluIndex==SeedList(CluL)))
                                Pos=find(MaxPos==SeedList(CluL));
                                if ~isempty(Pos)
                                    %add the seed cluster
                                    Pos=unique([SeedList(CluL),Pos]);
                                    %the first cluster (the seed
                                    %cluster) is automatically used
                                    Clusters{Round}{end+1,1}=Pos(1);
                                    PrevCluIndex=[PrevCluIndex;Pos(1)];
                                    CluPos=length(Clusters{Round});
                                    Cliques{Round}{CluPos,1}=Cliques{Round-1}{Pos(1)};
                                    if length(Pos)>1
                                        KeptClu=Pos(1);
                                        for PosL2=2:length(Pos)
                                            KeepIt=1;
                                            for PosL1=1:PosL2-1
                                                if ~isempty(find(KeptClu==Pos(PosL2)))
                                                    if Overlap(Pos(PosL1),Pos(PosL2))<OVERLAP_LIMIT|CurrCorr(Pos(PosL1),Pos(PosL2))<=FACTOR*CurrAnti(Pos(PosL1),Pos(PosL2))|CurrCorr(Pos(PosL1),Pos(PosL2))<CORR_LIMIT
                                                        KeepIt=0;
                                                        break
                                                    end
                                                end
                                            end
                                            if KeepIt
                                                KeptClu=[KeptClu;Pos(PosL2)];
                                                Clusters{Round}{CluPos,1}=[Clusters{Round}{CluPos,1};Pos(PosL2)];
                                                Cliques{Round}{CluPos,1}=[Cliques{Round}{CluPos,1};Cliques{Round-1}{Pos(PosL2)}];
                                                PrevCluIndex=[PrevCluIndex;Pos(PosL2)];
                                            end
                                        end
                                    end
                                    AddedCliques=[AddedCliques;Cliques{Round}{CluPos}];
                                    if length(AddedCliques)~=length(unique(AddedCliques))
                                        'stop'
                                        %                                     h=errordlg('cliques doublon');
                                        %                                     waitfor(h)
                                        %                                     error('process canceled')
                                    end
                                else
                                    PrevCluIndex=[PrevCluIndex;SeedList(CluL)];
                                    Clusters{Round}{end+1,1}=SeedList(CluL);
                                    Cliques{Round}{end+1,1}=Cliques{Round-1}{SeedList(CluL)};
                                    AddedCliques=[AddedCliques;Cliques{Round-1}{end}];
                                    if length(AddedCliques)~=length(unique(AddedCliques))
                                        'stop'
                                        %                                     h=errordlg('cliques doublon');
                                        %                                     waitfor(h)
                                        %                                     error('process canceled')
                                    end
                                end
                                if length(AddedCliques)~=length(unique(AddedCliques))
                                    'stop'
                                    %                                 h=errordlg('cliques doublon');
                                    %                                 waitfor(h)
                                    %                                 error('process canceled')
                                end
                            end
                            % end
                        end %if length(SeedList==0)
                    end %CluIndex>1
                end %CluIndex>0
            end %if LastStep
        end % while Continue2

        %complete Clusters with not used clusters
        %CluBindex=ones(CluNb,1);
        CluBindex(PrevCluIndex)=0;
        CluIndex=find(CluBindex);

        if ~isempty(CluIndex)
            for CluL=1:length(CluIndex)
                Clusters{Round}{end+1}=CluIndex(CluL);
                Cliques{Round}{end+1}=Cliques{Round-1}{CluIndex(CluL)};
            end
        end
        ChangeRound=1;
        if length(Clusters{Round})==length(Clusters{Round-1})
            %keep same round but change parameters
            if LastStep1==0
                'decrement limits'
                Percentile=Percentile-Delta;
                AntiPercentile=AntiPercentile-Delta;
                if Percentile<=0
                    LastStep1=1;
                end
                ChangeRound=0;
            else
                Continue1=0;
            end
        end
        if ChangeRound

            %calculate new overlap
            CluOverlap{Round}=zeros(length(Clusters{Round}));
            for CluL1=1:length(Clusters{Round})-1
                for CluL2=CluL1+1:length(Clusters{Round})
                    CluOverlap{Round}(CluL1,CluL2)=mean(mean(CluOverlap{Round-1}(Clusters{Round}{CluL1},Clusters{Round}{CluL2})));
                    CluOverlap{Round}(CluL2,CluL1)=mean(mean(CluOverlap{Round-1}(Clusters{Round}{CluL1},Clusters{Round}{CluL2})));
                end
            end
            %calculate new anti corr
            CluAnti{Round}=zeros(length(Clusters{Round}));
            CluCorr{Round}=zeros(length(Clusters{Round}));
            for CluL1=1:length(Clusters{Round})-1
                for CluL2=CluL1+1:length(Clusters{Round})
                    CluAnti{Round}(CluL1,CluL2)=mean(mean(CluAnti{Round-1}(Clusters{Round}{CluL1},Clusters{Round}{CluL2})));
                    CluAnti{Round}(CluL2,CluL1)=CluAnti{Round}(CluL1,CluL2);
                    CluCorr{Round}(CluL1,CluL2)=mean(mean(CluCorr{Round-1}(Clusters{Round}{CluL1},Clusters{Round}{CluL2})));
                    CluCorr{Round}(CluL2,CluL1)=CluCorr{Round}(CluL1,CluL2);
                end
            end
        else
            CluOverlap{Round}=CluOverlap{Round-1};
            CluCorr{Round}=CluCorr{Round-1};
            CluAnti{Round}=CluAnti{Round-1};
        end
        %recalculate limits
        Round
        %             CurrOverlap=CluOverlap{Round};
        %             CurrOverlap(find(CurrOverlap==0))=[];
        %             CurrOverlap=sort(CurrOverlap);
        %             OVERLAP_LIMIT=CurrOverlap(round(length(CurrOverlap)*0.95))
        %             CurrCorr=CluCorr{Round};
        %             CurrCorr(find(CurrCorr==0))=[];
        %             CurrCorr=sort(CurrCorr);
        %             CORR_LIMIT=CurrCorr(round(length(CurrCorr)*0.95))
        %             CurrAnti=CluAnti{Round};
        %             CurrAnti(find(CurrAnti==0))=[];
        %             CurrAnti=sort(CurrAnti);
        %             ANTI_LIMIT=CurrAnti(round(length(CurrAnti)*0.95))
        %             FACTOR=CORR_LIMIT/ANTI_LIMIT



        MaxOverlap=[];
        MaxAnti=[];
        MaxCorr=[];

        for CluL=1:length(CluOverlap{Round})
            Val=CluOverlap{Round}(CluL,:);
            if RefFlag
                Val(CurrClu)=0;
            end
            [MaxVal,MaxPos]=max(Val);
            MaxOverlap=[MaxOverlap;MaxVal];

            Val=CluCorr{Round}(CluL,:);
            if RefFlag
                Val(CurrClu)=0;
            end
            [MaxVal,MaxPos]=max(Val);
            MaxCorr=[MaxCorr;MaxVal];

            Val=CluAnti{Round}(CluL,:);
            if RefFlag
                Val(CurrClu)=0;
            end
            %[MaxVal,MaxPos]=max(Val);
            %MaxAnti=[MaxAnti;MaxVal];
            MaxAnti=[MaxAnti;Val(MaxPos)];
        end
        BinNb=histc(MaxOverlap,[0:100]);
        BinNb=cumsum(BinNb)/sum(BinNb);
        OVERLAP_LIMIT=find(BinNb>=Percentile);
        OVERLAP_LIMIT=OVERLAP_LIMIT(1);

        BinNb=histc(MaxCorr,[0:100]);
        BinNb=cumsum(BinNb)/sum(BinNb);
        CORR_LIMIT=find(BinNb>=Percentile);
        CORR_LIMIT=CORR_LIMIT(1);

        BinNb=histc(MaxAnti,[0:100]);
        BinNb=cumsum(BinNb)/sum(BinNb);
        ANTI_LIMIT=find(BinNb>=AntiPercentile);
        ANTI_LIMIT=ANTI_LIMIT(1);

        FACTOR=CORR_LIMIT/ANTI_LIMIT;
        Limits=[Limits;[Round,length(Clusters{Round}),Percentile,OVERLAP_LIMIT,CORR_LIMIT,FACTOR]]

        if ChangeRound
            if Round==54

                CurrCliList=[];
                CliTick=[];
                for CliL=1:length(Cliques{Round})
                    CurrCliList=[CurrCliList;Cliques{Round}{CliL}];
                    if CliL==1
                        CliTick=length(Cliques{Round}{CliL});
                    else
                        CliTick(CliL)=CliTick(CliL-1)+length(Cliques{Round}{CliL});
                    end
                end
                Overlap=RefOverlap(CurrCliList,CurrCliList);
                AntiCorr=RefAntiCorr(CurrCliList,CurrCliList);

                h1=figure;
                set(gcf,'color',[1,1,1])
                %subplot(2,2,1)
                set(h1,'name',sprintf('CLIQUES OVERLAP AT END OF ROUND %u',Round));
                image(Overlap)
                set(gca,'tickdir','out')
                set(gca,'xtick',CliTick)
                set(gca,'ytick',CliTick)
                set(gca,'xticklabel','')
                set(gca,'yticklabel','')
                title(sprintf('cliques overlap at end of round %u for m%un%u',Round,ModelRank,NetRank));
                xlabel('agregated quasi-cliques')
                ylabel('agregated quasi-cliques')
                %subplot(2,2,2)
                h2=figure;
                set(gcf,'color',[1,1,1])
                set(h2,'name',sprintf('CLIQUES ANTI/CORR AT END OF ROUND %u',Round));
                image(AntiCorr*3)
                set(gca,'tickdir','out')
                set(gca,'xtick',CliTick)
                set(gca,'ytick',CliTick)
                set(gca,'xticklabel','')
                set(gca,'yticklabel','')
                title(sprintf('cliques CORR at end of round %u for m%un%u',Round,ModelRank,NetRank));
                xlabel(sprintf('cliques ANTI at end of round %u for m%un%u',Round,ModelRank,NetRank))
                ylabel('agregated quasi-cliques')

                set(gca,'tickdir','out')
                set(gca,'xtick',CliTick)
                set(gca,'ytick',CliTick)
                set(gca,'xticklabel','')
                set(gca,'yticklabel','')
                title(sprintf('cliques anti/corr end of round %u',Round));
            end

            Round=Round+1;
        end

        %net_properties('clusterize quasi-cliques',8,1,54)

    end %while Continue1
end %if DoIt
%write real ranks of cliques
for RoundL=1:length(Cliques)
    for CliL=1:length(Cliques{RoundL})
        Cliques{RoundL}{CliL}=CliRanks(Cliques{RoundL}{CliL});
    end
end
cd(NetDir)
if RefFlag
    %complete with references cliques
    CurrCliques=Cliques{end};
    Cliques=RefCliques{end}';
    for CliL=1:length(CurrCliques)
        Cliques{end+1,1}=CliList(CurrCliques{CliL},1)+MaxRefCli;
    end
    
    eval(sprintf('save m%un%u_quasi-cliques_clusters_size%u_%02u Cliques',ModelRank,NetRank,CLIQUE_SIZE,CorrLimit))
else
eval(sprintf('save m%un%u_quasi-cliques_clusters_size%u_%02u CluOverlap Cliques Clusters',ModelRank,NetRank,CLIQUE_SIZE,CorrLimit))
end


if DISPLAY_IT
    
    CliList=[];
    CliTick=[];
    for CliL=1:length(Cliques)
        CliList=[CliList;Cliques{CliL}];
        if CliL==1
            CliTick=length(Cliques{CliL});
        else
            CliTick(CliL)=CliTick(CliL-1)+length(Cliques{CliL});
        end
    end
%     
%     for CliL=1:length(Cliques)
%         CliList=[CliList;Cliques{CliL}];
%         CliTick(CliL+RefRound)=CliTick(CliL+RefRound-1)+length(Cliques{Round}{CliL});
%     end
%     
    Overlap=MemCliOverlap(CliList,CliList);
    AntiCorr=MemCliAntiCorr(CliList,CliList);

    h1=figure;
    set(gcf,'color',[1,1,1])
    %subplot(2,2,1)
    set(h1,'name','CLIQUES OVERLAP AT END OF CLUSTERING');
    image(Overlap)
    set(gca,'tickdir','out')
    set(gca,'xtick',CliTick)
    set(gca,'ytick',CliTick)
    set(gca,'xticklabel','')
    set(gca,'yticklabel','')
    title('cliques overlap at end of clustering');
    xlabel('agregated quasi-cliques')
    ylabel('agregated quasi-cliques')
    %subplot(2,2,2)
    h2=figure;
    set(gcf,'color',[1,1,1])
    set(h2,'name','CLIQUES ANTI/CORR AT END OF CLUSTERING');
    image(AntiCorr*3)
    set(gca,'tickdir','out')
    set(gca,'xtick',CliTick)
    set(gca,'ytick',CliTick)
    set(gca,'xticklabel','')
    set(gca,'yticklabel','')
    title('cliques CORR at end clustering');
    xlabel('cliques ANTI at end of clustering')
    ylabel('agregated quasi-cliques')

    set(gca,'tickdir','out')
    set(gca,'xtick',CliTick)
    set(gca,'ytick',CliTick)
    set(gca,'xticklabel','')
    set(gca,'yticklabel','')
    title('cliques anti/corr at end of clustering');


    MODEL_RANK=11;
    NET_RANK=40;
    cd(fullfile(K.dir.net,sprintf('m%03u',MODEL_RANK),sprintf('n%05u',NET_RANK)))
    Round=length(Cliques);
    CluL1=1;
    CluL2=1;
        CurrCliques=Cliques{Round}{CluL1};
        %CurrCliques=Cliques{CluL1};
        PsRanks1=[];
        for CliL=1:length(CurrCliques)
            PsRanks1=[PsRanks1;find(Clu==CurrCliques(CliL))];
        end
        length(PsRanks1)
        CurrCliques=Cliques{Round}{CluL2};
        PsRanks2=[];
        for CliL=1:length(CurrCliques)
            PsRanks2=[PsRanks2;find(Clu==CurrCliques(CliL))];
        end
        length(PsRanks2)

        if length(PsRanks1)>=50
            C=load_data(sprintf('c_m%u_n%u.4mat',MODEL_RANK,NET_RANK),'./',PsNb,PsNb,'uint8','ieee-le',PsRanks1,PsRanks2);
            if CluL1==CluL2
                C=triu(C);
            end;
            A=load_data(sprintf('a_m%u_n%u.4mat',MODEL_RANK,NET_RANK),'./',PsNb,PsNb,'uint8','ieee-le',PsRanks1,PsRanks2);
            if CluL1==CluL2
                C=C+tril(A,-1);
                clear A
            end
            h=figure;
            set(h,'name',sprintf('CLU %u vs CLu %u',CluL1,CluL2));
            if CluL1==CluL2
                h=pcolor(flipud(double(C)));
                set(h,'linestyle','none')
                title(sprintf('m%un%u: CLU %u vs CLU %u',MODEL_RANK,NET_RANK,CluL1,CluL2))
            else
                subplot(1,2,1)
                h=pcolor(flipud(double(C)));
                set(h,'linestyle','none')
                title(sprintf('CORR m%un%u clu %u vs clu %u',MODEL_RANK,NET_RANK,CluL1,CluL2))
                subplot(1,2,2)
                h=pcolor(flipud(double(A)));
                title(sprintf('ANTI m%un%u clu %u vs clu %u',MODEL_RANK,NET_RANK,CluL1,CluL2))
                set(h,'linestyle','none')
            end
            
                
            set(gcf,'color',[1,1,1])
            
        end
    


    %compare overlap,anti corr

end

%% CLUSTERIZE QUASI-CLIQUES FIRST VERSION
function CLUSTERIZE_QUASI_CLIQUES_V1(ModelRank,NetRank,PsNb)
%net_properties('clusterize quasi-cliques',8,1,54)
global K
DISPLAY_IT=0;
%limit of significative overlap
OVERLAP_LIMIT=30;
CORR_LIMIT=0;
FACTOR=0;
CLIQUE_SIZE=5;
NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',NetRank));
cd(NetDir)
DoIt=1;
if exist(sprintf('m%un%u_quasi-clique_clusters_size%u.mat',ModelRank,NetRank,CLIQUE_SIZE),'file')
    load(sprintf('m%un%u_quasi-clique_clusters_size%u.mat',ModelRank,NetRank,CLIQUE_SIZE))
    DoIt=0;
else
    if exist(sprintf('m%un%u_quasi-cliques_size%u.mat',ModelRank,NetRank,CLIQUE_SIZE),'file')
        load(sprintf('m%un%u_quasi-cliques_size%u.mat',ModelRank,NetRank,CLIQUE_SIZE))
    else
        h=errordlg('Do QUASI CLIQUES');
        waitfor(h)
        error('process canceled')
    end
    if exist(sprintf('m%un%u_cliques_anticorr.mat',ModelRank,NetRank),'file')
        load(sprintf('m%un%u_cliques_anticorr.mat',ModelRank,NetRank))
    else
        h=errordlg('Do CALCULATE ANTI CORR');
        waitfor(h)
        error('process canceled')
    end
    if exist(sprintf('m%un%u_cliques_index.mat',ModelRank,NetRank),'file')
        load(sprintf('m%un%u_cliques_index.mat',ModelRank,NetRank))
    else
        h=errordlg('Save Clu in cliques_index');
        waitfor(h)
        error('process canceled')
    end

end

if DoIt
    %transform absolue rank of cliques into index
    %recover list of cliques
    CliList=[];
    for CliL=1:length(QClique)
        CliList=[CliList;QClique{CliL}];
    end
    CliList=sort(CliList);

    %reorder CluAntiCorr to match current cliques order
    CluAntiCorr=CluAntiCorr(CliList,CliList);
    %add index
    CliList=[CliList,[1:length(CliList)]'];
    CliIndex=zeros(max(CliList(:,1)),1);
    CliIndex(CliList(:,1))=CliList(:,2);

    % The merging of cliques is considered as the first fround or reordering
    % overlap between cluster at each round
    for CluL=1:length(QCliOverlap)
        QCliOverlap(CluL,CluL)=0;
    end
    CluOverlap{1}=QCliOverlap;

    %replace rank by indexes in Clusters and CLiques
    % composition of clusters indexed on clusters of the previous rank (apart
    % for the first round where it is indexed on cliques=
    Clusters{1}={};
    for CluL=1:length(QClique)
        Clusters{1}{CluL,1}=CliIndex(QClique{CluL});
    end
    % compositions of clusters indexed on cliques
    Cliques{1}=Clusters{1};

    %calculate  a matrix of anti.corr for  quasi-cliques
    %split anti and corr in separate matrix
    CluAnti=tril(CluAntiCorr)+tril(CluAntiCorr)';
    CluCorr=triu(CluAntiCorr)+triu(CluAntiCorr)';
    QCluAnti{1}=zeros(length(length(QClique)));
    QCluCorr{1}=zeros(length(length(QClique)));
    for CluL1=1:length(QClique)-1
        for CluL2=CluL1+1:length(QClique)
            QCluAnti{1}(CluL1,CluL2)=mean(mean(CluAnti(Clusters{1}{CluL1},Clusters{1}{CluL2})));
            QCluAnti{1}(CluL2,CluL1)=QCluAnti{1}(CluL1,CluL2);
            QCluCorr{1}(CluL1,CluL2)=mean(mean(CluCorr(Clusters{1}{CluL1},Clusters{1}{CluL2})));
            QCluCorr{1}(CluL2,CluL1)=QCluCorr{1}(CluL1,CluL2);
        end
    end
    SelfCorr=zeros(length(length(QClique)),1);
    SelfAnti=zeros(length(length(QClique)),1);
    for CluL=1:length(QClique)-1
        SelfCorr(CluL)=mean(mean(CluCorr(Clusters{1}{CluL},Clusters{1}{CluL})));
        SelfAnti(CluL)=mean(mean(CluAnti(Clusters{1}{CluL},Clusters{1}{CluL})));
    end

    %construct a reference clique overlap matrix
    %used to display the effect of rearrangement at each round
    CliList=[];
    for CliL=1:length(Cliques{1})
        CliList=[CliList;Cliques{1}{CliL}];
    end
    [Tmp,SortOrder]=sort(CliList);
    RefOverlap=CliOverlap(SortOrder,SortOrder);
    %construct a reference clique anti/corr matrix
    %used to display the effect of rearrangement at each round
    RefAntiCorr=CluAntiCorr(SortOrder,SortOrder);

    if DISPLAY_IT
        h1=figure;
        set(h1,'name',sprintf('M%uN%u - CLUSTERIZING QUASI-CLIQUES',ModelRank,NetRank))
        title(sprintf('M%uN%u - CLUSTERIZING QUASI-CLIQUES',ModelRank,NetRank))
        set(gcf,'color',[1,1,1])
    end
    Continue=1;
    RoundL=1;
    while Continue
        RoundL=RoundL+1;
        CurrOverlap=CluOverlap{RoundL-1};
        CurrAnti=QCluAnti{RoundL-1};
        CurrCorr=QCluCorr{RoundL-1};
        if RoundL==2
            %find cluster that have no many overlap with other and
            %deplace them towards the end

            MeanOverlap=zeros(length(CurrOverlap),1);
            StdOverlap=zeros(length(CurrOverlap),1);
            for CluL=1:length(CurrOverlap)
                Pos=find(CurrOverlap(CluL,:));
                if ~isempty(Pos)
                    MeanOverlap(CluL)=mean(CurrOverlap(CluL,Pos));
                    StdOverlap(CluL)=std(CurrOverlap(CluL,Pos));
                end
            end
            KeepIndex=find(MeanOverlap>=mean(MeanOverlap));
            MoveIndex=find(MeanOverlap<mean(MeanOverlap));
            NewIndex=[KeepIndex;MoveIndex]';
            CurrOverlap=CurrOverlap(NewIndex,NewIndex);
            CurrAnti=CurrAnti(NewIndex,NewIndex);
            CurrCorr=CurrCorr(NewIndex,NewIndex);
        else
            NewIndex=1:length(CurrOverlap);
        end
        if DISPLAY_IT
            h2=figure;
            set(h2,'name',sprintf('M%uN%u - BEFORE REARRANGEMENT ROUND µu',ModelRank,NetRank,RoundL))
            set(gcf,'color',[1,1,1])
            set(h2,'name',sprintf('clusters rearrangement (round %u)',RoundL))
            figure(h2)
            subplot(1,2,1)
            image(CurrOverlap)
            title(sprintf('clusters before rearrangement (round %u)',RoundL))
        end

        Overlap=triu(CurrOverlap,1);

        %¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹
        % GET RID OF NON SIGNIFICATIVE OVERLAP OR CORR
        Overlap(Overlap<OVERLAP_LIMIT|triu(CurrCorr)<FACTOR*tril(CurrAnti)'|triu(CurrCorr)<CORR_LIMIT)=0;
        %¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹

        [MaxVal,MaxPos]=max(Overlap);
        NullPos=find(MaxVal==0);
        MaxPos(NullPos)=NullPos;


        CurrClusters={};
        CurrCliques={};
        CluOrder=[];
        for CluL=1:length(CurrOverlap)
            Pos=find(MaxPos==CluL);
            if ~isempty(Pos)
                Pos=NewIndex(Pos);
                CluOrder=[CluOrder;Pos(1)];
                CurrClusters{end+1,1}=Pos(1);
                CliPos=length(CurrClusters);
                CurrCliques{CliPos,1}=Cliques{RoundL-1}{Pos(1)};
                if length(Pos)>1
                    for PosL1=2:length(Pos)
                        KeepIt=1;
                        for PosL2=1:PosL1-1
                            if Overlap(Overlap(Pos(PosL1),Pos(PosL2))<OVERLAP_LIMIT|triu(CurrCorr(Pos(PosL1),Pos(PosL2)))<FACTOR*tril(CurrAnti(Pos(PosL1),Pos(PosL2)))'|triu(CurrCorr(Pos(PosL1),Pos(PosL2)))<CORR_LIMIT)
                                KeepIt=0;
                                break
                            end
                        end
                        if KeepIt
                            CurrClusters{CliPos,1}=[CurrClusters{CliPos,1};Pos(PosL1)];
                            CurrCliques{CliPos,1}=[CurrCliques{CliPos,1};Cliques{RoundL-1}{Pos(PosL1)}];
                        end
                    end
                end
            end
        end

        if length(CurrClusters)<length(Clusters{RoundL-1})

            Clusters{RoundL}=CurrClusters;
            Cliques{RoundL}=CurrCliques;
            CurrOverlap=CurrOverlap(CluOrder,CluOrder);
            if DISPLAY_IT
                figure(h2)
                subplot(1,2,2)
                image(CurrOverlap)
                title(sprintf('clusters after rearrangement (round %u)',RoundL))
            end

            %calculate new overlap
            CluOverlap{RoundL}=zeros(length(Cliques{RoundL}));
            for CluL1=1:length(Clusters{RoundL})-1
                for CluL2=CluL1+1:length(Cliques{RoundL})
                    CluOverlap{RoundL}(CluL1,CluL2)=mean(mean(CluOverlap{RoundL-1}(Clusters{RoundL}{CluL1},Clusters{RoundL}{CluL2})));
                    CluOverlap{RoundL}(CluL2,CluL1)=mean(mean(CluOverlap{RoundL-1}(Clusters{RoundL}{CluL1},Clusters{RoundL}{CluL2})));
                end
            end
            %calculate new anti corr
            QCluAnti{RoundL}=zeros(length(Clusters{RoundL}));
            QCluCorr{RoundL}=zeros(length(Clusters{RoundL}));
            for CluL1=1:length(Clusters{RoundL})-1
                for CluL2=CluL1+1:length(Clusters{RoundL})
                    QCluAnti{RoundL}(CluL1,CluL2)=mean(mean(QCluAnti{RoundL-1}(Clusters{RoundL}{CluL1},Clusters{RoundL}{CluL2})));
                    QCluAnti{RoundL}(CluL2,CluL1)=QCluAnti{RoundL}(CluL1,CluL2);
                    QCluCorr{RoundL}(CluL1,CluL2)=mean(mean(QCluCorr{RoundL-1}(Clusters{RoundL}{CluL1},Clusters{RoundL}{CluL2})));
                    QCluCorr{RoundL}(CluL2,CluL1)=QCluCorr{RoundL}(CluL1,CluL2);
                end
            end

            SelfCorr=zeros(length(length(QClique{RoundL})),1);
            SelfAnti=zeros(length(length(QClique{RoundL})),1);
            for CluL=1:length(Clusters{RoundL})
                SelfCorr(CluL)=mean(mean(CluCorr(Clusters{1}{CluL},Clusters{RoundL}{CluL})));
                SelfAnti(CluL)=mean(mean(CluAnti(Clusters{1}{CluL},Clusters{RoundL}{CluL})));
            end
            if DISPLAY_IT
                h=figure;
                set(h,'name',sprintf('M%uN%u - ANTI CORR BETWEEN CLUSTERS ROUND %u',ModelRank,NetRank,RoundL))
                image(triu(QCluCorr{RoundL})+tril(QCluAnti{RoundL}))
                xlabel('ANTI')
                title('CORR')


                h=figure;
                set(h,'name',sprintf('M%uN%u - OVERLAP BETWEEN CLUSTERS ROUND %u',ModelRank,NetRank,RoundL))
                image(CluOverlap{RoundL})
                set(gcf,'color',[1,1,1])
                set(gca,'xlim',[0,200])
                set(gca,'ylim',[0,200])

            end

        else
            Continue=0;
        end

        CurrCliList=[];
        for CliL=1:length(CurrCliques)
            CurrCliList=[CurrCliList;CurrCliques{CliL}];
        end
        Overlap=RefOverlap(CurrCliList,CurrCliList);
        AntiCorr=RefAntiCorr(CurrCliList,CurrCliList);
        if DISPLAY_IT
            h=figure;
            set(h,'name',sprintf('CLIQUES OVERLAP AFTER REARRANGEMENT AT ROUND %u',RoundL));
            image(Overlap)
            title(sprintf('cliques overlap after rearrangement at round %u',RoundL));
            set(gcf,'color',[1,1,1])


            h=figure;
            set(h,'name',sprintf('CLIQUES ANTI/CORR AFTER REARRANGEMENT AT ROUND %u',RoundL));
            image(AntiCorr)
            title(sprintf('cliques anti/corr after rearrangement at round %u',RoundL));
            set(gcf,'color',[1,1,1])

            figure(h1)
            if RoundL<=10
                subplot(3,3,RoundL-1)
                image(Overlap)
                title(sprintf('cliques after rearrangement at round %u',RoundL));
            end
        end
    end
    cd(NetDir)
    eval(sprintf('save m%un%u_quasi-cliques_clusters_size%u CluOverlap Cliques Clusters',ModelRank,NetRank,CLIQUE_SIZE))
end

if DISPLAY_IT
    MODEL_RANK=8;
    NET_RANK=54;
    Round=length(Cliques);
    for CluL=1:length(Cliques{Round})
        CurrCliques=Cliques{Round}{CluL};
        PsRanks=[];
        for CluL=1:length(CurrCliques)
            PsRanks=[PsRanks;find(Clu==CurrCliques(CluL))];
        end
        length(PsRanks)
        if length(PsRanks>50)

            C=load_data(sprintf('c_m%u_n%u.4mat',MODEL_RANK,NET_RANK),'./',PsNb,PsNb,'uint8','ieee-le',PsRanks,PsRanks);
            C=triu(C);
            A=load_data(sprintf('a_m%u_n%u.4mat',MODEL_RANK,NET_RANK),'./',PsNb,PsNb,'uint8','ieee-le',PsRanks,PsRanks);
            C=C+tril(A,-1);
            clear A
            h=figure;
            set(h,'name',sprintf('CLU %u',CluL));
            image(C);
            set(gcf,'color',[1,1,1])
            title(sprintf('CLU %u',CluL))
        end
    end
end

%% FUNCTION FIND_LIMITS
function Limit=FIND_LIMITS(CluSize,InclusionLimits,ExclusionLimits,SizeTest,OutFlag)
FoundLimit=0;
for TestL=1:length(SizeTest)
    if eval(sprintf('CluSize%s',SizeTest{TestL}))
        InclusionLimit=InclusionLimits(TestL);
        ExclusionLimit=ExclusionLimits(TestL);
        FoundLimit=1;
        break
    end
end
if FoundLimit==0
    InclusionLimit=InclusionLimits(end);
    ExclusionLimit=ExclusionLimits(end);
end
if isequal(OutFlag,'IN')
    Limit=InclusionLimit;
else
    Limit=ExclusionLimit;
end
%% FUNCTION CLIQUES
function [Clu,CluSize,CliqueSize,IsSeed]=CLIQUES(ModelRank,NetRank,LoadedPsRanks,CurrConn,PsNb,CorrLimit,Clu,CluSize,CliqueSize,IsSeed,CliqueRank,MyAlgoFlag,DisplayFlag,varargin)
global K
Round=0;
CliNb=0;
Continue=1;
NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',NetRank));
if nargin==13
    FoundPsRanks=find(Clu);
else
    FoundPsRanks=varargin{1};
end

while Continue
    Round=Round+1;
    if Round>1
        %if Round==1, CurrConn has been passed by parameter
        %if Round>1 we construct CurrConn, the connectivity for the currently processed
        %sub-network. ps found in the previous step and the ps which connectivity is null
        %in CurrConn are deleted in order to recalculate the connectivity of the residual
        %network.
        
        FoundPsRanks=unique([FoundPsRanks;LoadedPsRanks(find(CurrConn==0))]);        
        LoadedPsRanks=[1:PsNb]';
        LoadedPsRanks(FoundPsRanks)=[];
        if ~isempty(LoadedPsRanks)
            CurrConn=zeros(length(LoadedPsRanks),1);
            Corr=load_data(sprintf('c_m%u_n%u.4mat',ModelRank,NetRank),NetDir,PsNb,PsNb,'uint8','ieee-le',LoadedPsRanks,LoadedPsRanks);
            for PsL=1:length(LoadedPsRanks)
                %self connectivity is not counted
                Corr(PsL,LoadedPsRanks(PsL))=0;                
                CurrConn(PsL)=length(find(Corr(PsL,:)>CorrLimit));
            end
        else
            CurrConn=0;           
        end
    end    
    if isempty(find(CurrConn>0))|(Round>1&max(Clu)==CliNb)
        %end of the process: either it does not exist ps with connectivity>0
        %or no new cliques has been created
        Continue=0;
        %integrate the last probe set in individual cliques
        SinglePsRanks=setdiff(LoadedPsRanks,FoundPsRanks);
        if ~isempty(SinglePsRanks)
            CliRank=max(Clu);
            for PsL=1:length(SinglePsRanks)
                CliRank=CliRank+1;
                Clu(SinglePsRanks(PsL))=CliRank;
                CluSize(CliRank)=1;
            end
        end
    else
        %keep in memory the current connectivity which is used if DisplayFlag=1
        MemConn=CurrConn;
        %keep in memory the rank of the last clique found in the previous step
        %in order to check at the end of the current round if new cliques have been found
        CliNb=max(Clu);
        %Connectivity values present in the network currently processed
        ConnVal=setdiff(unique(CurrConn),0);
        %indicates if the connectivity class has been effectively searched or
        %ignored because all the probe sets belonging to this class belong also
        %to a class of greated size and have been found in a previous step
        ConnSearched=zeros(length(ConnVal),1);
        %number of probe set found in each connectivity class
        ConnFound=zeros(length(ConnVal),1);
        %Total nb of edge linked to probe sets belonging to a given class of connectivity
        ConnTotal=zeros(length(ConnVal),1);
        for ConnL=1:length(ConnVal);
            ConnTotal(ConnL)=ConnVal(ConnL)*length(find(CurrConn==ConnVal(ConnL)));
        end
        %The real maximal connectivity of the processed subnetwork
        StartVal=zeros(length(ConnVal),1);
        %The size of the processed subnetwork
        CorrSize=zeros(length(ConnVal),1);
        %The size of the clique found on the processed subnetwotk by my algorithm
        MyVal=zeros(length(ConnVal),1);
        %The size of the clique found on the processed subnetwotk by cliquer algorithm
        CliquerVal=zeros(length(ConnVal),1);
        %The number of tested probe sets
        TestedNb=zeros(length(ConnVal),1);
        %The number of found ps in the previous step (ConnL-1)
        MemPsNb=0;
        for ConnL=length(ConnVal):-1:1           
            if ConnL<length(ConnVal)
                %fill the information for the previous ConnL step
                ConnFound(ConnL+1)=length(FoundPsRanks)-MemPsNb;
                TestedNb(ConnL+1)=TestedPsNb;
                sprintf('ConnL=%u: connectivity %u => %u found probesets among %u in previous step',ConnL+1,ConnVal(ConnL+1),ConnFound(ConnL+1),TestedPsNb)
                %update MemPsNb to find the number of found Ps in the next step
                MemPsNb=length(FoundPsRanks);
            end
            %process all the probe set which belong to the current connectivity class
            ConnPos=find(CurrConn==ConnVal(ConnL));
            TestedPsNb=0;
            %Ps with that connectivity could have been merged in a clique in the previous step, and the class could therefore be empty
            if ~isempty(ConnPos)
                Ok=1;
                %mark the current connectivity as searched
                ConnSearched(ConnL)=1;
                %several probe set may have the same connectivity
                while Ok
                    %ConnPos=find(CurrConn>100&CurrConn<300);
                    %uses the first probe set in the current connectivity class
                    CurrPsRank=LoadedPsRanks(ConnPos(1));
                    Corr=load_data(sprintf('c_m%u_n%u.4mat',ModelRank,NetRank),NetDir,PsNb,PsNb,'uint8','ieee-le',1:PsNb,CurrPsRank);
                    Anti=load_data(sprintf('a_m%u_n%u.4mat',ModelRank,NetRank),NetDir,PsNb,PsNb,'uint8','ieee-le',1:PsNb,CurrPsRank);
                    %forces current ps to be kept in final list
                    Corr(CurrPsRank)=100;
                    %don't keep Corr that are inferior to Anti
                    Corr(Corr<=Anti)=0;                                        
                    %do not process again found probesets => found clique or cluster will
                    %have his size less or equal to the connectivity of the current
                    %probeset
                    Corr(FoundPsRanks)=0;
                    %Load Corr
                    PsRanks=find(Corr>CorrLimit);                    
                    TestedPsNb=TestedPsNb+length(PsRanks);
                    if length(PsRanks)>1
                        PsPositions=[];
                        for PsL=1:length(PsRanks)
                            PsPositions(end+1,1)=find(LoadedPsRanks==PsRanks(PsL));
                        end
                        if length(PsRanks)==2
                            StartVal(ConnL,1)=1;
                            CorrSize(ConnL,1)=2;
                            Conn=1;
                        else
                            %load the sub network constituted of PsRanks
                            Corr=load_data(sprintf('c_m%u_n%u.4mat',ModelRank,NetRank),NetDir,PsNb,PsNb,'uint8','ieee-le',PsRanks,PsRanks);
                            Anti=load_data(sprintf('a_m%u_n%u.4mat',ModelRank,NetRank),NetDir,PsNb,PsNb,'uint8','ieee-le',PsRanks,PsRanks);
                            %don't keep Corr that are inferior to Anti
                            Corr(Corr<=Anti)=0;
                            clear Anti;
                            %make adjacency matrix
                            Corr=uint8(Corr>CorrLimit);
                            for i=1:length(Corr)
                                Corr(i,i)=0;
                            end                            
                            %sort by connectivity
                            Conn=sum(Corr);
                            StartVal(ConnL,1)=max(Conn);
                            CorrSize(ConnL,1)=length(Corr);
                            [Temp,SortIndex]=sort(Conn);
                        end
                        if max(Conn)==0
                            %no other probe set correlated to the current probe set
                            CurrConn(ConnPos(1))=0;
                            ConnPos=find(CurrConn==ConnVal(ConnL));
                            if isempty(ConnPos)|ConnVal(ConnL)==0
                                Ok=0;
                            end
                        else
                            if length(PsRanks)==2
                                CliquerVal(ConnL)=2;
                                CurrCluster=PsRanks;
                            else
                                %write ASCII file for cliquer program
                                fid=fopen(sprintf('m%un%u_cliquer_%u.txt',ModelRank,NetRank,CorrLimit),'w');
                                fprintf(fid,'c rank %u\n',CliqueRank+1);
                                fprintf(fid,'p edge %u %u\n',length(Corr),max(1,uint32((sum(sum(Corr))-length(Corr))/2)));
                                for PsL=1:length(Corr)
                                    Pos=find(Corr(SortIndex(PsL),:));
                                    Pos=Pos(find(Pos>SortIndex(PsL)));
                                    if ~isempty(Pos)
                                        for PosL=1:length(Pos)
                                            fprintf(fid,'e %u %u\n',SortIndex(PsL),Pos(PosL));
                                        end
                                    end
                                end
                                fclose(fid);
                                eval(sprintf('! /usr/local/cliquer/cl -su %s/m%03u/n%05u/m%un%u_cliquer_%u.txt>%s/m%03u/n%05u/m%un%u_cliquer_res_%u.txt;',K.dir.net,ModelRank,NetRank,ModelRank,NetRank,CorrLimit,K.dir.net,ModelRank,NetRank,ModelRank,NetRank,CorrLimit))
                                eval(sprintf('PsPos=load(''m%un%u_cliquer_res_%u.txt'');',ModelRank,NetRank,CorrLimit));
                                %The maximal clique found (can be inferior to the real greatest clique because calculus is stoped after 60"
                                Cliquer=PsRanks(PsPos);
                                if length(Corr)>500&MyAlgoFlag==1
                                    sprintf('using my algorithm')
                                    CurrLimit=length(Corr);
                                    while min(sum(Corr))<CurrLimit
                                        ClearIndex=find(sum(Corr)==min(sum(Corr)));
                                        PsRanks(ClearIndex)=[];
                                        Corr(ClearIndex,:)=[];
                                        Corr(:,ClearIndex)=[];
                                        CurrLimit=length(Corr);
                                        %sprintf('%u - %u',min(sum(Corr)),CurrLimit)
                                    end
                                    MyVal(ConnL)=length(PsRanks);
                                    CliquerVal(ConnL)=length(Cliquer);
                                    if length(Cliquer)>length(PsRanks)
                                        CurrCluster=Cliquer;
                                        CurrPsPos=PsPositions(PsPos);
                                    else
                                        sprintf('better than cliquer: %u instead of %u',length(PsRanks),length(Cliquer))
                                        CurrCluster=PsRanks;
                                        CurrPsPos=[];
                                        for PsL=1:length(PsRanks)
                                            CurrPsPos(end+1,1)=find(LoadedPsRanks==PsRanks(PsL));
                                        end
                                    end
                                else
                                    CliquerVal(ConnL)=length(Cliquer);
                                    CurrCluster=Cliquer;
                                    CurrPsPos=PsPositions(PsPos);
                                end
                            end
                            if ~isempty(CurrCluster)
                                %add a new clique
                                CliqueRank=CliqueRank+1;
                                Clu(CurrCluster)=CliqueRank;
                                CliqueSize(CurrCluster)=length(CurrCluster);
                                CluSize(CliqueRank)=length(CurrCluster);
                                FoundPsRanks=sort([FoundPsRanks;CurrCluster]);
                                if ~isempty(find( CurrCluster==CurrPsRank))
                                    IsSeed(CurrPsRank)=1;
                                else
                                    IsSeed(CurrPsRank)=2;
                                    %eliminate the CurrProbeSet to prevent
                                    %infinite loop
                                    CurrConn(ConnPos(1))=0;
                                end
                                %clear the probe sets belonging to the current
                                %clique
                                CurrConn(CurrPsPos)=0;
                                ConnPos=find(CurrConn==ConnVal(ConnL));
                                if isempty(ConnPos)||ConnVal(ConnL)==0
                                    Ok=0;
                                end
                            else %CurrCLuster is empty
                                h=errordlg(sprintf('no clique for probe set %u => error',CurrPsRank));
                                waitfor(h)
                                error('process canceled')
                            end
                        end
                    else
                        %the current ps is eliminated from the list
                        CurrConn(ConnPos(1))=0;
                        %to allow to proceed with the next one
                        ConnPos=find(CurrConn==ConnVal(ConnL));
                        if isempty(ConnPos)||ConnVal(ConnL)==0
                            Ok=0;
                        end
                    end
                end %while Ok
            end %if ~isempty(ConnPos)           
        end %ConnL
        ConnFound(1)=length(FoundPsRanks)-MemPsNb;
        TestedNb(1)=TestedPsNb;
        sprintf('ConnL=1: connectivity %u => %u found probesets among %u',ConnL,ConnVal(ConnL+1),ConnFound(ConnL+1),TestedPsNb)
        %update MemPsNb to find the number of found Ps in the next step
        MemPsNb=length(FoundPsRanks);
        if DisplayFlag
            if MyAlgoFlag
                PlotNb=3
            else
                PlotNb=2
            end
            h=figure;
            set(h,'name',sprintf('CLIQUES PROPERTIES ROUND %u',Round))

            subplot(1,PlotNb,1)
            plot(ConnVal,ConnFound,'r.')
            hold on
            plot(CorrSize,ConnFound,'g.')
            set(gca,'yscale','log')
            xlabel('Nominal (red) or real (g) connectivity')
            ylabel('Nb of probe sets agregated')
            title('NB OF PROBESET AGREGATED VS CONNECTIVITy')

            subplot(1,PlotNb,2)
            plot(ConnVal,CorrSize,'b.')
            xlabel('Nominal connectivity')
            ylabel('Connectivity when used')
            set(gcf,'color',[1,1,1])
            title('REAL VS NOMINAL CONNECTIVITy')
            
            if MyAlgoFlag
                subplot(1,PlotNb,3)
                Pos=find(MyVal>CliquerVal);
                Neg=find(MyVal<CliquerVal);
                Equ=find(MyVal==CliquerVal);
                plot(MyVal(Pos),CliquerVal(Pos),'b.')
                hold on
                plot(MyVal(Neg),CliquerVal(Neg),'c.')
                plot(MyVal(Equ),CliquerVal(Equ),'g.')
                xlabel('results MyCluster')
                ylabel('results Cliquer')
                title('PS FOUND BY MyCluster and by Cliquer')
            end


            h=figure;
            set(h,'name','CLIQUE SIZE VS CONNECTIVITy OF SEED PROBE SETS')
            plot(MemConn(IsSeed(LoadedPsRanks)==0),CliqueSize(IsSeed(LoadedPsRanks)==0),'b.')
            hold on
            plot(MemConn(IsSeed(LoadedPsRanks)==1),CliqueSize(IsSeed(LoadedPsRanks)==1),'ro')
            xlabel('connectivity')
            ylabel('clique size')
            title('red: seeding probe set, blue: agregated probe sets')
            set(gcf,'color',[1,1,1])


            h=figure;
            set(h,'name','NUMBER OF PS FOUND IN CLIQUES VS CONNECTIVITy')
            hold on
            plot(ConnVal,ConnTotal,'k.')
            plot(ConnVal(ConnSearched==1),ConnFound(ConnSearched==1),'b.')
            CumFound=cumsum(ConnFound);
            plot(ConnVal(ConnSearched==1),sum(ConnFound)-CumFound(ConnSearched==1)+ConnFound(end),'g.')
            set(gca,'box','on')
            set(gcf,'color',[1,1,1])
            xlabel('connectivity')
            ylabel('black: total number of linked ps, blue : ps found in clique, green ; cum sum')

            Sizes=unique(CluSize);
            CluNb=histc(CluSize,Sizes);
            h=figure;
            set(h,'name','CLIQUE SIZE DISTRIBUTION')
            plot(Sizes,CluNb,'ro')
            set(gca,'yscale','log')
            set(gcf,'color',[1,1,1])
            xlabel('clique size')
            ylabel('frequency')
            
            h=figure;
            set(h,'name','CLIQUE SIZE VS SIZE OF SEARCHED NETWORK')
            plot(TestedNb,ConnFound,'g.')
            set(gcf,'color',[1,1,1])
            ylabel('clique size')
            xlabel('searched network size')

        end %if DisplayFlag
    end %if isempty(find(CurrConn>0))|(Round>1&max(Clu)==CliNb) ... else
end %while Continue
