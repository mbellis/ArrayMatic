%==========================%
% FUNCTION NET_DISPLAYNTS  %
%==========================%

% NET_DISPLAYNTS display the reltionship between different clusters found by NetsTensor at different densities

%INPUT PARAMETERS
% 1    ChipRank: chip rank
% 2     PsList: list of probe set classes used
% 3    UsedPsNb: number of probeset to be exported (=0 => ALL)
% 4    NetRanks: networks used by TensorNet
% 5        Type: Type of calcul: either 1(C) or 2(C-A)
% 6   CorrLimit: inferior limit of C or C-A
% 7 GeneNbLimit: minimum number of gene in a TSN cluster
% 8  NetNbLimit: minimum of networks where a TSN cluster must be
% 9   Densities: densities used by tensornet

%OUTPUT FILE

%net_displaynts(8,11,[1:10],0,[212:226],2,1,2,10,5,5,[40:10:80]);
%net_displaynts(27,5,[1:10],0,[149:163],2,1,2,10,5,5,[40:10:80]);
%net_displaynts([8,27],{11,5},{[1:10],[1:10]},[0,0],{[212:226],[149:163]},[2,2],[1,1],[2,2],[10,10],[5,5],[5,5],{[40:10:80],[40:10:80]});
%net_displaynts([8],{11},{[1:10]},[0],{[212:226]},[2],[1],[2],[10],[5],[5],{[40:10:80]});
%net_displaynts([8,8],{11,11},{[1:10],[1:10]},[0,0],{[212:217,219],[220:226]},[2,2],[1,1],[2,2],[10,10],[5,5],[3,3],{[30:10:60],[30:10:60]});
%net_displaynts([27,27],{5,5},{[1:10],[1:10]},[0,0],{[149:155],[156:162]},[2,2],[1,1],[2,2],[10,10],[5,5],[3,3],{[30:10:60],[30:10:60]});
%net_displaynts([8,8],{11,11},{[1:10],[1:10]},[0,0],{[212:219],[220:226]},[2,2],[1,1],[2,2],[10,10],[5,5],[5,3],{[40:10:60],[40:10:60]});
%net_displaynts(27,{5},{[1:10]},0,{[149:155]},2,1,2,10,5,3,{[40:10:60]});
%net_displaynts(27,{5},{[1:10]},0,{[156:162]},2,1,2,10,5,3,{[40:10:60]});
%net_displaynts([8,8],{5,5},{0,0},[0,0],{[212:217],[218:223]},[2,2],[1,1],[2,2],[10,10],[5,5],[3,3],{[30,50],[30,50]});
%direct comparison between m8n228 and m27n164
%net_displaynts(8,{17},{0},0,{[228]},2,1,2,10,5,4,{[10:5:30]});

%net_displaynts(27,{{'crel_2_gene2psrank'}},{0},0,{[149:163]},0,0,2,10,3,5,{[30:5:70]});
%net_displaynts(8,{{'crel_2_gene2psrank'}},{0},0,{[212:226]},0,0,2,10,3,5,{[30:5:70]});
%net_displaynts([8,27],{{'crel_2_gene2psrank'},{'crel_2_gene2psrank'}},{0,0},[0,0],{[212:226],[149:163]},[0,0],[0,0],[2,2],[10,10],[3,3],[5,5],{[30:5:70],[30:5:70]});

%net_displaynts(27,{{'mouse_krebs_proteasome_mapk_gene2psrank'}},{0},0,{[149:163]},0,0,2,0,10,3,5,{[30:5:70]});
%net_displaynts(8,{{'mouse_krebs_proteasome_mapk_gene2psrank'}},{0},0,{[212:226]},0,0,2,0,10,3,5,{[30:5:70]});
%net_displaynts([8,27],{{'mouse_krebs_proteasome_mapk_gene2psrank'},{'mouse_krebs_proteasome_mapk_gene2psrank'}},{0,0},[0,0],{[212:226],[149:163]},[0,0],[0,0],[2,2],[10,10],[3,3],[5,5],{[30:5:70],[30:5:70]});

%test sur sous-ensemble de probesets ayant ou non
%net_displaynts(27,{[92,93]},{0},{[149:163]},2,0,5,10,{[30:5:70]});
%net_displaynts(8,{[185,186]},{0},{[212:226]},2,0,5,10,{[30:5:70]});
%net_displaynts(27,{[93]},{0},{[149:163]},2,0,5,10,{[30:5:70]});

%net_displaynts([91,91],{[1],[1]},{0,0},{[4:18],[19:34]},[2,2],[0,0],[5,5],[10,10],{[30:5:70],[30:5:70]});

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


function net_displaynts(ChipRank,PsList,UsedPsNb,NetRanks,Type,CorrLimit,GeneNbLimit,NetNbLimit,Densities)
global K

if length(ChipRank)>2
    h=errordlg('net_comparemcls need at most two chips');
    waitfor(h)
    error('process canceled')
end

%recover chip information
ChipNb=length(ChipRank);
ChipPos=find(K.chip.rank==ChipRank(1));
PsNb(1)=K.chip.probesetNb(ChipPos);
NetNb(1)=length(NetRanks{1});
if ChipNb==2
    ChipPos=find(K.chip.rank==ChipRank(2));
    PsNb(2)=K.chip.probesetNb(ChipPos);
    NetNb(2)=length(NetRanks{2});
end

%recover eventually correspondance between different chips
if ChipNb==2 & ChipRank(1)~=ChipRank(2)
    cd(K.dir.affyMetadata)
    FileName=sprintf('m%u_m%u_commonps.mat',min(ChipRank),max(ChipRank));
    if exist(FileName,'file')
        load(FileName)
        if ChipRank(1)>ChipRank(2)
            Temp=ComPsRank;
            ComPsRank(:,1)=ComPsRank(:,2);
            ComPsRank(:,2)=Temp(:,1);
            clear Temp
        end
        MatchPsRank=zeros(max(ComPsRank(:,1)),1);
        MatchPsRank(ComPsRank(:,1))=ComPsRank(:,2);
    else
        h=errordlg(sprintf('no correspondance file between m%y and m%u',ChipRank(1),ChipRank(2)));
        waitfor(h)
        error('process canceled')
    end
end



for ChipL=1:ChipNb
    %recover dir for TensorNet
    TensorDir=fullfile(K.dir.net,sprintf('m%u',ChipRank(ChipL)),'tsn','result');
    cd(TensorDir)
    if Type==2
        Prefix=sprintf('m%u_n%uton%u_c%u_TSN_diff',ChipRank(ChipL),NetRanks{ChipL}(1),NetRanks{ChipL}(end),CorrLimit(ChipL));
    else
        Prefix=sprintf('m%u_n%uton%u_c%u_TSN',ChipRank(ChipL),NetRanks{ChipL}(1),NetRanks{ChipL}(end),CorrLimit(ChipL));
    end
    %recover data
    for ListL=1:length(PsList{ChipL})
        if iscell(PsList{ChipL})
            CurrFile=PsList{ChipL}{ListL};
        else
            CurrClass=PsList{ChipL}(ListL);
        end
        if iscell(PsList{ChipL})
            cd(K.dir.list)
            load([CurrFile,'.mat'])
            CurrIndex=GenePos{ChipRank(ChipL)};
            %construct Index and GeneList
            Index1=[];
            GeneList{ChipL}{ListL}={};
            for GeneL=1:length(CurrIndex)
                Index1=[Index1;CurrIndex{GeneL}];
                for PsL=1:length(CurrIndex{GeneL})
                    GeneList{ChipL}{ListL}{end+1,1}=Gene{GeneL};
                end
            end
            [Index1,SortIndex]=unique(Index1);
            Index{ChipL}={};
            Index{ChipL}{ListL}{1}=Index1;
            GeneList{ChipL}{ListL}=GeneList{ChipL}{ListL}(SortIndex);
        end
        %if CurrFile exists, no clusters are used but the following loop is used once

        if ~iscell(PsList{ChipL})
            cd(fullfile(K.dir.net,sprintf('m%u',ChipRank(ChipL)),'list'));
            %recover ps lists
            ListRank=PsList{ChipL}(ListL);
            Index{ChipL}{ListL}=load_data(sprintf('m%u_pslist%u.u32',ChipRank(ChipL),ListRank),'./',0,0,'uint32','ieee-le');
            %restrict list if UsedPsNb>0
            if UsedPsNb{ChipL}>0&UsedPsNb{ChipL}<length(Index{ChipL}{ListL})
                Index{ChipL}{ListL}=Index{ChipL}{ListL}(1:UsedPsNb);
            end
        end
        cd(TensorDir)
        %recover used densities, used ps, and networks where clusters are found
        UsedDensities{ChipL}{ListL}=[];
        %UsedPs{ChipL}{ListL}=[];
        UsedNet{ChipL}{ListL}=zeros(1,NetNb(ChipL));
        for DensL=1:length(Densities{ChipL})
            if iscell(PsList{ChipL})
                ResultFile=sprintf('%s_%s_g%u_n%u_d%u.txt',Prefix,CurrFile,GeneNbLimit(ChipL),NetNbLimit(ChipL),Densities{ChipL}(DensL));
            else
                ResultFile=sprintf('%s_%u_g%u_n%u_d%u.txt',Prefix,ListRank,GeneNbLimit(ChipL),NetNbLimit(ChipL),Densities{ChipL}(DensL));
            end
            ResultFile
            exist(ResultFile)
            if exist(ResultFile)
                [NtsClusters{ChipL}{ListL}{DensL},NtsPsNbs{ChipL}{ListL}{DensL},NtsDensities{ChipL}{ListL}{DensL},NtsNetNbs{ChipL}{ListL}{DensL},NtsNetLists{ChipL}{ListL}{DensL}]=textread(ResultFile,'%s%u%.2f%u%s','delimiter','\t');
                UsedDensities{ChipL}{ListL}=[UsedDensities{ChipL}{ListL},DensL];
                CurrClu=NtsClusters{ChipL}{ListL}{DensL};
                for NtsL=1:length(CurrClu)
                    CurrClu{NtsL}=sort(eval(CurrClu{NtsL}));
                    if ChipL==2
                        %try
                        NtsClusters{ChipL}{ListL}{DensL}{NtsL}=Index{ChipL}{ListL}(CurrClu{NtsL});
                        %catch
                        %    NtsClusters{ChipL}{ListL}{DensL}{NtsL}=[];
                        %end
                    else
                        if ChipNb==2 & ChipRank(1)~=ChipRank(2)
                            %translation of probe set ranks of the first chip into the
                            %second one
                            Val=sort(MatchPsRank(Index{ChipL}{ListL}(CurrClu{NtsL})));
                            %some probe set of the first chip do not have a corresponding probe set in the
                            %second chip
                            Val(find(Val==0))=[];
                            NtsClusters{ChipL}{ListL}{DensL}{NtsL}=Val;
                        else
                            NtsClusters{ChipL}{ListL}{DensL}{NtsL}=Index{ChipL}{ListL}(CurrClu{NtsL});
                        end
                    end
                    %UsedPs{ChipL}{ListL}=union(UsedPs{ChipL}{ListL},NtsClusters{ChipL}{ListL}{DensL}{NtsL});
                end
                %recover networks positive for the current cluster
                CurrNet=NtsNetLists{ChipL}{ListL}{DensL};
                for NtsL=1:length(CurrClu)
                    CurrNetRanks=regexp([' ',CurrNet{NtsL}],' \d+(?=:)','match');
                    for NetL=1:length(CurrNetRanks)
                        CurrNetRank=str2num(CurrNetRanks{NetL});
                        try
                            UsedNet{ChipL}{ListL}(CurrNetRank)=UsedNet{ChipL}{ListL}(CurrNetRank)+1;
                        catch
                            'stop'
                        end
                    end
                end
            end
        end

    end

    %analyze data
    %frequency of net nb
    Colors=colors(colormap,length(Densities{ChipL}));
    h=figure;
    set(h,'name',sprintf('distribution of NetNb in m%u',ChipRank(ChipL)))
    set(gcf,'color',[1,1,1])
    for ListL=1:length(PsList{ChipL})
        ColNb=ceil(sqrt(length(PsList{ChipL})));
        RowNb=round(ColNb/length(PsList{ChipL}));
        subplot(RowNb,ColNb,ListL)
        hold on
        Legend={};
        for DensL=1:length(Densities{ChipL})
            CurrNetNb=[];
            try
                CurrNetNb=NtsNetNbs{ChipL}{ListL}{DensL};
            catch
            end
            if ~isempty(CurrNetNb)
                Legend{end+1}=['d',num2str(Densities{ChipL}(DensL))];
                Val=histc(CurrNetNb,[1:NetNb(ChipL)]);
                plot([1:NetNb(ChipL)],cumsum(fliplr(Val)),'color',Colors(DensL,:))
            end
        end
        set(gca,'box','on')
        if iscell(PsList{ChipL})
            title(sprintf('m%u %s nb of networks',ChipRank(ChipL),strrep(PsList{ChipL}{ListL},'_',' ')))
        else
            title(sprintf('m%u list%u',ChipRank(ChipL),PsList{ChipL}(ListL)))
        end
        xlabel('nb of networks positive')
        ylabel('nb of cluster')
        legend(Legend)
    end
end

%find correspondance between input clusters

'stop'


for ListL=1:length(PsList{1});
    Offset=0;
    DensNb=[0,0];
    NetNbs=[];
    %calculate the number of rows and columns for subplot
    CurrDens={};
    for ChipL=1:ChipNb     
        CurrDens{ChipL}=UsedDensities{ChipL}{ListL};
        DensNb(ChipL)=length(CurrDens{ChipL});
        for DensL=1:DensNb(ChipL)
            CurrNetNb=NtsNetNbs{ChipL}{ListL}{DensL};
            NetNbs=union(NetNbs,unique(CurrNetNb));
            if ~isempty(NetNbs)
            end
        end
    end
    MinNetNb=min(NetNbs);
    RowNb=length(NetNbs);
    if DensNb(1)>0
        Colors=colors(colormap,DensNb(1));
        h1=figure;
        set(gcf,'color',[1,1,1])
        if iscell(PsList{1})
            if ChipNb==2
                set(h1,'name',sprintf('m%u vs m%u display',ChipRank(1),ChipRank(2)))
            else
                set(h1,'name',sprintf('m%u display',ChipRank(1)))
            end
        else
            if ChipNb==2
                set(h1,'name',sprintf('m%u vs m%u NTS cluster display',ChipRank(1),ChipRank(2)))
            else
                set(h1,'name',sprintf('m%u NTS cluster display',ChipRank(1)))
            end
        end
        h2=figure;
        set(gcf,'color',[1,1,1])
        h3=figure;
        set(gcf,'color',[1,1,1])
        if ChipNb==2
            h4=figure;
            set(gcf,'color',[1,1,1])
        end
        figure(h2)
        set(h2,'name','ps redundancy')
        figure(h3)
        set(h3,'name','NTS cluster size')
        if ChipNb==2
            figure(h4)
            set(h4,'name','NTS cluster size')
        end
        CluSize{1}=[];
        CluSize{2}=[];
        for DensL=1:DensNb(1)
            DensPos=0;
            if ChipNb==2
                for DensL2=1:DensNb(2)
                    if Densities{1}(CurrDens{1}(DensL))==Densities{2}(CurrDens{2}(DensL2))+Offset
                        DensPos=DensL2;
                        break
                    end
                end
            end
            Clu=[];
            if ChipNb==1|(ChipNb==2&DensPos>0)
                CurrNetNb=[];
                CluNb=[0,0];
                %recover clusters of the first chip
                CurrClusters=NtsClusters{1}{ListL}{CurrDens{1}(DensL)};
                CluNb(1)=length(CurrClusters);
                CurrNetNb{1}=NtsNetNbs{1}{ListL}{CurrDens{1}(DensL)};
                if ChipNb==1
                    Clu=zeros(CluNb(1),PsNb(1));
                else
                    Clu=zeros(CluNb(1),PsNb(2));
                end
                %reorder Nts clusters{1} according to the number of positive networks
                [CurrNetNb{1},SortIndex]=sort(CurrNetNb{1},'descend');
                MemCurrNetNb{1}=CurrNetNb{1};
                CurrClusters=CurrClusters(SortIndex);
                for CluL=1:CluNb(1)
                    CluSize{1}(end+1)=length(CurrClusters{CluL});
                    Clu(CluL,CurrClusters{CluL})=round((64/RowNb)*(CurrNetNb{1}(CluL)-MinNetNb+1));
                end
                Clu=[Clu;zeros(1,size(Clu,2))];
                if ChipNb==2
                    CurrClusters=NtsClusters{2}{ListL}{CurrDens{2}(DensPos)};
                    CluNb(2)=length(CurrClusters);
                    CurrNetNb{2}=[NtsNetNbs{2}{ListL}{CurrDens{2}(DensPos)}];
                    CurrClu=zeros(CluNb(2),PsNb(2));
                    %reorder PredClusters{2} according to the number of positive networks
                    [CurrNetNb{2},SortIndex]=sort(CurrNetNb{2},'descend');
                    CurrClusters=CurrClusters(SortIndex);
                    for CluL=1:CluNb(2)
                        CluSize{2}(end+1)=length(CurrClusters{CluL});                     
                        CurrClu(CluL,CurrClusters{CluL})=round((64/RowNb)*(CurrNetNb{2}(CluL)-MinNetNb+1));
                    end
                    Clu=[Clu;CurrClu];
                    CurrNetNb{1}=[CurrNetNb{1};100;CurrNetNb{2}];
                end

                %make clusters displayed as continuous patch
                %eliminate not used probe sets
                Clu(:,find(sum(Clu,1)==0))=[];
                [temp,RowIndex]=sort(CurrNetNb{1},'descend');
                Index1=1:size(Clu,2);
                %process each cluster in descending order
                for CluL1=1:size(Clu,1)
                    %test if there still exist columns of Clu to be reordered
                    if length(Index1)>1
                        [temp SortIndex]=sort(Clu(RowIndex(CluL1),Index1),'descend');
                        Clu(:,Index1)=Clu(:,Index1(SortIndex));
                        Pos=find(Clu(RowIndex(CluL1),Index1));
                        Index1=Index1(1)+length(Pos):size(Clu,2);
                    end
                end

                NetNbs=setdiff(unique(CurrNetNb{1}),100);
                for NetNbL=1:length(NetNbs)
                    Pos=find(CurrNetNb{1}>=NetNbs(NetNbL));
                    CurrClu=Clu(Pos,:);
                    CurrClu(:,find(sum(CurrClu,1)==0))=[];
                    figure(h1)
                    subplot(RowNb,DensNb(1),(NetNbL-1)*DensNb(1)+DensL)
                    %Clu image
                    image(CurrClu)
                    if ChipNb==1
                        title(sprintf('d %u NetNb %u',Densities{1}(CurrDens{1}(DensL)),NetNbs(NetNbL)))
                    else
                        title(sprintf('d %u-%u NetNb %u',Densities{1}(CurrDens{1}(DensL)),Densities{2}(CurrDens{2}(DensPos)),NetNbs(NetNbL)))
                    end
                    %distribution of ps redondancy
                    CurrClu(find(CurrClu))=1;
                    SumClu=sum(CurrClu);
                    figure(h2)
                    subplot(RowNb,DensNb(1),(NetNbL-1)*DensNb(1)+DensL)
                    BinVal=histc(SumClu,[1:max(SumClu)]);
                    plot([1:max(SumClu)],BinVal,'color',Colors(DensL,:))
                    if ChipNb==1
                        title(sprintf('d %u NetNb %u',Densities{1}(CurrDens{1}(DensL)),NetNbs(NetNbL)))
                    else
                        title(sprintf('d %u-%u NetNb %u',Densities{1}(CurrDens{1}(DensL)),Densities{2}(CurrDens{2}(DensPos)),NetNbs(NetNbL)))
                    end
                    set(gca,'xlim',[0,15])
                    xlabel('network redundancy')
                    ylabel('probeset frequency')

                    figure(h3)
                    Pos1=find(MemCurrNetNb{1}>=NetNbs(NetNbL));
                    if ~isempty(Pos1)
                        subplot(RowNb,DensNb(1),(NetNbL-1)*DensNb(1)+DensL)
                        BinVal=histc(CluSize{1}(Pos1),[1:max(CluSize{1}(Pos1))]);
                        BinVal=BinVal';                        
                        bar([1:max(CluSize{1}(Pos1))],BinVal)                        
                        xlabel('cluster size')
                        ylabel('cluster frequency')
                        set(gca,'box','on')
                    end
                    if ChipNb==2
                        figure(h4)
                        subplot(RowNb,DensNb(1),(NetNbL-1)*DensNb(1)+DensL)
                        Pos2=find(CurrNetNb{2}>=NetNbs(NetNbL));
                        hold on
                        BinVal=histc(CluSize{2}(Pos2),[1:max(CluSize{2}(Pos2))]);
                        BinVal=BinVal';
                        bar([1:max(CluSize{2}(Pos2))],BinVal)
                        xlabel('cluster size')
                        ylabel('cluster frequency')
                        set(gca,'box','on')
                    end
                    if ChipNb==1
                        title(sprintf('d %u NetNb %u',Densities{1}(CurrDens{1}(DensL)),NetNbs(NetNbL)))
                    else
                        title(sprintf('d %u-%u NetNb %u',Densities{1}(CurrDens{1}(DensL)),Densities{2}(CurrDens{2}(DensPos)),NetNbs(NetNbL)))
                    end
                end
            end
        end
    end
end

'stop'

%display CORR ordored according to NTS clustering

dens=40;
netnb=12;
chip=1;
list=1;

NetRank=4;
%NetList{1}=[149:163];
%NetList{1}=134;
%NetList{1}=228;
%NetList{1}=1;
%NetList{2}=2;

%NetList{1}=[25:45];
%NetList{2}=[212:226];
%NetList{2}=196;
%ChipRank(2)=8;
%PsNb(2)=45101;
%NetList{2}=[149:163];
%NetList{1}=[212:226];
%ChipRank(2)=27;
%PsNb(2)=22690;
%NetList{2}=171;
%RowNb=ceil(sqrt(length(NetList{1})));
%ColNb=round(length(NetList{1})/RowNb);

if length(ChipRank)==2
    if ChipRank(1)~=ChipRank(2)
        cd(K.dir.affyMetadata)
        FileName=sprintf('m8_m27_commonps.mat');
        MatchPsRank={};
        if exist(FileName,'file')
            load(FileName)
            if ChipRank(1)<ChipRank(2)
                MatchPsRank{1}=zeros(max(ComPsRank(:,1)),1);
                MatchPsRank{1}(ComPsRank(:,1))=ComPsRank(:,2);
                MatchPsRank{2}=zeros(max(ComPsRank(:,2)),1);
                MatchPsRank{2}(ComPsRank(:,2))=ComPsRank(:,1);
            else
                MatchPsRank{2}=zeros(max(ComPsRank(:,1)),1);
                MatchPsRank{2}(ComPsRank(:,1))=ComPsRank(:,2);
                MatchPsRank{1}=zeros(max(ComPsRank(:,2)),1);
                MatchPsRank{1}(ComPsRank(:,2))=ComPsRank(:,1);
            end
        end
    end
end

DensPos=find(Densities{1}==dens);
if ~isempty(DensPos)
    Scale=[];
    MemPos=0.5;
    CurrClu=NtsClusters{chip}{list}{DensPos};
    CurrNetNb=NtsNetNbs{chip}{list}{DensPos};
    DiffRank={};
    DiffPos={};
    PsRank=[];
    CurrGene={};
    CluPsNb=0;
    for NetNbL=max(CurrNetNb):-1:netnb
        CluPos=find(CurrNetNb==NetNbL);
        if ~isempty(CluPos)
            for CluL=1:length(CluPos)
                CurrPsRank=setdiff(CurrClu{CluPos(CluL)},PsRank);
                if ~isempty(CurrPsRank)
                    DiffRank{end+1,1}=CurrPsRank;
                    DiffPos{end+1,1}=[1:length(CurrPsRank)]+CluPsNb;
                    CluPsNb=CluPsNb+length(CurrPsRank);
                    for PsL=1:length(CurrPsRank)
                        PsRank(end+1,1)=CurrPsRank(PsL);
                        PsPos=find(Index{chip}{list}==CurrPsRank(PsL));
                        if iscell(PsList{1})
                            CurrGene{end+1,1}=GeneList{chip}{list}{PsPos};
                        end
                    end
                    Scale(end+1,1)=MemPos+length(CurrPsRank);
                    MemPos=Scale(end,1);
                end
            end
        end
    end
    %    %complete with missing ps
    %    CurrPsRank=setdiff(Index{chip}{list},PsRank);
    %    for PsL=1:length(CurrPsRank)
    %        PsRank(end+1,1)=CurrPsRank(PsL);
    %        PsPos=find(Index{1}{c}{p}==CurrPsRank(PsL));
    %        CurrGene{end+1,1}=GeneList{chip}{list}{PsPos};
    %    end
end

%h=figure;
%set(h,'name',sprintf('m%u; NTS clusters of %s',ChipRank(1),PsList{1}{1}))
%set(gcf,'color',[1,1,1])
%for NetL=1:length(NetList{1})

cd(fullfile(K.dir.net,sprintf('m%u',ChipRank(chip)),sprintf('n%u',NetRank)))
C=load_data(sprintf('c_m%u_n%u.4mat',ChipRank(chip),NetRank),'./',PsNb(1),PsNb(1),'uint8','ieee-le',PsRank,PsRank);
C=double(C);
A=load_data(sprintf('a_m%u_n%u.4mat',ChipRank(chip),NetRank),'./',PsNb(1),PsNb(1),'uint8','ieee-le',PsRank,PsRank);
A=double(A);

if NetL==1
    C1=C;
    A1=A;
end
C=C-A;

if length(ChipRank)==2
    if ChipRank(1)~=ChipRank(2)
        if chip==2
            CurrPsRank=MatchPsRank{1}(PsRank);
            NullPos=find(CurrPsRank==0);
            C(NullPos,:)=[];
            C(:,NullPos)=[];
            A(NullPos,:)=[];
            A(:,NullPos)=[];
        end
    end
end

Cu=triu(C);
Cu(find(Cu<0))=0;
Cl=tril(C,-1);
Cl(find(Cl>0))=0;
C=Cu-Cl;
Val=setdiff(unique(C),100);
C(find(C>max(Val)))=max(Val);

h1=figure;
if iscell(PsList{1})
    set(h1,'name',sprintf('m%u: NTS clusters of %s (netw %u dens %u netwnb>=%u)',ChipRank(chip),strrep(PsList{chip}{list},'_','-'),NetRank,dens,netnb));
else
    set(h1,'name',sprintf('m%u: NTS clusters of list %u (netw %u dens %u netwnb>=%u)',ChipRank(chip),PsList{1}(list),NetRank,dens,netnb));
end
set(gcf,'color',[1,1,1])
image(C)
if iscell(PsList{1})
    set(gca,'ytick',[1:1:length(CurrGene)])
    set(gca,'yticklabel',CurrGene)
end
set(gca,'xtick',Scale)
set(gca,'tickdir','out')
set(gca,'xticklabel','')
if iscell(PsList{1})
    title(sprintf('m%u: NTS clusters of %s (netw %u dens%u netw>=%u)',ChipRank(chip),strrep(PsList{chip}{list},'_','-'),NetRank,dens,netnb));
else
    title(sprintf('m%u: NTS clusters of class %u (netw %u dens%u netw>=%u)',ChipRank(chip),PsList{chip}(list),NetRank,dens,netnb));
end




