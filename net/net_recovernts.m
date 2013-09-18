%==========================%
% FUNCTION NET_RECOVERNTS  %
%==========================%

% NET_RECOVERNTS recovers NetsTensor results

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

%net_recovernts(8,11,[1:10],0,[212:226],2,1,2,10,5,5,[40:10:80]);
%net_recovernts(27,5,[1:10],0,[149:163],2,1,2,10,5,5,[40:10:80]);
%net_recovernts([8,27],{11,5},{[1:10],[1:10]},[0,0],{[212:226],[149:163]},[2,2],[1,1],[2,2],[10,10],[5,5],[5,5],{[40:10:80],[40:10:80]});
%net_recovernts([8],{11},{[1:10]},[0],{[212:226]},[2],[1],[2],[10],[5],[5],{[40:10:80]});
%net_recovernts([8,8],{11,11},{[1:10],[1:10]},[0,0],{[212:217,219],[220:226]},[2,2],[1,1],[2,2],[10,10],[5,5],[3,3],{[30:10:60],[30:10:60]});
%net_recovernts([27,27],{5,5},{[1:10],[1:10]},[0,0],{[149:155],[156:162]},[2,2],[1,1],[2,2],[10,10],[5,5],[3,3],{[30:10:60],[30:10:60]});
%net_recovernts([8,8],{11,11},{[1:10],[1:10]},[0,0],{[212:219],[220:226]},[2,2],[1,1],[2,2],[10,10],[5,5],[5,3],{[40:10:60],[40:10:60]});
%net_recovernts(27,{5},{[1:10]},0,{[149:155]},2,1,2,10,5,3,{[40:10:60]});
%net_recovernts(27,{5},{[1:10]},0,{[156:162]},2,1,2,10,5,3,{[40:10:60]});
%net_recovernts([8,8],{5,5},{0,0},[0,0],{[212:217],[218:223]},[2,2],[1,1],[2,2],[10,10],[5,5],[3,3],{[30,50],[30,50]});
%direct comparison between m8n228 and m27n164
%net_recovernts(8,{17},{0},0,{[228]},2,1,2,10,5,4,{[10:5:30]});

%net_recovernts(27,{{'crel_2_gene2psrank'}},{0},0,{[149:163]},0,0,2,10,3,5,{[30:5:70]});
%net_recovernts(8,{{'crel_2_gene2psrank'}},{0},0,{[212:226]},0,0,2,10,3,5,{[30:5:70]});
%net_recovernts([8,27],{{'crel_2_gene2psrank'},{'crel_2_gene2psrank'}},{0,0},[0,0],{[212:226],[149:163]},[0,0],[0,0],[2,2],[10,10],[3,3],[5,5],{[30:5:70],[30:5:70]});

%net_recovernts(27,{{'mouse_krebs_proteasome_mapk_gene2psrank'}},{0},0,{[149:163]},0,0,2,0,10,3,5,{[30:5:70]});
%net_recovernts(8,{{'mouse_krebs_proteasome_mapk_gene2psrank'}},{0},0,{[212:226]},0,0,2,0,10,3,5,{[30:5:70]});
%net_recovernts([8,27],{{'mouse_krebs_proteasome_mapk_gene2psrank'},{'mouse_krebs_proteasome_mapk_gene2psrank'}},{0,0},[0,0],{[212:226],[149:163]},[0,0],[0,0],[2,2],[10,10],[3,3],[5,5],{[30:5:70],[30:5:70]});

%test sur sous-ensemble de probesets ayant ou non
%net_recovernts(27,{[92,93]},{0},{[149:163]},2,0,5,10,{[30:5:70]});
%net_recovernts(8,{[185,186]},{0},{[212:226]},2,0,5,10,{[30:5:70]});
%net_recovernts(27,{[93]},{0},{[149:163]},2,0,5,10,{[30:5:70]});

%net_recovernts(91,1,0,[4:18],2,0,5,10,[30:5:70]);
%net_recovernts(91,1,0,[4:18],2,0,5,10,[30:5:70]);
%net_recovernts(94,1,0,[1:3],2,0,5,3,70);
%net_recovernts(93,1,0,[1:3],2,0,5,3,[30:5:70]);

%net_recovernts(2,[24,26],0,[64:78],2,0,5,10,[30:5:55],64);
%net_recovernts(3,[11,13],0,[70:85],2,0,5,10,[30:5:45]);
%net_recovernts(5,[11,13],0,[107:122],2,0,5,10,[30:5:55]);
%net_recovernts(6,[11,13],0,[48:62],2,0,5,10,[30:5:60]);
%net_recovernts(8,[11,13],0,[212:226],2,0,5,10,[30:5:45]);
%net_recovernts(27,[11,13],0,[149:163],2,0,5,10,[30:5:50],164);



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


function net_recovernts(ChipRank,PsList,UsedPsNb,NetRank,Type,CorrLimit,GeneNbLimit,NetNbLimit,Densities,TestNetRank)
global K

if length(ChipRank)>1
    h=errordlg('net_recovernts needs at most one chip');
    waitfor(h)
    error('process canceled')
end
NetNb=length(NetRank);

%recover chip information
ChipPos=find(K.chip.rank==ChipRank);
PsNb=K.chip.probesetNb(ChipPos);

%recover dir for TensorNet
%TensorDir=fullfile(K.dir.net,sprintf('m%u',ChipRank),'tsn','result');
TensorDir=fullfile(K.dir.net,sprintf('m%u',ChipRank),'tsn');
cd(TensorDir)

%construct Prefix used to load ad to save data
if Type==2
    Prefix=sprintf('m%u_n%uton%u_c%u_TSN_diff',ChipRank,NetRank(1),NetRank(end),CorrLimit);
else
    Prefix=sprintf('m%u_n%uton%u_c%u_TSN',ChipRank,NetRank(1),NetRank(end),CorrLimit);
end

%recover data
for ListL=1:length(PsList)
    if iscell(PsList)
        CurrList=PsList{ListL};
        cd(K.dir.list)
        load([CurrList,'.mat'])
        %load GenePos (Ps corresponding to genes in each chip)
        % and Gene (list of gene name)
        CurrIndex=GenePos{ChipRank};
        %construct Index and GeneList for the current chip
        Index1=[];
        GeneList={};
        for GeneL=1:length(CurrIndex)
            Index1=[Index1;CurrIndex{GeneL}];
            for PsL=1:length(CurrIndex{GeneL})
                GeneList{end+1,1}=Gene{GeneL};
            end
        end
        [Index1,SortIndex]=unique(Index1);
        Index={};
        Index{1}=Index1;
        GeneList=GeneList(SortIndex);
    else
        cd(fullfile(K.dir.net,sprintf('m%u',ChipRank),'list'));
        %recover ps lists
        ListRank=PsList(ListL);
        Index=load_data(sprintf('m%u_pslist%u.u32',ChipRank,ListRank),'./',0,0,'uint32','ieee-le');
        %restrict list if UsedPsNb>0
        if UsedPsNb>0&UsedPsNb<length(Index)
            Index=Index(1:UsedPsNb);
        end
    end
    %recover used densities, used ps, and networks where clusters are found
    UsedDensities=[];
    UsedNet=zeros(1,NetNb);
    for DensL=1:length(Densities)
        if iscell(PsList)
            ResultFile=sprintf('%s_%s_g%u_n%u_d%u.txt',Prefix,CurrList,GeneNbLimit,NetNbLimit,Densities(DensL));
        else
            ResultFile=sprintf('%s_%u_g%u_n%u_d%u.txt',Prefix,ListRank,GeneNbLimit,NetNbLimit,Densities(DensL));
        end
        cd(TensorDir)
        if exist(ResultFile)
            [CurrClu,temp,temp,temp,CurrNetList]=textread(ResultFile,'%s%u%.2f%u%s','delimiter','\t');
            UsedDensities=[UsedDensities,DensL];
            %recover probe sets shared between several clusters
            SharedPs=sparse(zeros(length(CurrClu),length(Index)));
            for CluL=1:length(CurrClu)
                CurrClu{CluL}=eval(CurrClu{CluL});
                SharedPs(CluL,CurrClu{CluL})=1;
            end
            %             SharedFreq=sum(SharedPs);
            %             Pos=find(SharedFreq);
            %             figure
            %             set(gcf,'color',[1,1,1])
            %             plot(sort(SharedFreq(Pos)))
            %             set(gca,'box','on')
            %             title('number  of cluster containing a given probe set')
            %             ylabel('number of cluster')
            %             xlabel('ordered probe sets')
            %             FreqLimit=inputdlg('clear ps above freq limit','',1,{'0'});
            %             FreqLimit=str2num(FreqLimit{1});
            %             if FreqLimit>0
            %                 ClearPs=SharedFreq(Pos(find(SharedFreq(Pos)>=FreqLimit)));
            %                 ClearPs=Index(ClearPs);
            %             else
            %                 ClearPs=[];
            %             end
            ClearPs=[];


            CluRank=0;
            NtsClusters{DensL}={};
            for CluL=1:length(CurrClu)               
                PsCurrClu=Index(CurrClu{CluL});
                if ~isempty(ClearPs)
                    PsCurrClu=setdiff(PsCurrClu,ClearPs);
                end
                %                 %load CVM to verify homogeneity of cluster
                %                 [C,A]=load_cvm(ChipRank,NetRank(end),PsCurrClu,PsCurrClu,1,1,0);
                %                 C=C-A{1};
                %                 for i=1:length(C)
                %                     C(i,i)=0;
                %                 end
                %                 [temp,SortIndex]=sort(mean(C));
                %                 PsCurrClu=PsCurrClu(SortIndex);
                %                 C=C(SortIndex,:);
                %                 C=C(:,SortIndex);
                %
                %                 h=figure;
                %                 set(h,'name',sprintf('dens %u cluster %u upon %u',Densities(DensL),CluL,length(CurrClu)))
                %                 image(C)
                %                 set(gca,'ytick',[1:length(PsCurrClu)])
                %                 set(gca,'yticklabel',PsCurrClu)
                %                 set(gca,'tickdir','out')
                %                 SplitPos=inputdlg(sprintf('split at position\n(clear if <0)'),'',1,{'0'});
                %                 SplitPos=str2num(SplitPos{1});
                %                 if SplitPos>0
                %                     set(gca,'xtick',SplitPos+0.5)
                %                     set(gca,'xticklabel','S')
                %                     pause(1)
                %                 end
                %                 close(h)
                %                 if SplitPos>=0
                %                     CluRank=CluRank+1;
                %                     NtsClusters{DensL}{CluRank}=sort(PsCurrClu(SplitPos+1:end));
                %                 end
                NtsClusters{DensL}{CluL}=sort(PsCurrClu);
            end
            %recover networks positive for the current cluster
            for NtsL=1:length(CurrNetList)
                CurrNetRanks=regexp([' ',CurrNetList{NtsL}],' \d+(?=:)','match');
                for NetL=1:length(CurrNetRanks)
                    CurrNetRank=str2num(CurrNetRanks{NetL});
                    UsedNet(CurrNetRank)=UsedNet(CurrNetRank)+1;
                end
            end
            
            
            %% DISPLAY CLUSTERS

            %sort clusters to have in first position the one with the largest number of probe
            %sets           
            SharedPs=sparse(zeros(length(NtsClusters{DensL}),length(Index)));
            for CluL=1:length(NtsClusters{DensL})
                SharedPs(CluL,NtsClusters{DensL}{CluL})=1;
            end

%             CluSize=sum(full(SharedPs),2);
%             [CluSize,SortIndex]=sort(CluSize,'descend');
%             SharedPs=SharedPs(SortIndex,:);

            %process each cluster to allow a graphical representation (a PsNbxPsNb matrix showing the
            %relationships between biclusters
            Ps2Ps=zeros(length(Index),length(Index));
            IndexOrder=[];         
            % Probe sets already used
            UsedIndex=[];
            UsedIndexNb=0;
            PsOrder=zeros(length(Index),1);
            PsTick=[0];
            CluNb=size(SharedPs,1);
            for CluL=1:CluNb              
                CurrIndex=CurrClu{CluL};
                CurrNewIndex=setdiff(CurrIndex,UsedIndex);
                UsedIndex=union(UsedIndex,CurrNewIndex);
                %add a new tick if exist new probe sets
                if length(UsedIndex)>PsTick(end)
                    PsTick(end+1)=length(UsedIndex);
                end
                %indicates the position of new probe sets
                PsOrder(CurrNewIndex)=UsedIndexNb+1:UsedIndexNb+length(CurrNewIndex);
                IndexOrder=[IndexOrder,CurrNewIndex];
                %update UsedIndexNb
                UsedIndexNb=UsedIndexNb+length(CurrNewIndex);
                %indicates the redundancy of probe sets in clusters
                Ps2Ps(PsOrder(CurrIndex),PsOrder(CurrIndex))=Ps2Ps(PsOrder(CurrIndex),PsOrder(CurrIndex))+1;
            end

            %remove empty ps
            KeepPs=find(sum(Ps2Ps)>0);
            Ps2Ps=Ps2Ps(KeepPs,:);
            Ps2Ps=Ps2Ps(:,KeepPs);

            %reorder Ps
%             for PsL=1:length(PsTick)-1
%                 Interval=PsTick(PsL)+1:PsTick(PsL+1);
%                 [temp SortOrder]=sort(sum(Ps2Ps(Interval,:),2));
%                 Ps2Ps(Interval,:)=Ps2Ps(Interval(SortOrder),:);
%                 Ps2Ps(:,Interval)=Ps2Ps(:,Interval(SortOrder));
%                 IndexOrder(Interval)=IndexOrder(Interval(SortOrder));
%             end
%             %refine the order
%             UsedIndex=[];
%             MemPsL1=2;
%             %find the next Ps interval
%             for PsL1=MemPsL1:min(length(PsTick)-1,length(PsTick)-1)
%                 PsInterval=PsTick(PsL1)+1:PsTick(PsL1+1);
%                 %process the current cluster only if it has not probe sets used in
%                 %previous clusters
%                 if mean(mean(Ps2Ps(PsInterval,PsInterval)))>1                    
%                     %find the same probe sets in earlier clusters
%                     for PsL2=1:PsL1-1
%                         CurrPsInterval=PsTick(PsL2)+1:PsTick(PsL2+1);
%                         SamePsPos=find(min(Ps2Ps(CurrPsInterval,PsInterval),[],2)>0);
%                         %keep new ps
%                         SamePs=setdiff(CurrPsInterval(SamePsPos),UsedIndex);
%                         if ~isempty(SamePs)
%                             %reorder ps in the current cluster in last position
%                             FirstPs=setdiff(setdiff(CurrPsInterval,UsedIndex),SamePs);
%                             Ps2Ps(setdiff(CurrPsInterval,UsedIndex),:)=[Ps2Ps(FirstPs,:);Ps2Ps(SamePs,:)];
%                             Ps2Ps(:,setdiff(CurrPsInterval,UsedIndex))=[Ps2Ps(:,FirstPs,:),Ps2Ps(:,SamePs)];
%                             IndexOrder(setdiff(CurrPsInterval,UsedIndex))=[IndexOrder(FirstPs),IndexOrder(SamePs)];
%                             %update UsedIndex
%                             UsedIndex=union(UsedIndex,CurrPsInterval(length(FirstPs)+1:end));
%                         end
%                     end
% 
%                     break
%                 end
%             end

%             figure
%             set(gcf,'color',[1,1,1])
%             Val=Ps2Ps;
%             Val(find(Val>5))=5;
%             %Colors=colors(colormap,max(unique(Ps2Ps)));
%             Colors=colors(colormap,6);
%             colormap(Colors);
%             image(Val+1)
%             if iscell(PsList)
%                 title(sprintf('m%un%uton%u %u/%u ps in %u biclusters, dens=%u, list=%s',ChipRank,NetRank(1),NetRank(end),length(KeepPs),length(Index),size(SharedPs,1),Densities(DensL),CurrList));                
%             else
%                 title(sprintf('m%un%uton%u %u/%u ps in %u biclusters, dens=%u, list=%u',ChipRank,NetRank(1),NetRank(end),length(KeepPs),length(Index),size(SharedPs,1),Densities(DensL),ListRank));
%             end
%             set(gca,'xtick',PsTick)
%             set(gca,'ytick',PsTick)
%             set(gca,'tickdir','out')
%             set(gca,'xticklabel','')
%             set(gca,'yticklabel','')
%             xlabel('probe sets')
%             ylabel('probe sets')
            
            SortedPs{DensL}=Index(IndexOrder);
            if TestNetRank>0
                %load CVM to verify homogeneity of cluster
                [C,A]=load_cvm(ChipRank,TestNetRank,SortedPs{DensL},SortedPs{DensL},1,1,0);                C=triu(C)+tril(A{1});
                figure;
                image(C)
                set(gcf,'color',[1,1,1])
                if iscell(PsList)
                    title(sprintf('m%un%u %u/%u ps in %u biclusters, dens=%u, list=%s',ChipRank,TestNetRank,length(KeepPs),length(Index),size(SharedPs,1),Densities(DensL),CurrList));
                else
                    title(sprintf('m%un%u %u/%u ps in %u biclusters, dens=%u, list=%u',ChipRank,TestNetRank,length(KeepPs),length(Index),size(SharedPs,1),Densities(DensL),ListRank));
                end                
                set(gca,'xtick',PsTick)
                set(gca,'ytick',PsTick)
                set(gca,'tickdir','out')
                set(gca,'xticklabel','')
                set(gca,'yticklabel','')
            end                                   
        end
    end
    if ~isempty(UsedDensities)
        cd(TensorDir)
        if iscell(PsList)
            eval(sprintf('save %s_%s NtsClusters Densities UsedDensities UsedNet GeneList SortedPs',Prefix,CurrList))
        else
            eval(sprintf('save %s_l%u NtsClusters Densities UsedDensities UsedNet SortedPs',Prefix,ListRank))
        end
    end
end
