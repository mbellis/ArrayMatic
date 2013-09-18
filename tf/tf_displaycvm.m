%========================%
% FUNCTION TF_DISPLAYCVM %
%========================%

% TF_DISPLAYCVM displays for a given cluster of genes the CVM of genes and all transcription
% factors found in the promotor region of these genes

%INPUT PARAMETERS
% 1    ChipRank: chip rank
% 2     NtsFile: file containing NetsTensor results
% 3  TargetFile: file describing the merging process used to construct the chip (with rank
%                equal to ChipRank) from the first chip indicated in TargetFile
% 4     FromTsn: the first NtsFile used if the process is repeated across several chip with
%                successive mergin schemes ( human + mouse and then + rat for exemple)
% 5    InfoFlag: recover gene names and gene ids if equal to one

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

%tf_displaycvm(3,86,'m3_tsn_from_m93n1ton3_corr60',40,'m93n1ton3',[1])
%   tf_displaycvm([2,3,6,5,8,27],[80,86,63,123,228,164],{'m2_tsn_from_m93n1ton3_corr60',...
%   'm3_tsn_from_m93n1ton3_corr60','m6_tsn_from_m93n1ton3_corr60','m5_tsn_from_m93n1ton3_corr60',...
%   'm8_tsn_from_m93n1ton3_corr60','m27_tsn_from_m93n1ton3_corr60',},[70,70,70,70,70,70,70],...
%   {'m93n1ton3','m93n1ton3','m93n1ton3','m93n1ton3','m93n1ton3','m93n1ton3'},[1])

%   tf_displaycvm([3,6,5,8,27],[86,63,123,228,164],{'m3_tsn_from_m93n1ton3_corr60',...
%   'm6_tsn_from_m93n1ton3_corr60','m5_tsn_from_m93n1ton3_corr60',...
%   'm8_tsn_from_m93n1ton3_corr60','m27_tsn_from_m93n1ton3_corr60',},[70,70,70,70,70,70],...
%   {'m93n1ton3','m93n1ton3','m93n1ton3','m93n1ton3','m93n1ton3'},[1])

%   tf_displaycvm([5],[123],{'m5_tsn_from_m93n1ton3_corr60'},[70],{'m93n1ton3'},[1])

%   tf_displaycvm([2,3,5,8,27],[80,86,123,228,164],{'m2_tsn_from_m93n1ton3_corr60',...
%   'm3_tsn_from_m93n1ton3_corr60','m5_tsn_from_m93n1ton3_corr60',...
%   'm8_tsn_from_m93n1ton3_corr60','m27_tsn_from_m93n1ton3_corr60',},...
%   [30:5:70],...
%   ,{0,0,0,0,0},0)




function tf_displaycvm(ChipRank,NetRank,NtsFile,DensVal,CluRank,DisplayFlag)
global K

ChipNb=length(ChipRank);

%load tf definition and information
cd(K.dir.common)
load tflist

%% association between TF and clustered genes + display CVM

for DensL=1:length(DensVal)
    for ChipL=1:ChipNb
        %load chip information
        ChipPos=find(K.chip.rank==ChipRank(ChipL));
        Species=K.chip.species{ChipPos};


        %load cluster information
        TensorDir=fullfile(K.dir.net,sprintf('m%u',ChipRank(ChipL)),'tsn','result');
        cd(TensorDir)
        if exist(sprintf('%s.mat',NtsFile{ChipL}))
            eval(sprintf('load %s',NtsFile{ChipL}))
        else
            h=errordlg(sprintf('%s does \nnot exist',NtsFile{ChipL}));
            waitfor(h)
            error('process canceled');
        end


        DensPos=find(unique(Densities(UsedDensities))==DensVal(DensL));

        cd(K.dir.chip)
        ModifFlag=0;
        if exist(sprintf('m%u_ps2net.mat',ChipRank(ChipL)))
            load(sprintf('m%u_ps2net.mat',ChipRank(ChipL)))
            ModifFlag=1;
        end

        %load position of TF in genes of the current species
        if CluRank{ChipL}(1)==0
            CluRanks=1:length(NtsGeneNames{DensPos});
        else
            CluRanks=CluRank{ChipL};
        end
        for CluL=1:length(CluRanks)
            %recover EnsId
            CluGeneName{ChipL}=NtsGeneNames{DensPos}{CluRanks(CluL)};
            CluEnsId=NtsGeneIds{DensPos}{CluRanks(CluL)};
            EnsPos=strmatch('ENS',CluEnsId);
            ClearPos=setdiff(1:length(CluEnsId),EnsPos);
            CluEnsId(ClearPos)='';
            if ~isempty(CluEnsId)
                Prefix=regexp(CluEnsId{1},'[A-Z]+','match');
                CluEnsId=strrep(CluEnsId,Prefix,'');
                CluEnsId=str2num(cell2mat(CluEnsId));
                %number of each TF (WMP) in each gene
                TfNb=zeros(length(CluEnsId),length(Tf.name));
                cd(K.dir.wpm)
                BlocRank=1;
                UTfRank=[];
                CurrTfStart{ChipL}=[];
                CurrTfEnd{ChipL}=[];
                CurrTfStrand{ChipL}=[];
                CurrTfRank{ChipL}=[];
                CurrTfEnsId{ChipL}=[];
                CurrEnsId{ChipL}=CluEnsId;
                while exist(sprintf('%s_tfpos_%u.mat',Species,BlocRank),'file')
                    %fill TF information  for TF found in the genes of cluster
                    % load :
                    % EnsId : list of ensembl Id (numeric part) and for each position
                    % targeted by a TF
                    % TfRank refers to Tf order
                    % TfStart and TfEnd is the position of the TF binding site in the promoter region
                    % TfStrand is the strandedness of the TF binding site
                    load(sprintf('%s_tfpos_%u.mat',Species,BlocRank))
                    %find position of cluster Ensembl Id in EnsId (all genes with TF)
                    for EnsL=1:length(CluEnsId)
                        Pos=find(EnsId==CluEnsId(EnsL));
                        UTfRank=unique([UTfRank;TfRank(Pos)]);
                        CurrTfStat=histc(TfRank(Pos),[1:length(Tf.name)]);
                        TfNb(EnsL,:)=TfNb(EnsL,:)+CurrTfStat';
                        CurrTfEnsId{ChipL}=[CurrTfEnsId{ChipL};EnsId(Pos)];
                        CurrTfRank{ChipL}=[CurrTfRank{ChipL};TfRank(Pos)];
                        CurrTfStart{ChipL}=[CurrTfStart{ChipL};TfStart(Pos)];
                        CurrTfEnd{ChipL}=[CurrTfEnd{ChipL};TfEnd(Pos)];
                        CurrTfStrand{ChipL}=[CurrTfStrand{ChipL};TfStrand(Pos)];
                    end
                    BlocRank=BlocRank+1;
                end

                %presence of a particular Tf
                IsTf=zeros(size(TfNb));
                IsTf(find(TfNb))=1;
                %select Tf that are present in at elast 75% of genes belonging to the current cluster
                Pos=find(sum(IsTf)>=0.75*size(IsTf,1));
                %calculate a weigth for the presence
                TfRatio=zeros(1,size(IsTf,2));
                TfRatio(Pos)=sum(IsTf(:,Pos))./sum(TfNb(:,Pos));

                if DisplayFlag
                    h=figure;
                    set(gcf,'color',[1,1,1])
                    set(h,'name',sprintf('ordered stat on TF of cluster %u of m%u',CluRanks(CluL),ChipRank(ChipL)))
                    [temp,SortIndex]=sort(sum(IsTf),'descend');
                    [temp,SortIndex1]=sort(TfRatio(SortIndex),'descend');
                    SortIndex=SortIndex(SortIndex1);
                    subplot(5,1,1)
                    image(TfNb(:,SortIndex))
                    title(sprintf('Total nb of TF in each probe set of cluster %u of m%u',CluRanks(CluL),ChipRank(ChipL)))
                    subplot(5,1,2)
                    plot(sum(TfNb(:,SortIndex)))
                    set(gca,'xlim',[1,length(SortIndex)])
                    title('Total number of TF')
                    subplot(5,1,3)
                    plot(mean(TfNb(:,SortIndex)),'c')
                    set(gca,'xlim',[1,length(SortIndex)])
                    title('Mean number of TF')
                    subplot(5,1,4)
                    plot(sum(IsTf(:,SortIndex)),'m')
                    set(gca,'xlim',[1,length(SortIndex)])
                    title('Presence of TF')
                    subplot(5,1,5)
                    plot(TfRatio(SortIndex),'r')
                    set(gca,'xlim',[1,length(SortIndex)])
                    title('Weight of TF presence')
                end


                %display CVM

                %load probe set rank of the current cluster   to load CVM
                CluPsRank=NtsClusters{DensPos}{CluRanks(CluL)};
                if ModifFlag
                    CluPsRank=Ps2Net(CluPsRank);
                end
                CluSize=length(CluPsRank);
                [CClu]=load_cvm(ChipRank(ChipL),NetRank(ChipL),CluPsRank,CluPsRank,1,1,1);
                %load TF CVM
                TfPsRank=[];
                TfGeneName={};
                TfGeneRank=[];
                TfWeight=[];
                for TfL=1:length(UTfRank)
                    if ~isempty(Tf.psRank{ChipRank(ChipL)}{Tf.toUnique(UTfRank(TfL))})
                        TfPsRank=[TfPsRank;Tf.psRank{ChipRank(ChipL)}{Tf.toUnique(UTfRank(TfL))}];
                        TfGeneName=[TfGeneName;repmat(Tf.gene(UTfRank(TfL)),length(Tf.psRank{ChipRank(ChipL)}{Tf.toUnique(UTfRank(TfL))}),1)];
                        TfGeneRank=[TfGeneRank;repmat(Tf.toUnique(UTfRank(TfL)),length(Tf.psRank{ChipRank(ChipL)}{Tf.toUnique(UTfRank(TfL))}),1)];
                        TfWeight=[TfWeight;repmat(TfRatio(UTfRank(TfL)),length(Tf.psRank{ChipRank(ChipL)}{Tf.toUnique(UTfRank(TfL))}),1)];
                    end
                end
                TfWeight=TfWeight'*100;
                if ModifFlag
                    TfPsRank=Ps2Net(TfPsRank);
                end
                [CTf]=load_cvm(ChipRank(ChipL),NetRank(ChipL),TfPsRank,TfPsRank,1,1,1);

                %load Cluster x TF CVM
                [CCluxTf,ACluxTf]=load_cvm(ChipRank(ChipL),NetRank(ChipL),CluPsRank,TfPsRank,1,1,0);
                CCluxTf=CCluxTf-ACluxTf{1};
                CCluxTf(find(CCluxTf<0))=0;
                CTfxClu=CCluxTf';
                ATfxClu=ACluxTf{1}';
                ATfxClu=ATfxClu-CTfxClu;
                ATfxClu(find(ATfxClu<0))=0;

                %calculate mean correlation with the current cluster for each TF
                [MeanCorr,SortIndex]=sort(mean(CCluxTf),'descend');
                CorrCTf=CTf(:,SortIndex);
                CorrCTf=CorrCTf(SortIndex,:);

                STfGeneName{ChipL}{CluL}{1}=TfGeneName(SortIndex);
                STfGeneRank{ChipL}{CluL}{1}=TfGeneRank(SortIndex);
                SMeanCorr{ChipL}{CluL}{1}=MeanCorr;

                %combine CVM
                if DisplayFlag
                    h=figure;
                    set(gcf,'color',[1,1,1])
                    set(h,'name',sprintf('CVM cluster %u (dens %u) of m%u (TF ordered on mean(CORR))',CluRanks(CluL),DensVal(ChipL),ChipRank(ChipL)))
                    image([CClu,CCluxTf(:,SortIndex);[zeros(3,CluSize),repmat(MeanCorr,3,1)];...
                        [zeros(3,CluSize),repmat(TfWeight(SortIndex),3,1)];[ATfxClu(SortIndex,:),CorrCTf]])
                    set(gca,'xlim',[0,150])
                    set(gca,'ylim',[0,150])
                    set(gca,'ytick',[0:10:150])
                    set(gca,'yticklabel',flipud(STfGeneRank{ChipL}{CluL}{1}(1:16)))
                    set(gca,'xtick',[0:10:150])
                    set(gca,'xticklabel',STfGeneRank{ChipL}{CluL}{1}(17:32))
                    rotateXTickLabel
                    title(sprintf('CVM cluster %u (dens %u) of m%u (TF ordered on mean(CORR))',CluRanks(CluL),DensVal(ChipL),ChipRank(ChipL)))
                end


                %calculate mean correlation with the current cluster for each TF
                [MeanAnti,SortIndex]=sort(mean(ATfxClu'),'descend');
                AntiCTf=CTf(:,SortIndex);
                AntiCTf=AntiCTf(SortIndex,:);

                STfGeneName{ChipL}{CluL}{2}=TfGeneName(SortIndex);
                STfGeneRank{ChipL}{CluL}{2}=TfGeneRank(SortIndex);
                SMeanCorr{ChipL}{CluL}{2}=MeanAnti;

                %combine CVM
                if DisplayFlag
                    h=figure;
                    set(gcf,'color',[1,1,1])
                    set(h,'name',sprintf('CVM cluster %u (dens %u) of m%u (TF ordered on mean(ANTI))',CluRanks(CluL),DensVal(ChipL),ChipRank(ChipL)))
                    %image([CClu,CCluxTf(:,SortIndex);[zeros(3,CluSize),repmat(MeanCorr,3,1)];...
                    %    [zeros(3,CluSize),repmat(TfWeight(SortIndex),3,1)];[ATfxClu(SortIndex,:),AntiCTf]])
                    image([CClu,ATfxClu(SortIndex,:)';[zeros(3,CluSize),repmat(MeanAnti,3,1)];...
                        [zeros(3,CluSize),repmat(TfWeight(SortIndex),3,1)];[CCluxTf(:,SortIndex)',AntiCTf']])
                    set(gca,'xlim',[0,150])
                    set(gca,'ylim',[0,150])
                    set(gca,'ytick',[0:10:150])
                    set(gca,'yticklabel',flipud(STfGeneRank{ChipL}{CluL}{2}(1:16)))
                    set(gca,'xtick',[0:10:150])
                    set(gca,'xticklabel',STfGeneRank{ChipL}{CluL}{2}(17:32))
                    rotateXTickLabel
                    title(sprintf('CVM cluster %u (dens %u) of m%u (TF ordered on mean(ANTI))',CluRanks(CluL),DensVal(ChipL),ChipRank(ChipL)))
                end
            else
                STfGeneName{ChipL}{CluL}{1}={};
                STfGeneRank{ChipL}{CluL}{1}=[];
                SMeanCorr{ChipL}{CluL}{1}=[];
                STfGeneName{ChipL}{CluL}{2}={};
                STfGeneRank{ChipL}{CluL}{2}=[];
                SMeanCorr{ChipL}{CluL}{2}=[];
            end
        end % CluL
    end %ChipL


    %% Find common TF
    %PlotNb=ChipNb*(ChipNb-1)/2;
    %RowNb=ceil(sqrt(PlotNb));
    %ColNb=round(PlotNb/RowNb);
    RowNb=ChipNb;
    ColNb=ChipNb;
    PlotNb=ChipNb*ChipNb;
    PlotPos=reshape(1:PlotNb,ChipNb,ChipNb);
    PlotPos=PlotPos';
    Colors=colors(colormap,ChipNb);
    % TypeL = CORR and ANTI CVM

    %keep only the firs occurence of each TF
    for ChipL=1:ChipNb
        for TypeL=1:2
            for CluL=1:length(CluRanks)
                if ~isempty(STfGeneRank{ChipL}{CluL}{TypeL})
                    FTfGeneRank{ChipL}{CluL}{TypeL}=KEEP_FIRST(STfGeneRank{ChipL}{CluL}{TypeL});
                else
                    FTfGeneRank{ChipL}{CluL}{TypeL}=[];
                end
            end
        end
    end

    TOPNB=10;
    for TypeL=1:2
        for CluL=1:length(CluRanks)
            if ~isempty(FTfGeneRank{ChipL}{CluL}{TypeL})
                %calculate frequencies
                % number of chip which have at least one probe set targeting a TF
                TfNb=zeros(1,length(Tf.unique));
                for ChipL=1:ChipNb
                    TfNb(FTfGeneRank{ChipL}{CluL}{TypeL})=TfNb(FTfGeneRank{ChipL}{CluL}{TypeL})+1;
                end
                % list of TF that are in the TOPNB first positions in the ordered TF lists of all
                % chips (=> the TOPNB first TF)
                TopTfPos=[];
                % number of chip which have one of the TOPNB first TF in its own top list
                TopTfNb=zeros(1,length(Tf.unique));
                % for each of the TOPNB TF in the first chip top list, indicates the difference of
                % rank between the positions of the TF in the second and the first list
                TfDistance{CluL}{TypeL}=zeros(1,length(Tf.unique));
                %mean (on all chips) of mean correlations
                TfMeanCorr{CluL}{TypeL}=zeros(1,length(Tf.unique));
                if DisplayFlag
                    h=figure;
                    set(gcf,'color',[1,1,1])
                    if TypeL==1
                        set(h,'name',sprintf('Chip correspondance for Clu %u (CORR)',CluRanks(CluL)))
                    else
                        set(h,'name',sprintf('Chip correspondance for Clu %u (ANTI)',CluRanks(CluL)))
                    end
                end
                for ChipL1=1:ChipNb
                    for ChipL2=1:ChipNb
                        if ChipL1~=ChipL2
                            %keep only common Tf to the two current chips
                            [temp,Index1,Index2]=intersect(FTfGeneRank{ChipL1}{CluL}{TypeL},FTfGeneRank{ChipL2}{CluL}{TypeL});
                            %keep TF in the original order
                            Index1=sort(Index1);
                            FCommonTf{1}=FTfGeneRank{ChipL1}{CluL}{TypeL}(Index1);
                            TopNb=min(TOPNB,length(FCommonTf{1}));
                            Index2=sort(Index2);
                            FCommonTf{2}=FTfGeneRank{ChipL2}{CluL}{TypeL}(Index2);
                            %recover TOPNB TF of the first chip list
                            TopTfPos=[TopTfPos,FCommonTf{1}(1:TopNb)'];
                            TopTfNb(FCommonTf{1}(1:TopNb))=TopTfNb(FCommonTf{1}(1:TopNb))+1;
                            if DisplayFlag
                                if ChipL2~=ChipL1
                                    subplot(RowNb,ColNb,PlotPos(ChipL1,ChipL2))
                                    hold on
                                end
                            end
                            TfPos=zeros(TopNb,1);
                            for TfL=1:TopNb
                                TfPos(TfL)=find(FCommonTf{2}==FCommonTf{1}(TfL));
                                TfDistance{CluL}{TypeL}(FCommonTf{1}(TfL))=TfDistance{CluL}{TypeL}(FCommonTf{1}(TfL))+abs(TfPos(TfL)-TfL);
                                TfPos(TfL)=min(50,TfPos(TfL));
                                if DisplayFlag
                                    if ChipL2~=ChipL1
                                        plot(TfL,TfPos(TfL),'o','color',Colors(TfNb(FCommonTf{1}(TfL)),:))
                                    end
                                end
                            end
                            if DisplayFlag
                                set(gca,'xtick',[1:TopNb])
                                set(gca,'xticklabel',FCommonTf{1}(1:TopNb))
                                set(gca,'ytick',sort(unique(TfPos)))
                                set(gca,'yticklabel',sort(unique(TfPos)))
                                set(gca,'ylim',[0,50])
                                rotateXTickLabel
                                set(gca,'box','on')
                                title(sprintf('m%u vs m%u',ChipRank(ChipL2),ChipRank(ChipL1)))
                            end
                        end
                    end
                end

                TopTfPos=unique(TopTfPos);
                for ChipL1=1:ChipNb
                    for TfL=1:length(TopTfPos)
                        STfPos=find(STfGeneRank{ChipL}{CluL}{TypeL}==TopTfPos(TfL));
                        if ~isempty(STfPos)
                            TfMeanCorr{CluL}{TypeL}(TopTfPos(TfL))=TfMeanCorr{CluL}{TypeL}(TopTfPos(TfL))+SMeanCorr{ChipL}{CluL}{TypeL}(STfPos(1));
                        end
                    end
                end

                TfMeanCorr{CluL}{TypeL}=round(TfMeanCorr{CluL}{TypeL}./TfNb);
                TfNb=TfNb(TopTfPos);
                TopTfNb=TopTfNb(TopTfPos);
                TfDistance{CluL}{TypeL}=TfDistance{CluL}{TypeL}(TopTfPos);
                TfMeanCorr{CluL}{TypeL}=TfMeanCorr{CluL}{TypeL}(TopTfPos);
                SelPos=find(TopTfNb>(ChipNb-1)*ChipNb/4);
                TopTfPos=TopTfPos(SelPos);
                TfNb=TfNb(SelPos);
                TopTfNb=TopTfNb(SelPos);
                TfDistance{CluL}{TypeL}=TfDistance{CluL}{TypeL}(SelPos);
                TfMeanCorr{CluL}{TypeL}=TfMeanCorr{CluL}{TypeL}(SelPos);
                [temp,SortIndex]=sort(TfDistance{CluL}{TypeL}./TopTfNb);
                %[temp,SortIndex]=sort(TopTfNb,'descend');
                TopTfPos=TopTfPos(SortIndex);
                TfNb=TfNb(SortIndex);
                TopTfNb=TopTfNb(SortIndex);
                TfDistance{CluL}{TypeL}=TfDistance{CluL}{TypeL}(SortIndex);
                TfMeanCorr{CluL}{TypeL}=TfMeanCorr{CluL}{TypeL}(SortIndex);
                TfStat{DensL}{CluL}{TypeL}=[TfNb;TopTfNb;TfDistance{CluL}{TypeL};round(TfDistance{CluL}{TypeL}./TopTfNb);TfMeanCorr{CluL}{TypeL};TopTfPos];
            else
                TfStat{DensL}{CluL}{TypeL}=[];
            end
        end
    end


    %% DISPLAY TF POSITIONS
    if DisplayFlag
        %recover data
        CurrClu=1;
        SelTfRank{1}=TfStat{DensL}{CurrClu}{1}(6,[1:6]);
        SelTfRank{2}=TfStat{DensL}{CurrClu}{2}(6,[1:6]);
        SelTfName{1}=Tf.unique(SelTfRank{1});
        SelTfName{2}=Tf.unique(SelTfRank{2});
        Species={'human','mouse'};
        SpeciesNb=length(Species);
        SpeciesPos=[1,3];


        % GeneName and GeneId lists are out of register => reorder GeneId list
        cd(K.dir.chip)
        GeneName={};
        GeneId={};
        for SpeciesL=1:SpeciesNb
            [temp,CurrGeneId,CurrGeneName,temp,temp]=textread(sprintf('m%u_gene.txt',ChipRank(SpeciesPos(SpeciesL))),'%s%s%s%u%u','delimiter','\t');
            %reorder EnsId
            ClearPos=setdiff(1:length(CurrGeneId),strmatch('ENS',CurrGeneId));
            CurrGeneId(ClearPos)=[];
            CurrGeneName(ClearPos)=[];
            Prefix=regexp(CurrGeneId{1},'[A-Z]+','match');
            CurrGeneId=strrep(CurrGeneId,Prefix,'');
            CurrGeneId=str2num(cell2mat(CurrGeneId));

            GeneName{SpeciesL}=CluGeneName{SpeciesPos(SpeciesL)};
            GeneId{SpeciesL}=CurrEnsId{SpeciesPos(SpeciesL)};

            Order=zeros(length(GeneId{SpeciesL}),1);
            for IdL=1:length(GeneId{SpeciesL})
                IdPos=find(CurrGeneId==GeneId{SpeciesL}(IdL));
                GenePos=strmatch(CurrGeneName{IdPos(1)},GeneName{SpeciesL},'exact');
                Order(IdL)=GenePos;
            end
            GeneId{SpeciesL}=GeneId{SpeciesL}(Order);
        end

        DeltaGene=30;
        DeltaSpecies=20;
        SpeciesColor=colors(colormap,SpeciesNb);
        for TypeL=1:2
            TfName=Tf.unique(SelTfRank{TypeL});
            TfColors=colors(colormap,length(SelTfRank{TypeL}));
            h=figure;
            set(gcf,'color',[1,1,1])
            if TypeL==1
                set(h,'name','INDUCTIVE TF')
            else
                set(h,'name','REPRESSIVE TF')
            end
            LinePos=0;
            YTick=[];
            for GeneL=1:length(GeneName{1})
                if GeneL>1
                    LinePos=LinePos+DeltaGene;
                end
                YTick(end+1)=LinePos+DeltaSpecies/2;
                for SpeciesL=1:SpeciesNb
                    if SpeciesL>1
                        LinePos=LinePos+DeltaSpecies;
                    end
                    for TfL=1:length(SelTfRank{TypeL})
                        subplot(2,3,TfL)
                        line([-2100,50],[LinePos,LinePos],'color',SpeciesColor(SpeciesL,:),'linewidth',2)
                        hold on
                        TfPos=find(CurrTfRank{SpeciesPos(SpeciesL)}==SelTfRank{TypeL}(TfL)&CurrTfEnsId{SpeciesPos(SpeciesL)}==GeneId{SpeciesL}(GeneL));
                        for PosL=1:length(TfPos)
                            if CurrTfStrand{SpeciesPos(SpeciesL)}(TfPos(PosL))==1
                                line([-double(CurrTfStart{SpeciesPos(SpeciesL)}(TfPos(PosL))),-double(CurrTfStart{SpeciesPos(SpeciesL)}(TfPos(PosL)))],[LinePos,LinePos+10],'color',TfColors(TfL,:),'linewidth',4);
                            else
                                line([-double(CurrTfStart{SpeciesPos(SpeciesL)}(TfPos(PosL))),-double(CurrTfStart{SpeciesPos(SpeciesL)}(TfPos(PosL)))],[LinePos,LinePos-10],'color',TfColors(TfL,:),'linewidth',4);
                            end
                        end
                        if GeneL==length(GeneName{1})
                            set(gca,'box','on')
                            set(gca,'xlim',[-2100,100])
                            %set(gca,'xticklabel','')
                            set(gca,'ylim',[-20,-20+50*GeneL])
                            set(gca,'ytick',YTick)
                            if mod(TfL,3)==1
                                set(gca,'yticklabel',GeneName{SpeciesL})
                            else
                                set(gca,'yticklabel','')
                            end
                            TfPos=find(TfStat{DensL}{CurrClu}{TypeL}(6,:)==SelTfRank{TypeL}(TfL));
                            title(sprintf('%s (tf%u %u %u %u %u %u)',SelTfName{TypeL}{TfL},TfStat{DensL}{CurrClu}{TypeL}(6,TfPos),...
                                TfStat{DensL}{CurrClu}{TypeL}(1,TfPos),TfStat{DensL}{CurrClu}{TypeL}(2,TfPos),TfStat{DensL}{CurrClu}{TypeL}(3,TfPos),...
                                TfStat{DensL}{CurrClu}{TypeL}(4,TfPos),TfStat{DensL}{CurrClu}{TypeL}(5,TfPos)))
                        end
                    end
                end
            end
        end
    end %DisplayFlag
end %DensL

'stop'
DensL=1;
for TypeL=1:2
    TfMat{TypeL}=zeros(length(TfStat{DensL}),length(Tf.unique));
    for CluL=1:length(TfStat{DensL})
        if ~isempty(TfStat{DensL}{CluL}{TypeL})
            TfMat{TypeL}(CluL,TfStat{DensL}{CluL}{TypeL}(6,:))=[1:length(TfStat{DensL}{CluL}{TypeL}(6,:))];
        end
    end
    for RepL=1
        TfMatBin{TypeL,RepL}=zeros(size(TfMat{TypeL}));
        TfPos=find(TfMat{TypeL}<=3&TfMat{TypeL}>0);
        TfMatBin{TypeL,RepL}(TfPos)=TfMat{TypeL,RepL}(TfPos)*20+10;
    end
end

figure
subplot(2,2,1)
image(TfMat{1}*10)
for RepL=1:3
subplot(2,2,RepL+1)
image(TfMatBin{1,RepL}*30)
end

figure
image(TfMat{2}*10)


function List=KEEP_FIRST(List)

ClearPos=zeros(length(List),1);
for PosL=1:length(List)
    if ClearPos(PosL)==0
        Pos=find(List==List(PosL));
        if length(Pos)>1
            ClearPos(Pos(2:end))=1;
        end
    end
end
List(find(ClearPos))=[];

