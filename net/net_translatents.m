%===========================%
% FUNCTION NET_TRANSLATENTS %
%===========================%

% NET_TRANSLATENTS recover original probesets ranks of NTS clusters found in composit chip
% (merge of different chip model belonging to the same species or to different species)

%INPUT PARAMETERS
% 1    ChipRank: chip rank
% 2     NtsFile: file containing NetsTensor results
% 3  TargetFile: file describing the merging process used to construct the chip (with rank
%                equal to ChipRank) from the first chip indicated in TargetFile
% 4     FromNts: the first NtsFile used if the process is repeated across several chip with
%                successive mergin schemes ( human + mouse and then + rat for exemple)
% 5    InfoFlag: recover gene names and gene ids if equal to one
% 6  PsRankFlag: find correspondance with probe set rank


%net_translatents(94,'m94_n1ton3_c0_TSN_diff_l1','m93n3_m6n63_combinedps_corr60','m94_n1ton3',0,0)
%net_translatents(93,'m93_tsn_from_m94_n1ton3_corr60','m91n3_m92n4_combinedps_corr60','m94_n1ton3',0,0)
%net_translatents(93,'m93_tsn_from_m94_n1ton3_corr60','m92n4_m91n3_combinedps_corr60','m94_n1ton3',0,0)
%net_translatents(91,'m91_tsn_from_m94_n1ton3_corr60','m2_combinedps_corr60','m94_n1ton3',1,0)
%net_translatents(91,'m91_tsn_from_m94_n1ton3_corr60','m3_combinedps_corr60','m94_n1ton3',1,0)
%net_translatents(92,'m92_tsn_from_m94_n1ton3_corr60','m5_combinedps_corr60','m94_n1ton3',1,0)
%net_translatents(92,'m92_tsn_from_m94_n1ton3_corr60','m8_combinedps_corr60','m94_n1ton3',1,0)
%net_translatents(92,'m92_tsn_from_m94_n1ton3_corr60','m27_combinedps_corr60','m94_n1ton3',1,0)

%net_translatents(93,'m93_n1ton3_c0_TSN_diff_l1','m91n3_m92n4_combinedps_corr60','m93n1ton3',0,0)
%net_translatents(93,'m93_n1ton3_c0_TSN_diff_l1','m92n4_m91n3_combinedps_corr60','m93n1ton3',0,0)

%net_translatents(91,'m91_tsn_from_m93n1ton3_corr60','m2n80_m3n86_combinedps_corr60','m93n1ton3',1,0)
% net_translatents(91,'m91_tsn_from_m93n1ton3_corr60','m3n86_m2n80_combinedps_corr60','m93n1ton3',1,0)
% net_translatents(92,'m92_tsn_from_m93n1ton3_corr60','m5n123_m27n164_combinedps_corr60','m93n1ton3',1,0)
% net_translatents(92,'m92_tsn_from_m93n1ton3_corr60','m8_combinedps_corr60','m93n1ton3',1,0)
% net_translatents(92,'m92_tsn_from_m93n1ton3_corr60','m27n164_m5n123_combinedps_corr60','m93n1ton3',1,0)


% net_translatents(27,'m27_tsn_from_m93n1ton3_corr60','m6n63_m27n164_combinedps_corr60','m93n1ton3',1,1)







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


function net_translatents(ChipRank,NtsFile,TargetFile,FromNts,InfoFlag,PsRankFlag)
global K


%load NetTensor file
TensorDir=fullfile(K.dir.net,sprintf('m%u',ChipRank),'tsn','result');
cd(TensorDir)
if exist(sprintf('%s.mat',NtsFile))
    eval(sprintf('load %s',NtsFile))
else
  h=errordlg(sprintf('%s does \nnot exist',NtsFile));
  waitfor(h)
  error('process canceled');  
end


%load file allowing convertion of ps ranks
MemNtsClusters=NtsClusters;
NtsClusters={};
%load merging scheme
cd(K.dir.chip)
eval(sprintf('load %s',TargetFile))

%find chips used in combined file
CurrChipRank=regexp(TargetFile,'(?<=m)\d+(?=n)','match');
if isempty(CurrChipRank)
    CurrChipRank=regexp(TargetFile,'(?<=m)\d+(?=_)','match');
end
ChipRanks=[];
for ChipL=1:length(CurrChipRank)
    ChipRanks(ChipL)=str2num(CurrChipRank{ChipL});
end

%find info on first chip
ChipPos=find(K.chip.rank==ChipRanks(1));
MaxPsNb=K.chip.probesetNb(ChipPos);

%load chip information for the first chip
if InfoFlag
    cd(K.dir.chip)
    [ProbeSetId,GeneId,GeneName,Class,PNb]=textread(sprintf('m%u_gene.txt',ChipRanks(1)),'%s%s%s%u%u','delimiter','\t');
end
%translate rank for the first chip
for DensL=1:length(MemNtsClusters)
    NtsClusters{DensL,1}={};
    if InfoFlag
        NtsGeneNames{DensL,1}={};
        NtsGeneIds{DensL,1}={};
    end
    for NtsL=1:length(MemNtsClusters{DensL})
        CurrPsRank=[];
        for PsL=1:length(MemNtsClusters{DensL}{NtsL})
            if PsRankFlag
                PsPos=find(PsRanks(:,2)==MemNtsClusters{DensL}{NtsL}(PsL));
                if ~isempty(PsPos)
                    if max(PsRank{PsPos})>MaxPsNb
                        h=errordlg('error exist rank>%u in PsL==%u',MaxPsNb,PsL);
                        waitfor(h)
                        error('process canceled')
                    else
                        CurrPsRank=[CurrPsRank;PsRank{PsPos}];
                    end
                end
            else
                if max(PsRank{MemNtsClusters{DensL}{NtsL}(PsL)})>MaxPsNb                  
                    h=errordlg('error exist rank>%u in PsL==%u',MaxPsNb,PsL);
                    waitfor(h)
                    error('process canceled')
                else
                    CurrPsRank=[CurrPsRank;PsRank{MemNtsClusters{DensL}{NtsL}(PsL)}];
                end
            end
        end
        CurrPsRank=sort(CurrPsRank);
        NtsClusters{DensL,1}{NtsL,1}=CurrPsRank;
        if InfoFlag
            MemGeneNames=GeneName(CurrPsRank);
            MemGeneIds=GeneId(CurrPsRank);            
            GeneClearPos=strmatch('',MemGeneNames,'exact');            
            GeneClearPos=union(GeneClearPos,strmatch('-',MemGeneNames,'exact'));
            IdClearPos=strmatch('',MemGeneIds,'exact');
            IdClearPos=union(IdClearPos,strmatch('-',MemGeneIds,'exact'));
            ClearPos=intersect(GeneClearPos,IdClearPos);
            GeneNames=MemGeneNames;
            GeneIds=MemGeneIds;
            GeneNames(ClearPos)=[];
            GeneIds(ClearPos)=[];        
            GeneIdName=cell(size(GeneNames));
            for GeneL=1:length(GeneNames)
                GeneIdName{GeneL}=[GeneIds{GeneL},'@',GeneNames{GeneL}];
            end
            GeneIdName=unique(GeneIdName);
            NtsGeneNames{DensL,1}{NtsL,1}=cell(size(GeneIdName));
            NtsGeneIds{DensL,1}{NtsL,1}=cell(size(GeneIdName));
            for GeneL=1:length(GeneIdName)
                CurGeneId=regexp(GeneIdName{GeneL},'.+(?=@)','match');
                NtsGeneIds{DensL}{NtsL}{GeneL,1}=CurGeneId{1};
                CurrGeneName=regexp(GeneIdName{GeneL},'(?<=@).+','match');
                NtsGeneNames{DensL}{NtsL}{GeneL,1}=CurrGeneName{1};
            end
            NtsGenePos{DensL,1}{NtsL,1}=zeros(length(CurrPsRank),1);
            for PsL=1:length(CurrPsRank)
                GenePos=strmatch(MemGeneIds{PsL},NtsGeneIds{DensL}{NtsL},'exact');
                if isempty(GenePos)
                    GenePos=strmatch(MemGeneNames{PsL},NtsGeneNames{DensL}{NtsL},'exact');
                end
                if~isempty(GenePos)                    
                    NtsGenePos{DensL,1}{NtsL,1}(PsL)=GenePos;
                end
            end            
        end        
    end
end


for ChipL=1
    TensorDir=fullfile(K.dir.net,sprintf('m%u',ChipRanks(ChipL)),'tsn','result');
    CorrVal=regexp(TargetFile,'(?<=_corr)\d+','match');
    CorrVal=str2num(CorrVal{1});
    if ~exist(TensorDir)
        mkdir(TensorDir)
    end
    cd(TensorDir)
    if InfoFlag
        eval(sprintf('save m%u_tsn_from_%s_corr%u NtsClusters NtsGeneNames NtsGeneIds NtsGenePos Densities UsedDensities ',ChipRanks(ChipL),FromNts,CorrVal))
    else
        eval(sprintf('save m%u_tsn_from_%s_corr%u NtsClusters Densities UsedDensities',ChipRanks(ChipL),FromNts,CorrVal))
    end
end



