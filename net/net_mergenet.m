% NET_MERGENET makes a synthesis of several networks constructed from
% different combinations of biological conditions.
% One option is to take the intersection of all the merged networks
% (synthetic values of CORR and ANTI are then for a particular pair of probe sets the mean or the median of all the non zero values, if a
% given percent of networks have non zero values (otherwise the synthetic
% values are null).
%
% INPUT
% ChipRank: chiprank(s)
% NetRanks: list of existing networks to be processed
% MergeType: type of merging : mean, intersection, union, fitting
%            fitting: correct values of one chip to normalise with another chip
%                     needs m%u_to_m%u_corr file
%            merge: merge probe sets that are similar in a given percent of network
%            combine: extract a subset of a chip following a given pattern
% QlimitFlag: indiquates if qlimit has been used in processed networks
% VARARGIN
% if MergeType=='merge'
% 1       MinNetNb: the minimal number of networks which must have significative
%                   values (default = NetNb) for integrate corr and anti intersection network
% 2       NetRanks: ranks of networks usedfor ps assignation (psawn)
% 3   ProbeNbLimit: minimum number of probes targeting a gene
% 4     PvCorrRank: pv(overlap) is calculated for corr limit >[0,40,50,60]. PvCorrRank
%                    indicates the corr limit to be used, by giving its index in
%                    the corr list([0,40,50,60])
% 5 NetFrequencies: list of net frequencies that mus be considered ([1,25,50,75,100],
%                    recommended)
% 6         Suffix: suffix for construct result directory name
% 7        Species: chip species
% 8  MaxFlag: take max(merged values) otherwise mean(merged values)
% if MergeType=='combine'
% 1 NewChipRank
% 2 NewNetRank
% 3 NormFile




%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤%
%                          c) Michel Bellis                                                %
%                          michel.bellis@crbm.cnrs.fr                                      %
%            Affiliation:  CNRS (Centre National de la Recherche Scientifique - France)    %
%  Bioinformatic Project:  ARRAYMATIC => http://code.google.com/p/arraymatic               %
%        Code Repository:  GITHUB => http://github.com/mbellis                             %
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤%

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
%  THIS CODE IS DISTRIBUTED UNDER THE CeCILL LICENSE, WHICH IS COMPATIBLE WITH       %
%  THE GNU GENERAL PUBLIC LICENCE AND IN ACCORDANCE WITH THE EUROPEAN LEGISLATION.   %
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

%net_mergenet(27,[2:22],'union',1)
%net_mergenet(8,[7:21],'union',1)
%net_mergenet(27,[88:108],'mean',0)
%net_mergenet(8,[39:53],'mean',0)
%net_mergenet(8,[126:140],'mean',0)
%net_mergenet(8,[196:210],'mean',0)
%net_mergenet(8,[212:226],'mean',0)
%net_mergenet(27,[119:133],'mean',0)
%net_mergenet(27,[134:148],'mean',0)
%net_mergenet(27,[149:163],'mean',0)
%net_mergenet(8,[212:217,219],'mean',0)
%net_mergenet(8,[220:226],'mean',0)
%net_mergenet(27,[134:148],'mean',0)
%net_mergenet(27,[149:163],'mean',0)
%net_mergenet(8,[212:217],'mean',0)
%net_mergenet(8,[218:223],'mean',0)
%net_mergenet(27,[149:154],'mean',0)
%net_mergenet(27,[155:160],'mean',0)
%net_mergenet([27,8],[164],'fitting',0)
%net_mergenet(5,[89:104],'mean',0)
%net_mergenet(2,[45,46,48:59],'mean',0)
%net_mergenet(3,[54:69],'mean',0)
%net_mergenet(6,[33:47],'mean',0)
%net_mergenet(5,[107:122],'mean',0)
%net_mergenet(2,[64:78],'mean',0)
%net_mergenet(3,[70:85],'mean',0)
%net_mergenet(6,[48:62],'mean',0)
%net_mergenet(6,[33:47],'mean',0)
%net_mergenet(2,80,'combine',0,91,1,60,'')
%net_mergenet(3,86,'combine',0,91,2,60,'m3tom2_normcurve.mat')
%net_mergenet(5,123,'combine',0,92,1,60,'')
%net_mergenet(8,228,'combine',0,92,2,60,'m8tom5_normcurve.mat')
%net_mergenet(27,164,'combine',0,92,3,60,'m27tom5_normcurve.mat')
%net_mergenet(91,[1:2],'mean',0)
%net_mergenet(92,[1:3],'mean',0)
%net_mergenet(91,3,'combine',0,93,1,60,'')
%net_mergenet(92,4,'combine',0,93,2,60,'m92tom91_normcurve.mat')
%net_mergenet(93,[1:2],'mean',0)
%  net_mergenet(93,1,'combine',0,94,1,60,'');
%  net_mergenet(93,2,'combine',0,94,2,60,'');
%  net_mergenet(6,63,'combine',0,94,3,60,'m6tom93_normcurve.mat');
%  net_mergenet(94,[1:3],'mean',0);
% Pos=4;
% for i=107:122
%     Pos=Pos+1;
%     net_mergenet(5,i,'combine',0,92,Pos,60,'')
% end
% for i=212:226
%     Pos=Pos+1;
%     net_mergenet(8,i,'combine',0,92,Pos,60,'')
% end
% for i=149:163
%     Pos=Pos+1;
%     net_mergenet(27,i,'combine',0,92,Pos,60,'')
% end
    
% net_mergenet(2,64,'combine',0,91,4,60,'')
% net_mergenet(2,65,'combine',0,91,5,60,'')
% net_mergenet(2,66,'combine',0,91,6,60,'')
% net_mergenet(2,67,'combine',0,91,7,60,'')
% net_mergenet(2,68,'combine',0,91,8,60,'')
% net_mergenet(2,69,'combine',0,91,9,60,'')
% net_mergenet(2,70,'combine',0,91,10,60,'')
% net_mergenet(2,71,'combine',0,91,11,60,'')
% net_mergenet(2,72,'combine',0,91,12,60,'')
% net_mergenet(2,73,'combine',0,91,13,60,'')
% net_mergenet(2,74,'combine',0,91,14,60,'')
% net_mergenet(2,75,'combine',0,91,15,60,'')
% net_mergenet(2,76,'combine',0,91,16,60,'')
% net_mergenet(2,77,'combine',0,91,17,60,'')
% net_mergenet(2,77,'combine',0,91,18,60,'')

% net_mergenet(3,70,'combine',0,91,19,60,'m3tom2_normcurve.mat')
% net_mergenet(3,71,'combine',0,91,20,60,'m3tom2_normcurve.mat')
% net_mergenet(3,72,'combine',0,91,21,60,'m3tom2_normcurve.mat')
% net_mergenet(3,73,'combine',0,91,22,60,'m3tom2_normcurve.mat')
% net_mergenet(3,74,'combine',0,91,23,60,'m3tom2_normcurve.mat')
% net_mergenet(3,75,'combine',0,91,24,60,'m3tom2_normcurve.mat')
% net_mergenet(3,76,'combine',0,91,25,60,'m3tom2_normcurve.mat')
% net_mergenet(3,77,'combine',0,91,26,60,'m3tom2_normcurve.mat')
% net_mergenet(3,78,'combine',0,91,27,60,'m3tom2_normcurve.mat')
% net_mergenet(3,79,'combine',0,91,28,60,'m3tom2_normcurve.mat')
% net_mergenet(3,80,'combine',0,91,29,60,'m3tom2_normcurve.mat')
% net_mergenet(3,81,'combine',0,91,30,60,'m3tom2_normcurve.mat')
% net_mergenet(3,82,'combine',0,91,31,60,'m3tom2_normcurve.mat')
% net_mergenet(3,83,'combine',0,91,32,60,'m3tom2_normcurve.mat')
% net_mergenet(3,84,'combine',0,91,33,60,'m3tom2_normcurve.mat')
% net_mergenet(3,85,'combine',0,91,34,60,'m3tom2_normcurve.mat')

% net_mergenet([92,6],[1,4,63],'combine',0,95,1,60,'');
%  net_mergenet([92,6],[2,4,63],'combine',0,95,2,60,'');
%  net_mergenet([92,6],[3,4,63],'combine',0,95,3,60,'');
%  net_mergenet([6,92],[63,63,4],'combine',0,95,4,60,'m6tom92_normcurve.mat');



%net_mergenet(8,55,'merge',0,0,[7:21],1,1,[1,25,50,75,100],'_1p','mouse',1)
%net_mergenet(27,24,'merge',0,0,[2:22],1,1,[1,25,50,75,100],'_2to22_1p','mouse',1)
%net_mergenet(8,228,'merge',0,0,[7:21],1,1,[1,25,50,75,100],'_1p','mouse',0)
%net_mergenet(27,164,'merge',0,0,[2:22],1,1,[1,25,50,75,100],'_2to22_1p','mouse',0)
% 
% 
% net_mergenet(27,149,'merge',0,0,[2:22],1,1,[1,25,50,75,100],'_2to22_1p','mouse',1)
% net_mergenet(27,150,'merge',0,0,[2:22],1,1,[1,25,50,75,100],'_2to22_1p','mouse',1)
% net_mergenet(27,151,'merge',0,0,[2:22],1,1,[1,25,50,75,100],'_2to22_1p','mouse',1)
% net_mergenet(27,152,'merge',0,0,[2:22],1,1,[1,25,50,75,100],'_2to22_1p','mouse',1)
% net_mergenet(27,153,'merge',0,0,[2:22],1,1,[1,25,50,75,100],'_2to22_1p','mouse',1)
% net_mergenet(27,154,'merge',0,0,[2:22],1,1,[1,25,50,75,100],'_2to22_1p','mouse',1)
% net_mergenet(27,155,'merge',0,0,[2:22],1,1,[1,25,50,75,100],'_2to22_1p','mouse',1)
% net_mergenet(27,156,'merge',0,0,[2:22],1,1,[1,25,50,75,100],'_2to22_1p','mouse',1)
% net_mergenet(27,157,'merge',0,0,[2:22],1,1,[1,25,50,75,100],'_2to22_1p','mouse',1)
% net_mergenet(27,158,'merge',0,0,[2:22],1,1,[1,25,50,75,100],'_2to22_1p','mouse',1)
% net_mergenet(27,159,'merge',0,0,[2:22],1,1,[1,25,50,75,100],'_2to22_1p','mouse',1)
% net_mergenet(27,160,'merge',0,0,[2:22],1,1,[1,25,50,75,100],'_2to22_1p','mouse',1)
% net_mergenet(27,161,'merge',0,0,[2:22],1,1,[1,25,50,75,100],'_2to22_1p','mouse',1)
% net_mergenet(27,162,'merge',0,0,[2:22],1,1,[1,25,50,75,100],'_2to22_1p','mouse',1)
% net_mergenet(27,163,'merge',0,0,[2:22],1,1,[1,25,50,75,100],'_2to22_1p','mouse',1)



function net_mergenet(ChipRank,NetRanks,MergeType,QLimitFlag,varargin)
global K

Continue=1;
if nargin<4
    h=warndlg('net_mergenet needs at least four parameters (ChipRank,NetRanks,MergeType & QLimitFlag');
    waitfor(h)
    Continue=0;
elseif nargin==5
    MinNetNb=varargin{1};
elseif nargin==5
    ChipRank(1)=ChipRank;
    NetRanks(1)=NetRanks;
    ChipRank(2)=varargin{1};
    NetRanks(2)=varargin{2};
elseif nargin==8
    NetRank=NetRanks;
    NewChipRank=varargin{1};
    NewNetRank=varargin{2};
    CorrLimit=varargin{3};    
    NormFile=varargin{4};
elseif nargin==12
    NetRank=NetRanks;
    NetRanks=varargin{2};
    ProbeNbLimit=varargin{3};
    PvCorrRank=varargin{4};
    NetFrequencies=varargin{5};
    Suffix=varargin{6};
    Species=varargin{7};
    MaxFlag=varargin{8};
end
ChipPos=find(K.chip.rank==ChipRank(1));
ProbesetNb=K.chip.probesetNb(ChipPos);
NetNb=length(NetRanks);
ChipDir=fullfile(K.dir.net,sprintf('m%u',ChipRank(1)));
if ~isfield(K.dir,'mldata')
    K.dir.mldata='/home/mbellis/sosma/data/psawn/mldata';
    K.dir.pydata='/home/mbellis/sosma/data/psawn/pydata';
    K.dir.rawdata='/home/mbellis/sosma/data/psawn/rawdata';
end

if Continue
    %% merge
    if isequal(MergeType,'combine')
        tic
        cd(K.dir.chip)
        eval(sprintf('load m%un%u_m%un%u_combinedps_corr%u.mat',ChipRank(1),NetRank(2),ChipRank(2),NetRank(3),CorrLimit))
        if ~isempty(NormFile)
            load(NormFile)
            NormFlag=1;
        else
            NormFlag=0;            
        end
        if exist(sprintf('m%u_ps2net.mat',ChipRank(1)))
            ModifFlag=1;
            load(sprintf('m%u_ps2net.mat',ChipRank(1)))
            NetPs=Ps2Net;
        else
            ModifFlag=0;
        end

        %make new directory for merged network
        NewNetDir=fullfile(K.dir.net,sprintf('m%u',NewChipRank),sprintf('n%u',NewNetRank));
        mkdir(NewNetDir)
        if isempty(find(K.chip.rank==NewChipRank))
            NewChipFlag=1;
        else
            NewChipFlag=0;
        end

        NetDir=fullfile(K.dir.net,sprintf('m%u',ChipRank(1)),sprintf('n%u',NetRank(1)));
        cd(NewNetDir)
        fidc=fopen(sprintf('c_m%u_n%u.4mat',NewChipRank,NewNetRank),'w','ieee-le');
        fida=fopen(sprintf('a_m%u_n%u.4mat',NewChipRank,NewNetRank),'w','ieee-le');


        %WRITE COMPLETE LINES IN A FIRST ROUND        
        cd(NetDir)
        PsPos=0;
        BlocNb=ceil(length(PsRank)/500);        
        for BlocL=1:BlocNb
            if BlocL<BlocNb
                ValC=uint8(zeros(ProbesetNb,500));
                ValA=uint8(zeros(ProbesetNb,500));
            else
                ValC=uint8(zeros(ProbesetNb,length(PsRank)-(500*(BlocNb-1))));
                ValA=uint8(zeros(ProbesetNb,length(PsRank)-(500*(BlocNb-1))));
            end
            BlocPsPos=0;
            for PsL=1:size(ValC,2)          
                PsPos=PsPos+1;
                BlocPsPos=BlocPsPos+1;
                CurrPsRank=PsRank{PsPos};
                if ModifFlag
                    CurrPsRank=NetPs(CurrPsRank);
                end
                Val=load_data(sprintf('c_m%u_n%u.4mat',ChipRank(1),NetRank(1)),'./',ProbesetNb,ProbesetNb,'uint8','ieee-le',1:ProbesetNb,CurrPsRank);
                if length(CurrPsRank)>1
                    Val=uint8(round(mean(double(Val),2)));
                end;
                ValC(:,BlocPsPos)=Val;

                Val=load_data(sprintf('a_m%u_n%u.4mat',ChipRank(1),NetRank(1)),'./',ProbesetNb,ProbesetNb,'uint8','ieee-le',1:ProbesetNb,CurrPsRank);
                if length(CurrPsRank)>1
                    Val=uint8(round(mean(double(Val),2)));
                end;
                ValA(:,BlocPsPos)=Val;
            end
            %WRITE CORRECT LINES IN A SECOND ROUND
            NewVal=zeros(length(PsRank),size(ValC,2));
            for PsL=1:length(PsRank)
                CurrPsRank=PsRank{PsL};
                if ModifFlag
                    CurrPsRank=NetPs(CurrPsRank);
                end
                if length(CurrPsRank)==1
                    NewVal(PsL,:)=ValC(CurrPsRank,:);
                else
                    try
                    NewVal(PsL,:)=uint8(mean(double(ValC(CurrPsRank,:))));
                    catch
                        'stop'
                    end
                end
            end
            if NormFlag
                %normalize Corr{2}
                CurrVal=NewVal;
                MaxC=max(max(NewVal));           
                MaxPos=find(OldVal<=MaxC);                
                for ValL=1:min(MaxPos(end),99)
                    Pos=find(CurrVal>=OldVal(ValL)&CurrVal<OldVal(ValL+1));
                    NewVal(Pos)=ValL;
                end
                if MaxPos(end)<99
                    Pos=find(CurrVal>=OldVal(MaxPos(end)+1)&CurrVal<100);
                    NewVal(Pos)=MaxPos(end)+1;
                end
            end
            fwrite(fidc,NewVal,'uint8',0,'ieee-le')

            NewVal=zeros(length(PsRank),size(ValA,2));         
            for PsL=1:length(PsRank)
                CurrPsRank=PsRank{PsL};
                if ModifFlag
                    CurrPsRank=NetPs(CurrPsRank);
                end
                if length(CurrPsRank)==1
                    NewVal(PsL,:)=ValA(CurrPsRank,:);
                else
                    NewVal(PsL,:)=uint8(mean(double(ValA(CurrPsRank,:))));
                end
            end
            if NormFlag
                %normalize Corr{2}
                CurrVal=NewVal;
                MaxC=max(max(NewVal));
                MaxPos=find(OldVal<=MaxC);
                for ValL=1:min(MaxPos(end),99)
                    Pos=find(CurrVal>=OldVal(ValL)&CurrVal<OldVal(ValL+1));
                    NewVal(Pos)=ValL;
                end
                if MaxPos(end)<99
                    Pos=find(CurrVal>=OldVal(MaxPos(end)+1)&CurrVal<100);
                    NewVal(Pos)=MaxPos(end)+1;
                end
            end
            fwrite(fida,NewVal,'uint8',0,'ieee-le')
        end
        fclose(fidc)
        fclose(fida)
        toc
        if NewChipFlag
            cd(K.dir.common)
            load chiplist
            Pos=length(Tempo.rank)+1;
            Tempo.rank(Pos,1)=NewChipRank;
            Tempo.myName{Pos,1}=sprintf('m%u',NewChipRank);
            Tempo.name{Pos,1}=sprintf('%s reduced',Tempo.name{ChipPos});
            Tempo.shortName{Pos,1}=sprintf('%s reduced',Tempo.shortName{ChipPos});
            Tempo.species{Pos,1}=Tempo.species{ChipPos};
            Tempo.probesetNb(Pos,1)=length(PsRank);
            Tempo.probeNb(Pos,1)=Tempo.probeNb(ChipPos);
            Tempo.compName{Pos,1}=Tempo.compName{ChipPos};
            Tempo.geoName{Pos,1}='';
            Tempo.type{Pos,1}=Tempo.type{ChipPos};
            cd(K.dir.common)
            save chiplist Tempo
            K.chip=Tempo;
        end

        %register new networks

        cd(K.dir.common)
        load netlist
        if NewChipFlag
            Pos=1;
        else
            if length(Tempo)>=NewChipRank
                if isempty(Tempo{NewChipRank})
                    Pos=1;
                else
                    Pos=length(Tempo{NewChipRank}.rank)+1;
                end
            else
                Pos=1;
            end
        end
        NetPos=find(Tempo{ChipRank(1)}.rank==NetRank(1));
        Tempo{NewChipRank}.name{Pos,1}=sprintf('m%un%u: %s',ChipRank(1),NetRank(1),Tempo{ChipRank(1)}.name{NetPos});
        Tempo{NewChipRank}.rank(Pos,1)=NewNetRank;
        Tempo{NewChipRank}.biolRank{Pos,1}=Tempo{ChipRank(1)}.biolRank{NetPos};
        Tempo{NewChipRank}.compNb(Pos,1)=Tempo{ChipRank(1)}.compNb(NetPos);
        Tempo{NewChipRank}.fdr(Pos,:)=Tempo{ChipRank(1)}.fdr(NetPos,:);
        Tempo{NewChipRank}.s(Pos,:)=Tempo{ChipRank(1)}.s(NetPos,:);
        Tempo{NewChipRank}.blocNb(Pos,1)=BlocNb;
        Tempo{NewChipRank}.blocSize(Pos,1)=500;
        Tempo{NewChipRank}.comment{Pos,1}=Tempo{ChipRank(1)}.comment{NetPos};
        Tempo{NewChipRank}.netMade(Pos,1)=1;
        Tempo{NewChipRank}.nb=length(Tempo{NewChipRank}.name);
        cd(K.dir.common)
        save netlist Tempo
        K.net=Tempo;
        
        Continue=0;


    elseif isequal(MergeType,'merge')
        Continue=0;
        BlocSize=5000;
        BlocNb=ceil(ProbesetNb/BlocSize);

        tic
        for FreqL=1:length(NetFrequencies)
            FreqL
            sprintf('merging net %u freq %u',NetRank,FreqL)
            %load PsMatrix corresponding to the current net frequency
            cd(fullfile(K.dir.mldata,Species,sprintf('m%u%s',ChipRank,Suffix)))
            StatFileName=sprintf('m%u_n%u_netnb%u_probenb%u_newps_stat_netprc%03u_pvcorr%u',ChipRank,NetRanks(1),NetNb,ProbeNbLimit,NetFrequencies(FreqL),PvCorrRank);
            eval(sprintf('load %s',StatFileName))
            if FreqL==1
                RoundNb=2;
            else
                RoundNb=1;
            end
            for RoundL=1:RoundNb
                tic
                %make new directory for merged network
                NewNetPos=length(K.net{ChipRank(1)}.name)+1;
                NewNetRank=setdiff([1:NewNetPos],K.net{ChipRank(1)}.rank);
                NewNetRank=NewNetRank(1);
                NewNetDir=fullfile(K.dir.net,sprintf('m%u',ChipRank(1)),sprintf('n%u',NewNetRank));
                mkdir(NewNetDir)


                %construct new networks by merging nodes that really target the same group of transcripts
                %MergeInfo contains
                %MergedPs: a cell array containing list of sorted old probeset ranks
                %the reference probeset is in first position (lower rank).
                %MergedPsNames:a cell array containing list of probeset names sorted as their ranks are in MergedPs
                %the reference probeset name is in first position (lower rank).

                %AllPsRanks: two column matrix with first column indicating the old probeset rank
                %and the new column indicating the new probeset rank (merged probesets have their
                %new rank set to the rank of their reference probeset)
                %
                NetDir=fullfile(K.dir.net,sprintf('m%u',ChipRank),sprintf('n%u',NetRank));

                %merge info and clear index for columns
                ClearIndex=[];
                %rank of merging ps in the original probe sets list
                MergingPsRank=[];
                %rank of merged ps in the original probe sets list
                MergedPsRank=[];
                if FreqL==1&RoundL==1
                    PsMade=zeros(size(PsMatrix,1),1);
                    for PsL=1:size(PsMatrix,1)
                        %for PsL=1:100
                        if PsMade(PsL)==0
                            %find all ps related to the current ps
                            CurrPs=find(PsMatrix(:,1)==PsMatrix(PsL,1));
                            PsMade(CurrPs)=1;
                            %process multiple ps
                            if length(CurrPs)>1
                                %take the ps with the smallest rank as the reference ps
                                MergingPsRank(end+1,1)=min(CurrPs);
                                %all other ps are merged with the reference ps
                                MergedPsRank{end+1,1}=setdiff(CurrPs,MergingPsRank(end));
                            end
                        end
                    end
                else

                    %don't consider group of transcripts of pivot probe set
                    %PsMatrix(find(PsMatrix(:,12)==1),10)=0;
                    %process each probe set once
                    PsMade=zeros(size(PsMatrix,1),1);
                    %PivotPb=[];
                    for PsL=1:size(PsMatrix,1)
                        %for PsL=1:100
                        if PsMade(PsL)==0
                            if PsMatrix(PsL,1)==0
                                PsMade(PsL)=1;
                            else
                                %find all ps related to the current ps at condition that the ps
                                %is assigned to something
                                CurrPs=find(PsMatrix(:,1)==PsMatrix(PsL,1));
                                PsMade(CurrPs)=1;
                                %process multiple ps
                                if length(CurrPs)>1
                                    %recover group of transcripts (GOT) targeted by the current ps
                                    CurrGrp=setdiff(unique(PsMatrix(CurrPs,10)),0);
                                    %verify if exist a pivot
                                    %                             PivotPos=find(PsMatrix(CurrPs,12)==1);
                                    %                             if ~isempty(PivotPos)
                                    %                                 %find pivots
                                    %                                 PivotNb=length(PivotPos);
                                    %                                 PivotGrp=zeros(PivotNb,length(CurrGrp));
                                    %                                 for PivotL=1:PivotNb
                                    %                                     PivotGrp(PivotL,find(dec2bin(PsMatrix(CurrPs(PivotPos(PivotL)),11))=='1'))=1;
                                    %                                 end
                                    %                                 %process pivots
                                    %                                 UsedGrp=[];
                                    %                                 for PivotL=1:length(PivotPos)
                                    %                                     CurrPivotGrp=find(PivotGrp(PivotL,:));
                                    %                                     for GrpL=1:length(CurrPivotGrp)
                                    %                                         try
                                    %                                             GrpPos=find(PsMatrix(CurrPs,10)==CurrGrp(CurrPivotGrp(GrpL)));
                                    %                                             %take the ps with the smallest rank as the reference ps
                                    %                                             if isempty(intersect(MergingPsRank,min(CurrPs(GrpPos))))
                                    %                                                 MergingPsRank(end+1,1)=min(CurrPs(GrpPos));
                                    %                                                 %all other ps are merged with the reference ps
                                    %                                                 MergedPsRank{end+1,1}=setdiff(CurrPs(GrpPos),MergingPsRank(end));
                                    %                                                 %add current pivot
                                    %                                                 MergedPsRank{end,1}=[MergedPsRank{end,1},CurrPs(PivotPos(PivotL))];
                                    %                                                 %plus eventually other pivots
                                    %                                                 OtherPivot=find(PivotGrp(PivotL+1:end,CurrPivotGrp(GrpL)));
                                    %                                                 if ~isempty(OtherPivot)
                                    %                                                     MergedPsRank{end,1}=sort([MergedPsRank{end,1},CurrPs(PivotPos(OtherPivot))']);
                                    %                                                 end
                                    %                                             end
                                    %                                             %eliminate these grp from PivotGrp
                                    %                                             AllPivot=find(PivotGrp(:,CurrPivotGrp(GrpL)));
                                    %                                             PivotGrp(AllPivot,CurrPivotGrp(GrpL))=0;
                                    %                                             %keep track of processed GOT
                                    %                                             UsedGrp(end+1,1)=CurrGrp(CurrPivotGrp(GrpL));
                                    %                                         catch
                                    %                                             PivotPb(end+1,1)=PsL;
                                    %                                         end
                                    %                                     end
                                    %                                 end
                                    %                                 %update CurrGrp to process probe set not involved in pivot
                                    %                                 %relationships
                                    %                                 CurrGrp=setdiff(CurrGrp,UsedGrp);
                                    %                             end
                                    if ~isempty(CurrGrp)
                                        %process each GOT
                                        for GrpL=1:length(CurrGrp)
                                            GrpPos=find(PsMatrix(CurrPs,10)==CurrGrp(GrpL));
                                            if length(GrpPos)>1 & isempty(intersect(MergingPsRank,min(CurrPs(GrpPos))))
                                                %take the ps with the smallest rank as the reference ps
                                                MergingPsRank(end+1,1)=min(CurrPs(GrpPos));
                                                %all other ps are merged with the reference ps
                                                MergedPsRank{end+1,1}=setdiff(CurrPs(GrpPos),MergingPsRank(end));
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

                MergedNb=0;
                MergingNb=length(MergingPsRank);
                for PsL=1:MergingNb
                    MergedNb=MergedNb+length(MergedPsRank{PsL});
                end

                %list of merged ps to be cleared
                for PsL=1:length(MergedPsRank)
                    PsRanks=MergedPsRank{PsL};
                    if length(PsRanks)==0
                        h=errordlg('%u MergedPsRank is empty',PsL);
                        waitfor(h)
                        error('process canceled')
                    else
                        try
                            ClearIndex=[ClearIndex;PsRanks];
                        catch
                            ClearIndex=[ClearIndex;PsRanks'];
                        end
                    end
                end
                ClearIndex=unique(ClearIndex);

                %construct correspondance between old ans new ranks
                NewPsNb=ProbesetNb-length(ClearIndex);
                NewBlocNb=ceil(NewPsNb/BlocSize);
                Trans=zeros(ProbesetNb,4);
                Trans(:,1)=1:ProbesetNb;
                Rank=0;
                for PsL=1:ProbesetNb
                    if isempty(find(ClearIndex==PsL))
                        Rank=Rank+1;
                        Trans(PsL,2)=Rank;
                        Pos=find(MergingPsRank==PsL);
                        if ~isempty(Pos)
                            Trans(PsL,3)=1;
                            for PsL1=1:length(MergedPsRank{Pos})
                                Trans(MergedPsRank{Pos}(PsL1),2)=Rank;
                                Trans(MergedPsRank{Pos}(PsL1),4)=1;
                            end
                        end
                    end
                end

                %CLEAR MERGED LINES IN A FIRST ROUND
                for BlocL=1:BlocNb
                    for NetL=1:2
                        if NetL==1
                            NetFile=sprintf('c_m%u_n%u.4mat',ChipRank,NetRank);
                            NewNetFile=sprintf('cc_m%u_n%u.4mat',ChipRank,NewNetRank);
                        else
                            NetFile=sprintf('a_m%u_n%u.4mat',ChipRank,NetRank);
                            NewNetFile=sprintf('aa_m%u_n%u.4mat',ChipRank,NewNetRank);
                        end
                        cd(NetDir)
                        Val=load_data(NetFile,'./',ProbesetNb,ProbesetNb,'uint8','ieee-le',1:ProbesetNb,((BlocL-1)*BlocSize)+1:min(BlocL*BlocSize,ProbesetNb));
                        for PsL=1:MergingNb
                            PsRanks=[MergingPsRank(PsL);MergedPsRank{PsL}];
                            if MaxFlag
                                %take the maximum of all the group
                                Val(MergingPsRank(PsL),:)=max(Val(PsRanks,:));
                            else
                                %take the mean of all the group
                                Val(MergingPsRank(PsL),:)=mean(Val(PsRanks,:));
                            end
                        end
                        %clear merged lines
                        Val(ClearIndex,:)=[];
                        cd(NewNetDir)
                        if NetL==1
                            C=Val;
                            clear Val
                            save_data(C,NewNetFile,'./','a+','uint8','ieee-le')
                        else
                            A=Val;
                            clear Val
                            save_data(A,NewNetFile,'./','a+','uint8','ieee-le')
                        end
                    end
                end
                %CLEAR MERGED COLUMNS IN A SND ROUND
                cd(NewNetDir)
                for BlocL=1:NewBlocNb
                    for NetL=1:2
                        if NetL==1
                            NetFile=sprintf('cc_m%u_n%u.4mat',ChipRank,NewNetRank);
                            NewNetFile=sprintf('c_m%u_n%u.4mat',ChipRank,NewNetRank);
                        else
                            NetFile=sprintf('aa_m%u_n%u.4mat',ChipRank,NewNetRank);
                            NewNetFile=sprintf('a_m%u_n%u.4mat',ChipRank,NewNetRank);
                        end
                        Val=load_data(NetFile,'./',NewPsNb,ProbesetNb,'uint8','ieee-le',((BlocL-1)*BlocSize)+1:min(BlocL*BlocSize,NewPsNb),1:ProbesetNb);
                        for PsL=1:MergingNb
                            PsRanks=[MergingPsRank(PsL);MergedPsRank{PsL}];
                            if MaxFlag
                                %take the maximum of all the group
                                Val(:,MergingPsRank(PsL))=max(Val(:,PsRanks),[],2);
                            else
                                %take the mean of all the group
                                Val(:,MergingPsRank(PsL))=mean(Val(:,PsRanks),2);
                            end
                        end
                        %clear merged lines
                        Val(:,ClearIndex)=[];
                        if NetL==1
                            C=Val';
                            clear Val;
                            save_data(C,NewNetFile,'./','a+','uint8','ieee-le')
                        else
                            A=Val';
                            clear Val;
                            save_data(A,NewNetFile,'./','a+','uint8','ieee-le')
                        end
                    end
                end
                eval(sprintf('save t_m%u_n%u Trans',ChipRank,NewNetRank))

                %del intermediate files
                delete(sprintf('cc_m%u_n%u.4mat',ChipRank,NewNetRank))
                delete(sprintf('aa_m%u_n%u.4mat',ChipRank,NewNetRank))

                %register new networks
                cd(K.dir.common)
                load netlist
                RefPos=find(Tempo{ChipRank(1)}.rank==NetRank);
                %intersect of networks
                if FreqL==1&RoundL==1
                    Tempo{ChipRank(1)}.name{NewNetPos,1}=sprintf('%s of networks %u NetPrc0',MergeType,NetRank);
                else
                    Tempo{ChipRank(1)}.name{NewNetPos,1}=sprintf('%s of networks %u NetPrc%u',MergeType,NetRank,NetFrequencies(FreqL));
                end
                Tempo{ChipRank(1)}.rank(NewNetPos,1)=NewNetRank;
                Tempo{ChipRank(1)}.biolRank{NewNetPos,1}='';
                Tempo{ChipRank(1)}.compNb(NewNetPos,1)=Tempo{ChipRank}.compNb(RefPos,1);
                Tempo{ChipRank(1)}.fdr(NewNetPos,:)=Tempo{ChipRank}.fdr(RefPos,:);
                Tempo{ChipRank(1)}.s(NewNetPos,:)=Tempo{ChipRank}.s(RefPos,:);
                Tempo{ChipRank(1)}.blocNb(NewNetPos,1)=NewBlocNb;
                Tempo{ChipRank(1)}.blocSize(NewNetPos,1)=5000;
                Tempo{ChipRank(1)}.comment{NewNetPos,1}=Tempo{ChipRank}.comment{RefPos,1};
                Tempo{ChipRank(1)}.netMade(NewNetPos,1)=1;
                Tempo{ChipRank(1)}.nb=length(Tempo{ChipRank}.name);
                cd(K.dir.common)
                save netlist Tempo
                K.net=Tempo;
                toc
            end %RoundL
        end %freqL
        toc
        %% fitting
    elseif isequal(MergeType,'fitting')

        if length(ChipRank)~=2
            h=warndlg('In case of fitting You must select a second chip');
            waitfor(h)
            Continue=0;
        else
            cd(fullfile(K.dir.net,sprintf('m%u',ChipRank(1))))
            eval(sprintf('load m%utom%u_normcurve',ChipRank(1),ChipRank(2)))
        end
        for NetL=1:NetNb
            ChipDir=fullfile(K.dir.net,sprintf('m%u',ChipRank(1)),sprintf('n%u',NetRanks(NetL)));
            Pos=length(K.net{ChipRank(1)}.name)+1;
            NetRank=setdiff([1:Pos],K.net{ChipRank(1)}.rank);
            NetRank=NetRank(1);
            NetDir=fullfile(K.dir.net,sprintf('m%u',ChipRank(1)),sprintf('n%u',NetRank));
            mkdir(NetDir)
            %copy f file
            try
                copyfile(fullfile(ChipDir,sprintf('f_m%u_n%u.4mat',ChipRank(1),NetRanks(NetL))),...
                    fullfile(NetDir,sprintf('f_m%u_n%u.4mat',ChipRank(1),NetRank)))
            catch
                h=warndlg(sprintf('f_m%u_n%u.4mat not copied',ChipRank(1),NetRanks(NetL)));
                waitfor(h)
            end
            BlocNb=ceil(ProbesetNb/5000);
            for BlocL=1:BlocNb
                BlocL
                cd(ChipDir)
                tic
                CurrC=load_data(sprintf('c_m%u_n%u.4mat',ChipRank(1),NetRanks(NetL)),'./',ProbesetNb,ProbesetNb,'uint8','ieee-le',1:ProbesetNb,((BlocL-1)*5000)+1:min(BlocL*5000,ProbesetNb));
                toc
                C=CurrC;
                MaxC=max(max(C));
                MaxPos=find(OldVal<=MaxC);
                tic
                for ValL=1:min(MaxPos(end),99)
                    Pos=find(CurrC>=OldVal(ValL)&CurrC<OldVal(ValL+1));
                    C(Pos)=ValL;
                end
                if MaxPos(end)<99
                    Pos=find(CurrC>=OldVal(MaxPos(end)+1)&CurrC<100);
                    C(Pos)=MaxPos(end)+1;
                end
                toc
                clear CurrC
                cd(NetDir)
                save_data(C,sprintf('c_m%u_n%u.4mat',ChipRank(1),NetRank),'./','a+','uint8','ieee-le');
                clear C

                cd(ChipDir)
                CurrA=load_data(sprintf('a_m%u_n%u.4mat',ChipRank(1),NetRanks(NetL)),'./',ProbesetNb,ProbesetNb,'uint8','ieee-le',1:ProbesetNb,((BlocL-1)*5000)+1:min(BlocL*5000,ProbesetNb));
                A=CurrA;
                MaxA=max(max(A));
                MaxPos=find(OldVal<=MaxA);
                tic
                for ValL=1:MaxPos(end)
                    Pos=find(CurrA>=OldVal(ValL)&CurrA<OldVal(ValL+1));
                    A(Pos)=ValL;
                end
                if MaxPos(end)<99
                    Pos=find(CurrA>=OldVal(MaxPos(end)+1)&CurrA<100);
                    A(Pos)=MaxPos(end)+1;
                end
                toc
                clear CurrA
                cd(NetDir)
                save_data(A,sprintf('a_m%u_n%u.4mat',ChipRank(1),NetRank),'./','a+','uint8','ieee-le');
                clear A
            end

            cd(K.dir.common)
            load netlist
            Pos=length(K.net{ChipRank(1)}.name)+1;
            RefPos=find(Tempo{ChipRank(1)}.rank==NetRanks(NetL));
            %intersect of networks
            Tempo{ChipRank(1)}.name{Pos,1}=sprintf('%s of networks %u towards m%u',MergeType,NetRanks(NetL),ChipRank(2));
            Tempo{ChipRank(1)}.rank(Pos,1)=NetRank;
            Tempo{ChipRank(1)}.biolRank{Pos,1}='';
            Tempo{ChipRank(1)}.compNb(Pos,1)=Tempo{ChipRank(1)}.compNb(RefPos,1);
            Tempo{ChipRank(1)}.fdr(Pos,:)=Tempo{ChipRank(1)}.fdr(RefPos,:);
            Tempo{ChipRank(1)}.s(Pos,:)=Tempo{ChipRank(1)}.s(RefPos,:);
            Tempo{ChipRank(1)}.blocNb(Pos,1)=Tempo{ChipRank(1)}.blocNb(RefPos,1);
            Tempo{ChipRank(1)}.blocSize(Pos,1)=Tempo{ChipRank(1)}.blocSize(RefPos,1);
            Tempo{ChipRank(1)}.comment{Pos,1}=Tempo{ChipRank(1)}.comment{RefPos,1};
            Tempo{ChipRank(1)}.netMade(Pos,1)=1;
            Tempo{ChipRank(1)}.nb=length(Tempo{ChipRank(1)}.name);
            cd(K.dir.common)
            save netlist Tempo
            K.net=Tempo;
        end
    else
        Pos=length(K.net{ChipRank(1)}.name)+1;
        NetRank=setdiff([1:Pos],K.net{ChipRank(1)}.rank);
        NetRank=NetRank(1);
        NetDir=fullfile(K.dir.net,sprintf('m%u',ChipRank(1)),sprintf('n%u',NetRank));
        mkdir(NetDir)
        cd(NetDir)
        Continue=1;
        if nargin==4
            MinNetNb=NetNb;
        end
        if NetNb<MinNetNb
            h=warndlg(sprintf('You must select at least %u networks.',MinNetNb));
            waitfor(h)
            Continue=0;
        elseif NetNb<2
            h=warndlg('You must select at least two networks.');
            waitfor(h)
            Continue=0;
        end
    end
    if Continue
        %% union, intersect, mean
        %         if QLimitFlag
        %             HistRepNbC=zeros(NetNb,1);
        %             HistIfExistC=zeros(100,1);
        %             HistIfExistA=zeros(100,1);
        %             HistIfSupeC=zeros(100,1);
        %             HistIfSupeA=zeros(100,1);
        %         end

        %create NetFile
        NetA1=sprintf('a_m%u_n%u.4mat',ChipRank(1),NetRank);
        NetC1=sprintf('c_m%u_n%u.4mat',ChipRank(1),NetRank);
        NetF1=sprintf('f_m%u_n%u.4mat',ChipRank(1),NetRank);
        %NetRepC=sprintf('c_m%u_n%u_rep.4mat',ChipRank(1),NetRank);
        cd(NetDir)
        NetA1Fid=fopen(NetA1,'w','ieee-le');
        NetC1Fid=fopen(NetC1,'w','ieee-le');
        NetF1Fid=fopen(NetF1,'w','ieee-le');
        %NetRepCFid=fopen(NetRepC,'w','ieee-le');
        HReadC=zeros(NetNb,1);
        HReadA=zeros(NetNb,1);
        HReadF=zeros(NetNb,1);

        for NetL=1:NetNb
            CurrNet=NetRanks(NetL);
            if NetL==1
                NetList=num2str(CurrNet);
            else
                NetList=sprintf('%s %u',NetList,CurrNet);
            end
            NetDir=fullfile(ChipDir,sprintf('n%u',CurrNet));
            cd(NetDir)
            NetFile=sprintf('a_m%u_n%u.4mat',ChipRank(1),CurrNet);
            HReadA(NetL)=fopen(NetFile,'r','ieee-le');
            NetFile=sprintf('c_m%u_n%u.4mat',ChipRank(1),CurrNet);
            HReadC(NetL)=fopen(NetFile,'r','ieee-le');
            NetFile=sprintf('f_m%u_n%u.4mat',ChipRank(1),CurrNet);
            HReadF(NetL)=fopen(NetFile,'r','ieee-le');
        end


        %process each probe set
        tic
        t=0;
        PsNb=0;
        for PsL=1:ProbesetNb
            PsNb=PsNb+1;
            %recover CORR and ANTI values for all the networks
            CurrC=uint8(ones(ProbesetNb,NetNb)*NaN);
            CurrA=uint8(ones(ProbesetNb,NetNb)*NaN);
            %CurrF=uint8(ones(ProbesetNb,NetNb)*NaN);
            CORR=uint8(zeros(ProbesetNb,1));
            ANTI=uint8(zeros(ProbesetNb,1));
            FREQ=uint8(zeros(ProbesetNb,1));
            for NetL=1:NetNb
                cd(fullfile(ChipDir,sprintf('n%u',NetRanks(NetL))));
                CurrA(:,NetL)=load_data(HReadA(NetL),'',ProbesetNb,ProbesetNb,'uint8','ieee-le',1:ProbesetNb,PsL);
                CurrC(:,NetL)=load_data(HReadC(NetL),'',ProbesetNb,ProbesetNb,'uint8','ieee-le',1:ProbesetNb,PsL);
                try
                CurrF(:,NetL)=load_data(HReadF(NetL),'',ProbesetNb,ProbesetNb,'uint8','ieee-le',1:ProbesetNb,PsL);
                catch
                end
            end
            switch MergeType
                case {'intersection','union'}
                    %number of significative values
                    RepInfo=max([sum(CurrC>0,2),sum(CurrA>0,2)],[],2);
                    %write RepInfo
                    %save_data(RepInfo,NetRepCFid,'','a','uint8','ieee-le',ProbesetNb,[],PsL);
                    %recover data if sup or equal to MinNetNb
                    SelIndex=RepInfo>=MinNetNb;
                    if isequal(MergeType,'intersection')
                        CORR(SelIndex)=uint8(ceil(mean(CurrC(SelIndex,:),2)));
                        ANTI(SelIndex)=uint8(ceil(mean(CurrA(SelIndex,:),2)));
                        RAWCORR(SelIndex)=mean((double(CurrC(SelIndex,:)).*double(CurrF(SelIndex,:)))/100,2);
                        FREQ(SelIndex)=uint8((double(CORR(SelIndex))./RAWCORR(SelIndex))*100);
                        %remove 100
                        CORR(PsL)=0;
                        %frequency of values
                        HistIfSupeC=HistIfSupeC+histc(CORR,1:100);
                        HistIfSupeA=HistIfSupeA+histc(ANTI,1:100);
                    else
                        %recover data if exist at least one value
                        CurrC=single(CurrC);
                        CurrA=single(CurrA);
                        CurrC(find(CurrC==0))=nan;
                        CurrA(find(CurrA==0))=nan;
                        SelIndex=RepInfo>=1;
                        CORR(SelIndex)=uint8(ceil(nanmean(CurrC(SelIndex,:)')'));
                        ANTI(SelIndex)=uint8(ceil(nanmean(CurrA(SelIndex,:)')'));
                        %RAWCORR(SelIndex)=nanmean((double(CurrC(SelIndex,:)).*double(CurrF(SelIndex,:)))/100,2);
                        %FREQ(SelIndex)=uint8((double(CORR(SelIndex))./RAWCORR(SelIndex))*100);
                        %         %write ANTI & CORR
                        %save_data(CORR,NetC2Fid,'','a','uint8','ieee-le',ProbesetNb,[],PsL);
                        %save_data(ANTI,NetA2Fid,'','a','uint8','ieee-le',ProbesetNb,[],PsL);
                        %remove 100
                        CORR(PsL)=0;
                        %frequency of values
                        %HistIfExistC=HistIfExistC+histc(CORR,1:100);
                        %HistIfExistA=HistIfExistA+histc(ANTI,1:100);


                        %frequency of repetition
                        %HistRepNbC=HistRepNbC+histc(RepInfo,1:NetNb);
                        %if PsNb==100
                        %    t=t+toc;
                        %    sprintf('%u Ps processed in %u min : estimates that rests %u min',PsL,round(t/60),round((ProbesetNb-PsL)*t/(PsL*60)))
                        %    tic
                        %    PsNb=0;
                        %end

                    end
                case 'mean'
                    CORR=uint8(ceil(mean(CurrC,2)));
                    ANTI=uint8(ceil(mean(CurrA,2)));
                    try
                    FREQ=uint8(ceil(mean(CurrF,2)));
                    catch
                    end
            end
            %write ANTI & CORR
            save_data(CORR,NetC1Fid,'','a','uint8','ieee-le',ProbesetNb,[],PsL);
            save_data(ANTI,NetA1Fid,'','a','uint8','ieee-le',ProbesetNb,[],PsL);
            try
            save_data(FREQ,NetF1Fid,'','a','uint8','ieee-le',ProbesetNb,[],PsL);
            catch
            end
        end

        %         if QLimitFlag
        %             %remove 100 on corr diagonal
        %             HistRepNbC(end)=HistRepNbC(end)-ProbesetNb;
        %             eval(sprintf('save m%u_n%u_hist Hist*',ChipRank(1),NetRank))
        %             h=figure;
        %             set(h,'color',[1,1,1])
        %             subplot(1,2,1)
        %             plot(1:100,HistIfExistC/sum(HistIfExistC),'m')
        %             hold on
        %             plot(1:100,HistIfExistA/sum(HistIfExistA),'c')
        %             plot(1:100,HistIfSupeC/sum(HistIfSupeC),'r')
        %             plot(1:100,HistIfSupeA/sum(HistIfSupeA),'b')
        %             title('CORR and ANTI distributions')
        %             legend({'Union CORR','Union ANTI','Intersect CORR','Intersect ANTI'})
        %             subplot(1,2,2)
        %             set(h,'color',[1,1,1])
        %             plot(1:NetNb,HistRepNbC/sum(HistRepNbC),'K')
        %             title('reproducibility (nb of networks)')
        %
        %             fclose(NetC1Fid);
        %             fclose(NetA1Fid);
        %             fclose(NetC2Fid);
        %             fclose(NetA2Fid);
        %             fclose(NetRepCFid);
        for NetL=1:NetNb
            fclose(HReadA(NetL));
            fclose(HReadC(NetL));
        end
        %  end

        %register new networks
        cd(K.dir.common)
        load netlist
        RefPos=find(Tempo{ChipRank(1)}.rank==NetRanks(1));
        %intersect of networks
        Tempo{ChipRank(1)}.name{Pos,1}=sprintf('%s of networks %s',MergeType,NetList);
        Tempo{ChipRank(1)}.rank(Pos,1)=NetRank;
        Tempo{ChipRank(1)}.biolRank{Pos,1}='';
        Tempo{ChipRank(1)}.compNb(Pos,1)=Tempo{ChipRank(1)}.compNb(RefPos,1);
        Tempo{ChipRank(1)}.fdr(Pos,:)=Tempo{ChipRank(1)}.fdr(RefPos,:);
        Tempo{ChipRank(1)}.s(Pos,:)=Tempo{ChipRank(1)}.s(RefPos,:);
        Tempo{ChipRank(1)}.blocNb(Pos,1)=Tempo{ChipRank(1)}.blocNb(RefPos,1);
        Tempo{ChipRank(1)}.blocSize(Pos,1)=Tempo{ChipRank(1)}.blocSize(RefPos,1);
        Tempo{ChipRank(1)}.comment{Pos,1}=Tempo{ChipRank(1)}.comment{RefPos,1};
        Tempo{ChipRank(1)}.netMade(Pos,1)=1;
        Tempo{ChipRank(1)}.nb=length(Tempo{ChipRank(1)}.name);
        cd(K.dir.common)
        save netlist Tempo
        K.net=Tempo;
    end
end

