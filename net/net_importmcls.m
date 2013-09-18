%=========================%
% FUNCTION NET_IMPORTMCLS %
%=========================%

%NET_IMPORTMCLS: import several results of MCL clustering
%INPUT PARAMETERS
% 1  ChipRank: chip rank
% 2  NetRanks: list of net ranks
% 3  Inflation: MCL paramter
% 4  ListRank: rank of probe set list used
% 7  CorrLimits: CORR limits
% 8  DiffLimits: CORR-ANTI limits
% 9  RawLImits RAW CORR limits
% 10 DiffRawLimits: RAW CORR - ROW ANTI limits
% 11 ImportFlag: if 1 import data, if 0 control that all results exist (if not write command
%    files used to complete the results
% 12 FileName: Result file name

% net_importmcls(2,8,[24:38],[1:12],[20:10:70],[0:10:60],[10:10:60],[0:10:50],'n24ton38',1)
% net_importmcls(2,8,[230],[201,205],[],[10:2:20],[],[],'n230A',1)
% 
% net_importmcls(2,27,[25:45],2,[1:6],1,1,1,[20:10:70],[0:10:60],[10:10:60],[0:10:50],1,'n25ton45')
% net_importmcls(2,27,[88:108],2,[1:6],1,1,1,[20:10:70],[0:10:60],[10:10:60],[0:10:50],1,'n88ton108')
% net_importmcls(27,[118],2,[1:6],1,1,1,[20:10:70],[10:10:60],[10:10:60],[10:10:50],1,'n118')

% net_importmcls(8,[125],2,[1:12],1,1,1,[20:10:70],[0:10:60],[10:10:60],[0:10:50],3,0)
% net_importmcls(8,[125],2,[1:12],1,1,1,[20:10:70],[10:10:60],[10:10:60],[10:10:50],3,1)
% net_importmcls(8,[125],2,[7:12],1,1,1,[20:10:70],[20:10:60],[20:10:60],[10:10:50],4,1)

% net_importmcls(27,[164],2,[1:6],1,1,1,[20:10:70],[10:10:60],[10:10:60],[10:10:50],1,'n164')
% net_importmcls(27,[165],2,[1:6],1,1,1,[20:10:70],[10:10:60],[10:10:60],[10:10:50],1,'n165')
% net_importmcls(27,[166],2,[1:6],1,1,1,[20:10:70],[10:10:60],[10:10:60],[10:10:50],1,'n166')
% net_importmcls(8,[228],2,[7:12],1,1,1,[20:2:30],[20:2:30],[],[],1,'n228A')
% net_importmcls(27,[164],2,[1:6],1,1,1,[20:2:30],[20:2:30],[],[],1,'n164A')
% net_importmcls(8,[228],2,[7:12],1,1,1,[12:2:18],[12:2:18],[],[],1,'n228B')
% net_importmcls(2,27,[164],[1:6],[12:2:18],[12:2:18],[],[],'n164B',1)
% net_importmcls(27,[169],2,[201,205],[],[10:2:20],[],[],1,'n169A')
% net_importmcls(2,27,[149],[201,205],[],[14],[],[],'n149A',1)
% net_importmcls(27,[170],2,[201,205],[],[10:2:20],[],[],1,'n170A')
% net_importmcls(27,[172],2,[12:22],[],[14],[],[],'172A',1)
% net_importmcls(2,27,[24,164],[1:13],[],[14],[],[],'n24n164A',1)
% net_importmcls(2,8,[55,228],[14:26],[],[14],[],[],'n55n228A',1)
% net_importmcls(2,27,[172],[15:26],[],[14],[],[],'172A',1)
% net_importmcls(2,27,[173],[28:39],[],[14],[],[],'173A',1)
% net_importmcls(2,27,[174],[41:52],[],[14],[],[],'174A',1)
% net_importmcls(2,27,[175],[54:65],[],[14],[],[],'175A',1)
% net_importmcls(2,27,[176],[67:78],[],[14],[],[],'176A',1)
% net_importmcls(2,27,[177],[80:91],[],[14],[],[],'177A',1)
% net_importmcls(2,27,[178],[15:26],[],[14],[],[],'178A',1)
% net_importmcls(2,27,[179],[28:39],[],[14],[],[],'179A',1)
% net_importmcls(2,27,[180],[41:52],[],[14],[],[],'180A',1)
% net_importmcls(2,27,[181],[54:65],[],[14],[],[],'181A',1)
% net_importmcls(2,27,[182],[67:78],[],[14],[],[],'182A',1)
% net_importmcls(2,27,[183],[80:91],[],[14],[],[],'183A',1)
% net_importmcls(2,27,[24,164],[1:13],[],[14],[],[],'n24n164A',1)
% net_importmcls(2,27,[2:22],[5,6],[],[14],[],[],'n2n22B',1)
% net_importmcls(2,8,[234],[41:52],[],[14],[],[],'n234A',1)
% net_importmcls(2,8,[235],[67:78],[],[14],[],[],'n235A',1)
% net_importmcls(2,8,[236],[93:104],[],[14],[],[],'n236A',1)
% net_importmcls(2,8,[237],[119:130],[],[14],[],[],'n237A',1)
% net_importmcls(2,8,[238],[145:156],[],[14],[],[],'n238A',1)
% net_importmcls(2,8,[239],[171:182],[],[14],[],[],'n239A',1)
% net_importmcls(2,8,[240],[41:52],[],[14],[],[],'n240A',1)
% net_importmcls(2,8,[241],[67:78],[],[14],[],[],'n241A',1)
% net_importmcls(2,8,[242],[93:104],[],[14],[],[],'n242A',1)
% net_importmcls(2,8,[243],[119:130],[],[14],[],[],'n243A',1)
% net_importmcls(2,8,[244],[145:156],[],[14],[],[],'n244A',1)
% net_importmcls(2,8,[245],[171:182],[],[14],[],[],'n245A',1)
%net_importmcls(2,91,[1,2,3],1,[],[10],[],[],'m91a',1)
%net_importmcls(2,92,[1,2,3,4],1,[],[10],[],[],'m92a',1)
%net_importmcls(2,93,[1,2,3],1,[],[10],[],[],'m93a',1)
%net_importmcls(2,94,[1,2,3,4],1,[],[10],[],[],'m94a',1)
%net_importmcls(3,27,[24],[1:13],[],[14],[],[],'n24A',1)
%net_importmcls(4,27,[24],[1:13],[],[14],[],[],'n24B',1)
%net_importmcls(5,27,[24],[1:13],[],[14],[],[],'n24C',1)
% net_importmcls(3,8,[55],[14:26],[],[14],[],[],'n55A',1)
% net_importmcls(4,8,[55],[14:26],[],[14],[],[],'n55B',1)
% net_importmcls(5,8,[55],[14:26],[],[14],[],[],'n55C',1)
% net_importmcls(2,2,[80],[11,13],[],[0:10:70],[],[],'n80A',1)
% net_importmcls(2,3,[86],[11,13],[],[0:10:70],[],[],'n86A',1)
% net_importmcls(2,5,[123],[11,13],[],[0:10:70],[],[],'n123A',1)
% net_importmcls(2,6,[63],[11,13],[],[0:10:70],[],[],'n63A',1)
% net_importmcls(2,8,[228],[11,13],[],[0:10:70],[],[],'n228F',1)
% net_importmcls(2,27,[164],[11,13],[],[0:10:70],[],[],'n164F',1)
% net_importmcls(2,2,[80],[14,26],[],[0:10:70],[],[],'n80A',1)
%
% net_importmcls(2,2,[64:78],[24,26],[],[0:10:70],[],[],'n64n78A',1)
% net_importmcls(2,3,[70:85],[11,13],[],[0:10:70],[],[],'n70n85A',1)
% net_importmcls(2,5,[107:122],[11,13],[],[0:10:70],[],[],'n107n122A',1)
% net_importmcls(2,6,[48:62],[11,13],[],[0:10:70],[],[],'n48n62A',1)
% net_importmcls(2,8,[212:226],[11,13],[],[0:10:70],[],[],'n212n226A',1)
% net_importmcls(2,27,[149:163],[11,13],[],[10:10:70],[],[],'n149n163A',1)





%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%
%                         c) Michel Bellis                                                %
%                         michel.bellis@crbm.cnrs.fr                                      %
%            Affiliation:  CNRS (Centre National de la Recherche Scientifique - France)    %
%  Bioinformatic Project:  ARRAYMATIC => http://code.google.com/p/arraymatic               %
%        Code Repository:  GITHUB => http://github.com/mbellis                             %
%          Personal Page:  http://bns.crbm.cnrs.fr                                         %
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
%  THIS CODE IS DISTRIBUTED UNDER THE CeCILL LICENSE, WHICH IS COMPATIBLE WITH       %
%  THE GNU GENERAL PUBLIC LICENCE AND IN ACCORDANCE WITH THE EUROPEAN LEGISLATION.   %
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%


function net_importmcls(Inflation,ChipRank,NetRanks,ListRank,CorrLimits,DiffLimits,RawLimits,DiffRawLimits,FileName,ImportFlag)
global K
SList=sort(ListRank);
cd(K.dir.common)
ResNb=length(NetRanks)*length(ListRank)*(length(CorrLimits)+length(DiffLimits)+length(RawLimits)+length(DiffRawLimits));
load mcljury
CorrFlag=0;
DiffFlag=0;
RawFlag=0;
DiffRawFlag=0;
if ~isempty(CorrLimits)
    CorrFlag=1;
end
if ~isempty(DiffLimits)
    DiffFlag=1;
end
if ~isempty(RawLimits)
    RawFlag=1;
end
if ~isempty(DiffRawLimits)
    DiffRawFlag=1;
end

cd(K.dir.chip)
try
    eval(sprintf('load m%u_ps2net',ChipRank))
    %replace 0 by their current position (do not change)
    NetZeroPos=find(Net2Ps==0);
    if ~isempty(NetZeroPos)        
        Net2Ps(NetZeroPos)=NetZeroPos;
    end
    ModifFlag=1;
catch
    ModifFlag=0;
end
ListNb=length(ListRank);
ChipPos=find(K.chip.rank==ChipRank);
PsNb=K.chip.probesetNb(ChipPos);
CmlDir=fullfile(K.dir.net,sprintf('m%u',ChipRank),'mcl');
cd(CmlDir)


MissingNb=0;
Jury={};
Clu={};
CluNb={};
%Exists for measures of inter node correlations
for CorrTypeL=1:4
    DoIt=0;
    switch CorrTypeL
        case 1
            %CORR
            DoIt=CorrFlag;
            Suffix='';
            CurrLimits=CorrLimits;
        case 2
            %CORR-ANTI
            DoIt=DiffFlag;
            Suffix='_diff';
            CurrLimits=DiffLimits;

        case 3
            % RAW CORR
            DoIt=RawFlag;
            Suffix='_raw';
            CurrLimits=RawLimits;
        case 4
            % RAW CORR - RAW ANTI
            DoIt=DiffRawFlag;
            Suffix='_diff_raw';
            CurrLimits=DiffRawLimits;
    end
    Missing{CorrTypeL}=[];

    if DoIt
        ErrFile=sprintf('mcl_%s.mclerr',FileName)
        if exist(ErrFile,'file')
            fid=fopen(ErrFile,'r');            
            while 1
                CurrLine=fgetl(fid);
                if ~ischar(CurrLine)
                    break
                end
                CurrJury=regexp(CurrLine,'(?<=\[mcl\] jury pruning synopsis: <).*(?=>)','match');
                if ~isempty(CurrJury)
                    JuryVal=regexp(CurrJury{1},'\d+.\d+','match');
                    JuryVal=str2num(JuryVal{1});
                    JuryPhrase=regexp(CurrJury{1},'(?<=\d+ or ).+','match');
                    if isempty(strmatch(JuryPhrase{1},MclJury.name,'exact'))
                        MclJury.name{end+1,1}=JuryPhrase{1};
                        MclJury.val(end+1,1)=JuryVal;
                    end
                    CurrList={};
                    while isempty(CurrList)
                        CurrLine=fgetl(fid);
                        if ~ischar(CurrLine)
                            break
                        else
                            CurrList=regexp(CurrLine,'(?<=.*_c\d+_)\d+(?=_)','match');
                        end
                    end
                    if ischar(CurrLine)
                        CurrList=find(SList==str2num(CurrList{1}));
                        CurrCorr=regexp(CurrLine,'(?<=.*_c)\d+(?=_.*)','match');
                        CurrCorr=find(CurrLimits==str2num(CurrCorr{1}));
                        CurrNet=regexp(CurrLine,'(?<=.*_n)\d+(?=_c.*)','match');
                        CurrNet=find(NetRanks==str2num(CurrNet{1}));                         
                        Jury{CorrTypeL}{CurrList}{CurrNet}(1,CurrCorr)=JuryVal;                        
                    else
                        break
                    end
                end
            end
            fclose(fid)
        end


        Clu{CorrTypeL}={};
        CluNb{CorrTypeL}={};
        %process all probe set lists
        for ListL=1:ListNb
            ListPos=find(SList==ListRank(ListL));            
            Clu{CorrTypeL}{ListPos}={};
            CluNb{CorrTypeL}{ListPos}={};
            UsedLimit{CorrTypeL}{ListPos}={};
            %process all ranks
            for NetL=1:length(NetRanks)
                cd(fullfile(K.dir.net,sprintf('m%u',ChipRank),sprintf('n%u',NetRanks(NetL))))
                try
                    eval(sprintf('load t_m%u_n%u',ChipRank,NetRanks(NetL)))
                    %the current network has been made from another network by merging probe
                    %sets
                    MergeFlag=1;
                catch
                    MergeFlag=0;
                end
                cd(CmlDir)
                Clu{CorrTypeL}{ListPos}{NetL}=[];
                CluNb{CorrTypeL}{ListPos}{NetL}={};
                UsedLimit{CorrTypeL}{ListPos}{NetL}=[];
                for LimitL=1:length(CurrLimits)                    
                    MclFile=sprintf('mcl_i%u_m%u_n%u_c%u_%u%s.txt',Inflation,ChipRank,NetRanks(NetL),CurrLimits(LimitL),ListRank(ListL),Suffix);
                    if exist(MclFile,'file')
                        if ImportFlag
                            ErrorFlag=0;
                            CurrCluNb=[];
                            CurrClu=zeros(PsNb,1);
                            Rank=0;
                            fid=fopen(MclFile,'r');
                            while 1
                                Rank=Rank+1;
                                CurrPsRank=fgetl(fid);
                                if ~ischar(CurrPsRank), break, end
                                CurrPsRank=str2num(CurrPsRank);
                                if max(CurrPsRank)>PsNb
                                   % ErrorFlag=1;
                                    sprintf('error in %s at rank %u: exist PsRank>PsNb',MclFile,Rank)
                                   
                                else
                                    if MergeFlag
                                        for PsL=1:length(CurrPsRank)
                                            %put all probe sets eventually merged with
                                            %the current probe set in the current cluster
                                            PsRank=find(Trans(:,2)==CurrPsRank(PsL));
                                            CurrClu(PsRank)=Rank;
                                        end
                                    else
                                        CurrClu(CurrPsRank)=Rank;
                                    end
                                    CurrCluNb(Rank,1)=length(CurrPsRank);                                    
                                end
                            end
                            fclose(fid);
                            if ErrorFlag
                                Missing{CorrTypeL}(end+1,:)=[NetRanks(NetL),CurrLimits(LimitL),ListRank(ListL)];
                                MissingNb=MissingNb+1;
                            else
                                if ModifFlag                                    
                                    CurrClu=CurrClu(Net2Ps);
                                    if ~isempty(NetZeroPos)   
                                        CurrClu(NetZeroPos)=0;
                                    end
                                end
                                Clu{CorrTypeL}{ListPos}{NetL}(:,LimitL)=CurrClu;
                                CluNb{CorrTypeL}{ListPos}{NetL}{1,LimitL}=CurrCluNb;
                            end
                        end
                    else
                        sprintf('missing %s',MclFile)
                        Missing{CorrTypeL}(end+1,:)=[NetRanks(NetL),CurrLimits(LimitL),ListRank(ListL)];
                        MissingNb=MissingNb+1;
                    end
                end
            end
        end
    end
end

if MissingNb>0
    cd(CmlDir)
    h=warndlg(sprintf('miss %u files among %u',MissingNb,ResNb));
    waitfor(h)
    fid=fopen(sprintf('m%u_mcl_i2_%s.sh',ChipRank,FileName),'w');
    Suffixes={'','_diff','_raw','_diff_raw'};
    for CorrTypeL=1:4
        if ~isempty(Missing{CorrTypeL})
            for FileL=1:size(Missing{CorrTypeL},1)
                fprintf(fid,'/home/cinbell/mywork/MCL/bin/mcl m%u_n%u_c%u_MCL_%u%s.txt -I 2.0 --abc -o mcl_i2_m%u_n%u_c%u_%u%s.txt\n',ChipRank,Missing{CorrTypeL}(FileL,:),Suffixes{CorrTypeL},ChipRank,Missing{CorrTypeL}(FileL,:),Suffixes{CorrTypeL});
            end
        end
    end
    fclose(fid)
    fid=fopen(sprintf('exp_%s.sh',FileName),'w');
    for CorrTypeL=1:4
        if ~isempty(Missing{CorrTypeL})
            for FileL=1:size(Missing{CorrTypeL},1)
                switch CorrTypeL
                    case 1
                        fprintf(fid,'/home/cinbell/mywork/multiple/exportnet %u %u %u MCL %u 1.0 0 %u 0\n',ChipRank,Missing{CorrTypeL}(FileL,1),PsNb,Missing{CorrTypeL}(FileL,2),Missing{CorrTypeL}(FileL,3));
                    case 2
                        fprintf(fid,'/home/cinbell/mywork/multiple/exportnet %u %u %u MCL %u 1.0 1 %u 0\n',ChipRank,Missing{CorrTypeL}(FileL,1),PsNb,Missing{CorrTypeL}(FileL,2),Missing{CorrTypeL}(FileL,3));
                    case 3
                        fprintf(fid,'/home/cinbell/mywork/multiple/exportnet %u %u %u MCL %u 1.0 0 %u 1\n',ChipRank,Missing{CorrTypeL}(FileL,1),PsNb,Missing{CorrTypeL}(FileL,2),Missing{CorrTypeL}(FileL,3));
                    case 4
                        fprintf(fid,'/home/cinbell/mywork/multiple/exportnet %u %u %u MCL %u 1.0 1 %u 1\n',ChipRank,Missing{CorrTypeL}(FileL,1),PsNb,Missing{CorrTypeL}(FileL,2),Missing{CorrTypeL}(FileL,3));
                end
            end
        end
    end
    fclose(fid)
end


if ImportFlag
    cd(CmlDir)
    Info.netRanks=NetRanks;
    Info.listNb=ListNb;
    Info.listRank=ListRank;
    Info.type=[CorrFlag,DiffFlag,RawFlag,DiffRawFlag];
    Info.limits{1}=CorrLimits;
    Info.limits{2}=DiffLimits;
    Info.limits{3}=RawLimits;
    Info.limits{4}=DiffRawLimits;
    MatFile=sprintf('m%u_mcl_%s.mat',ChipRank,FileName);
    eval(sprintf('save %s Clu CluNb Info Jury',MatFile))
end
[temp SortIndex]=sort(MclJury.val);
MclJury.name=MclJury.name(SortIndex);
MclJury.val=MclJury.val(SortIndex);
cd(K.dir.common)
save mcljury MclJury
