%=====================%
% FUNCTION NET_MCLBAT %
%=====================%

%NET_MCLBAT: write files containing commands to run MCL 
%INPUT PARAMETERS
% 1 ChipRank: chip rank
% 2 NetRanks: list of net ranks
% 3  ListRanks: rank of list used
% 4  Limits: CORR limits
% 5 DiffLimits: CORR-ANTI limits
% 6 RawLimits: raw(CORR) limits
% 7 DiffRawLimits: raw(CORR)-raw(ANTI) limits
% 8 FileName: used in .bat, .sh and .pbs files
% 9 DistalFlag: 0=>write commands for local station 
%               1=> write for CINES station


% net_mclbat(27,[88:108],[1],[30:10:70],[20:10:60],[20:10:60],[10:10:50])
% net_mclbat(8,[39:53],[7:12],[20:10:70],[0:10:60],[10:10:60],[0:10:50])
% net_mclbat(27,[118],[1:6],[20:10:70],[0:10:60],[10:10:60],[0:10:50])
% net_mclbat(27,[118],[1],[20:10:70],[0:10:60],[10:10:60],[0:10:50])
% net_mclbat(8,[227],[7:12],[20:10:70],[10:10:60],[10:10:60],[0:10:50],'n227')
% net_mclbat(27,[164],[1,3:6],[20:10:70],[10:10:60],[10:10:60],[10:10:50],'n164')
% net_mclbat(8,[228],[7:12],[20:2:30],[20:2:30],[],[],'n228A')
% net_mclbat(27,[164],[1:6],[20:2:30],[20:2:30],[],[],'n164A')
% net_mclbat(8,[228],[7:12],[12:2:18],[12:2:18],[],[],'n228B')
% net_mclbat(27,[164],[1:6],[12:2:18],[12:2:18],[],[],'n164B')
% net_mclbat(27,[164],[1:6],[12:2:18],[12:2:18],[],[],'n164B')

% net_mclbat(8,[230],[201,205],[],[10:2:20],[],[],'n230A',1)
% net_mclbat(8,[232],[201,205],[],[10:2:20],[],[],'n232A',1)
% net_mclbat(27,[169],[201,205],[],[10:2:20],[],[],'n169A',1)
% net_mclbat(27,[170],[201,205],[],[10:2:20],[],[],'n170A',1)


%  net_mclbat(8,[230],[207,211],[],[14],[],[],'n230B',0)
%  net_mclbat(8,[232],[207,211],[],[14],[],[],'n232B',0)
%  net_mclbat(8,[230],[201,205],[],[20],[],[],'n230C',0)
%  net_mclbat(8,[232],[201,205],[],[20],[],[],'n232C',0)
%  net_mclbat(8,[212],[207,211],[],[14],[],[],'n212A',0)
%  net_mclbat(8,[213],[207,211],[],[14],[],[],'n213A',0)
%  net_mclbat(27,[149],[201,205],[],[14],[],[],'n149A',0)
%  net_mclbat(27,[150],[201,205],[],[14],[],[],'n150A',0)

% 2012/11/02 : new ps lists
% net_mclbat(8,[55,228],[14:26],[],[14],[],[],'n55n228A',1)
% 
% net_mclbat(8,[234],[41:52],[],[14],[],[],'n234A',1)
% net_mclbat(8,[235],[67:78],[],[14],[],[],'n235A',1)
% net_mclbat(8,[236],[93:104],[],[14],[],[],'n236A',1)
% net_mclbat(8,[237],[119:130],[],[14],[],[],'n237A',1)
% net_mclbat(8,[238],[145:156],[],[14],[],[],'n238A',1)
% net_mclbat(8,[239],[171:182],[],[14],[],[],'n239A',1)
% net_mclbat(8,[240],[41:52],[],[14],[],[],'n240A',1)
% net_mclbat(8,[241],[67:78],[],[14],[],[],'n241A',1)
% net_mclbat(8,[242],[93:104],[],[14],[],[],'n242A',1)
% net_mclbat(8,[243],[119:130],[],[14],[],[],'n243A',1)
% net_mclbat(8,[244],[145:156],[],[14],[],[],'n244A',1)
% net_mclbat(8,[245],[171:182],[],[14],[],[],'n245A',1)
% 
% net_mclbat(27,[24,164],[1:13],[],[14],[],[],'n24n164A',1)
% 
% net_mclbat(27,[172],[15:26],[],[14],[],[],'172A',1)
% net_mclbat(27,[173],[28:39],[],[14],[],[],'173A',1)
% net_mclbat(27,[174],[41:52],[],[14],[],[],'174A',1)
% net_mclbat(27,[175],[54:65],[],[14],[],[],'175A',1)
% net_mclbat(27,[176],[67:78],[],[14],[],[],'176A',1)
% net_mclbat(27,[177],[80:91],[],[14],[],[],'177A',1)
% net_mclbat(27,[178],[15:26],[],[14],[],[],'178A',1)
% net_mclbat(27,[179],[28:39],[],[14],[],[],'179A',1)
% net_mclbat(27,[180],[41:52],[],[14],[],[],'180A',1)
% net_mclbat(27,[181],[54:65],[],[14],[],[],'181A',1)
% net_mclbat(27,[182],[67:78],[],[14],[],[],'182A',1)
% net_mclbat(27,[183],[80:91],[],[14],[],[],'183A',1)
%



% for i=7:21
% eval(sprintf('net_mclbat(8,[%u],[9,10,20,21,3,5,11,1,14,6,22,16,7,8,17,12,18,19,2,13,4,15],[],[14],[],[],''n%uA'',1)',i,i))
% end
% for i=212:226
% eval(sprintf('net_mclbat(8,[%u],[9,10,20,21,3,5,11,1,14,6,22,16,7,8,17,12,18,19,2,13,4,15],[],[14],[],[],''n%uA'',1)',i,i))
% end
% for i=2:22
% eval(sprintf('net_mclbat(27,[%u],[9,10,3,11,5,1,7,6,2,8,4],[],[14],[],[],''n%uA'',1)',i,i))
% end
% for i=149:163
% eval(sprintf('net_mclbat(27,[%u],[9,10,3,11,5,1,7,6,2,8,4],[],[14],[],[],''n%uA'',1)',i,i))
% end

% net_mclbat(8,[7:21],[9,10,20,21,3,5,11,1,14,6,22,16,7,8,17,12,18,19,2,13,4,15],[],[14],[],[],'n7n21',1)
% net_mclbat(8,[212:226],[9,10,20,21,3,5,11,1,14,6,22,16,7,8,17,12,18,19,2,13,4,15],[],[14],[],[],'n212n226',1)
% net_mclbat(27,[2:22],[9,10,3,11,5,1,7,6,2,8,4],[],[14],[],[],'n2n22',1)
% net_mclbat(27,[149:163],[9,10,3,11,5,1,7,6,2,8,4],[],[14],[],[],'n149n163',1)

%correction liste (exp)

% net_mclbat(8,[7:21],[5,6,16,17],[],[14],[],[],'n7n21B',1)
% net_mclbat(8,[212:226],[5,6,16,17],[],[14],[],[],'n212n226B',1)
% net_mclbat(27,[2:22],[5,6],[],[14],[],[],'n2n22B',1)
% net_mclbat(27,[149:163],[5,6],[],[14],[],[],'n149n163B',1)

% ajout inflation
% net_mclbat(27,[24],[1:13],[],[14],[],[],[3,4,5],'24B',1)
% net_mclbat(8,[55],[14:26],[],[14],[],[],[3,4,5],'55B',1)
% net_mclbat(27,[164],[1:13],[],[14],[],[],[3,4,5],'164B',1)
% net_mclbat(8,[228],[14:26],[],[14],[],[],[3,4,5],'228B',1)
% net_mclbat(8,[55],[1:13],[],[14],[],[],[2,3,4,5],'55C',1)
% net_mclbat(8,[228],[1:13],[],[14],[],[],[2,3,4,5],'228C',1)

%suppression différent type corr , ajout SupLimit
% net_mclbat(27,[164],[1:13],[14],[60],[2],'164D',1)
% net_mclbat(8,[228],[14:26],[14],[60],[2],'228D',1)
% net_mclbat(8,[228],[1:13],[14],[60],[2],'228E',1)

%   net_mclbat(91,[1,2,3],1,10,2,'m91a',1)
%   net_mclbat(92,[1,2,3,4],1,10,2,'m92a',1)
%   net_mclbat(93,[1,2,3],1,10,2,'m93a',1)
%   net_mclbat(94,[1,2,3,4],1,10,2,'m94a',1)

% net_mclbat(2,[80],[24,26],[0:10:70],2,'n80A',1)
% net_mclbat(2,[80],[24],[0:10:70],2,'n80B',1)
% net_mclbat(3,[86],[11,13],[0:10:70],2,'n86A',1)
% net_mclbat(5,[123],[11,13],[0:10:70],2,'n123A',1)
% net_mclbat(6,[63],[11,13],[0:10:70],2,'n63A',1)
% net_mclbat(8,[228],[11,13],[0:10:70],2,'n228F',1)
% net_mclbat(27,[164],[11,13],[0:10:70],2,'n164F',1)

% net_mclbat(2,[64:78],[24,26],[0:10:70],2,'n64n78A',1)
% net_mclbat(3,[70:85],[11,13],[0:10:70],2,'n70n85A',1)
% net_mclbat(5,[107:122],[11,13],[0:10:70],2,'n107n122A',1)
% net_mclbat(6,[48:62],[11,13],[0:10:70],2,'n48n62A',1)
% net_mclbat(8,[212:226],[11,13],[0:10:70],2,'n212n226A',1)
% net_mclbat(27,[149:163],[11,13],[0:10:70],2,'n149n163A',1)
% net_mclbat(8,[221],[13],[10],2,'n221A',0)





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
function net_mclbat(ChipRank,NetRanks,ListRanks,InfLimits,Inflation,FileName,DistalFlag)
global K

NetNb=length(NetRanks);
ListNb=length(ListRanks);
ChipPsNb=K.chip.probesetNb(ChipRank);
DiffFlag=1;
RawFlag=0;
Suffix='diff';



CluDir=fullfile(K.dir.net,sprintf('m%u',ChipRank),'mcl');
try
    cd(CluDir)
catch
    mkdir(CluDir)
    cd(CluDir)
end
if DistalFlag
    ExportProgDir='/home/cinbell/mywork/multiple/exportnet';
    ProgDir='/home/cinbell/mywork/MCL/bin/mcl';    
else
    ExportProgDir='/home/mbellis/sosma/projet/C/exportnet/bin/Debug/exportnet';    
    ProgDir='/usr/bin/mcl';    
end

fid2=zeros(length(Inflation),1);
for InfL=1:length(Inflation)
    fid2(InfL)=fopen(sprintf('m%u_mcl_i%u_%s.sh',ChipRank,round(Inflation(InfL)),FileName),'w');
    if DistalFlag
        fid3=fopen(sprintf('m%u_mcl_i%u_%s.pbs',ChipRank,round(Inflation(InfL)),FileName),'w');
        fprintf(fid3,'#PBS -S /bin/bash\n');
        fprintf(fid3,'#PBS -N mcl_i%u_%s\n',round(Inflation(InfL)),FileName);
        fprintf(fid3,'#PBS -e mcl_i%u_%s.mclerr\n',round(Inflation(InfL)),FileName);
        fprintf(fid3,'#PBS -o mcl_i%u_%s.mcllog\n',round(Inflation(InfL)),FileName);
        fprintf(fid3,'#PBS -l walltime=2:00:00\n');
        fprintf(fid3,'#PBS -l select=1:ncpus=8:mpiprocs=9\n');
        fprintf(fid3,'module unload intel\n');
        fprintf(fid3,'module load intel/11.1.072\n');
        fprintf(fid3,'module load numpy/1.5.1\n');
        fprintf(fid3,'module load scipy/0.8.0\n');
        fprintf(fid3,'module load pserie/9.3.26\n');
        %fprintf(fid3,'export LD_PRELOAD=/work/SGI/misc/getcwd_wrapper.so\n');
        fprintf(fid3,'export MPI_DSM_DISTRIBUTE=ON\n');
        fprintf(fid3,'cd /scratch/cinbell/mcl\n');
        fprintf(fid3,'mpiexec pserie_load_balanced < m%u_mcl_i%u_%s.sh\n',ChipRank,round(Inflation(InfL)),FileName);
        fclose(fid3);
    end
end


fid1=fopen(sprintf('m%u_exp_%s.sh',ChipRank,FileName),'w');



if ~isempty(InfLimits)
    LimitNb=length(InfLimits);
    for ListL=1:ListNb
        CurrList=ListRanks(ListL);
        for LimitL=1:LimitNb
            CurrInfLimit=InfLimits(LimitL);
%             if ~isempty(SupLimits)
%                 CurrSupLimit=SupLimits(LimitL);
%             else
%                 CurrSupLimit=100;
%             end
            for NetL=1:NetNb
                CurrNet=NetRanks(NetL);
                cd(fullfile(K.dir.net,sprintf('m%u',ChipRank),sprintf('n%u',CurrNet)))
                try
                    eval(sprintf('load t_m%u_n%u',ChipRank,CurrNet))
                    PsNb=max(Trans(:,2));
                catch
                    PsNb=ChipPsNb;
                end
                fprintf(fid1,'%s %u %u %u MCL %u 1.0 %u %u %u %u\n',ExportProgDir,ChipRank,CurrNet,PsNb,CurrInfLimit,DiffFlag,CurrList,RawFlag,DistalFlag);
                for InfL=1:length(Inflation)
                    fprintf(fid2(InfL),'%s m%u_n%u_c%u_MCL_%u_%s.txt -I %.1f --abc -o mcl_i%u_m%u_n%u_c%u_%u_%s.txt\n',ProgDir,ChipRank,CurrNet,CurrInfLimit,CurrList,Suffix,Inflation(InfL),Inflation(InfL),ChipRank,CurrNet,CurrInfLimit,CurrList,Suffix);
                end
            end
        end
    end
end


fclose(fid1);
for InfL=1:length(Inflation)
    fclose(fid2(InfL));
end




