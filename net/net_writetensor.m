%==========================%
% FUNCTION NET_WRITETENSOR %
%==========================%

%NET_WRITETENSOR: write files used by TensorNet program
%INPUT PARAMETERS
% 1    ChipRank: chip rank
% 2        List:  rank of list of probe set (uint32 file) used
%                 or mat file name pointing on a file
%                 containing either probe set ranks (Index variable)
%                 or a cell structure such that GenePos{ChipRank}{GenePos} contains all the
%                 the probe set rank targeting this gene
%                 in this last case ExportFlag=0 (exportation by this program) and DistalFlag =0
% 3    UsedPsNb: number of probeset to be exported (=0 => ALL)
% 4    NetRanks: networks used by TensorNet
% 5        Type: Type of calcul: either 1(C) or 2(C-A)
% 6       Limit: inferior limit of C or C-A
% 7    GeneNbs: minimum number of gene in a TSN cluster
% 8     NetNbs: minimum of networks where a TSN cluster must be
% 9  Densities: densities used by tensornet
% 10 ExportFlag: if ==0 => exportation by Matlab (very slow), else exportation by ExportNet(c++)
% 11 DistalFlag: uses another station with different folder


%all probeset targeting one gene (no MCL clusters) (do not restrict m8 to common m27 probesets)
%NetsTensor on three networks
% net_writetensor(8,5,0,[212:214],2,0,5,3,[35:5:60],1,1);
% net_writetensor(8,5,0,[215:217],2,0,5,3,[35:5:60],1,1);
% net_writetensor(8,5,0,[218:221],2,0,5,3,[35:5:60],1,1);
% net_writetensor(8,5,0,[222:223],2,0,5,3,[35:5:60],1,1);
% net_writetensor(8,5,0,[224:226],2,0,5,3,[35:5:60],1,1);
% net_writetensor(27,5,0,[149:151],2,0,5,3,[35:5:60],1,1);
% net_writetensor(27,5,0,[152:154],2,0,5,3,[35:5:60],1,1);
% net_writetensor(27,5,0,[155:157],2,0,5,3,[35:5:60],1,1);
% net_writetensor(27,5,0,[158:160],2,0,5,3,[35:5:60],1,1);
% net_writetensor(27,5,0,[161:163],2,0,5,3,[35:5:60],1,1);

% %NetsTensor on six networks
% net_writetensor(8,5,0,[212:217],2,0,5,3,[30:10:60],1,1);
% net_writetensor(8,5,0,[218:223],2,0,5,3,[30:10:60],1,1);
% net_writetensor(27,5,0,[149:154],2,0,5,3,[30:10:60],1,1);
% net_writetensor(27,5,0,[155:162],2,0,5,3,[30:10:60],1,1);
% %NetsTensor on 15 networks
% net_writetensor(8,5,0,[212:226],2,0,5,5,[30:10:60],1,1);
% net_writetensor(27,5,0,[149:163],2,0,5,5,[30:10:60],1,1);
% see effect of merging multiple ps
% net_writetensor(27,13,0,[149:163],2,0,5,5,[30],1,1);
% net_writetensor(27,26,0,[184:6:268],2,0,5,5,[30],1,1);
% net_writetensor(27,39,0,[185:6:269],2,0,5,5,[30],1,1);
% net_writetensor(27,52,0,[186:6:270],2,0,5,5,[30],1,1);
% net_writetensor(27,65,0,[187:6:271],2,0,5,5,[30],1,1);
% net_writetensor(27,78,0,[188:6:272],2,0,5,5,[30],1,1);
% net_writetensor(27,91,0,[189:6:273],2,0,5,5,[30],1,1);



%comparaison m8/m27 : sortir m8 dans le même ordre que m27
%net_writetensor(8,17,0,[228],2,0,5,4,[10:5:50],1,1);
%net_writetensor(27,5,0,[164],2,0,5,4,[10:5:50],1,1);
%fitted m27 n164
%net_writetensor(27,5,0,[171],2,0,5,4,[10:5:50],1,1);

%Etude petite listes
%net_writetensor(8,{'crel_2_gene2psrank'},0,[212:226],2,0,3,5,[30:5:70],0,0);
%net_writetensor(27,{'crel_2_gene2psrank'},0,[149:163],2,0,3,5,[30:5:70],0,0);
%net_writetensor(8,{'mouse_krebs_proteasome_mapk_gene2psrank'},0,[212:226],2,0,3,5,[30:5:70],0,0);
%net_writetensor(27,{'mouse_krebs_proteasome_mapk_gene2psrank'},0,[149:163],2,0,3,5,[30:5:70],0,0);
%net_writetensor(27,{'mouse_krebs_proteasome_mapk_gene2psrank'},0,[134:148],2,0,3,5,[30:5:70],0,0);

%net_writetensor(27,{'m27_freqrank50_sup','m27_freqrank50_inf'},0,[149:163],2,0,3,5,[30:5:70],0,0);
%net_writetensor(27,[92,93],0,[149:163],2,0,5,10,[30:5:70],1,0);
%net_writetensor(8,[185,186],0,[212:226],2,0,5,10,[30:5:70],1,0);

% Etude réseaux fusionnés
%   net_writetensor(91,1,0,[1:3],2,0,5,3,[30:5:70],1,1);
%   net_writetensor(92,1,0,[1:3],2,0,5,3,[30:5:70],1,1);
%   net_writetensor(92,2,0,[1:3],2,0,5,3,[30:5:70],1,1);
%   net_writetensor(93,1,0,[1:3],2,0,5,3,[30:5:70],1,1);
%   net_writetensor(94,1,0,[1:3],2,0,5,3,[30:5:70],1,1);

% net_writetensor(91,1,0,[4:18],2,0,5,10,[30:5:70],1,1);
% net_writetensor(91,1,0,[19:34],2,0,5,10,[30:5:70],1,1);
% net_writetensor(92,1,0,[1:3],2,0,5,3,[30:5:70],1,1);
% net_writetensor(92,1,0,[5:20],2,0,5,10,[30:5:70],1,1);
% net_writetensor(92,1,0,[21:35],2,0,5,10,[30:5:70],1,1);
% net_writetensor(92,1,0,[36:50],2,0,5,10,[30:5:70],1,1);
% net_writetensor(93,1,0,[1:3],2,0,5,3,[20,25],1,1);
% net_writetensor(94,1,0,[1:3],2,0,5,3,[15:5:25],1,1);
% net_writetensor(92,3,0,[3],2,0,5,10,[30:5:70],1,1);
% net_writetensor(92,2,0,[2],2,0,5,10,[30:5:70],1,0);

% net_writetensor(95,1,0,[1:4],2,0,5,4,[30:5:70],1,1);

% %NetsTensor on 15 networks 
% net_writetensor(2,[24,26],0,[64:78],2,0,5,10,[30:5:60],1,1);
% net_writetensor(3,[11,13],0,[70:85],2,0,5,10,[30:5:60],1,1);
% net_writetensor(5,[11,13],0,[107:122],2,0,5,10,[30:5:60],1,1);
% net_writetensor(6,[11,13],0,[48:62],2,0,5,10,[30:5:60],1,1);
% net_writetensor(8,[11,13],0,[212:226],2,0,5,10,[30:5:60],1,1);
% net_writetensor(27,[11,13],0,[149:163],2,0,5,10,[30:5:60],1,1);
% net_writetensor(27,[11],0,[149:163],2,0,5,10,[35:5:60],1,1);


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

function net_writetensor(ChipRank,List,UsedPsNb,NetRanks,Type,Limit,GeneNbs,NetNbs,Densities,ExportFlag,DistalFlag)
global K

ChipPos=find(K.chip.rank==ChipRank);
ChipPsNb=K.chip.probesetNb(ChipPos);
NetNb=length(NetRanks);
%make dir for TensorNet
TensorDir=fullfile(K.dir.net,sprintf('m%u',ChipRank),'tsn');
mkdir(TensorDir)
%write files used by tensor net
%gene ids
cd(TensorDir)
fid=fopen('gene_ids','w');
for PsL=1:ChipPsNb
    fprintf(fid,'%u\n',PsL);
end
fclose(fid)

if Type==2
    Prefix=sprintf('m%u_n%uton%u_c%u_TSN_diff',ChipRank,NetRanks(1),NetRanks(end),Limit);
else
    Prefix=sprintf('m%u_n%uton%u_c%u_TSN',ChipRank,NetRanks(1),NetRanks(end),Limit);
end
%write tensortnet command
CommandFile=sprintf('%s.sh',Prefix);
fid1=fopen(CommandFile,'w');
%write exportnet commands
if ExportFlag
    ExportFile=sprintf('%s_export.sh',Prefix);
    fid3=fopen(ExportFile,'w');
end
for ListL=1:length(List)
    if iscell(List)
        CurrFile=List{ListL};
        cd(K.dir.list)
        load([CurrFile,'.mat'])
        if  exist('GenePos')==1
            CurrIndex=GenePos{ChipRank};
            %construct Index
            Index=[];
            %GeneList={};
            for GeneL=1:length(CurrIndex)
                Index=[Index;CurrIndex{GeneL}];
            end
            [Index,SortIndex]=unique(Index);
        else
            %index loaded directly
        end
    else
        ListRank=List(ListL);
        if ListRank==0
            Index=1:ChipPsNb;
        else
            cd(fullfile(K.dir.net,sprintf('m%u',ChipRank),'list'));
            Index=load_data(sprintf('m%u_pslist%u.u32',ChipRank,ListRank),'./',0,0,'uint32','ieee-le');
        end
    end
    if UsedPsNb>0&UsedPsNb<length(Index)
        Index=Index(1:UsedPsNb);
    end
    %load data if ExportFlag==0
    if ExportFlag==0
        NodeNb=0;
        ColNb=K.chip.probesetNb(ChipPos);
        LineNb=ColNb;
        for NetL=1:length(NetRanks)
            NetRank=NetRanks(NetL);
            cd(fullfile(K.dir.net,sprintf('m%u',ChipRank),sprintf('n%u',NetRank)));
            c=load_data(sprintf('c_m%u_n%u.4mat',ChipRank,NetRank),'./',LineNb,ColNb,'uint8','ieee-le',Index,Index);
            if Type==2
                a=load_data(sprintf('a_m%u_n%u.4mat',ChipRank,NetRank),'./',LineNb,ColNb,'uint8','ieee-le',Index,Index);
                c=c-a;
            end
            i=find(c<=Limit);
            c(i)=0;
            NodeNb=NodeNb+length(find(c));
            cd(TensorDir)
            if iscell(List)
                if Type==2
                    fid=fopen(sprintf('m%u_n%u_c%u_TSN_%s_diff.txt',ChipRank,NetRank,Limit,CurrFile),'w');
                else
                    fid=fopen(sprintf('m%u_n%u_c%u_TSN_%s.txt',ChipRank,NetRank,Limit,CurrFile),'w');
                end
            else
                if Type==2
                    fid=fopen(sprintf('m%u_n%u_c%u_TSN_%u_diff.txt',ChipRank,NetRank,Limit,ListRank),'w');
                else
                    fid=fopen(sprintf('m%u_n%u_c%u_TSN_%u.txt',ChipRank,NetRank,Limit,ListRank),'w');
                end
            end
            for LineL=1:length(c)
                Pos=find(c(LineL,LineL+1:end));
                if ~isempty(Pos)
                    for PosL=1:length(Pos)
                        %index
                        fprintf(fid,'%u\t%u\t%.2f\n',LineL,Pos(PosL)+LineL,single(c(LineL,LineL+Pos(PosL)))/100);
                    end
                end
            end
            fclose(fid)
        end
    end
    cd(TensorDir)
    %write list of networks
    if iscell(List)
        NetFile=sprintf('%s_%s.txt',Prefix,CurrFile);
    else
        NetFile=sprintf('%s_%u.txt',Prefix,ListRank);
    end
    fid2=fopen(NetFile,'w');
    %write list of networks
    for NetL=1:NetNb
        CurrNet=NetRanks(NetL);
        if iscell(List)
            if Type==1
                fprintf(fid2,'m%u_n%u_c%u_TSN_%s\n',ChipRank,CurrNet,Limit,CurrFile);
            else
                fprintf(fid2,'m%u_n%u_c%u_TSN_%s_diff\n',ChipRank,CurrNet,Limit,CurrFile);
            end

        else
            if Type==1
                fprintf(fid2,'m%u_n%u_c%u_TSN_%u\n',ChipRank,CurrNet,Limit,ListRank);
            else
                fprintf(fid2,'m%u_n%u_c%u_TSN_%u_diff\n',ChipRank,CurrNet,Limit,ListRank);
            end
        end
    end
    fclose(fid2)
    %write tensortnet command
    for GeneNbL=1:length(GeneNbs)
        for NetNbL=1:length(NetNbs)
            if iscell(List)
                for DensL=1:length(Densities)
                    if DistalFlag==0
                        fprintf(fid1,'/usr/local/netstensor/netstensor --DatasetsListFile="%s" --ResultFile="%s_%s_g%u_n%u_d%u.txt" --networksPath="." --resultPath="./result" --minGene=%u --minNet=%u --minDensity=%.2f --suffixdataFile=".txt" --overlapPatternChoose="PATTERN_WITH_MORE_GENES"\n',...
                            NetFile,Prefix,CurrFile,GeneNbs(GeneNbL),NetNbs(NetNbL),Densities(DensL),GeneNbs(GeneNbL),NetNbs(NetNbL),Densities(DensL)/100);
                    else
                        fprintf(fid1,'/home/cinbell/mywork/TSN/NetsTensor --DatasetsListFile="%s" --ResultFile="%s_%s_g%u_n%u_d%u.txt" --networksPath="." --resultPath="./result" --minGene=%u --minNet=%u --minDensity=%.2f --suffixdataFile=".txt" --overlapPatternChoose="PATTERN_WITH_MORE_GENES"\n',...
                            NetFile,Prefix,CurrFile,GeneNbs(GeneNbL),NetNbs(NetNbL),Densities(DensL),GeneNbs(GeneNbL),NetNbs(NetNbL),Densities(DensL)/100);
                    end
                end

            else
                for DensL=1:length(Densities)
                    if DistalFlag==0
                        fprintf(fid1,'/usr/local/netstensor/netstensor --DatasetsListFile="%s" --ResultFile="%s_%u_g%u_n%u_d%u.txt" --networksPath="." --resultPath="./result" --minGene=%u --minNet=%u --minDensity=%.2f --suffixdataFile=".txt" --overlapPatternChoose="PATTERN_WITH_MORE_GENES"\n',...
                            NetFile,Prefix,ListRank,GeneNbs(GeneNbL),NetNbs(NetNbL),Densities(DensL),GeneNbs(GeneNbL),NetNbs(NetNbL),Densities(DensL)/100);
                    else
                        fprintf(fid1,'/home/cinbell/mywork/TSN/NetsTensor --DatasetsListFile="%s" --ResultFile="%s_%u_g%u_n%u_d%u.txt" --networksPath="." --resultPath="./result" --minGene=%u --minNet=%u --minDensity=%.2f --suffixdataFile=".txt" --overlapPatternChoose="PATTERN_WITH_MORE_GENES"\n',...
                            NetFile,Prefix,ListRank,GeneNbs(GeneNbL),NetNbs(NetNbL),Densities(DensL),GeneNbs(GeneNbL),NetNbs(NetNbL),Densities(DensL)/100);
                    end
                end
            end
        end
    end
    if iscell(List)==0
        if ExportFlag
            %write exportnet commands
            for NetL=1:NetNb
                CurrNet=NetRanks(NetL);
                cd(fullfile(K.dir.net,sprintf('m%u',ChipRank),sprintf('n%u',CurrNet)))
                try
                    eval(sprintf('load t_m%u_n%u',ChipRank,CurrNet))
                    PsNb=max(Trans(:,2));
                catch
                    PsNb=ChipPsNb;
                end
                if Type==1
                    % exportnet needs modelRank netRank psNb  exportType cLimit cRatio  diffFlag listRank rawFlag distalFlag
                    if DistalFlag
                        fprintf(fid3,'/home/cinbell/mywork/multiple/exportnet %u %u %u TSN %u 1.0 0 %u %u %u\n',ChipRank,CurrNet,PsNb,Limit,ListRank,0,1);
                    else
                        fprintf(fid3,'/home/mbellis/sosma/projet/C/exportnet/bin/Debug/exportnet %u %u %u TSN %u 1.0 0 %u %u %u\n',ChipRank,CurrNet,PsNb,Limit,ListRank,0,0);
                    end
                else
                    if DistalFlag
                        fprintf(fid3,'/home/cinbell/mywork/multiple/exportnet %u %u %u TSN %u 1.0 1 %u %u %u\n',ChipRank,CurrNet,PsNb,Limit,ListRank,0,1);
                    else
                        fprintf(fid3,'/home/mbellis/sosma/projet/C/exportnet/bin/Debug/exportnet %u %u %u TSN %u 1.0 1 %u %u %u\n',ChipRank,CurrNet,PsNb,Limit,ListRank,0,0);
                    end
                end
            end
        end
    end

end
fclose(fid1)
if ExportFlag
    fclose(fid3);
end
%write pbs command file
if DistalFlag
    cd(TensorDir)

    fid=fopen(sprintf('%s_export.pbs',Prefix),'w');
    fprintf(fid,'#PBS -S /bin/bash\n');
    fprintf(fid,'#PBS -N tsn\n');
    fprintf(fid,sprintf('#PBS -e %s.eerr\n',Prefix));
    fprintf(fid,sprintf('#PBS -o %s.elog\n',Prefix));
    fprintf(fid,'#PBS -l walltime=2:00:00\n');
    fprintf(fid,'#PBS -l select=1:ncpus=8:mpiprocs=9\n');
    fprintf(fid,'module unload intel\n');
    fprintf(fid,'module load intel/11.1.072\n');
    fprintf(fid,'module load pserie/9.3.26\n');
    %fprintf(fid,'export LD_PRELOAD=/work/SGI/misc/getcwd_wrapper.so\n');
    fprintf(fid,'export MPI_DSM_DISTRIBUTE=ON\n');
    fprintf(fid,'cd /scratch/cinbell/tsn\n',ChipRank);
    fprintf(fid,'mpiexec pserie_load_balanced < %s\n',ExportFile);
    fclose(fid);





    fid=fopen(sprintf('%s.pbs',Prefix),'w');
    fprintf(fid,'#PBS -S /bin/bash\n');
    fprintf(fid,'#PBS -N tsn\n');
    fprintf(fid,sprintf('#PBS -e %s.terr\n',Prefix));
    fprintf(fid,sprintf('#PBS -o %s.tlog\n',Prefix));
    fprintf(fid,'#PBS -l walltime=24:00:00\n');
    fprintf(fid,'#PBS -l select=1:ncpus=8:mpiprocs=9\n');
    fprintf(fid,'module unload intel\n');
    fprintf(fid,'module load intel/11.1.072\n');
    fprintf(fid,'module load pserie/9.3.26\n');
    %fprintf(fid,'export LD_PRELOAD=/work/SGI/misc/getcwd_wrapper.so\n');
    fprintf(fid,'export MPI_DSM_DISTRIBUTE=ON\n');
    fprintf(fid,'cd /scratch/cinbell/tsn\n',ChipRank);
    fprintf(fid,'mpiexec pserie_load_balanced < %s\n',CommandFile);
    fclose(fid);
end


