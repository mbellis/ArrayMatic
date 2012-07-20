%=====================%
% FUNCTION NET_MCLBAT %
%=====================%

%NET_MCLBAT: write files containing commands to run MCL 
%INPUT PARAMETERS
% 1 CluType: either 'MCL' or 'TSN'
% 2 ChipRank: chip rank
% 3 NetRanks: list of net ranks
% 4  ListRank: rank of list used
% 5  Limits: CORR limits
% 6 DiffLimits: CORR-ANTI limits
% 7 RawLImits: raw(CORR) limits
% 8 DiffRawLimits: raw(CORR)-raw(ANTI) limits


% net_mclbat('MCL',27,[88:108],[1],[30:10:70],[20:10:60],[20:10:60],[10:10:50])
% net_mclbat('MCL',8,[39:53],[7:12],[20:10:70],[0:10:60],[10:10:60],[0:10:50])
% net_mclbat('MCL',27,[118],[1:6],[20:10:70],[0:10:60],[10:10:60],[0:10:50])
% net_mclbat('MCL',27,[118],[1],[20:10:70],[0:10:60],[10:10:60],[0:10:50])
% net_mclbat('MCL',8,[125],[7:12],[20:10:70],[0:10:60],[10:10:60],[0:10:50])
% net_mclbat('MCL',27,[118],[2],[20:10:70],[10:10:60],[10:10:60],[10:10:50])

% net_mclbat('TSN',27,[119:133],[1:6],[30:10:60],[20:10:50],[],[])
 
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
function net_mclbat(CluType,ChipRank,NetRanks,ListRank,Limits,DiffLimits,RawLimits,DiffRawLimits)
global K

NetNb=length(NetRanks);
ListNb=length(ListRank);
PsNb=K.chip.probesetNb(ChipRank);



if isequal(CluType,'MCL')
    CluDir=fullfile(K.dir.net,sprintf('m%u',ChipRank),'mcl');
elseif isequal(CluType,'TSN')
    CluDir=fullfile(K.dir.net,sprintf('m%u',ChipRank),'tsn');
end
try
    cd(CluDir)
catch
    mkdir(CluDir)
    cd(CluDir)
end
if isequal(CluType,'MCL')
    fid2=fopen('mcl.bat','w');
elseif isequal(CluType,'TSN')
    %write list of gene_ids
    cd(CluDir)
    if ~exist('gene_ids1','file')
        fid=fopen('gene_ids1','w');
        for PsL=1:PsNb
            fprintf(fid,'%u\n',PsL);
        end
        fclose(fid)
    end
end


fid1=fopen('exp.bat','w');

Suffixes={'','_diff','_raw','_diff_raw'};

for TypeL=1:4
    switch TypeL
        case 1
            CurrLimits=Limits;
            DiffFlag=0;
            RawFlag=0;
        case 2
            CurrLimits=DiffLimits;
            DiffFlag=1;
            RawFlag=0;
        case 3
            CurrLimits=RawLimits;
            DiffFlag=0;
            RawFlag=1;
        case 4
            CurrLimits=DiffRawLimits;
            DiffFlag=1;
            RawFlag=1;
    end
    if ~isempty(CurrLimits)
        LimitNb=length(CurrLimits);
        for ListL=1:ListNb
            CurrList=ListRank(ListL);            
            for LimitL=1:LimitNb
                CurrLimit=CurrLimits(LimitL);
                if isequal(CluType,'TSN')
                    fid2=fopen(sprintf('tsn_l%u_c%u%s.bat',CurrList,CurrLimit,Suffixes{TypeL}),'w');
                    fprintf(fid2,'/work/cinbell/TSN/NetsTensor --DatasetsListFile="netlist_l%u_c%u%s.txt" --ResultFile="m%u_n%uto%u_l%u_c%u_g10_n%u_d%u%s" --networksPath="." --resultPath="./result" --minGene=5 --minNet=%u --minDensity=%.2f --overlapPatternChoose="PATTERN_WITH_MORE_GENES"\n',...
                        CurrList,CurrLimit,Suffixes{TypeL},...
                        ChipRank,NetRanks(1),NetRanks(2),CurrList,CurrLimit,NetNb,CurrLimit+10,Suffixes{TypeL},...
                        NetNb,(CurrLimit+10)/100);
                    fclose(fid2);
                    fid3=fopen(sprintf('netlist_l%u_c%u%s.txt',CurrList,CurrLimit,Suffixes{TypeL}),'w');
                end
                for NetL=1:NetNb
                    CurrNet=NetRanks(NetL);
                    if isequal(CluType,'MCL')
                        fprintf(fid1,'/work/cinbell/multiple/exportnet %u %u %u MCL %u 1.0 %u %u %u\n',ChipRank,CurrNet,PsNb,CurrLimit,DiffFlag,CurrList,RawFlag);
                        fprintf(fid2,'/work/cinbell/MCL/bin/mcl m%u_n%u_c%u_MCL_%u%s.txt -I 2.0 --abc -o mcl_i2_m%u_n%u_c%u_%u%s.txt\n',ChipRank,CurrNet,CurrLimit,CurrList,Suffixes{TypeL},ChipRank,CurrNet,CurrLimit,CurrList,Suffixes{TypeL});
                    elseif isequal(CluType,'TSN')
                        fprintf(fid1,'/work/cinbell/multiple/exportnet %u %u %u TSN %u 1.0 %u %u %u\n',ChipRank,CurrNet,PsNb,CurrLimit,DiffFlag,CurrList,RawFlag);
                        fprintf(fid3,'m%u_n%u_c%u_TSN_%u%s\n',ChipRank,CurrNet,CurrLimit,CurrList,Suffixes{TypeL});
                    end
                end
                if isequal(CluType,'TSN')

                    fclose(fid3);
                end
            end
        end
    end
end
fclose(fid1);
if isequal(CluType,'MCL')
    fclose(fid2);
end




