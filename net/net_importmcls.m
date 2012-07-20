%=========================%
% FUNCTION NET_IMPORTMCLS %
%=========================%

%NET_IMPORTMCLS: import several results of MCL clustering
%INPUT PARAMETERS
% 1 ChipRank: chip rank
% 2  NetRanks: list of net ranks
% 3 Inflation: MCL paramter
% 4  ListRank: rank of list used
% 5  DiffFlag: corr-anti values used in MCL
% 6   RawFlag: raw values (corr*freq, anti*frq) used in MCL
% 7 DiffRawFlag: raw(corr _ raw(anti)
% 8  Limits:
% 9 DiffLimits
% 10 RawLImits
% 11 DiffRawLimits
% 12 OutRank

% net_importmcls(8,[24:38],2,[1:6],1,1,1,[20:10:70],[0:10:60],[10:10:60],[0:10:50],2)
% net_importmcls(8,[24:38],2,[1:12],1,1,1,[20:10:70],[0:10:60],[10:10:60],[0:10:50],1)
% net_importmcls(27,[25:45],2,[1:6],1,1,1,[20:10:70],[0:10:60],[10:10:60],[0:10:50],1)
% net_importmcls(27,[88:108],2,[1:6],1,1,1,[20:10:70],[0:10:60],[10:10:60],[0:10:50],2)
% net_importmcls(27,[88:108],2,[1:6],1,1,1,[20:10:70],[0:10:60],[10:10:60],[0:10:50],2,0)
% net_importmcls(27,[88:108],2,[1:6],1,1,1,[20:10:70],[0:10:60],[10:10:60],[0:10:50],2,1)
% net_importmcls(8,[125],2,[1:12],1,1,1,[20:10:70],[0:10:60],[10:10:60],[0:10:50],3,0)
% net_importmcls(8,[125],2,[1:12],1,1,1,[20:10:70],[10:10:60],[10:10:60],[10:10:50],3,1)
% net_importmcls(8,[125],2,[7:12],1,1,1,[20:10:70],[20:10:60],[20:10:60],[10:10:50],4,1)
% net_importmcls(27,[118],2,[1:6],1,1,1,[20:10:70],[10:10:60],[10:10:60],[10:10:50],3,1)

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


function net_importmcls(ChipRank,NetRanks,Inflation,ListRank,DiffFlag,RawFlag,DiffRawFlag,Limits,DiffLimits,RawLimits,DiffRawLimits,OutRank,ImportFlag)
global K

ListNb=length(ListRank);
PsNb=K.chip.probesetNb(ChipRank);
CmlDir=fullfile(K.dir.net,sprintf('m%u',ChipRank),'mcl');
cd(CmlDir)

ResNb=0;
for TypeL=1:4
    DoIt=0;
    switch TypeL
        case 1
            CurrLimits=Limits;

        case 2
            CurrLimits=DiffLimits;

        case 3
            CurrLimits=RawLimits;
        case 4
            CurrLimits=DiffRawLimits;
    end
    ResNb=ResNb+length(CurrLimits)*length(NetRanks)*ListNb;
end

MissingNb=0;
for TypeL=1:4
    DoIt=0;
    switch TypeL
        case 1
            DoIt=1;
            Suffix='';
            CurrLimits=Limits;

        case 2
            DoIt=DiffFlag;
            Suffix='_diff';
            CurrLimits=DiffLimits;

        case 3
            DoIt=RawFlag;
            Suffix='_raw';
            CurrLimits=RawLimits;
        case 4
            DoIt=DiffRawFlag;
            Suffix='_diff_raw';
            CurrLimits=DiffRawLimits;
    end
    Missing{TypeL}=[];

    if DoIt
        Clu{TypeL}={};
        CluNb{TypeL}={};
        for ListL=1:ListNb
            Clu{TypeL}{ListL}={};
            CluNb{TypeL}{ListL}={};
            UsedLimit{TypeL}{ListL}={};
            for NetL=1:length(NetRanks)
                Clu{TypeL}{ListL}{NetL}=[];
                CluNb{TypeL}{ListL}{NetL}={};
                UsedLimit{TypeL}{ListL}{NetL}=[];
                for LimitL=1:length(CurrLimits)
                    MclFile=sprintf('mcl_i%u_m%u_n%u_c%u_%u%s.txt',Inflation,ChipRank,NetRanks(NetL),CurrLimits(LimitL),ListRank(ListL),Suffix);
                    if exist(MclFile,'file')
                        if ImportFlag
                            CurrCluNb=[];
                            CurrClu=zeros(PsNb,1);
                            Rank=0;
                            fid=fopen(MclFile,'r');
                            while 1
                                Rank=Rank+1;
                                CurrPsRank=fgetl(fid);
                                if ~ischar(CurrPsRank), break, end
                                CurrPsRank=str2num(CurrPsRank);
                                CurrClu(CurrPsRank)=Rank;
                                CurrCluNb(Rank,1)=length(CurrPsRank);
                            end
                            fclose(fid);
                            Clu{TypeL}{ListL}{NetL}(:,LimitL)=CurrClu;
                            CluNb{TypeL}{ListL}{NetL}{1,LimitL}=CurrCluNb;
                        end
                    else
                        Missing{TypeL}(end+1,:)=[NetRanks(NetL),CurrLimits(LimitL),ListRank(ListL)];
                        MissingNb=MissingNb+1;
                    end
                end
            end
        end
    end
end

if MissingNb>0
    h=warndlg(sprintf('miss %u files among %u',MissingNb,ResNb));
    waitfor(h)
    fid=fopen('mcl.bat','w');
    Suffixes={'','_diff','_raw','_diff_raw'};
    for TypeL=1:4
        if ~isempty(Missing{TypeL})
            for FileL=1:size(Missing{TypeL},1)
                fprintf(fid,'/work/cinbell/MCL/bin/mcl m%u_n%u_c%u_MCL_%u%s.txt -I 2.0 --abc -o mcl_i2_m%u_n%u_c%u_%u%s.txt\n',ChipRank,Missing{TypeL}(FileL,:),Suffixes{TypeL},ChipRank,Missing{TypeL}(FileL,:),Suffixes{TypeL});
            end
        end
    end
    fclose(fid)
    fid=fopen('exp.bat','w');
    for TypeL=1:4
        if ~isempty(Missing{TypeL})
            for FileL=1:size(Missing{TypeL},1)
                switch TypeL
                    case 1
                        fprintf(fid,'/work/cinbell/multiple/exportnet %u %u %u MCL %u 1.0 0 %u 0\n',ChipRank,Missing{TypeL}(FileL,1),PsNb,Missing{TypeL}(FileL,2),Missing{TypeL}(FileL,3));
                    case 2
                        fprintf(fid,'/work/cinbell/multiple/exportnet %u %u %u MCL %u 1.0 1 %u 0\n',ChipRank,Missing{TypeL}(FileL,1),PsNb,Missing{TypeL}(FileL,2),Missing{TypeL}(FileL,3));
                    case 3
                        fprintf(fid,'/work/cinbell/multiple/exportnet %u %u %u MCL %u 1.0 0 %u 1\n',ChipRank,Missing{TypeL}(FileL,1),PsNb,Missing{TypeL}(FileL,2),Missing{TypeL}(FileL,3));
                    case 4
                        fprintf(fid,'/work/cinbell/multiple/exportnet %u %u %u MCL %u 1.0 1 %u 1\n',ChipRank,Missing{TypeL}(FileL,1),PsNb,Missing{TypeL}(FileL,2),Missing{TypeL}(FileL,3));
                end
            end
        end
    end
    fclose(fid)
end

if ImportFlag
    Info.netRanks=NetRanks;
    Info.listNb=ListNb;
    Info.listRank=ListRank;
    Info.diffFlag=DiffFlag;
    Info.rawFlag=RawFlag;
    Info.diffRawFlag=DiffRawFlag;
    Info.limits{1}=Limits;
    Info.limits{2}=DiffLimits;
    Info.limits{3}=RawLimits;
    Info.limits{4}=DiffRawLimits;

    MatFile=sprintf('m%u_mcl%u.mat',ChipRank,OutRank);
    eval(sprintf('save %s Clu CluNb Info',MatFile))
end