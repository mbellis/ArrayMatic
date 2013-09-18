%=========================%
% FUNCTION NET_DISPLAYMCL %
%=========================%

%NET_DISPLAYMCL: displays a CVM structured by a MCL result and eventually by a NTS result
%INPUT PARAMETERS
% 1 ChipRank: chip rank
% 2  NetRank: network(s) to be displayed 
% 3 MclFile: MCL result
% 4 MclType: rank of type in {'C','C-A','rawC','raw(C-A)'};
% 5 MclRank: rank of MCL results (MCL may have been conducted on several networks)
% varargin:
% 4 NtsFile: NTS result
% 5 Density: NTS clusters to be used are specified by their density

% net_displaymcl(93,[1,2],'m93a',2,3,'m93_tsn_from_m93n1ton3_corr60',30)

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


function net_displaymcl(ChipRank,NetRank,MclFile,MclType,MclRank,varargin)
global K

TypeName={'C','C-A','rawC','raw(C-A)'};

if length(ChipRank)>1
    h=errordlg('net_displaymcl uses only one chip');
    waitfor(h)
    error('process canceled')
end

NtsFlag=0;
if nargin==7
    NtsFile=varargin{1};
    Density=varargin{2};
    NtsFlag=1;
    NtsChipRank=regexp(NtsFile,'(?<=m)\d+(?=_)','match');
    NtsChipRank=str2num(NtsChipRank{1});
end


ChipPos=find(K.chip.rank==ChipRank);
PsNb=K.chip.probesetNb(ChipPos);

%load MCL results
CmlDir=fullfile(K.dir.net,sprintf('m%u',ChipRank),'mcl');
cd(CmlDir)
eval(sprintf('load m%u_mcl_%s.mat',ChipRank,MclFile))


%load eventually NTS results
if NtsFlag
    NtsDir=fullfile(K.dir.net,sprintf('m%u',NtsChipRank),'tsn','result');
    cd(NtsDir)
    eval(sprintf('load %s',NtsFile))
    DensPos=find(Densities(UsedDensities)==Density);
    NtsClu=NtsClusters{DensPos};
    NtsCluNb=length(NtsClu);
end




%load ps lists
ListNb=Info.listNb;
    NetNb=length(Info.netRanks);
    ListSize=zeros(ListNb,1);
    PsRankList=cell(ListNb,1);
    %cd(K.dir.chip)
    cd(fullfile(K.dir.net,sprintf('m%u',ChipRank),'list'))
    for ListL=1:ListNb
        fid=fopen(sprintf('m%u_pslist%u.u32',ChipRank,Info.listRank(ListL)),'r','ieee-le');
        PsRankList{ListL}=fread(fid,inf,'uint32');
        ListSize(ListL)=length(PsRankList{ListL});
        fclose(fid);
    end

%translate MCL results
Clu=Clu{MclType}{1}{MclRank};
MclCluNb=max(Clu);

%reorder Mcl according to NtsResults and annealing clustering
if NtsFlag
    NtsMcl=zeros(NtsCluNb,MclCluNb);
    for MclL=1:MclCluNb
        MclClu{MclL}=find(Clu==MclL);
        NewMclClu{MclL}=[];
        FullMclClu{MclL}=[];
    end

    for NtsL=1:NtsCluNb
        for MclL=1:MclCluNb
            Ps=intersect(MclClu{MclL},NtsClu{NtsL});
            if ~isempty(Ps)
                NtsMcl(NtsL,MclL)=length(Ps);
                FullMclClu{MclL}=[FullMclClu{MclL};Ps];
                Ps=setdiff(Ps,NewMclClu{MclL});
                if ~isempty(Ps)
                    NewMclClu{MclL}=[NewMclClu{MclL};Ps];
                end
            end
        end        
    end
    for MclL=1:MclCluNb
        Ps=setdiff(MclClu{MclL},NewMclClu{MclL});
        if ~isempty(Ps)
            C=load_cvm(ChipRank,NetRank(1),Ps,Ps,1,1,1);
            [PsOrder,temp,temp]=annealclust(single(C)/100);
            NewMclClu{MclL}=[NewMclClu{MclL};Ps(PsOrder)];
            FullMclClu{MclL}=[FullMclClu{MclL};Ps(PsOrder)];
        end
    end
else
    for MclL=1:MclCluNb
        MclClu=find(Clu==MclL);
        if length(MclClu)>50
            C=load_cvm(ChipRank,NetRank(1),MclClu,MclClu,1,1,1);
            [NewMclClu{MclL},temp,temp]=annealclust(single(C)/100);
        else
            NewMclClu{MclL}=MclClu;
        end
    end    
end
    
%construct final PsOrder
PsOrder=[];
for MclL=1:MclCluNb
    PsOrder=[PsOrder;NewMclClu{MclL}];
end
PsOrder=[PsOrder;setdiff([1:PsNb]',PsOrder)];
if NtsFlag
    FullPsOrder=[];
    for MclL=1:MclCluNb
        FullPsOrder=[FullPsOrder;FullMclClu{MclL}];
    end
end
FullPsOrder=[FullPsOrder;setdiff([1:PsNb]',FullPsOrder)];

%% DISPLAY CVM

if NtsFlag
    OutNb=2;
else
    OutNb=1;
end
Factor=1;
for OutL=1:OutNb
    if OutL==1
        BlocNb=ceil(PsNb/500);
    else
        BlocNb=ceil(length(FullPsOrder)/500);
    end
    for NetL=1:length(NetRank)
        if OutL==1
            Output='';
        else
            Output='full';
        end        
        net_displaycvm(ChipRank,NetRank(NetL),PsOrder,Output,100,3000,1)                        
    end        
end