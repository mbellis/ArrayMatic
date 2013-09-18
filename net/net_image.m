%====================%
% FUNCTION NET_IMAGE %
%====================%

%NET_IMAGE: display an image of CVM
%INPUT PARAMETERS
% 1 ChipRank: chip rank
% 2 NetRank: net rank
% 3 Cluster: list(s) of probe set ranks (each list corresponding for exemple to a cluster)
%             or name fo cluster file
% 4 MaxVal: used to determine which colormap use to export image
% 5 AnnealFlag: do annealing clustering before displaying result
% 6 Output: variable part of output file name (chip and net ranks added)

% varargin
% if Cluster is the a file name
% 7 ListPos
% 8 ChipPos
% 9 CorrType (C,C-A,...)
% 10 Column (corr value in MCL)


% net_image(27,164,c,30,1,'mcl_c14_l3')
% net_image(27,183,t,30,1,'mcl_c20_l69')
% net_image(27,183,q,30,1,'mcl_c20_l69b')
% net_image(27,183,r,30,1,'mcl_c20_l69c')
%net_image(27,164,'m27_mcl_n24n164A',30,1,'mcl_c14_164b',8,2,2,1)
%net_image(27,178,'m27_mcl_178A',30,1,'mcl_c14_178a',8,1,2,1)
%net_image(27,183,'m27_mcl_183A',30,1,'mcl_c14_183a',8,1,2,1)

% net_image(27,164,'m27_mcl_n164B',30,1,'mcl_c14_164ss_corr',2,1,1,2)
% net_image(27,164,'m27_mcl_n164B',30,1,'mcl_c14_164ss_diff',2,1,2,2)
%net_image(27,164,'m27_mcl_n164B',30,1,'mcl_c14_164ms_corr',3,1,1,2)
%net_image(27,164,'m27_mcl_n164B',30,1,'mcl_c14_164ms_diff',3,1,2,2)
%net_image(27,164,'m27_mcl_n164B',30,1,'mcl_c14_164sx_corr',4,1,1,2)
%net_image(27,164,'m27_mcl_n164B',30,1,'mcl_c14_164sx_diff',4,1,2,2)



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


function net_image(ChipRank,NetRank,Cluster,MaxVal,AnnealFlag,Output,varargin)
global K
if nargin==10
    Index1=varargin{1};
    Index2=varargin{2};
    CorrType=varargin{3};
    Column=varargin{4};
end
ChipPos=find(K.chip.rank==ChipRank);
cd(fullfile(K.dir.net,sprintf('m%u',ChipRank),sprintf('n%u',NetRank)))
try
    eval(sprintf('load t_m%u_n%u',ChipRank,NetRank))
    ProbesetNb=max(Trans(:,2));
    TransFlag=1;
catch
    ProbesetNb=K.chip.probesetNb(ChipPos);
    TransFlag=0;
end

if ischar(Cluster)
    cd(fullfile(K.dir.net,sprintf('m%u',ChipRank),'mcl'))
    eval(sprintf('load %s',Cluster))
    CluNb=10;
    Cluster={};
    for CluL=1:CluNb
        Cluster{CluL}=find(Clu{CorrType}{Index1}{Index2}(:,Column)==CluL);
        if TransFlag
            Cluster{CluL}=unique(Trans(Cluster{CluL},2));            
        end
    end
else
    CluNb=length(Cluster);
end
Colors=colors(colormap,MaxVal);
Output=sprintf('%s_m%u_n%u',Output,ChipRank,NetRank);


%recover all probeset ranks
PsRank=[];
Limit=[];
EndPos=0;
for CluL=1:CluNb
    try
        PsRank=[PsRank;Cluster{CluL}];
    catch
        PsRank=[PsRank,Cluster{CluL}];        
    end
    Limit(end+1,:)=[EndPos+1,EndPos+length(Cluster{CluL})];
    EndPos=Limit(end,2);
end

cd(fullfile(K.dir.net,sprintf('m%u',ChipRank),sprintf('n%u',NetRank)))

%recover anneal clustering order
PsOrder={};

for CluL=1:CluNb
    if AnnealFlag
        c=load_data(sprintf('c_m%u_n%u.4mat',ChipRank,NetRank),'.',...
            ProbesetNb,ProbesetNb,'uint8','ieee-le',...
            Cluster{CluL},Cluster{CluL});
        [PsOrder{CluL},temp,temp]=annealclust(single(c)/100);
%         h=figure;
%         set(h,'name',sprintf('%s clu%u',Output,CluL))
%         image(c(PsOrder{CluL},PsOrder{CluL}))
    else
        PsOrder{CluL}=[1:length(Cluster{CluL})];
    end
end

for CluL1=1:CluNb
    CluL1
    c=load_data(sprintf('c_m%u_n%u.4mat',ChipRank,NetRank),'.',...
        ProbesetNb,ProbesetNb,'uint8','ieee-le',...
        PsRank,Cluster{CluL1});
    a=load_data(sprintf('a_m%u_n%u.4mat',ChipRank,NetRank),'.',...
        ProbesetNb,ProbesetNb,'uint8','ieee-le',...
        PsRank,Cluster{CluL1});
    if CluL1==1
        if MaxVal==0            
            MaxVal=max(max(a));            
        end
        Colors=colors(colormap,MaxVal+1);
    end
    %make a c horizontal
    c=c';
    a=a';
    %order ps lines
    c=c(PsOrder{CluL1},:);
    a=a(PsOrder{CluL1},:);
    %order ps column
    for CluL2=1:CluNb
        ColOrder=[Limit(CluL2,1):Limit(CluL2,2)];
        c(:,ColOrder)=c(:,ColOrder(PsOrder{CluL2}));
        a(:,ColOrder)=a(:,ColOrder(PsOrder{CluL2}));
    end

    %keep c horizontal
    %make a vertical
    a=a';
    %process all clu and save bmp files
    for CluL2=CluL1:CluNb
        CurrC=c(:,Limit(CluL2,1):Limit(CluL2,2));
        CurrA=a(Limit(CluL2,1):Limit(CluL2,2),:);
        if CluL2==CluL1
            %make composit bloc (C and A values)
            CurrC=triu(CurrC,1)+tril(CurrA,-1);
        end
        FileName=sprintf('%s_clu%u_vs_clu%u.bmp',Output,CluL1,CluL2);
        imwrite(CurrC,Colors,FileName,'bmp')
        if CluL2~=CluL1
            FileName=sprintf('%s_clu%u_vs_clu%u.bmp',Output,CluL2,CluL1); 
            imwrite(CurrA,Colors,FileName,'bmp')
        end
    end
    
end

%write commands for montage (image magick)
FileNameList={};
for CluL1=1:CluNb
    for CluL2=CluL1:CluNb
        FileName=sprintf('%s_clu%u_vs_clu%u.bmp',Output,CluL1,CluL2);
        FileNameList{end+1,1}=FileName;
        if CluL2~=CluL1
            FileName=sprintf('%s_clu%u_vs_clu%u.bmp',Output,CluL2,CluL1);
            FileNameList{end+1,1}=FileName;
        end
    end
end


FileOrder=zeros(CluNb,CluNb);
Rank=0;
for CluL1=1:CluNb
    for CluL2=CluL1:CluNb
        Rank=Rank+1;
        FileOrder(CluL1,CluL2)=Rank;
        if CluL2~=CluL1
            Rank=Rank+1;
            FileOrder(CluL2,CluL1)=Rank;
        end
    end
end
FileOrder=FileOrder';




cd(fullfile(K.dir.net,sprintf('m%u',ChipRank),sprintf('n%u',NetRank)))
fid=fopen(sprintf('%s_imagelist.txt',Output),'w');
CurrFileOrder=FileOrder(1:CluNb,1:CluNb);
CurrFileNameList=FileNameList(CurrFileOrder(:));
for NameL=1:length(CurrFileNameList)
    fprintf(fid,'%s\n',CurrFileNameList{NameL});
end
fclose(fid)
fid=fopen(sprintf('%s.bat',Output),'w');
fprintf(fid,'montage @%s_imagelist.txt -tile %ux%u -geometry "10x10<" -mode Concatenate  %s.bmp',Output,CluNb,CluNb,Output)
fclose(fid)


