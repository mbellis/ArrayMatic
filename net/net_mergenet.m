% NET_MERGENET makes a synthesis of several networks constructed from
% different combinations of biological conditions.
% One option is to take the intersection of all the merged networks
% (synthetic values of CORR and ANTI are then for a particular pair of probe sets the mean or the median of all the non zero values, if a
% given percent of networks have non zero values (otherwise the synthetic
% values are null).
%
% INPUT
% MergeType: type of merging : mean, intersection, union
% NewNetRank: the rank of the first merged net (intersection or mean)
% QLimitFlag: indicates if QLimit has been used for constructing processed networks
% VARARGIN
% MinNetNb: the minimal number of networks which must have significative
% values (default = NetNb) for integrate corr and anti intersection network

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

%net_mergenet(27,[88:108],'mean',0)
%net_mergenet(8,[39:53],'mean',0)

function net_mergenet(ChipRank,NetRanks,MergeType,QLimitFlag,varargin)
global K
Continue=1;
if nargin<4
    h=warndlg('net_mergenet needs at least two parameters (MergeType & QLimitFlag');
    waitfor(h)
    Continue=0;
elseif nargin==5
    MinNetNb=varargin{1};
end
ChipPos=find(K.chip.rank==ChipRank);
ProbesetNb=K.chip.probesetNb(ChipPos);
if Continue
    Pos=length(K.net{ChipRank}.name)+1;
    NetRank=setdiff([1:Pos],K.net{ChipRank}.rank);
    NetRank=NetRank(1);    
    NetDir=fullfile(K.dir.net,sprintf('m%u',ChipRank),sprintf('n%u',NetRank));
    mkdir(NetDir)
    cd(NetDir)
    %[ChipRank,NetRank,NetPos]=select_net('multiple');
    NetNb=length(NetRanks);
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

% if Continue
%     ProbesetNb=K.chip.probeSetNbs{ChipRank}(1);
%     %CONTROL THAT THE FIRST MERGED NETWORK DOES NOT EXIST (CREATED LATTER IN K.net)
%     if ~isempty(find(K.net{ChipRank}.rank==NetRank))
%         Suggest=setdiff(1:max(K.net{ChipRank}.rank),K.net{ChipRank}.rank);
%         if isempty(Suggest)
%             Suggest1=max(K.net{ChipRank}.rank)+1;
%         else
%             Suggest1=Suggest(1);
%         end
%         if ~isempty(find(K.net{ChipRank}.rank==NetRank))
%             h=warndlg(sprintf('You must indicate not used ranks\n for the merged networks\n (e.g. %u)',Suggest1));
%         elseif ~isempty(find(K.net{ChipRank}.rank==NetRank))
%             h=warndlg(sprintf('%u is already used. Select another\n rank for the merged networks\n (e.g. %u)',NetRank,Suggest1));
%             waitfor(h)
%             Continue=0;
%         end
%     end
% end

if Continue
    ChipSetDir=fullfile(K.dir.net,sprintf('m%u',ChipRank));
    if QLimitFlag
        HistRepNbC=zeros(NetNb,1);
        HistIfExistC=zeros(100,1);
        HistIfExistA=zeros(100,1);
        HistIfSupeC=zeros(100,1);
        HistIfSupeA=zeros(100,1);
    end

    %create NetFile
    NetA1=sprintf('a_m%u_n%u.4mat',ChipRank,NetRank);
    NetC1=sprintf('c_m%u_n%u.4mat',ChipRank,NetRank);
    NetF1=sprintf('f_m%u_n%u.4mat',ChipRank,NetRank);
    NetRepC=sprintf('c_m%u_n%u_rep.4mat',ChipRank,NetRank);
    cd(NetDir)
    NetA1Fid=fopen(NetA1,'w','ieee-le');
    NetC1Fid=fopen(NetC1,'w','ieee-le');
    NetF1Fid=fopen(NetF1,'w','ieee-le');
    NetRepCFid=fopen(NetRepC,'w','ieee-le');
    HReadC=zeros(NetNb,1);
    HReadA=zeros(NetNb,1);

    for NetL=1:NetNb
        CurrNet=NetRanks(NetL);
        if NetL==1
            NetList=num2str(CurrNet);
        else
            NetList=sprintf('%s %u',NetList,CurrNet);
        end
        NetDir=fullfile(ChipSetDir,sprintf('n%u',CurrNet));
        cd(NetDir)
        NetFile=sprintf('a_m%u_n%u.4mat',ChipRank,CurrNet);
        HReadA(NetL)=fopen(NetFile,'r','ieee-le');
        NetFile=sprintf('c_m%u_n%u.4mat',ChipRank,CurrNet);
        HReadC(NetL)=fopen(NetFile,'r','ieee-le');
        NetFile=sprintf('f_m%u_n%u.4mat',ChipRank,CurrNet);
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
        CurrF=uint8(ones(ProbesetNb,NetNb)*NaN);
        CORR=uint8(zeros(ProbesetNb,1));
        ANTI=uint8(zeros(ProbesetNb,1));
        RAWCORR=uint8(zeros(ProbesetNb,1));
        FREQ=uint8(zeros(ProbesetNb,1));
        for NetL=1:NetNb
            CurrA(:,NetL)=load_data(HReadA(NetL),'',ProbesetNb,ProbesetNb,'uint8','ieee-le',1:ProbesetNb,PsL);
            CurrC(:,NetL)=load_data(HReadC(NetL),'',ProbesetNb,ProbesetNb,'uint8','ieee-le',1:ProbesetNb,PsL);
            CurrF(:,NetL)=load_data(HReadF(NetL),'',ProbesetNb,ProbesetNb,'uint8','ieee-le',1:ProbesetNb,PsL);
        end
        if ~isequal(MergeType,'mean')
            %number of significative values
            RepInfo=max([sum(CurrC>0,2),sum(CurrA>0,2)],[],2);
            %write RepInfo
            save_data(RepInfo,NetRepCFid,'','a','uint8','ieee-le',ProbesetNb,[],PsL);
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
                RAWCORR(SelIndex)=nanmean((double(CurrC(SelIndex,:)).*double(CurrF(SelIndex,:)))/100,2);
                FREQ(SelIndex)=uint8((double(CORR(SelIndex))./RAWCORR(SelIndex))*100);
                %         %write ANTI & CORR
                save_data(CORR,NetC2Fid,'','a','uint8','ieee-le',ProbesetNb,[],PsL);
                save_data(ANTI,NetA2Fid,'','a','uint8','ieee-le',ProbesetNb,[],PsL);
                %remove 100
                CORR(PsL)=0;
                %frequency of values
                HistIfExistC=HistIfExistC+histc(CORR,1:100);
                HistIfExistA=HistIfExistA+histc(ANTI,1:100);


                %frequency of repetition
                HistRepNbC=HistRepNbC+histc(RepInfo,1:NetNb);
                if PsNb==100
                    t=t+toc;
                    sprintf('%u Ps processed in %u min : estimates that rests %u min',PsL,round(t/60),round((ProbesetNb-PsL)*t/(PsL*60)))
                    tic
                    PsNb=0;
                end

            end
        else
            CORR=uint8(ceil(mean(CurrC,2)));
            ANTI=uint8(ceil(mean(CurrA,2)));
            RAWCORR=mean((double(CurrC).*double(CurrF))/100,2);
            FREQ=uint8((double(CORR)./RAWCORR)*100);
        end
        %write ANTI & CORR
        save_data(CORR,NetC1Fid,'','a','uint8','ieee-le',ProbesetNb,[],PsL);
        save_data(ANTI,NetA1Fid,'','a','uint8','ieee-le',ProbesetNb,[],PsL);
        save_data(FREQ,NetF1Fid,'','a','uint8','ieee-le',ProbesetNb,[],PsL);

    end

    if QLimitFlag
        %remove 100 on corr diagonal
        HistRepNbC(end)=HistRepNbC(end)-ProbesetNb;
        eval(sprintf('save m%u_n%u_hist Hist*',ChipRank,NetRank))
        h=figure;
        set(h,'color',[1,1,1])
        subplot(1,2,1)
        plot(1:100,HistIfExistC/sum(HistIfExistC),'m')
        hold on
        plot(1:100,HistIfExistA/sum(HistIfExistA),'c')
        plot(1:100,HistIfSupeC/sum(HistIfSupeC),'r')
        plot(1:100,HistIfSupeA/sum(HistIfSupeA),'b')
        title('CORR and ANTI distributions')
        legend({'Union CORR','Union ANTI','Intersect CORR','Intersect ANTI'})
        subplot(1,2,2)
        set(h,'color',[1,1,1])
        plot(1:NetNb,HistRepNbC/sum(HistRepNbC),'K')
        title('reproducibility (nb of networks)')

        fclose(NetC1Fid);
        fclose(NetA1Fid);
        fclose(NetC2Fid);
        fclose(NetA2Fid);
        fclose(NetRepCFid);
        for NetL=1:NetNb
            fclose(HReadA(NetL));
            fclose(HReadC(NetL));
        end
    end

    %register new networks
    Tempo=K.net;
    RefPos=find(Tempo{ChipRank}.rank==NetRanks(1));
    %intersect of networks
    Tempo{ChipRank}.name{Pos,1}=sprintf('%s of networks %s',MergeType,NetList);
    Tempo{ChipRank}.rank(Pos,1)=NetRank;
    Tempo{ChipRank}.biolRank{Pos,1}='';
    Tempo{ChipRank}.compNb(Pos,1)=Tempo{ChipRank}.compNb(RefPos,1);    
    Tempo{ChipRank}.fdr(Pos,:)=Tempo{ChipRank}.fdr(RefPos,:);
    Tempo{ChipRank}.s(Pos,:)=Tempo{ChipRank}.s(RefPos,:);    
    Tempo{ChipRank}.blocNb(Pos,1)=Tempo{ChipRank}.blocNb(RefPos,1);
    Tempo{ChipRank}.blocSize(Pos,1)=Tempo{ChipRank}.blocSize(RefPos,1);
    Tempo{ChipRank}.comment{Pos,1}=Tempo{ChipRank}.comment{RefPos,1};        
    Tempo{ChipRank}.netMade(Pos,1)=1;    
    Tempo{ChipRank}.nb=length(Tempo{ChipRank}.name);
    cd(K.dir.common)
    save netlist Tempo
    K.net=Tempo;
end
