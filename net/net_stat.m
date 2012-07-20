%===================%
% FUNCTION NET_STAT %
%===================%

%NET_STAT calculates the mean corr, anti for each probe set,
% and the mean connectivity for several combination of corr
% and anti values

%INPUT PARAMETERS

%EXTERNAL FILES
% 1 ModelRank: chip model rank
% 2   NetRanks: list of net ranks

%OUTPUT FILES
% sprintf('save MeanCorr_m%u_n%u MeanCorr
% sprintf('save Connect_C%02u_A%02u_m%u_n%u Connect',CLimit,ALimit,ModelRank,NetRank))    

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
%calculate mean corr, anti and connectivity

function net_stat(ModelRank,NetRanks)


global K
for NetL=1:length(NetRanks)
    CurrNetRank=NetRanks(NetL);
    tic
    NetPos=find(K.net{ModelRank}.rank==CurrNetRank);
    if length(NetPos)~=1
        h=errordlg(sprintf('%u net with rank equal to %u instead of one and only one. Process canceled',length(NetPos),CurrNetRank));
        waitfor(h)
        error('process canceled')
    end

    %BlocSize=K.net{ModelRank}.blocSize(NetPos);
    PsNb=K.chip.probesetNb(ModelRank);
    BlocSize=5000;
    BlocNb=ceil(PsNb/BlocSize);
    %LastBlocSize=mod(PsNb,BlocSize);
    if BlocNb==0
        h=errordlg(sprintf('K.net{%u}.blocnb(%u)=0 ! Process canceled',ModelRank,NetPos));
        waitfor(h)
        error('process canceled')
    end
    NetDir=fullfile(K.dir.net,sprintf('m%u',ModelRank),sprintf('n%u',CurrNetRank));
    StatDir=fullfile(NetDir,'stat');

    ValRange=0:10:60;
    cd(NetDir)

    for i=1:length(ValRange)
        Conn{i}=zeros(PsNb,2);
    end
    MeanCorr=zeros(PsNb,2);

    AFile=sprintf('a_m%u_n%u.4mat',ModelRank,CurrNetRank);
    CFile=sprintf('c_m%u_n%u.4mat',ModelRank,CurrNetRank);
    for BlocL=1:BlocNb
        BlocL
        if BlocL<BlocNb
            Range=BlocSize*(BlocL-1)+1:BlocSize*BlocL;
        else
            Range=BlocSize*(BlocL-1)+1:PsNb;
        end
        for ValPos=1:length(ValRange)
            ALimit=ValRange(ValPos);
            CLimit=ALimit;
            if ValPos==1
                if BlocL<BlocNb
                    A=load_data(AFile,'./',PsNb,PsNb,'uint8','ieee-le',1:PsNb,(BlocL-1)*BlocSize+1:BlocL*BlocSize);
                else
                    A=load_data(AFile,'./',PsNb,PsNb,'uint8','ieee-le',1:PsNb,(BlocL-1)*BlocSize+1:PsNb);
                end
                A=single(A);
                % replaces 0 by NaN
                A(A==0)=NaN;
                MeanCorr(Range,2)=nanmean(A)';
                %NaN are are replaced by 0
                A=uint8(A');
                if BlocL<BlocNb
                    C=load_data(CFile,'./',PsNb,PsNb,'uint8','ieee-le',1:PsNb,(BlocL-1)*BlocSize+1:BlocL*BlocSize);
                else
                    C=load_data(CFile,'./',PsNb,PsNb,'uint8','ieee-le',1:PsNb,(BlocL-1)*BlocSize+1:PsNb);
                end
                C=single(C);
                % replaces 0 by NaN
                C(C==0)=NaN;
                %replace 100 by NaN on the main diagonal
                Offset=BlocSize*(BlocL-1);
                for i=1:size(C,2)
                    C(Offset+i,i)=NaN;
                end
                MeanCorr(Range,1)=nanmean(C)';
                %NaN are are replaced by 0
                C=uint8(C');
            end
            % replaces positive values by one for calculating connectivity
            % at the current CORR and ANTI limits combination
            R=uint8(zeros(size(C)));
            R(C>CLimit)=1;
            Conn{ValPos}(Range,1)=sum(R,2);

            R(R==1)=0;
            R(A>ALimit)=1;
            Conn{ValPos}(Range,2)=sum(R,2);

            R(R==1)=0;
            R(C>CLimit&A>ALimit)=1;
            Conn{ValPos}(Range,3)=sum(R,2);

            R(R==1)=0;
            R(C>CLimit&A==0)=1;
            Conn{ValPos}(Range,4)=sum(R,2);

            R(R==1)=0;
            R(A>ALimit&C==0)=1;
            Conn{ValPos}(Range,5)=sum(R,2);
        end
        clear A C
    end
    try
        cd(StatDir)
    catch
        cd(NetDir)
        mkdir('stat')
        cd(StatDir)
    end
    eval(sprintf('save MeanCorr_m%u_n%u MeanCorr',ModelRank,CurrNetRank))
    for ValPos=1:length(ValRange)
        Connect=Conn{ValPos};
        Connect(isnan(Connect))=0;
        CLimit=ValRange(ValPos);
        ALimit=CLimit;
        eval(sprintf('save Connect_C%02u_A%02u_m%u_n%u Connect',CLimit,ALimit,ModelRank,CurrNetRank))
    end
    toc
end