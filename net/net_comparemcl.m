%=========================%
% FUNCTION NET_COMPAREMCL %
%=========================%

%NET_COMPAREMCL compares MCL results obtained with CORR-ANTI values (Type=2)
%INPUT PARAMETERS
% 1        ChipRank: list of two chip ranks
% 2         TypePos: type position (unique) to be used in each chip
% 3         ListPos: list position (unique) to be used in each chip
% 4          NetPos: network position (unique) to be used in each chip
% 5        NetRanks: networks existing in each chip
% 6         Postfix: list of postfix using to load MCL results
% 7       LimitFlag: display MCL results with identical CORR limit (if = 1)
% varargin:
% 8           MclPos: MCL result to be used
% 9          TwinPos: correspondance between clusters between the two chips
% 10 CommonPsFile{1}: correspondance bewteen two chips for first chip
% 11 CommonPsFile{2}: correspondance bewteen two chips for second chip


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


function net_comparemcl(ChipRank,TypePos,ListPos,NetPos,NetRanks,Postfix,LimitFlag,varargin)
%net_comparemcl([8,27],[125,118])
%net_comparemcl([8,27],[55,24])
%net_comparemcl([8,27],[54,23])
%net_comparemcl([8,27],[54,23],[2,4],{[1;2;5;6;7],[2;7;3;10;9]})
    %net_comparemcl([8,42],[55,2])
%net_comparemcl([27,42],[24,2])
%net_comparemcl([8,27],[55,24],[3,4],{[4,7;5,5;3,3;2,2;1,1],[1:5]'})
%net_comparemcl([8,42],[55,2],[1,4],{[4,4;1,1;6,9;3,3;5,5;7,7;10,10;8,8],[1:8]'})
%net_comparemcl([27,42],[24,2],[4,5],{[1,4,3,5]',[4,4,4,4;1,2,5,6;8,8,8,8;3,9,10,10]})
%net_comparemcl([8,8],[7,8])
%net_comparemcl([8,8],[7,8],[4,4],{[1;2;3;4;5;6;7;9],[1;3;4;2;6;7;8;5]})
%net_comparemcl([27,27],[2,3])
%net_comparemcl([8,27],[7,2])
% net_comparemcl([8,27],[55,24],[4,4],{{[4,7,24,30,33,34,36,45,47,48,49,51,53,55,63,65,71,79,85,92,95,96,99,100],...
% [5,12,14,17,23,31,50,62,82],[3,15,43,76,80,87,90],[2,11,13,21,25,28,32,41,42,44,46,54,56,66,75,78,89,93],[1,6]},...
% {[1],[2],[3],[4],[5]}})
% net_comparemcl([8,27],[55,24],'b',[3,4],{{[4,7,24,30,33,34,36,45,47,48,49,51,53,55,63,65,71,79,85,92,95,96,99,100],...
% [5,12,14,17,23,31,50,62,82],[3,15,76,80,87,90],[2,11,13,21,25,28,32,41,42,44,46,54,56,66,75,78,89,93],[1,64]},...
% {[1],[2],[3],[4],[5]}})
%net_comparemcl([8,27],[55,24],{'n55','n24'},[3,4],{{[3,4],[5],[3],[2],[1],[9]},{[1],[2],[3],[4],[5],[6]}},'m8n228_m27n164_combinedps_corr60','m27n164_m8n228_combinedps_corr60')
% verification chip2chip is correct => diagonal
%net_comparemcl([5,27],[5,24],{'2004n5','m5_2004n5'},[],{},'m5n123_m27n164_combinedps_corr60','m27n164_m5n123_combinedps_corr60')
%net_comparemcl([5,8],[5,55],{'2004n5','m5_2004n5'},[],{},'m5n123_m8n228_combinedps_corr60','m8n228_m5n123_combinedps_corr60')
%assignation region
%net_comparemcl([5,8],[5,228],{'2006n5','n228'},[],{},'m5n123_m8n228_combinedps_corr60','m8n228_m5n123_combinedps_corr60')
%net_comparemcl([2,3],[80,86],{'2006n5','2006n7'},[],{},'m2n80_m3n86_combinedps_corr60','m3n86_m2n80_combinedps_corr60')
%net_comparemcl([5,5],[2,2],[1,2],[1,1],{[5],[123]},{'2006n5','n123A'},0)
%net_comparemcl([5,5],[2,2],[1,2],[1,1],{[123],[123]},{'n123A','n123A'},1)
%net_comparemcl([5,5],[2,2],[1,1],[1,1],{[5],[123]},{'2006n5','n123B'},0)
%net_comparemcl([2,3],[2,2],[1,1],[1,1],[80,86],{'n80B','n86B'},1,[],{},'m2n80_m3n86_combinedps_corr60','m3n86_m2n80_combinedps_corr60')
%net_comparemcl([3,5],[2,2],[1,1],[1,1],[86,123],{'n86B','n123B'},1,[],{},'m3n86_m5n123_combinedps_corr60','m5n123_m3n86_combinedps_corr60')
%net_comparemcl([8,27],[2,2],[1,1],[1,1],{[228],[164]},{'n228F','n164F'},1,[],{},'m8n228_m27n164_combinedps_corr60','m27n164_m8n228_combinedps_corr60')
%net_comparemcl([8,27],[2,2],[1,1],[1,1],{[228],[164]},{'n228G','n164G'},1,[],{},'m8n228_m27n164_combinedps_corr60','m27n164_m8n228_combinedps_corr60')
%net_comparemcl([8,27],[2,2],[2,2],[1,1],{[228],[164]},{'n228F','n164F'},1,[],{},'m8n228_m27n164_combinedps_corr60','m27n164_m8n228_combinedps_corr60')
%net_comparemcl([8,27],[2,2],[2,2],[1,1],{[228],[164]},{'n228G','n164G'},1,[],{},'m8n228_m27n164_combinedps_corr60','m27n164_m8n228_combinedps_corr60')
%net_comparemcl([5,6],[2,2],[1,1],[1,1],{[123],[63]},{'n123B','n63B'},1,[],{},'m5n123_m6n63_combinedps_corr60','m6n63_m5n123_combinedps_corr60')
%net_comparemcl([5,6],[2,2],[2,2],[1,1],{[123],[63]},{'n123B','n63B'},1,[],{},'m5n123_m6n63_combinedps_corr60','m6n63_m5n123_combinedps_corr60')
global K

%% VERIFICATION OF PARAMETERS
if nargin==7
    TwinFlag=0;
elseif nargin==9
    TwinFlag=1;
    MclPos=varargin{1};
    TwinPos=varargin{2};
elseif nargin==11
    if ~isempty(varargin{1})
        TwinFlag=1;
        MclPos=varargin{1};
        TwinPos=varargin{2};    
    else
        TwinFlag=0;
    end
    CommonPsFile{1}=varargin{3};
    CommonPsFile{2}=varargin{4};
else
    h=errordlg('needs 7 or 9 or 11 parameters');
    waitfor(h)
    error('process canceled')
end


if length(ChipRank)>2 |length(NetRanks)>2
    h=errordlg('at most 2 networks');
    waitfor(h)
    error('process canceled')
end

if length(ChipRank)==1
    ChipRank=[ChipRank,ChipRank];
    ChipPos=find(K.chip.rank==ChipRank);
    ChipPos=[ChipPos,ChipPos];    
else
    ChipPos=find(K.chip.rank==ChipRank(1));
    ChipPos=[ChipPos,find(K.chip.rank==ChipRank(2))];
end
PsNb=[K.chip.probesetNb(ChipPos(1)),K.chip.probesetNb(ChipPos(2))];
if length(NetRanks)==1
    NetRanks=[NetRanks,NetRanks];
end

%% LOAD CORESPONDANCE BETWEEN CHIPS IF NECESSARY

if ChipRank(1)==8&ChipRank(2)==27
    cd(K.dir.chip)
    load m8_m27_commonps
    ComPsRanks=ComPsRank;
    ComPsRank={};
    ComPsRank{1}=ComPsRanks(:,1);
    ComPsRank{2}=ComPsRanks(:,2);
else
    if ChipRank(1)~=ChipRank(2)
        cd(K.dir.chip)
        for ChipL=1:2
            eval(sprintf('load %s;',CommonPsFile{ChipL}))
            ComPsRank{ChipL}=zeros(length(PsRank),1);
            %construct correspondance between the current chip and the common order of ps
            %keep one probe set if exist alternative probe sets (similarity >=75%)
            for PsL=1:length(PsRank)
                ComPsRank{ChipL}(PsL)=PsRank{PsL}(1);
            end
        end
        clear PsRank
    else
        for ChipL=1:2
            ComPsRank{ChipL}=[1:PsNb(1)]';
        end
    end
end




% if ChipRank(1)~=ChipRank(2)
%     Postfix=sprintf('m%u_m%u_commonps.mat',min(ChipRank),max(ChipRank));
%     if exist(Postfix,'file')
%         load(Postfix)
%         if ChipRank(1)>ChipRank(2)
%             Temp=ComPsRank;
%             ComPsRank(:,1)=ComPsRank(:,2);
%             ComPsRank(:,2)=Temp(:,1);
%             clear Temp
%         end
%     else
%         h=errordlg(sprintf('no correspondance file between m%y and m%u',ChipRank(1),ChipRank(2)));
%         waitfor(h)
%         error('process canceled')
%     end
% else
%     ComPsRank=repmat([1:K.chip.probeSetNbs{ChipRank(1)}(1)]',1,2);
% end

%% LOAD MCL CLUSTERS
%recover CorrLimit
CorrLimits=[];
for ChipL=1:2
    NetDir=fullfile(K.dir.net,sprintf('m%u',ChipRank(ChipL)),'mcl');
    cd(NetDir)
    load(sprintf('m%u_mcl_%s.mat',ChipRank(ChipL),Postfix{ChipL}))
    if exist('CorrLimit')
        CorrLimits=union(CorrLimits,CorrLimit);
        CLimit{ChipL}=CorrLimit;
        clear CorrLimit
    else
        CorrLimits=union(CorrLimits,Info.limits{2});
        CLimit{ChipL}=Info.limits{2};        
    end
end
CluSize=cell(1,2);
for ChipL=1:2
    NetDir=fullfile(K.dir.net,sprintf('m%u',ChipRank(ChipL)),'mcl');
    cd(NetDir)
    load(sprintf('m%u_mcl_%s.mat',ChipRank(ChipL),Postfix{ChipL}))       
    %eliminate empty results
    ClearPos=find(sum(Clu{TypePos(ChipL)}{ListPos(ChipL)}{NetPos(ChipL)})==0);
    if ~isempty(ClearPos)
        Clu{TypePos(ChipL)}{ListPos(ChipL)}{NetPos(ChipL)}(:,ClearPos)=[];
        CLimit{ChipL}(ClearPos)=[];
    end
    CluSize{ChipL}=cell(1,2);
    for LimitL=1:size(Clu{TypePos(ChipL)}{ListPos(ChipL)}{NetPos(ChipL)},2)
        LimitPos=find(CorrLimits==CLimit{ChipL}(LimitL));
        for CluL=1:10
            CluSize{ChipL}{1}(CluL,LimitPos)=length(find(Clu{TypePos(ChipL)}{ListPos(ChipL)}{NetPos(ChipL)}(:,LimitL)==CluL));
        end
    end
    if ChipRank(1)~=ChipRank(2)
        %ComClu{ChipL}=Clu{2}{1}{1}(ComPsRank(:,ChipL),:);        
        ComClu{ChipL}=Clu{TypePos(ChipL)}{ListPos(ChipL)}{NetPos(ChipL)}(ComPsRank{ChipL},:);        
        for LimitL=1:size(Clu{TypePos(ChipL)}{ListPos(ChipL)}{NetPos(ChipL)},2)
            LimitPos=find(CorrLimits==CLimit{ChipL}(LimitL));
            for CluL=1:10
                CluSize{ChipL}{2}(CluL,LimitPos)=length(find(ComClu{ChipL}(:,LimitL)==CluL));
            end
        end
    else
        ComClu{ChipL}=Clu{TypePos(ChipL)}{ListPos(ChipL)}{NetPos(ChipL)};
        CluSize{ChipL}{2}= CluSize{ChipL}{1};
    end
    clear Clu
end

%display the percentage of ps eliminated in the chip common to both chips
round((CluSize{1}{1}-CluSize{1}{2})*100./CluSize{1}{1})
round((CluSize{2}{1}-CluSize{2}{2})*100./CluSize{2}{1})


%% CALCULATE AND DISPLAY SIMILARITY BETWEEN MCL REGIONS

%calculate similarity on a small grid
[CommonLimit,CommonPos1,CommonPos2]=intersect(CLimit{1},CLimit{2});
Pos=0;
Sim={};
NbSup=[];
ComNb={};
LimitSup=50;
if LimitFlag
    EndLimit1=length(CommonLimit);
else
    EndLimit1=size(ComClu{1},2);
end
for LimitL1=1:EndLimit1
%for LimitL1=1:4
    Continue=1;
    if LimitFlag
        if isempty(intersect(LimitL1,CommonPos1))
            Continue=0;
        else
            StartLimit2=CommonPos2(LimitL1);
            EndLimit2=CommonPos2(LimitL1);
        end
    else
        StartLimit2=1;
        EndLimit2=size(ComClu{2},2);
    end
    if Continue
        CurrCluSize=histc(ComClu{1}(:,LimitL1),unique(ComClu{1}(:,LimitL1)));
        LIMIT1=find(CurrCluSize>=LimitSup);
        LIMIT1=LIMIT1(end)+1;
        %LIMIT1=min(max(ComClu{1}(:,LimitL1))+1,LimitSup);
    end
    %
    for LimitL2=StartLimit2:EndLimit2
    %for LimitL2=1:4
        Pos=Pos+1;      
        CurrCluSize=histc(ComClu{2}(:,LimitL2),unique(ComClu{2}(:,LimitL2)));
        LIMIT2=find(CurrCluSize>=LimitSup);
        LIMIT2=LIMIT2(end)+1;
        %LIMIT2=min(max(ComClu{2}(:,LimitL2))+1,LimitSup);
        Sim{Pos}=zeros(LIMIT1,LIMIT2);
        ComNb{Pos}=zeros(LIMIT1,LIMIT2);
        for CluL1=1:LIMIT1-1
            for CluL2=1:LIMIT2-1
                ComNb{Pos}(CluL1,CluL2)=length(intersect(find(ComClu{1}(:,LimitL1)==CluL1),find(ComClu{2}(:,LimitL2)==CluL2)));
                Sim{Pos}(CluL1,CluL2)=round(ComNb{Pos}(CluL1,CluL2)*100/min(length(find(ComClu{1}(:,LimitL1)==CluL1)),length(find(ComClu{2}(:,LimitL2)==CluL2))));
            end
        end
        NbSup=[NbSup;[LimitL1,LimitL2,length(find(Sim{Pos}>=50)),length(find(Sim{Pos}>=75))]];
    end
end



%display similarity between MCL regions on a small grid
CompNb=length(Sim);
ColNb=ceil(sqrt(CompNb));
RawNb=floor(CompNb/ColNb);
if ColNb*RawNb<CompNb
    RawNb=RawNb+1;
end
for RoundL=1:3
    h=figure;
    set(gcf,'color',[1,1,1])
    for CompL=1:CompNb
        if CompL
            subplot(RawNb,ColNb,CompL)
            try
            if RoundL==1
                CurrSim=Sim{CompL};
                CurrSim(find(CurrSim<70))=0;
                pcolor(CurrSim);
            elseif RoundL==2
                pcolor(Sim{CompL});
            else
                pcolor(ComNb{CompL});
            end
            if CompL==1
                ylabel(sprintf('%u - m%u',CLimit{1}(NbSup(CompL,1)),ChipRank(1)))
                xlabel(sprintf('%u - m%u',CLimit{2}(NbSup(CompL,2)),ChipRank(2)))
            else
                ylabel(sprintf('%u',CLimit{1}(NbSup(CompL,1))))
                xlabel(sprintf('%u',CLimit{2}(NbSup(CompL,2))))
            end
            catch
            end
        end
    end
end



%calculate similarity on a large grid

Pos=0;
Sim={};
NbSup=[];
ComNb={};
Limit=[];

%for LimitL=1:length(CommonLimit)
for LimitL=1:length(CLimit{2})
    %LimitPos1=find(CLimit{1}==CommonLimit(LimitL));
    LimitPos1=1;
    %LimitPos2=find(CLimit{2}==CommonLimit(LimitL));
    LimitPos2=LimitL;
    LIMIT1=min(max(ComClu{1}(:,LimitPos1))+1,100);
    LIMIT2=min(max(ComClu{2}(:,LimitPos2))+1,100);
    Pos=Pos+1;
    Sim{Pos}=zeros(LIMIT1,LIMIT2);
    ComNb{Pos}=zeros(LIMIT1,LIMIT2);
    for CluL1=1:LIMIT1-1
        for CluL2=1:LIMIT2-1
            ComNb{Pos}(CluL1,CluL2)=length(intersect(find(ComClu{1}(:,LimitPos1)==CluL1),find(ComClu{2}(:,LimitPos2)==CluL2)));
            Sim{Pos}(CluL1,CluL2)=round(ComNb{Pos}(CluL1,CluL2)*100/min(length(find(ComClu{1}(:,LimitPos1)==CluL1)),length(find(ComClu{2}(:,LimitPos2)==CluL2))));
        end
    end
    Limit(Pos,:)=[CLimit{1}(LimitPos1),CLimit{2}(LimitPos2)];
    %NbSup=[NbSup;[CommonLimit(LimitL),length(find(Sim{Pos}>=50)),length(find(Sim{Pos}>=75))]];
end

%display similarity
CompNb=length(Sim);
ColNb=ceil(sqrt(CompNb));
RawNb=floor(CompNb/ColNb);
for RoundL=1:2
    for CompL=1:CompNb
        h=figure;
        set(gcf,'color',[1,1,1])
        if RoundL==1
            %CurrSim=Sim{CompL}(1:LIMIT1,1:LIMIT2);
            CurrSim=Sim{CompL};
            CurrSim(find(CurrSim<70))=0;
            pcolor(CurrSim);
        elseif RoundL==2
            %pcolor(Sim{CompL}(1:LIMIT1,1:LIMIT2));
            pcolor(Sim{CompL});
        else
            %pcolor(ComNb{CompL}(1:LIMIT1,1:LIMIT2));
            pcolor(ComNb{CompL});
        end      
        title(sprintf('m%un%u %u vs m%un%u %u',ChipRank(1),NetRanks{1},Limit(CompL,1),ChipRank(2),NetRanks{2},Limit(CompL,2)))
        %ylabel(sprintf('%u - m%u',NbSup(CompL,1),ChipRank(1)))
        %xlabel(sprintf('%u - m%u',NbSup(CompL,1),ChipRank(2)))
    end
end



%% DISPLAY CORR AND ANTI PLOTS BETWEEN CORRESPONDING REGIONS
if TwinFlag
    %%display corr and anti CVM of corresponding regions
    %recover ps ranks common to corresponding regions
    TwinNb=length(TwinPos{1});
    for TwinL=1:TwinNb
        for ChipL=1:2
            CurrMcl=TwinPos{ChipL}{TwinL};
            PsPos{ChipL}=[];
            for MclL=1:length(CurrMcl)
                PsPos{ChipL}=union(PsPos{ChipL},find(ComClu{ChipL}(:,MclPos(ChipL))==CurrMcl(MclL)));
            end
        end
        PsRank{1}{TwinL}=ComPsRank{1}(intersect(PsPos{1},PsPos{2}));
        PsRank{2}{TwinL}=ComPsRank{2}(intersect(PsPos{1},PsPos{2}));
    end

    %load corr and anti values
    %display heatmap
    NetDir=[];
    for ChipL=1:2
        NetDir{ChipL}=fullfile(K.dir.net,sprintf('m%u',ChipRank(ChipL)),sprintf('n%u',NetRanks(ChipL)));
        CFile{ChipL}=sprintf('c_m%u_n%u.4mat',ChipRank(ChipL),NetRanks(ChipL));
        AFile{ChipL}=sprintf('a_m%u_n%u.4mat',ChipRank(ChipL),NetRanks(ChipL));
    end
    CorrNb=zeros(TwinNb,3);
    Colors=colors(colormap,TwinNb);
    ColNb=ceil(sqrt(TwinNb));
    RawNb=floor(TwinNb/ColNb);
    if RawNb*ColNb<TwinNb
        RawNb=RawNb+1;
    end
    hc=figure;
    set(hc,'name','FIG1 - CORR');
    set(gcf,'color',[1,1,1])
    ha=figure;
    set(ha,'name','FIG1 - ANTI');
    set(gcf,'color',[1,1,1])
    for TwinL=1:TwinNb
        for ChipL=1:2
            C{ChipL}=load_data(CFile{ChipL},NetDir{ChipL},PsNb(ChipL),PsNb(ChipL),'uint8','ieee-le',PsRank{ChipL}{TwinL},PsRank{ChipL}{TwinL});
            A{ChipL}=load_data(AFile{ChipL},NetDir{ChipL},PsNb(ChipL),PsNb(ChipL),'uint8','ieee-le',PsRank{ChipL}{TwinL},PsRank{ChipL}{TwinL});
        end
        CorrNb(TwinL,1)=length(PsRank{1}{TwinL});
        CorrNb(TwinL,2)=round(length(find(C{1}))*100/length(C{1})^2);
        CorrNb(TwinL,3)=round(length(find(C{2}))*100/length(C{2})^2);
        CorrNb(TwinL,4)=round(length(find(C{1}&C{2}))*100/length(find(C{1})));
        figure(hc)
        subplot(RawNb,ColNb,TwinL)
        plot(C{1},C{2},'.','color',Colors(TwinL,:))
        figure(ha)
        subplot(RawNb,ColNb,TwinL)
        plot(A{1},A{2},'.','color',Colors(TwinL,:))

    end
    %display for each cluster
    % the number of common Ps (col 1)
    % the percentage of CORR >0 in each chip (col 2 and 3)
    % the percentage of common CORR>0 in the two chips relatively to the number of comon Ps (col4)
    CorrNb
end
