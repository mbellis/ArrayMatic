function net_comparemcl(ModelRank,NetRank,MclStart,varargin)
%net_comparemcl([8,27],[125,118],[1,1])
%net_comparemcl([8,27],[55,24],[1,1])
%net_comparemcl([8,27],[54,23],[1,1])
%net_comparemcl([8,27],[54,23],[1,1],[2,4],{[1;2;5;6;7],[2;7;3;10;9]})
%net_comparemcl([8,42],[55,2],[1,1])
%net_comparemcl([27,42],[24,2],[1,1])
%net_comparemcl([8,27],[55,24],[1,1],[3,4],{[4,7;5,5;3,3;2,2;1,1],[1:5]'})
%net_comparemcl([8,42],[55,2],[1,1],[1,4],{[4,4;1,1;6,9;3,3;5,5;7,7;10,10;8,8],[1:8]'})
%net_comparemcl([27,42],[24,2],[1,1],[4,5],{[1,4,3,5]',[4,4,4,4;1,2,5,6;8,8,8,8;3,9,10,10]})
%net_comparemcl([8,8],[7,8],[1,1])
%net_comparemcl([8,8],[7,8],[1,1],[4,4],{[1;2;3;4;5;6;7;9],[1;3;4;2;6;7;8;5]})
%net_comparemcl([27,27],[2,3],[1,1])
%net_comparemcl([8,27],[7,2],[1,1])
% net_comparemcl([8,27],[55,24],[1,1],[4,4],{{[4,7,24,30,33,34,36,45,47,48,49,51,53,55,63,65,71,79,85,92,95,96,99,100],...
% [5,12,14,17,23,31,50,62,82],[3,15,43,76,80,87,90],[2,11,13,21,25,28,32,41,42,44,46,54,56,66,75,78,89,93],[1,6]},...
% {[1],[2],[3],[4],[5]}})
% net_comparemcl([8,27],[55,24],[1,1],'b',[3,4],{{[4,7,24,30,33,34,36,45,47,48,49,51,53,55,63,65,71,79,85,92,95,96,99,100],...
% [5,12,14,17,23,31,50,62,82],[3,15,76,80,87,90],[2,11,13,21,25,28,32,41,42,44,46,54,56,66,75,78,89,93],[1,64]},...
% {[1],[2],[3],[4],[5]}})
%net_comparemcl([8,27],[55,24],[1,1],[3,4],{{[4,7],[5],[3],[2],[1]},{[1],[2],[3],[4],[5]}})

global K
if nargin==3
    TwinFlag=0;
elseif nargin==5
    TwinFlag=1;
    MclPos=varargin{1};
    TwinPos=varargin{2};
else
    h=errordlg('needs 3 or 5 parameters');
    waitfor(h)
    error('process canceled')
end

%load probe set correspondance between two chip model
cd(K.dir.affyMetadata)
if length(ModelRank)>2 |length(NetRank)>2
    h=errordlg('at most 2 networks');
    waitfor(h)
    error('process canceled')
end
if length(ModelRank)==1
    ModelRank=[ModelRank,ModelRank];
end
if length(NetRank)==1
    NetRank=[NetRank,NetRank];
end
if ModelRank(1)~=ModelRank(2)
    FileName=sprintf('m%u_m%u_commonps.mat',min(ModelRank),max(ModelRank));
    if exist(FileName,'file')
        load(FileName)
        if ModelRank(1)>ModelRank(2)
            Temp=ComPsRank;
            ComPsRank(:,1)=ComPsRank(:,2);
            ComPsRank(:,2)=Temp(:,1);
            clear Temp
        end
    else
        h=errordlg(sprintf('no correspondance file between m%y and m%u',ModelRank(1),ModelRank(2)));
        waitfor(h)
        error('process canceled')
    end
else
    ComPsRank=repmat([1:K.chipSet.probeSetNbs{ModelRank(1)}(1)]',1,2);
end

%load MCL clusters
CorrLimit=[0:10:70];
CluSize=cell(1,2);
for ModelL=1:2
    NetDir=fullfile(K.dir.net,sprintf('m%u',ModelRank(ModelL)),sprintf('n%u',NetRank(ModelL)),'mcl');
    cd(NetDir)
    load(sprintf('m%un%u_mcl.mat',ModelRank(ModelL),NetRank(ModelL)))
    CLimit{ModelL}=CorrLimits;
    CluSize{ModelL}=cell(1,2);
    for MclL=MclStart(ModelL):size(Clu,2)
        Pos=find(CorrLimit==CorrLimits(MclL));
        for CluL=1:10
            CluSize{ModelL}{1}(CluL,Pos)=length(find(Clu(:,MclL)==CluL));
        end
    end
    if ModelRank(1)~=ModelRank(2)
        ComClu{ModelL}=Clu(ComPsRank(:,ModelL),:);        
        for MclL=MclStart(ModelL):size(Clu,2)
            Pos=find(CorrLimit==CorrLimits(MclL));
            for CluL=1:10
                CluSize{ModelL}{2}(CluL,Pos)=length(find(ComClu{ModelL}(:,MclL)==CluL));
            end
        end
    else
        ComClu{ModelL}=Clu;
        CluSize{ModelL}{2}= CluSize{ModelL}{1};
    end
    clear Clu
end
round((CluSize{1}{1}-CluSize{1}{2})*100./CluSize{1}{1})
round((CluSize{2}{1}-CluSize{2}{2})*100./CluSize{2}{1})


%calculate similarity between MCL regions
Pos=0;
Sim={};
NbSup=[];
ComNb={};
for MclL1=MclStart(1):size(ComClu{1},2)
    for MclL2=MclStart(2):size(ComClu{2},2)
        Pos=Pos+1;
        Sim{Pos}=zeros(100,100);
        ComNb{Pos}=zeros(100,100);
        for CluL1=1:10
            for CluL2=1:10
                ComNb{Pos}(CluL1,CluL2)=length(intersect(find(ComClu{1}(:,MclL1)==CluL1),find(ComClu{2}(:,MclL2)==CluL2)));
                Sim{Pos}(CluL1,CluL2)=round(ComNb{Pos}(CluL1,CluL2)*100/min(length(find(ComClu{1}(:,MclL1)==CluL1)),length(find(ComClu{2}(:,MclL2)==CluL2))));
            end
        end
        NbSup=[NbSup;[MclL1,MclL2,length(find(Sim{Pos}>=50)),length(find(Sim{Pos}>=75))]];
    end
end


if TwinFlag==0
    %display similarity between MCL regions
    CompNb=length(Sim);
    ColNb=ceil(sqrt(CompNb));
    RawNb=floor(CompNb/ColNb);
    LIMIT=11;
    if ColNb*RawNb<CompNb
        RawNb=RawNb+1;
    end
    for RoundL=1:3
        h=figure;
        set(gcf,'color',[1,1,1])
        for CompL=1:CompNb
            subplot(RawNb,ColNb,CompL)
            if RoundL==1
                CurrSim=Sim{CompL}(1:LIMIT,1:LIMIT);
                CurrSim(find(CurrSim<70))=0;
                pcolor(CurrSim);
            elseif RoundL==2
                pcolor(Sim{CompL}(1:LIMIT,1:LIMIT));
            else
                pcolor(ComNb{CompL}(1:LIMIT,1:LIMIT));
            end
            ylabel(sprintf('%u',CLimit{1}(NbSup(CompL,1))))
            xlabel(sprintf('%u',CLimit{2}(NbSup(CompL,2))))
        end
    end
end



%
% Pos=0;
% Sim={};
% NbSup=[];
% ComNb={};
% for MclL1=MclStart(1):size(ComClu{1},2)
%     if size(ComClu{1},2)>=MclL1+1
%     for MclL2=MclL1+1
%         Pos=Pos+1;
%         Sim{Pos}=zeros(100,100);
%         ComNb{Pos}=zeros(100,100);
%         for CluL1=1:100
%             for CluL2=1:100
%                 ComNb{Pos}(CluL1,CluL2)=length(intersect(find(ComClu{1}(:,MclL1)==CluL1),find(ComClu{2}(:,MclL2)==CluL2)));
%                 Sim{Pos}(CluL1,CluL2)=round(ComNb{Pos}(CluL1,CluL2)*100/min(length(find(ComClu{1}(:,MclL1)==CluL1)),length(find(ComClu{2}(:,MclL2)==CluL2))));
%             end
%         end
%         NbSup=[NbSup;[MclL1,MclL2,length(find(Sim{Pos}>=50)),length(find(Sim{Pos}>=75))]];
%     end
%     end
% end
% % 
% % % if TwinFlag==0
% %     %display similarity between MCL regions
%     CompNb=length(Sim);
%     ColNb=ceil(sqrt(CompNb));
%     RawNb=floor(CompNb/ColNb);
%     LIMIT=100;
%     if ColNb*RawNb<CompNb
%         RawNb=RawNb+1;
%     end
%     for RoundL=1:2        
%         for CompL=1:CompNb
%             h=figure;
%             set(gcf,'color',[1,1,1])
%             %subplot(RawNb,ColNb,CompL)
%             if RoundL==1
%                 CurrSim=Sim{CompL}(1:LIMIT,1:LIMIT);
%                 CurrSim(find(CurrSim<70))=0;
%                 pcolor(CurrSim);
%             elseif RoundL==2
%                 pcolor(Sim{CompL}(1:LIMIT,1:LIMIT));
%             else
%                 pcolor(ComNb{CompL}(1:LIMIT,1:LIMIT));
%             end
%             ylabel(sprintf('%u',CLimit{1}(NbSup(CompL,1))))
%             xlabel(sprintf('%u',CLimit{2}(NbSup(CompL,2))))
%         end
%     end
% end
% 

if TwinFlag
    %%display corr and anti CVM of corresponding regions
    %recover ps ranks common to corresponding regions
    TwinNb=length(TwinPos{1});
    for TwinL=1:TwinNb
        for ModelL=1:2
            CurrMcl=TwinPos{ModelL}{TwinL};
            PsPos{ModelL}=[];
            for MclL=1:length(CurrMcl)
                PsPos{ModelL}=union(PsPos{ModelL},find(ComClu{ModelL}(:,MclPos(ModelL))==CurrMcl(MclL)));
            end
        end
        PsRank{1}{TwinL}=ComPsRank(intersect(PsPos{1},PsPos{2}),1);
        PsRank{2}{TwinL}=ComPsRank(intersect(PsPos{1},PsPos{2}),2);
    end

    %load corr and anti values
    %display heatmap
    PsNb=[];
    NetDir=[];
    for ModelL=1:2
        PsNb(ModelL)=K.chipSet.probeSetNbs{ModelRank(ModelL)}(1);
        NetDir{ModelL}=fullfile(K.dir.net,sprintf('m%03u',ModelRank(ModelL)),sprintf('n%05u',NetRank(ModelL)));
        CFile{ModelL}=sprintf('c_m%u_n%u.4mat',ModelRank(ModelL),NetRank(ModelL));
        AFile{ModelL}=sprintf('a_m%u_n%u.4mat',ModelRank(ModelL),NetRank(ModelL));
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
    %for TwinL=4
        for ModelL=1:2
            C{ModelL}=load_data(CFile{ModelL},NetDir{ModelL},PsNb(ModelL),PsNb(ModelL),'uint8','ieee-le',PsRank{ModelL}{TwinL},PsRank{ModelL}{TwinL});
            A{ModelL}=load_data(AFile{ModelL},NetDir{ModelL},PsNb(ModelL),PsNb(ModelL),'uint8','ieee-le',PsRank{ModelL}{TwinL},PsRank{ModelL}{TwinL});
        end
        %A{ModelL}=load_data(AFile{ModelL},NetDir{ModelL},PsNb(ModelL),PsNb(ModelL),'uint8','ieee-le',PsRank{ModelL}{TwinL},PsRank{ModelL}{TwinL});
        CorrNb(TwinL,1)=round(length(find(C{1}))*100/length(C{1})^2);
        CorrNb(TwinL,2)=round(length(find(C{2}))*100/length(C{2})^2);
        CorrNb(TwinL,3)=round(length(find(C{1}&C{2}))*100/min(length(find(C{1})),length(find(C{2}))));
%         %eliminate uncommon corr
%         Uncommon=find(C{1}.*C{2}==0|A{1}.*A{2}==0);
%         C{1}(Uncommon)=0;
%         C{2}(Uncommon)=0;
%         A{1}(Uncommon)=0;
%         A{2}(Uncommon)=0;
        figure(hc)
        subplot(RawNb,ColNb,TwinL)
        plot(C{1},C{2},'.','color',Colors(TwinL,:))
        figure(ha)
        subplot(RawNb,ColNb,TwinL)
        plot(A{1},A{2},'.','color',Colors(TwinL,:))
        %write data
        for ModelL=1:2
        %for ModelL=2
            cd(fullfile(NetDir{ModelL},'mcl'))
            CurrC=single(C{ModelL})/100;
            CurrA=single(A{ModelL})/100;

            %print WGCNA input files
            Diff=CurrC-CurrA;
            for PsL=1:length(Diff)
                Diff(PsL,PsL)=0;
            end
            FileName=sprintf('wgcna_m%un%uc%umcl%u_alldiff.csv',ModelRank(ModelL),NetRank(ModelL),MclPos(ModelL),TwinL);
            FileName=strrep(FileName,'__','_');
            fid=fopen(FileName,'w');
                        for PsL=1:length(Diff)
                fprintf(fid,'%.2f,',Diff(PsL,1:end-1));
                fprintf(fid,'%.2f\n',Diff(PsL,end));
            end
            fclose(fid)
            
            %eliminate uncommon corr
            Uncommon=find(C{1}.*C{2}==0|A{1}.*A{2}==0);
            CurrC(Uncommon)=0;
            CurrA(Uncommon)=0;
            Diff=CurrC-CurrA;            
            for PsL=1:length(Diff)
                Diff(PsL,PsL)=0;
            end
            FileName=sprintf('wgcna_m%un%uc%umcl%u_alldiffcom.csv',ModelRank(ModelL),NetRank(ModelL),MclPos(ModelL),TwinL);
            FileName=strrep(FileName,'__','_');
            fid=fopen(FileName,'w');
                        for PsL=1:length(Diff)
                fprintf(fid,'%.2f,',Diff(PsL,1:end-1));
                fprintf(fid,'%.2f\n',Diff(PsL,end));
            end
            fclose(fid)

            Diff=CurrC-CurrA;
            for PsL=1:length(Diff)
                Diff(PsL,PsL)=0;
            end
            Diff(find(Diff<0))=0;
            FileName=sprintf('wgcna_m%un%uc%umcl%u_posdiff.csv',ModelRank(ModelL),NetRank(ModelL),MclPos(ModelL),TwinL);
            FileName=strrep(FileName,'__','_');
            fid=fopen(FileName,'w');
                        for PsL=1:length(Diff)
                fprintf(fid,'%.2f,',Diff(PsL,1:end-1));
                fprintf(fid,'%.2f\n',Diff(PsL,end));
            end
            fclose(fid)        
            
            Diff=CurrC-CurrA;
            for PsL=1:length(Diff)
                Diff(PsL,PsL)=0;
            end
            Diff(find(Diff<0))=0;            
            Diff(Uncommon)=0;
            FileName=sprintf('wgcna_m%un%uc%umcl%u_posdiffcom.csv',ModelRank(ModelL),NetRank(ModelL),MclPos(ModelL),TwinL);
            FileName=strrep(FileName,'__','_');
            fid=fopen(FileName,'w');
                        for PsL=1:length(Diff)
                fprintf(fid,'%.2f,',Diff(PsL,1:end-1));
                fprintf(fid,'%.2f\n',Diff(PsL,end));
            end
            fclose(fid)                                            
            %print GAPinput files
            
            Diff=CurrC-CurrA;
            for PsL=1:length(Diff)
                Diff(PsL,PsL)=0;
            end
            FileName=sprintf('gap_m%un%uc%umcl%u_alldiff.csv',ModelRank(ModelL),NetRank(ModelL),MclPos(ModelL),TwinL);
            FileName=strrep(FileName,'__','_');
            fid=fopen(FileName,'w');
            %header
            fprintf(fid,'ID')
            fprintf(fid,'\tV%u,',1:length(Diff))
            fprintf(fid,'\n')
            for PsL=1:length(Diff)
                fprintf(fid,'%u',PsL)
                fprintf(fid,'\t%.2f',Diff(PsL,1:end));
                fprintf(fid,'\n')
            end
            fclose(fid)


            %eliminate uncommon corr
            Uncommon=find(C{1}.*C{2}==0|A{1}.*A{2}==0);
            CurrC(Uncommon)=0;
            CurrA(Uncommon)=0;
            Diff=CurrC-CurrA;            
            for PsL=1:length(Diff)
                Diff(PsL,PsL)=0;
            end
            FileName=sprintf('gap_m%un%uc%umcl%u_alldiffcom.csv',ModelRank(ModelL),NetRank(ModelL),MclPos(ModelL),TwinL);
            FileName=strrep(FileName,'__','_');
            fid=fopen(FileName,'w');
            %header
            fprintf(fid,'ID')
            fprintf(fid,'\tV%u,',1:length(Diff))
            fprintf(fid,'\n')
            for PsL=1:length(Diff)
                fprintf(fid,'%u',PsL)
                fprintf(fid,'\t%.2f',Diff(PsL,1:end));
                fprintf(fid,'\n')
            end
            fclose(fid)

        end
    end
    CorrNb
end
