function net_compareclique(ModelRank,NetRank,Corr)
%net_compareclique([8,27],[55,24],[0,0])
%net_compareclique([8,27],[55,24],[50,50])
global K
if nargin==3
else
    h=errordlg('needs 3');
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

%load cliques
for ModelL=1:2
    NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank(ModelL)),sprintf('n%05u',NetRank(ModelL)));
    cd(NetDir)
    load(sprintf('m%un%u_cliques_%02u.mat',ModelRank(ModelL),NetRank(ModelL),Corr(ModelL)))
    if ModelRank(1)~=ModelRank(2)
        ComClu{ModelL}=Clu(ComPsRank(:,ModelL),:);
    else
        ComClu{ModelL}=Clu;
    end    
    clear Clu
end

%load MCL regions
for ModelL=1:2
    NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank(ModelL)),sprintf('n%05u',NetRank(ModelL)),'mcl');
    cd(NetDir)
    load(sprintf('m%un%u_mcl.mat',ModelRank(ModelL),NetRank(ModelL)))
    if ModelRank(1)~=ModelRank(2)
        ComReg{ModelL}=Clu(ComPsRank(:,ModelL),:);
    else
        ComReg{ModelL}=Clu;
    end
    CorrLimit{ModelL}=CorrLimits;
    clear Clu
end

%calculate similarity between cliques
NbSup=[];
ComNb={};
Sim=zeros(100,100);
ComNb=zeros(100,100);
for CluL1=1:100
    for CluL2=1:100
        ComNb(CluL1,CluL2)=length(intersect(find(ComClu{1}==CluL1),find(ComClu{2}==CluL2)));
        Sim(CluL1,CluL2)=round(ComNb(CluL1,CluL2)*100/min(length(find(ComClu{1}==CluL1)),length(find(ComClu{2}==CluL2))));    
        NbSup=[NbSup;[CluL1,CluL2,length(find(Sim>=50)),length(find(Sim>=75))]];
    end
end





    %display similarity between MCL regions   
    for RoundL=1:3
        h=figure;
        set(gcf,'color',[1,1,1])        
            if RoundL==1                
                CurrSim=Sim;
                CurrSim(find(CurrSim<70))=0;
                pcolor(CurrSim);
                
            elseif RoundL==2
                pcolor(Sim);
                title(sprintf('m%un%u vs m%un%u - overlap',ModelRank(1),ModelRank(1),NetRank(1),NetRank(1)))
            else
                pcolor(ComNb);
                title(sprintf('m%un%u vs m%un%u - ps nb',ModelRank(1),ModelRank(1),NetRank(1),NetRank(1)))
            end
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

    for ModelL=1:2
        Pos=find(CorrLimit{ModelL}==Corr(ModelL));
        if isempty(Pos)
            CorrPos(ModelL)=1;
        else
            CorrPos(ModelL)=Pos;
        end
    end
    
    if ModelRank(1)~=ModelRank(2)|Netrank(1)~=NetRank(2)
    
    PsPos=find(ComClu{1}<=50&ComClu{1}>0);
    [temp,SortIndex]=sort(ComClu{1}(PsPos));
    PsRank{1}=ComPsRank(PsPos(SortIndex),1);
    PsRank{2}=ComPsRank(PsPos(SortIndex),2);
    [Region,SortIndex]=sort(ComReg{1}(PsPos(SortIndex),CorrPos(1)));
    RegionRank=unique(Region);
    RegionRank=setdiff(RegionRank,0);
    RegLabel=cell(length(RegionRank),1);
    RegTick=ones(length(RegionRank),1);
    %RegPos=zeros(length(RegionRank),1);
    for RegL=1:length(RegionRank)
        RegLabel{RegL}=sprintf('%u',RegionRank(RegL));
        CurrRegPos=find(Region==RegionRank(RegL));
        if RegL<length(RegionRank)
            RegTick(RegL+1)=CurrRegPos(end);
        end
    end
    PsRank{1}=PsRank{1}(SortIndex);
    PsRank{2}=PsRank{2}(SortIndex);

    for ModelL=1:2
        C{ModelL}=load_data(CFile{ModelL},NetDir{ModelL},PsNb(ModelL),PsNb(ModelL),'uint8','ieee-le',PsRank{ModelL},PsRank{ModelL});
        A{ModelL}=load_data(AFile{ModelL},NetDir{ModelL},PsNb(ModelL),PsNb(ModelL),'uint8','ieee-le',PsRank{ModelL},PsRank{ModelL});
    end
    
    hc=figure;
    set(hc,'name','FIG1 - CORR');
    set(gcf,'color',[1,1,1])
    ha=figure;
    set(ha,'name','FIG1 - ANTI');
    set(gcf,'color',[1,1,1])
    figure(hc)
    subplot(1,2,1)
    image(C{1})
    set(gca,'xtick',RegTick)
    set(gca,'xticklabel',RegLabel)
    set(gca,'ytick',RegTick)
    set(gca,'yticklabel',RegLabel)
    set(gca,'tickdir','out')
    subplot(1,2,2)
    image(C{2})
    set(gca,'xtick',RegTick)
    set(gca,'xticklabel',RegLabel)
    set(gca,'ytick',RegTick)
    set(gca,'yticklabel',RegLabel)
    set(gca,'tickdir','out')
    set(hc,'position',[5 365 1266 528])
    figure(ha)
    subplot(1,2,1)
    image(A{1})
    set(gca,'xtick',RegTick)
    set(gca,'xticklabel',RegLabel)
    set(gca,'ytick',RegTick)
    set(gca,'yticklabel',RegLabel)
    set(gca,'tickdir','out')
    subplot(1,2,2)
    image(A{2})
    set(gca,'xtick',RegTick)
    set(gca,'xticklabel',RegLabel)
    set(gca,'ytick',RegTick)
    set(gca,'yticklabel',RegLabel)
    set(gca,'tickdir','out')
    set(ha,'position',[5 365 1266 528])


    PsPos=find(ComClu{2}<=50&ComClu{2}>0);
    [temp,SortIndex]=sort(ComClu{2}(PsPos));
    PsRank{1}=ComPsRank(PsPos(SortIndex),1);
    PsRank{2}=ComPsRank(PsPos(SortIndex),2);
    [Region,SortIndex]=sort(ComReg{2}(PsPos(SortIndex),CorrPos(1)));
    PsRank{1}=PsRank{1}(SortIndex);
    PsRank{2}=PsRank{2}(SortIndex);
    RegionRank=unique(Region);
    RegionRank=setdiff(RegionRank,0);
    RegLabel=cell(length(RegionRank),1);
    RegTick=ones(length(RegionRank),1);
    %RegPos=zeros(length(RegionRank),1);
    for RegL=1:length(RegionRank)
        RegLabel{RegL}=sprintf('%u',RegionRank(RegL));
        CurrRegPos=find(Region==RegionRank(RegL));
        if RegL<length(RegionRank)
            RegTick(RegL+1)=CurrRegPos(end);
        end
    end


    for ModelL=1:2
        C{ModelL}=load_data(CFile{ModelL},NetDir{ModelL},PsNb(ModelL),PsNb(ModelL),'uint8','ieee-le',PsRank{ModelL},PsRank{ModelL});
        A{ModelL}=load_data(AFile{ModelL},NetDir{ModelL},PsNb(ModelL),PsNb(ModelL),'uint8','ieee-le',PsRank{ModelL},PsRank{ModelL});
    end
    hc=figure;
    set(hc,'name','FIG2 - CORR');
    set(gcf,'color',[1,1,1])
    ha=figure;
    set(ha,'name','FIG2 - ANTI');
    set(gcf,'color',[1,1,1])
    figure(hc)
    subplot(1,2,1)
    image(C{1})
    set(gca,'xtick',RegTick)
    set(gca,'xticklabel',RegLabel)
    set(gca,'ytick',RegTick)
    set(gca,'yticklabel',RegLabel)
    set(gca,'tickdir','out')
    subplot(1,2,2)
    image(C{2})
    set(gca,'xtick',RegTick)
    set(gca,'xticklabel',RegLabel)
    set(gca,'ytick',RegTick)
    set(gca,'yticklabel',RegLabel)
    set(gca,'tickdir','out')
    set(hc,'position',[5 365 1266 528])
    figure(ha)
    subplot(1,2,1)
    image(A{1})
    set(gca,'xtick',RegTick)
    set(gca,'xticklabel',RegLabel)
    set(gca,'ytick',RegTick)
    set(gca,'yticklabel',RegLabel)
    set(gca,'tickdir','out')
    subplot(1,2,2)
    image(A{2})
    set(gca,'xtick',RegTick)
    set(gca,'xticklabel',RegLabel)
    set(gca,'ytick',RegTick)
    set(gca,'yticklabel',RegLabel)
    set(gca,'tickdir','out')
    set(ha,'position',[5 365 1266 528])
    else
        RegTick=[];
        RegionRank=[];
        for ModelL=1:2
            PsPos=find(ComClu{ModelL}<=50&ComClu{ModelL}>0);
            [temp,SortIndex]=sort(ComClu{ModelL}(PsPos));
            PsRank{ModelL}=ComPsRank(PsPos(SortIndex),ModelL);
            [Region,SortIndex]=sort(ComReg{ModelL}(PsPos(SortIndex),CorrPos(ModelL)));
            RegionRank{ModelL}=unique(Region);
            RegionRank{ModelL}=setdiff(RegionRank{ModelL},0);
            RegLabel{ModelL}=cell(length(RegionRank{ModelL}),1);
            RegTick{ModelL}=ones(length(RegionRank{ModelL}),1);
            %RegPos=zeros(length(RegionRank),1);
            for RegL=1:length(RegionRank{ModelL})
                RegLabel{ModelL}{RegL}=sprintf('%u',RegionRank{ModelL}(RegL));
                CurrRegPos=find(Region==RegionRank{ModelL}(RegL));
                if RegL<length(RegionRank{ModelL})
                    RegTick{ModelL}(RegL+1)=CurrRegPos(end);
                end
            end
            PsRank{ModelL}=PsRank{ModelL}(SortIndex);
        end
    
   


        C=load_data(CFile{1},NetDir{1},PsNb(1),PsNb(1),'uint8','ieee-le',PsRank{1},PsRank{2});
        A=load_data(AFile{1},NetDir{1},PsNb(1),PsNb(1),'uint8','ieee-le',PsRank{1},PsRank{2});

    
    hc=figure;
    set(hc,'name','FIG1 - CORR & ANTI');
    set(gcf,'color',[1,1,1])
    figure(hc)
    subplot(1,2,1)
    image(C)
    set(gca,'xtick',RegTick{2})
    set(gca,'xticklabel',RegLabel{2})
    set(gca,'ytick',RegTick{1})
    set(gca,'yticklabel',RegLabel{1})
    set(gca,'tickdir','out')
    subplot(1,2,2)
    image(A)
    set(gca,'xtick',RegTick{2})
    set(gca,'xticklabel',RegLabel{2})
    set(gca,'ytick',RegTick{1})
    set(gca,'yticklabel',RegLabel{1})
    set(gca,'tickdir','out')
    set(hc,'position',[5 365 1266 528])

    end
'stop'