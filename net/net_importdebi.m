%=========================%
% FUNCTION NET_IMPORTDEBI %
%=========================%

%NET_IMPORTDEBI: import several results of DEBI bi-clustering
%INPUT PARAMETERS
% 1  ChipRank: chip rank
% 2  DebiType: either 'exp' (intra experimental comparisons)
%              or 'nr' (not redundant biological conditions)
%              or 'biol' (all biological conditions)
% 3  NetRanks: list of net ranks
% 4  Fdr: Fdr limit used to discretize the fdr values (0 or 1)
% 5 NrRank : non redundant used (if DebiType='nr')

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

%net_importdebi(2,'exp',[80],1000,0)
%net_importdebi(2,'exp',[80],100,0)
%net_importdebi(2,'nr',[1],1000,1)
%net_importdebi(2,'nr',[1],100,1)
%net_importdebi(6,'nr',[64],100,1)
%net_importdebi(27,'nr',[164],1000,51)

function net_importdebi(ChipRank,DebiType,NetRank,Fdr,NrRank)
global K
switch ChipRank
    case 2
        ChipName='affy_hs54k_01';
    case 3
        ChipName='affy_hs39k_01';
    case 5
        ChipName='affy_md36k_01';
    case 6
        ChipName='affy_rn27k_01';
    case 8
        ChipName='affy_md34k_01';
    case 27
        ChipName='affy_md45k_01';
end
cd('/home/mbellis/sosma/raydata/mlf/point/')
%TableData=fullfile('/home/mbellis/sosma/raydata/mlf/table',ChipName,'sosma');
TableData=fullfile('/home/mbellis/array1/sosma/raydata/mlf/table',ChipName,'sosma');
eval(sprintf('load m%u_point',ChipRank))

ChipPos=find(K.chip.rank==ChipRank);
PsNb=K.chip.probesetNb(ChipPos);
if isempty(NetRank)
    LoopNb=1;
else
    LoopNb=length(NetRank);
end

%% RECOVER CLUSTERS
mkdir(fullfile('/home/mbellis/array1/sosma/net',sprintf('m%u',ChipRank),'debi'))
for LoopL=1:LoopNb
    cd('/home/mbellis/array1/sosma/debi')
    if isequal(DebiType,'biol')
        %construct correspondance between comparisons and column
        BiolRank1=K.net{ChipRank}.biolRank{NetRank(LoopL)}{1};
        BiolRank2=K.net{ChipRank}.biolRank{NetRank(LoopL)}{2};
        CommonBiol=intersect(BiolRank1,BiolRank2);
        if ~isempty(CommonBiol)
            BiolRank2=setdiff(BiolRank2,CommonBiol);
        end
        ColPos=zeros(length(BiolRank1),length(BiolRank2));
        CompNb=0;
        for BiolL1=1:length(BiolRank1)
            for BiolL2=1:length(BiolRank2)
                CompNb=CompNb+1;
                ColPos(BiolL1,BiolL2)=CompNb;
            end
        end
    end
    %count the number of clusters and the number of comparisons
    %the first line count the number of clusters and the number of comparisons
    switch DebiType
        case 'biol'
            fid=fopen(sprintf('m%un%u_fdr%u_debi.txt',ChipRank,NetRank(LoopL),Fdr(LoopL)),'r');
            OutMat=sprintf('m%un%u_fdr%u_debi',ChipRank,NetRank(LoopL),Fdr(LoopL));
        case 'exp'
            fid=fopen(sprintf('m%uexp%u_fdr%u_debi.txt',ChipRank,NetRank(LoopL),Fdr(LoopL)),'r');
            OutMat=sprintf('m%uexp%u_fdr%u_debi',ChipRank,NetRank(LoopL),Fdr(LoopL));
        case 'nr'
            fid=fopen(sprintf('m%unr%u_fdr%u_debi.txt',ChipRank,NrRank(LoopL),Fdr(LoopL)),'r');
            OutMat=sprintf('m%unr%u_fdr%u_debi',ChipRank,NrRank(LoopL),Fdr(LoopL));
    end
    LineNb=0;
    if isequal(DebiType,'biol')
        while 1
            CurrLine=fgetl(fid);
            if ~ischar(CurrLine)
                break
            end
            LineNb=LineNb+1;

        end
    else
        CompRank=[];
        while 1
            CurrLine=fgetl(fid);
            if ~ischar(CurrLine)
                break
            end
            LineNb=LineNb+1;
            if mod(LineNb,3)==0
                CurrComp=str2num(CurrLine);                
                CompRank=union(CompRank,CurrComp);
            end
        end
        fclose(fid)
        CompNb=length(CompRank);
    end
    CluNb=LineNb/3;
    
    Debi{LoopL}.clu2comp=uint8(zeros(CluNb,CompNb));
    Debi{LoopL}.clu2ps=uint8(zeros(CluNb,PsNb));
    
    switch DebiType
        case 'biol'
            fid=fopen(sprintf('m%un%u_fdr%u_debi.txt',ChipRank,NetRank(LoopL),Fdr(LoopL)),'r');
        case 'exp'
            fid=fopen(sprintf('m%uexp%u_fdr%u_debi.txt',ChipRank,NetRank(LoopL),Fdr(LoopL)),'r');
        case 'nr'
            fid=fopen(sprintf('m%unr%u_fdr%u_debi.txt',ChipRank,NrRank(LoopL),Fdr(LoopL)),'r');
    end

    CluPos=0;
    LinePos=0;
    while 1        
        CurrLine=fgetl(fid);
        if ~ischar(CurrLine)
            break
        end
        LinePos=LinePos+1;
        switch(mod(LinePos,3))
            case 0
                if isequal(DebiType,'biol')
                    %recover comparisons
                    CurrCol=regexp(CurrLine,'\w+','match');
                    CurrComp=zeros(1,length(CurrCol));
                    for CompL=1:length(CurrComp)
                        Biol1=regexp(CurrCol{CompL},'\d+(?=_)','match');
                        Biol2=regexp(CurrCol{CompL},'(?<=_)\d+','match');
                        Biol1=str2num(Biol1{1});
                        Biol2=str2num(Biol2{1});
                        CurrComp(CompL)=ColPos(Biol1,Biol2);
                    end                    
                else
                    CurrComp=str2num(CurrLine);
                    for CompL=1:length(CurrComp)
                        CurrComp(CompL)=find(CompRank==CurrComp(CompL));
                    end
                end
                %comparisons where where cluster exist (nb of clusters x nb of comparisons)
                Debi{LoopL}.clu2comp(CluPos,CurrComp)=1;
            case 1
                %recover information                                
                CluPos=CluPos+1;
            case 2
                %recover probe sets rank
                %CurrLine=regexp(CurrLine,'\d+','match'); 
                CurrLine=str2num(CurrLine); 
                CurrPs=zeros(length(CurrLine),1);
                for PsL=1:length(CurrLine)
                    CurrPs(PsL)=CurrLine(PsL);
                end
                %
                Debi{LoopL}.clu2ps(CluPos,CurrPs)=1;
        end
    end
    fclose(fid)    
end

%% DISPLAY CLUSTERS

%sort clusters to have in first position the one with the largest number of comparisons
for LoopL=1:LoopNb
    CompSize=sum(Debi{LoopL}.clu2comp,2);
    [CompSize,SortIndex]=sort(CompSize,'descend');
    Debi{LoopL}.clu2comp=Debi{LoopL}.clu2comp(SortIndex,:);
    Debi{LoopL}.clu2ps=Debi{LoopL}.clu2ps(SortIndex,:);    
    PsSize=sum(Debi{LoopL}.clu2ps,2);
    [PsSize,SortIndex]=sort(PsSize,'descend');
    Debi{LoopL}.clu2comp=Debi{LoopL}.clu2comp(SortIndex,:);
    Debi{LoopL}.clu2ps=Debi{LoopL}.clu2ps(SortIndex,:);            
end



%process each cluster to allow a graphical representation (a PsNbxCompNb matrix showing the
%relationships between biclusters
for LoopL=1:LoopNb
    Debi{LoopL}.ps2comp=zeros(PsNb,CompNb);
    Debi{LoopL}.psOrder=[];
    Debi{LoopL}.compOrder=[];    
    % Comparisons already used
    UsedComp=[];
    % Probe sets already used
    UsedPs=[];
    UsedCompNb=0;
    UsedPsNb=0;
    CompOrder=zeros(1,CompNb);
    PsOrder=zeros(PsNb,1);
    PsTick=[0];
    CompTick=[0];
    CluNb=size(Debi{LoopL}.clu2comp,1);
    for CluL=1:CluNb
        CurrPs=find(Debi{LoopL}.clu2ps(CluL,:));
        CurrComp=find(Debi{LoopL}.clu2comp(CluL,:));
        CurrNewPs=setdiff(CurrPs,UsedPs);
        UsedPs=union(UsedPs,CurrNewPs);
        %add a new tick if exist new probe sets
        if length(UsedPs)>PsTick(end)
            PsTick(end+1)=length(UsedPs);
        end
        CurrNewComp=setdiff(CurrComp,UsedComp);
        UsedComp=union(UsedComp,CurrNewComp);
        %add new tick if exist new comparisons
        if length(UsedComp)>CompTick(end)
            CompTick(end+1)=length(UsedComp);
        end
        %indicates the position of new comparisons
        CompOrder(CurrNewComp)=UsedCompNb+1:UsedCompNb+length(CurrNewComp);
        Debi{LoopL}.compOrder=[Debi{LoopL}.compOrder,CurrNewComp];
        %update UsedCompNb
        UsedCompNb=UsedCompNb+length(CurrNewComp);
        %indicates the position of new probe sets
        PsOrder(CurrNewPs)=UsedPsNb+1:UsedPsNb+length(CurrNewPs);
        Debi{LoopL}.psOrder=[Debi{LoopL}.psOrder,CurrNewPs];
        %update UsedPsNb
        UsedPsNb=UsedPsNb+length(CurrNewPs);
        %indicates the redundancy of probe sets in clusters
        Debi{LoopL}.ps2comp(PsOrder(CurrPs),CompOrder(CurrComp))=Debi{LoopL}.ps2comp(PsOrder(CurrPs),CompOrder(CurrComp))+1;
    end
    
    %remove empty ps
    KeepPs=find(sum(Debi{LoopL}.ps2comp,2)>0);
    Debi{LoopL}.ps2comp=Debi{LoopL}.ps2comp(KeepPs,:);
    
    %reorder Ps and Comp
    for PsL=1:length(PsTick)-1
        Interval=PsTick(PsL)+1:PsTick(PsL+1);
        [temp SortOrder]=sort(sum(Debi{LoopL}.ps2comp(Interval,:),2));
        Debi{LoopL}.ps2comp(Interval,:)=Debi{LoopL}.ps2comp(Interval(SortOrder),:);
        Debi{LoopL}.psOrder(Interval)=Debi{LoopL}.psOrder(Interval(SortOrder));
    end
    for CompL=1:length(CompTick)-1
        Interval=CompTick(CompL)+1:CompTick(CompL+1);
        [temp SortOrder]=sort(sum(Debi{LoopL}.ps2comp(:,Interval)));
        Debi{LoopL}.ps2comp(:,Interval)=Debi{LoopL}.ps2comp(:,Interval(SortOrder));
        Debi{LoopL}.compOrder(Interval)=Debi{LoopL}.compOrder(Interval(SortOrder));
    end
    
    %refine the order
    UsedPs=[];
    UsedComp=[];
    MemPsL1=2;
    %process the clusters that are round the diagonal
    for CompL1=2:min(length(PsTick)-1,length(CompTick)-1)    
        CompInterval=CompTick(CompL1)+1:CompTick(CompL1+1);
        %find the next Ps interval
        for PsL1=MemPsL1:min(length(PsTick)-1,length(CompTick)-1)        
            PsInterval=PsTick(PsL1)+1:PsTick(PsL1+1);
            %process the current cluster only if it has not probe sets and comparisons used in
            %previous clusters
            if mean(mean(Debi{LoopL}.ps2comp(PsInterval,CompInterval)))>1
                %update MemPsL
                MemPsL1=PsL1+1;

                %find the same probe sets in earlier clusters
                for PsL2=1:PsL1-1                    
                    CurrPsInterval=PsTick(PsL2)+1:PsTick(PsL2+1);
                    SamePsPos=find(min(Debi{LoopL}.ps2comp(CurrPsInterval,CompInterval),[],2)>0);
                    %keep new ps
                    SamePs=setdiff(CurrPsInterval(SamePsPos),UsedPs);
                    if ~isempty(SamePs)                        
                        %reorder ps in the current cluster in last position
                        FirstPs=setdiff(setdiff(CurrPsInterval,UsedPs),SamePs);
                        Debi{LoopL}.ps2comp(setdiff(CurrPsInterval,UsedPs),:)=[Debi{LoopL}.ps2comp(FirstPs,:);Debi{LoopL}.ps2comp(SamePs,:)];                 
                        Debi{LoopL}.psOrder(setdiff(CurrPsInterval,UsedPs))=[Debi{LoopL}.psOrder(FirstPs),Debi{LoopL}.psOrder(SamePs)];
                        %update UsedPs
                        UsedPs=union(UsedPs,CurrPsInterval(length(FirstPs)+1:end));                        
                    end
                end

                %find the same comparisons in earlier clusters
                for CompL2=1:CompL1-1
                    CurrCompInterval=CompTick(CompL2)+1:CompTick(CompL2+1);
                    SameCompPos=find(min(Debi{LoopL}.ps2comp(PsInterval,CurrCompInterval))>0);
                    %keep new comp
                    SameComp=setdiff(CurrCompInterval(SameCompPos),UsedComp);
                    if ~isempty(SameComp)
                        %reorder comp in the current cluster in last position
                        FirstComp=setdiff(setdiff(CurrCompInterval,UsedComp),SameComp);
                        Debi{LoopL}.ps2comp(:,setdiff(CurrCompInterval,UsedComp))=[Debi{LoopL}.ps2comp(:,FirstComp),Debi{LoopL}.ps2comp(:,SameComp)];
                        Debi{LoopL}.compOrder(setdiff(CurrCompInterval,UsedComp))=[Debi{LoopL}.compOrder(FirstComp),Debi{LoopL}.compOrder(SameComp)];
                        %update UsedComp
                        UsedComp=union(UsedComp,CurrCompInterval(length(FirstComp)+1:end));
                    end
                end
                break
            end
        end
    end

 
    figure
    set(gcf,'color',[1,1,1])
    Val=Debi{LoopL}.ps2comp;
    %Val(find(Val>5))=5;
    Colors=colors(colormap,6);
    CurrColormap=colormap;
    colormap(Colors);
    image(Val)       
    title(sprintf('m%un%u fdr%u %u/%u ps %u/%u comp',ChipRank,NetRank(LoopL),Fdr(LoopL),length(KeepPs),PsNb,size(Debi{LoopL}.ps2comp,2),CompNb));
    set(gca,'xtick',CompTick)
    set(gca,'ytick',PsTick)
    set(gca,'tickdir','out')
    set(gca,'xticklabel','')
    set(gca,'yticklabel','')
    xlabel('comparisons')
    ylabel('probe sets')
    for LineL=1:20
        line([0,CompNb],[PsTick(LineL),PsTick(LineL)])
    end
    for LineL=1:20
        line([CompTick(LineL),CompTick(LineL)],[0,length(KeepPs)])
    end

    %image of CVM ordered according to biclustering
    %without repetition of probe sets
%     PsOrder=Debi{LoopL}.psOrder;    
%     PsOrder=PsOrder(find(PsOrder>0));        
%     switch DebiType
%         case 'biol'
%             net_displaycvm(ChipRank,NetRank(LoopL),PsOrder,sprintf('m%un%u_fdr%u_debi',ChipRank,NetRank(LoopL),Fdr(LoopL)),100,1000,1)
%         case 'exp'
%             net_displaycvm(ChipRank,NetRank(LoopL),PsOrder,sprintf('m%uexp%u_fdr%u_debi',ChipRank,NetRank(LoopL),Fdr(LoopL)),100,1000,1)            
%         case 'nr'
%             net_displaycvm(ChipRank,NetRank(LoopL),PsOrder,sprintf('m%unr%u_fdr%u_debi',ChipRank,NrRank(LoopL),Fdr(LoopL)),100,1000,1)                        
%     end
       
end
cd(fullfile('/home/mbellis/array1/sosma/net',sprintf('m%u',ChipRank),'debi'))
eval(sprintf('save %s Debi',OutMat));

%% CALCULATE DISTANCES BETWEEN CLUSTERS IF EXIST SEVERAL NETWORKS
%distances between clusters of a control network and all other test networks
if LoopNb>1
    ControlNet=1;
    ControlCluNb=size(Debi{ControlNet}.clu2comp,1);
    for LoopL=1:length(NetRank)
        if LoopL~=NetRank
            %Distances{LoopL}=zeros(size(Debi{LoopL}.clu2comp,1)*ControlCluNb/2,1);
            Distances=zeros(50);
            DistPos=0;
            %for CluL1=1:ControlCluNb
            for CluL1=1:50
                %for CluL2=1:size(Debi{LoopL}.clu2comp,1)
                for CluL2=1:50
                    DistPos=DistPos+1;
                    %Distances{LoopL}(DistPos)=1-length(intersect(Debi{ControlNet}.clu2ps(CluL1,:),Debi{LoopL}.clu2ps(CluL2,:)))/length(union(Debi{ControlNet}.clu2ps(CluL1,:),Debi{LoopL}.clu2ps(CluL2,:)));
                    Distances(CluL1,CluL2)=1-length(intersect(find(Debi{ControlNet}.clu2ps(CluL1,:)),find(Debi{LoopL}.clu2ps(CluL2,:))))/length(union(find(Debi{ControlNet}.clu2ps(CluL1,:)),find(Debi{LoopL}.clu2ps(CluL2,:))));
                    Distances(CluL2,CluL1)=Distances(CluL1,CluL2);
                end
            end
        end
    end
end










