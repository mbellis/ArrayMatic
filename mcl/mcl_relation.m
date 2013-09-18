%=======================%
% FUNCTION MCL_RELATION %
%=======================%


% MCL_RELATION find the relationship between clusters defined at different limit

%INPUT PARAMETERS

% 1   ChipRank : chip rank
% 2    Postfix : postfix used to construct MCL file name to be loaded
% 3    MinSize : minmum size of clusters to be considered for searching parent cluster
% 4     NewMcl : postfix for writing the new clustering
% varargin : used if variable Info does not exists
% 5 CorrLimits : cell (1x4) containing for each type the correlation limit(s) used for MCL
% 6    NetRank : network(s) used
% 7   ListRank : probe set list(s) used


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

%mcl_relation(5,'n123A',1,'n123B')
%mcl_relation(2,'n80A',1,'n80B')
%mcl_relation(3,'n86A',1,'n86B')
%mcl_relation(6,'n63A',1,'n63B')
%mcl_relation(8,'n228F',1,'n228G')
%mcl_relation(27,'n164F',1,'n164G')
%mcl_relation(27,'n88ton108',1,'n88ton108A')
%mcl_relation(27,'n118',1,'n118A')

%mcl_relation(2,'n64n78A',1,'n64n78B')
%mcl_relation(3,'n70n85A',1,'n70n85B')
%mcl_relation(6,'n48n62A',1,'n48n62B')
%mcl_relation(8,'n212n226A',1,'n212n226B')
%mcl_relation(27,'n149n163A',1,'n149n163B')

%mcl_relation(27,'n164F',50,'n164H')
%mcl_relation(27,'n149n163A',50,'n149n163C')
% mcl_relation(27,'n164F',100,'n164I')
% mcl_relation(27,'n149n163A',100,'n149n163D')
% mcl_relation(27,'n164F',200,'n164J')
% mcl_relation(27,'n149n163A',200,'n149n163E')
% mcl_relation(27,'n164F',300,'n164K')
% mcl_relation(27,'n149n163A',300,'n149n163F')
% mcl_relation(27,'n164F',500,'n164L')
% mcl_relation(27,'n149n163A',500,'n149n163G')






function mcl_relation(ChipRank,Postfix,MinSize,NewMcl,varargin)
global K

%ChipPos=find(K.chip.rank==ChipRank);
%PsNb=K.chip.probesetNb(ChipPos);

NetDir=fullfile(K.dir.net,sprintf('m%u',ChipRank),'mcl');
cd(NetDir)
load(sprintf('m%u_mcl_%s.mat',ChipRank,Postfix)) 

if nargin==4
    CorrLimits=Info.limits;
    NetRank=Info.netRanks;
    ListRank=Info.listRank;
else
    if nargin==7
        CorrLimits=varargin{1};
        NetRank=varargin{2};
        ListRank=varargin{3};
    else
        h=errordlg('mcl_relation needs 4 or 7 parameters');
        waitfor(h)
        error('process canceled')
    end
end


MclClu=Clu;
if MinSize>1
    for TypeL=1:length(MclClu)
        CorrLimit=CorrLimits{TypeL};
        if ~isempty(CorrLimit)
            %process each probe set list
            for ListL=1:length(MclClu{TypeL})
                %process each network
                for NetL=1:length(MclClu{TypeL}{ListL})
                    for ColL=1:size(MclClu{TypeL}{ListL}{NetL},2)
                        CluVal=unique(MclClu{TypeL}{ListL}{NetL}(:,ColL));
                        CluNb=histc(MclClu{TypeL}{ListL}{NetL}(:,ColL),CluVal);
                        ClearPos=find(CluNb<MinSize);
                        if ~isempty(ClearPos)
                            for CluL=1:length(ClearPos)
                                PsPos=find(MclClu{TypeL}{ListL}{NetL}(:,ColL)==CluVal(ClearPos(CluL)));
                                MclClu{TypeL}{ListL}{NetL}(PsPos,ColL)=0;
                            end
                        end
                    end
                end
            end
        end
    end
end



%% PARENT CLUSTERS
%output txt file
fid=fopen(sprintf('m%u_mcl_%s_min%u_parent.txt',ChipRank,Postfix,MinSize),'w');
fidstat=fopen(sprintf('m%u_mcl_%s_min%u_stat.txt',ChipRank,Postfix,MinSize),'w');
fprintf(fidstat,'TypeL\tListL\tnetwork\tLimitL\tparent clunb\tclunb\tps nb\n');
%newcluster
OutClu={};
for TypeL=1:length(Clu)
    CorrLimit=CorrLimits{TypeL};
    if ~isempty(CorrLimit)
        %process each probe set list
        for ListL=1:length(Clu{TypeL})
            PsList=load_pslist(ChipRank,ListRank(ListL));
            PsNb=length(PsList);
            %process each network
            for NetL=1:length(Clu{TypeL}{ListL})
                for LoopL=1:2
                    if LoopL==1
                        Cluster=Clu{TypeL}{ListL}{NetL};
                    else
                        Cluster=NewCluster;
                        OutClu{TypeL}{ListL}{NetL}=NewCluster;
                    end
                    %eliminate column without clustering (corr limit not used)
                    ClearPos=find(sum(Cluster)==0);
                    if ~isempty(ClearPos)
                        Cluster(:,ClearPos)=[];
                    end
                    %number of corr limit effectively used
                    LimitNb=size(Cluster,2);
                    %process the current MCL results
                    if LimitNb>1
                        %search for the the parent cluster that is the cluster in the previous
                        %limit run that contains most of the probe sets of
                        %the current cluster in the current limit run
                        Link=cell(1,LimitNb-1);
                        %percentage of probe sets that comes from the parent cluster
                        LinkPrc=cell(1,LimitNb-1);
                        %size of the current cluster
                        LinkSize=cell(1,LimitNb-1);
                        for LimitL=LimitNb-1:-1:1
                            %process all the clusters of the current limit run
                            Link{LimitL}=zeros(max(Cluster(:,LimitL+1)),1);
                            LinkPrc{LimitL}=zeros(max(Cluster(:,LimitL+1)),1);
                            for CluL=1:length(Link{LimitL})
                                %process the current cluster
                                Pos=find(Cluster(:,LimitL+1)==CluL);
                                %find all the putative parent clusters
                                CluVal=setdiff(unique(Cluster(Pos,LimitL)),0);
                                if ~isempty(CluVal)
                                    %recover the number of probe sets coming from these clusters
                                    CluNb=histc(Cluster(Pos,LimitL),CluVal);
                                    %select the one with the maximal number of probe sets
                                    [MaxCluNb,MaxPos]=max(CluNb);
                                    Link{LimitL}(CluL,1)=CluVal(MaxPos);
                                    LinkPrc{LimitL}(CluL,1)=round(MaxCluNb*100/length(Pos));
                                end
                                LinkSize{LimitL}(CluL,1)=length(Pos);
                            end
                        end

                        %put the results in tabular form
                        %fill the first two columns with the MCL result run with the lower CORR
                        %limit
                        Res=[];
                        Res(:,1)=setdiff(unique(Cluster(:,1)),0);
                        Res(:,2)=histc(Cluster(:,1),Res(:,1));
                        %process all other MCL results with higher CORR limit
                        for LimitL=1:LimitNb-1
                            Res(:,LimitL*3+1)=zeros(size(Res,1),1);
                            Res(:,LimitL*3+2)=zeros(size(Res,1),1);
                            Res(:,LimitL*3+3)=zeros(size(Res,1),1);
                            %process only clusters with a size greater than MinSize
                            %Pos contains the cluster rank of the currently processed MCL result
                            Pos=find(LinkSize{LimitL}>=MinSize);
                            %sort according to the parent clusters
                            [temp,SortIndex]=sort(Link{LimitL}(Pos));
                            Pos=Pos(SortIndex);
                            CurrLink=Link{LimitL}(Pos);
                            CluVal=setdiff(unique(CurrLink),0);
                            if ~isempty(CluVal)
                                CluNb=histc(CurrLink,[1:max(CluVal)]);
                                CluNb=CluNb(CluVal);
                                %process each parent cluster
                                for CluL1=1:length(CluVal)
                                    CurrPos=find(CurrLink==CluVal(CluL1));
                                    CurrClu=Pos(CurrPos);
                                    CurrLinkSize=LinkSize{LimitL}(Pos(CurrPos));
                                    CurrLinkPrc=LinkPrc{LimitL}(Pos(CurrPos));
                                    %find the position in Res of the parent cluster and eventually add
                                    %empty lines if CluNb is greater
                                    PrevCluPos=find(Res(:,(LimitL-1)*3+1)==CluVal(CluL1));
                                    %the parent cluster might miss in the previous result if its size is
                                    %smaller than MinSize
                                    if ~isempty(PrevCluPos)
                                        if length(PrevCluPos)<CluNb(CluL1)
                                            Res=[Res(1:PrevCluPos-1+length(PrevCluPos),:);zeros(CluNb(CluL1)-length(PrevCluPos),size(Res,2));Res(PrevCluPos(end)+1:end,:)];
                                        end
                                        for CluL2=1:CluNb(CluL1)
                                            Res(PrevCluPos(1)-1+CluL2,LimitL*3+1)=CurrClu(CluL2);
                                            Res(PrevCluPos(1)-1+CluL2,LimitL*3+2)=CurrLinkSize(CluL2);
                                            Res(PrevCluPos(1)-1+CluL2,LimitL*3+3)=CurrLinkPrc(CluL2);
                                        end
                                    end
                                end
                            end
                        end
                        if LoopL==2
                            %recover relation between the cluster and original cluster in the MCL
                            %run with the lowest limit
                            Reg=zeros(size(Res,1),max(Clu{TypeL}{ListL}{NetL}(:,1)));
                            for ResL=1:size(Res,1)
                                Pos=find(Res(ResL,:));
                                Pos=Pos(1);
                                Limit=round(Pos/3)+1;
                                ParentClu=Clu{TypeL}{ListL}{NetL}(find(Cluster(:,Limit)==Res(ResL,Pos)),1);
                                for CluL=1:size(Reg,2)
                                    Reg(ResL,CluL)=round(length(find(ParentClu==CluL))*100/length(ParentClu));
                                end
                            end
                            %add some zeros to have same length with Res
                            Reg=[Reg;zeros(4,size(Reg,2))];
                        end
                        %eliminate third empty column
                        Res(:,3)=[];
                        %line of separator (O)
                        Res=[Res;zeros(1,size(Res,2))];
                        %process first corr limit
                        CluVal=setdiff(unique(Cluster(:,1)),0);
                        CluSize=histc(Cluster(:,1),[1:max(CluVal)]);
                        CluSize=CluSize(CluVal);
                        Pos=find(CluSize>=MinSize);
                        %number of probe sets in clusters greater than MinSize
                        Res(end+1,2)=sum(CluSize(Pos));
                        %number of probe sets in clusters smaller than MinSize
                        Res(end+1,2)=sum(CluSize)-Res(end,2);
                        %number of probe sets not clusterized
                        Res(end+1,2)=PsNb-Res(end,2)-Res(end-1,2);
                        %process other corr limits
                        for LimitL=1:LimitNb-1
                            %list of cluster in the parent
                            %do not use unique since some cluster may be absent due to Ps2Net
                            %utilisation
                            CluVal=[1:max(Cluster(:,LimitL+1))];
                            CluSize=histc(Cluster(:,LimitL+1),[1:max(CluVal)]);
                            CluSize=CluSize(CluVal);
                            Pos=find(LinkSize{LimitL}>=MinSize&Link{LimitL}>0);
                            Res(end-2,4+(LimitL-1)*3)=sum(CluSize(Pos));
                            Res(end-1,4+(LimitL-1)*3)=sum(CluSize)-Res(end-2,4+(LimitL-1)*3);
                            Res(end,4+(LimitL-1)*3)=PsNb-Res(end-1,4+(LimitL-1)*3)-Res(end-2,4+(LimitL-1)*3);
                        end


                        %output results
                        if LoopL==2
                            fprintf(fid,'NEWCLUSTER\n');
                            Res=[Reg,Res];
                        end
                        fprintf(fid,'Type %u\n',TypeL);
                        fprintf(fid,'List %u\n',ListRank(ListL));
                        fprintf(fid,'Network %u\n',NetRank(NetL));
                        if LoopL==2
                            for RegL=1:size(Reg,2)
                                fprintf(fid,'R%u\t',RegL)
                            end
                        end
                        fprintf(fid,'C%u\t\t',CorrLimit(1));
                        for LimitL=2:length(CorrLimit)
                            fprintf(fid,'C%u\t\t\t',CorrLimit(LimitL));
                        end
                        fprintf(fid,'\n');
                        Format=[repmat('%u\t',1,size(Res,2)-1),'%u\n'];
                        for ResL=1:size(Res,1)
                            fprintf(fid,Format,Res(ResL,:));
                        end
                        fprintf(fid,'\n');
                        fprintf(fid,'\n');

                        if LoopL==1
                            NewCluster=zeros(size(Cluster));
                            for LimitL=2:size(Cluster,2)                                
                                CluVal=setdiff(unique(Res(:,3+(LimitL-2)*3)),0);
                                if ~isempty(CluVal)
                                    for CluL=1:length(CluVal)
                                        PsPos=find(Cluster(:,LimitL)==CluVal(CluL));
                                        NewCluster(PsPos,LimitL)=CluVal(CluL);
                                    end
                                end
                            end
                            'stop'
%                             TmpCluster=Cluster;
%                             NewCluster=zeros(size(TmpCluster));
%                             %start with lowest corr limit
%                             for LimitL=2:size(TmpCluster,2)                                
%                                 CluVal=setdiff(unique(TmpCluster(:,LimitL)),0);
%                                 %greatest cluster in current corr limit used to add new
%                                 %clusters
%                                 NewCluRank=max(TmpCluster(:,LimitL));                                
%                                 %process each cluster of the current corr limit
%                                 for CluL1=1:length(CluVal)
%                                     Pos=find(TmpCluster(:,LimitL)==CluVal(CluL1));
%                                     %find the parent clusters in the previous corr limit
%                                     PrevCluVal=unique(TmpCluster(Pos,LimitL-1));
%                                     if length(PrevCluVal)>1
%                                         for CluL2=2:length(PrevCluVal)
%                                             %assign a new rank for each of the cluster
%                                             %splitted according to their parent clusters
%                                             %keep the same rank for the first cluster
%                                             NewCluRank=NewCluRank+1;
%                                             PrevPos=find(TmpCluster(Pos,LimitL-1)==PrevCluVal(CluL2));
%                                             TmpCluster(Pos(PrevPos),LimitL)=NewCluRank;
%                                         end
%                                     end
%                                 end
%                                 %reassign new ranks to clusters
%                             
%                             
%                             
%                                 %sort all clusters belonging to the same parent cluster according to their size
%                                 PrevCluVal=setdiff(unique(TmpCluster(:,LimitL-1)),0);
%                                 NewCluRank=0;
%                                 for CluL1=1:length(PrevCluVal)
%                                     PrevCluPos=find(TmpCluster(:,LimitL-1)==PrevCluVal(CluL1));
%                                     CurrCluVal=setdiff(unique(TmpCluster(PrevCluPos,LimitL)),0);
%                                     CurrCluSize=histc(TmpCluster(PrevCluPos,LimitL),[1:max(CurrCluVal)]);
%                                     CurrCluSize=CurrCluSize(CurrCluVal);
%                                     [tmp,SortIndex]=sort(CurrCluSize,'descend');
%                                     CurrCluVal=CurrCluVal(SortIndex);
%                                     for CluL2=1:length(CurrCluVal)
%                                         NewCluRank=NewCluRank+1;
%                                         CluPos=find(TmpCluster(PrevCluPos,LimitL)==CurrCluVal(CluL2));
%                                         NewCluster(PrevCluPos(CluPos),LimitL)=NewCluRank;
%                                     end
%                                 end
% 
%                                 %relationship between old, new rank and size of clusters
%                                 NewCluVal=[0:max(NewCluster(:,LimitL))]';
%                                 NewCluSize=histc(NewCluster(:,LimitL),NewCluVal);
%                                 TmpCluVal=zeros(size(NewCluVal));
%                                 PrevCluVal=zeros(size(NewCluVal));
%                                 FirstCluVal=zeros(size(NewCluVal));
%                                 for CluL=1:length(NewCluVal)
%                                     CluPos=find(NewCluster(:,LimitL)==NewCluVal(CluL));
%                                     if length(unique(TmpCluster(CluPos(1),LimitL)))==1
%                                     TmpCluVal(CluL,1)=TmpCluster(CluPos(1),LimitL);
%                                     else
%                                         'stop'
%                                     end
%                                     if length(unique(TmpCluster(CluPos(1),LimitL-1)))==1
%                                         PrevCluVal(CluL,1)=TmpCluster(CluPos(1),LimitL-1);
%                                         FirstCluVal(CluL,1)=TmpCluster(CluPos(1),1);
%                                     else
%                                         'stop'
%                                     end
%                                 end
%                                 %[FirstCluVal PrevCluVal TmpCluVal,NewCluVal,NewCluSize]
%                                 fprintf(fidstat,'%u\t%u\t%u\t%u\t%u\t%u\t%u\n',TypeL,ListL,NetRank(NetL),LimitL,PrevCluVal(end),NewCluVal(end),sum(NewCluSize)-NewCluSize(1));
% 
%                                 %update TmpCluster with new cluster ranks
%                                 TmpCluster(:,LimitL)=NewCluster(:,LimitL);
%                             end                        
                        end %if LoopL==1
                    end %if LimitNb>1
                end % for LoopL
            end % for NetL
        end % for ListL
    end % if ~isempty(CorrLimit)
end % TypeL
fclose(fid)
fclose(fidstat)

cd(NetDir)
Clu=OutClu;
if nargin==7
    Info.limits=CorrLimits;
    Info.netRanks=NetRank;
    Info.listRank=ListRank;
end
eval(sprintf('save m%u_mcl_%s.mat MclClu Clu Info',ChipRank,NewMcl))
