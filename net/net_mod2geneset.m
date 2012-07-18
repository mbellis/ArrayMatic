%write a list of modules in the format used to conduct gsea analysis
%parameters
%ModelRank = rank of chipset model
%NetRank = rank of network where modules are constructed
%ModSize = minimal size of a module to be considered in gsea analysis
%module2geneset(11,42,5,0)

function net_mod2geneset(ModelRank,NetRank,ModSize,CorrLimit)
global K
NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',NetRank));
cd(NetDir)
load(sprintf('m%un%u_cliques_%02u',ModelRank,NetRank,CorrLimit))
load(sprintf('m%un%u_quasi-cliques_clusters_size%u_%02u',ModelRank,NetRank,ModSize,CorrLimit))
Gs.name={};
Gs.index={};
Gs.region=[];
Gs.cliquePos=[];
Gs.cliqueRank=[];


for RegL=1:length(Cliques)
    ModPos=0;
    for ModL=1:length(Cliques{RegL})
        ModPos=ModPos+1;
        ModRank=Cliques{RegL}(ModL);
        Gs.name{end+1,1}=sprintf('r%u-mp%u_mr%u',RegL,ModPos,ModRank);
        Gs.index{end+1,1}=find(Clu==ModRank);            
        Gs.region(end+1,1)=RegL;
        Gs.cliquePos(end+1,1)=ModPos;
        Gs.cliqueRank(end+1,1)=ModRank;
    end
end

Ps=[];
for i=1:length(Gs.name)
    Ps=[Ps;Gs.index{i}];
end
if length(Ps)==length(unique(Ps))
    cd(K.dir.cliques)
eval(sprintf('save m%u_cliques Gs',ModelRank))
else
    h=warndlg(sprintf('exist %u doublons in modules',length(Ps)-length(unique(Ps))));
    waitfor(h)
end
    
