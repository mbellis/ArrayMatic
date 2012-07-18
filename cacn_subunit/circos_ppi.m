%draw CIRCOS image for protein-protein interactions found with PsicQuick
cd('/home/mbellis/sosma/data/cacn')
load psicquic

SubUnitNb=length(FileList);
SubUnits=FileList;
for FileL=1:length(FileList)
    SubUnits{FileL}=FileList{FileL}(1:length(FileList{FileL})-4);
end
%correct first 10 subunits
SubUnits{1}='Cav1.1';
SubUnits{2}='Cav1.2';
SubUnits{3}='Cav1.3';
SubUnits{4}='Cav1.4';
SubUnits{5}='Cav2.1';
SubUnits{6}='Cav2.2';
SubUnits{7}='Cav2.3';
SubUnits{8}='Cav3.1';
SubUnits{9}='Cav3.2';
SubUnits{10}='Cav3.3';
%find differrent CACN subunits (alpha 1,2, beta and gamma)
CacnTypes={'Cav','CACNA2','CACNB','CACNG'};
CacnType=zeros(SubUnitNb,1);
for CacnL=1:length(CacnTypes)
    CacnType(strmatch(CacnTypes{CacnL},SubUnits))=CacnL;
end
SubUnitShortName={'','alpha2','beta','gamma'};


%make matrix
PsiMat=zeros(length(AllAlias),SubUnitNb);
for FileL=1:length(FileList)
    CurrAliases=unique(MainAlias{FileL});
    for AliasL=1:length(CurrAliases)
        PsiMat(strmatch(CurrAliases{AliasL},AllAlias,'exact'),FileL)=1;
    end
end
%codage according to repetition
for i=1:length(PsiMat)
    PsiMat(i,find(PsiMat(i,:)))=sum(PsiMat(i,:));
end

%reorder CACNA1 according to accepted nomenclature
[SortOrder]=[10,3,4,6,1,2,5,7,8,9];
PsiMat(:,1:10)=PsiMat(:,SortOrder);

%max nb of subunits having the same interactor
%used to make color map
PpiNb=max(max(PsiMat));
LinkColorMap=colors(colormap,PpiNb);


MemAllAlias=AllAlias;
MemPsiMat=PsiMat;
%reorder PPi according to the number of subunits intrracting with them
MaxSubUnit=max(PsiMat,[],2);
[MaxSubUnit SortOrder]=sort(MaxSubUnit);
PsiMat=PsiMat(SortOrder,:);
AllAlias=AllAlias(SortOrder);
%in each class of interaction, reorder according to the subunits
for NbL=1:PpiNb
    Pos=find(MaxSubUnit==NbL);
    SubPos=zeros(length(Pos),1);
    for PosL=1:size(Pos)
        BinString=repmat('0',1,SubUnitNb);
        BinString(find(PsiMat(Pos(PosL),:)==NbL))='1';
        SubPos(PosL)=bin2dec(BinString);
    end
    [temp,SortOrder]=sort(SubPos);
    PsiMat(Pos,:)=PsiMat(Pos(SortOrder),:);
    AllAlias(Pos)=AllAlias(Pos(SortOrder));
end
AliasNb=length(AllAlias);


%% MCL regions
% %LOAD MCL clusters on m27n23
% cd('//home/mbellis/array1/sosma/net/m027/n00023/mcl')
% load m27n23_mcl
% %LOAD PsMatrix
% ModelRank=27;NetRanks=[2:22];
% NetNb=21;
% cd('/home/mbellis/sosma/raydataraw/2011_CACN')
% eval(sprintf('load cacn_m%u_net%u_netnb%u',ModelRank,NetRanks(1),NetNb))
% cd('/home/mbellis/sosma/data/psawn/mldata/mouse/m27_2to22')
% eval(sprintf('load m%u_n%u_netnb%u_probenb1_newps_stat_netprc100_pvcorr1',ModelRank,NetRanks(1),NetNb));
% clear NewPs PsBy Stat LinkedPs
% %write ENSEMBL ID FOR six first MCL clusters
% cd('/home/mbellis/sosma/data/cacn/')
% for CluL=1:6
%     fid=fopen(sprintf('m27n23_mlc_clu%u_ensids.txt',CluL),'w');
%     CurrPs=find(Clu(:,1)==CluL);
%     CurrGenePos=PsMatrix(CurrPs,1);
%     CurrGenePos(find(CurrGenePos==0))=[];
%     CurrGenes=Genes.name(CurrGenePos);
%     KeepPos=strmatch('ENS',CurrGenes);
%     CurrGenes=CurrGenes(KeepPos);
%     for GeneL=1:length(CurrGenes)
%         fprintf(fid,'%s\n',CurrGenes{GeneL});
%     end
%     fclose(fid)
% end

%LOAD MCL clusters on m27n24
cd('//home/mbellis/array1/sosma/net/m027/n00024/mcl')
load m27n24_mcl
%LOAD PsMatrix
ModelRank=27;NetRanks=[2:22];
NetNb=21;
cd('/home/mbellis/sosma/raydataraw/2011_CACN')
eval(sprintf('load cacn_m%u_net%u_netnb%u',ModelRank,NetRanks(1),NetNb))
cd('/home/mbellis/sosma/data/psawn/mldata/mouse/m27_2to22')
eval(sprintf('load m%u_n%u_netnb%u_probenb1_newps_stat_netprc100_pvcorr1',ModelRank,NetRanks(1),NetNb));
clear NewPs PsBy Stat LinkedPs
%write ENSEMBL ID FOR six first MCL clusters
cd('/home/mbellis/sosma/data/cacn/')
for CluL=1:6
    fid=fopen(sprintf('m27n24_mlc_clu%u_ensids.txt',CluL),'w');
    CurrPs=find(Clu(:,1)==CluL);
    CurrGenePos=PsMatrix(CurrPs,1);
    CurrGenePos(find(CurrGenePos==0))=[];
    CurrGenes=Genes.name(CurrGenePos);
    KeepPos=strmatch('ENS',CurrGenes);
    CurrGenes=CurrGenes(KeepPos);
    for GeneL=1:length(CurrGenes)
        fprintf(fid,'%s\n',CurrGenes{GeneL});
    end
    fclose(fid)
end

%load correspondance between m27 Ensembl Ids and Genes Ids
cd('/home/mbellis/sosma/data/cacn')
[EnsIds,GeneIds]=textread('EnsId_GeneId_m27.txt','%s%s','delimiter','\t');
%find correspondance between MCL region and genes ids 
Regions=zeros(AliasNb,6);
for RegL=1:6
    CurrGenePos=PsMatrix(find(Clu(:,1)==RegL),1);
    CurrGenePos(find(CurrGenePos==0))=[];
    CurrIds=unique(Genes.name(CurrGenePos));
    [temp AliasPos,temp]=intersect(EnsIds,CurrIds);
    [temp AliasPos temp]=intersect(AllAlias,upper(GeneIds(AliasPos)));
    Regions(AliasPos,RegL)=1;
end

RegionOrder=[5,1,2,3,6,4];


%reorder PsiMat and AllAlias according to regions
MemPsiMat=PsiMat;
MemAllAlias=AllAlias;
MaxSubUnit=max(PsiMat,[],2);
CurrRegions=Regions;
%PROCESS SINGLE INTERACTIONS
RegLimits1={};
for SubL=1:SubUnitNb
    Pos=find(PsiMat(:,SubL)==1);
    SortOrder=[];
    RegLimits1{SubL}=zeros(7,2);
    Index=min(Pos);
    for RegL=1:6
        CurrPos=find(CurrRegions(Pos,RegionOrder(RegL))==1);
        if ~isempty(CurrPos)
            SortOrder=[SortOrder;CurrPos];
            RegLimits1{SubL}(RegL,:)=[Index,Index+length(CurrPos)-1];
            Index=Index+length(CurrPos);
            %don't count other regions of the selected genes
            CurrRegions(Pos(CurrPos),:)=zeros(length(CurrPos),6);
        else
            RegLimits1{SubL}(RegL,:)=[0,0];
        end
    end
    CurrPos=find(sum(Regions(Pos,:),2)==0);
    if ~isempty(CurrPos)
        SortOrder=[SortOrder;CurrPos];
        RegLimits1{SubL}(7,:)=[Index,Index+length(CurrPos)-1];
    else
        RegLimits1{SubL}(7,:)=[0,0];
    end
    if~isempty(Pos)
        PsiMat(Pos,:)=PsiMat(Pos(SortOrder),:);
        AllAlias(Pos)=AllAlias(Pos(SortOrder));
    end
end

%PROCESS MULTIPLE INTERACTIONS
RegLimits2={};
for NbL=2:PpiNb
    Pos=find(MaxSubUnit==NbL);
    RegLimits2{NbL-1}=zeros(7,2);
    SortOrder=[];
    Index=min(Pos);
    for RegL=1:6
        CurrPos=find(CurrRegions(Pos,RegionOrder(RegL))==1);
        if ~isempty(CurrPos)
            SortOrder=[SortOrder;CurrPos];
            RegLimits2{NbL-1}(RegL,:)=[Index,Index+length(CurrPos)-1];
            Index=Index+length(CurrPos);
            %don't count other regions of the selected genes
            CurrRegions(Pos(CurrPos),:)=zeros(length(CurrPos),6);
        else
            RegLimits2{NbL-1}(RegL,:)=[0,0];
        end
    end
    CurrPos=find(sum(Regions(Pos,:),2)==0);
    if ~isempty(CurrPos)
        SortOrder=[SortOrder;CurrPos];
        RegLimits2{NbL-1}(7,:)=[Index,Index+length(CurrPos)-1];
    else
        RegLimits2{NbL-1}(7,:)=[0,0];
    end
    if~isempty(Pos)
        PsiMat(Pos,:)=PsiMat(Pos(SortOrder),:);
        AllAlias(Pos)=AllAlias(Pos(SortOrder));
    end
end



%write ENSEMBL ID for whole chip
%     fid=fopen(sprintf('m27_ensids.txt',CluL),'w');
%     CurrGenePos=find(PsMatrix(:,1));
%     CurrGenes=Genes.name(CurrGenePos);
%     KeepPos=strmatch('ENS',CurrGenes);
%     CurrGenes=CurrGenes(KeepPos);
%     for GeneL=1:length(CurrGenes)
%         fprintf(fid,'%s\n',CurrGenes{GeneL});
%     end
%     fclose(fid)


%LOAD GO TERMS FOUND BY GoTermFinder on PPI
GoBranch={'bp','mf','cc'};
cd('/home/mbellis/sosma/data/cacn/')
for TypeL=1:3
    sprintf('gotermfinder_%s_mm_only.txt',GoBranch{TypeL})
    [Go{TypeL},Def{TypeL},Pv1{TypeL},Pv2{TypeL},Found{TypeL},ListSize{TypeL},TotalFound{TypeL},TotalSize{TypeL},Fdr{TypeL},FalsePos{TypeL},GeneList{TypeL}]=textread(sprintf('gotermfinder_%s_mm_only.txt',GoBranch{TypeL}),'%s%s%s%s%u%u%u%u%s%s%s','delimiter','\t','headerlines',12,'bufsize',30000);
end
SearchedGO{1}={'nervous system development';...
'neurogenesis';...
'generation of neurons';...
'neuron differentiation';...
'synaptic transmission';...
'transmission of nerve impulse';...
'neuron development';...
'brain development';...
'regulation of neurogenesis';...
'response to light stimulus';...
'dendrite development';...
'neurotransmitter secretion'};
SearchedGO{3}={'synapse';...
    'dendrite';...
    'neuron projection';...
    'neuronal cell body';...
    'axon';...
    'synaptic membrane';...
    'synaptosome'};
%immmunity
SearchedGO{3}{1}={'GO:0002376';'GO:0006955';'GO:0050896';'GO:0006952';'GO:0045321'};
SearchedGO{3}{2}={'GO:0004871';'GO:0060089';'GO:0005515';'GO:0019955';'GO:0003823'};
SearchedGO{3}{3}={'GO:0000323';'GO:0005764';'GO:0042611';'GO:0005773';'GO:0001772'};
%nucleus
SearchedGO{2}{1}={'GO:0006139';'GO:0043170';'GO:0044238';'GO:0044237';'GO:0006396'};
SearchedGO{2}{2}={'GO:0003723';'GO:0003676';'GO:0000166';'GO:0003735';'GO:0005515'};
SearchedGO{2}{3}={'GO:0005634';'GO:0044424';'GO:0005622';'GO:0043226';'GO:0043229'};
%development
SearchedGO{5}{1}={'GO:0048856';'GO:0007155';'GO:0022610';'GO:0009653';'GO:0032502'};
SearchedGO{5}{2}={'GO:0005201';'GO:0005515';'GO:0008092';'GO:0003779';'GO:0005509'};
SearchedGO{5}{3}={'GO:0005578';'GO:0031012';'GO:0005576';'GO:0044421';'GO:0044420'};
%metabolism I
SearchedGO{4}{1}={'GO:0006091';'GO:0006936';'GO:0006119';'GO:0015980';'GO:0042773'};
SearchedGO{4}{2}={'GO:0003954';'GO:0008137';'GO:0050136';'GO:0016655';'GO:0016651'};
SearchedGO{4}{3}={'GO:0005739';'GO:0044429';'GO:0044444';'GO:0005740';'GO:0005743'};
%nervous system
SearchedGO{1}{1}={'GO:0019226';'GO:0007399';'GO:0007268';'GO:0007267';'GO:0050877'};
SearchedGO{1}{2}={'GO:0015075';'GO:0022891';'GO:0005216';'GO:0022838';'GO:0005215'};
SearchedGO{1}{3}={'GO:0043005';'GO:0045202';'GO:0030424';'GO:0042995';'GO:0044456'};
%metabolism II
SearchedGO{6}{1}={'GO:0019752';'GO:0006082';'GO:0006807';'GO:0006629';'GO:0009308'};
SearchedGO{6}{2}={'GO:0016491';'GO:0003824';'GO:0004497';'GO:0016705';'GO:0020037'};
SearchedGO{6}{3}={'GO:0005783';'GO:0005792';'GO:0042598';'GO:0044444';'GO:0005737'};





for RegL=1:6
    GoGeneList{RegL}={};
    for BranchL=[1,3]
        for i=1:length(SearchedGO{RegL}{BranchL})
            Pos=strmatch(SearchedGO{RegL}{BranchL}{i},Go{BranchL},'exact');
            if ~isempty(Pos)
                CurrGeneList={};
                CurrGenes=GeneList{BranchL}{Pos};
                CommaPos=findstr(',',CurrGenes);
                if ~isempty(CommaPos)
                    CurrGeneList{1}=CurrGenes(1:CommaPos(1)-1);
                    for j=1:length(CommaPos)-1
                        CurrGeneList{end+1}=CurrGenes(CommaPos(j)+2:CommaPos(j+1)-1);
                    end
                    CurrGeneList{end+1}=CurrGenes(CommaPos(end)+2:end);
                else
                    CurrGeneList{1}=CurrGenes;
                end
                GoGeneList{RegL}=union(GoGeneList{RegL},CurrGeneList);
            end
        end
    end
end

%% write histogram files
cd('/home/mbellis/sosma/data/cacn/circos/ppi')
for RegL1=1:6
    fid=fopen(sprintf('cacn_ppi_histogram_%u.txt',RegL1),'w');
    for SubL=1:SubUnitNb
        for RegL=1:7
            if RegLimits1{SubL}(RegL,2)-RegLimits1{SubL}(RegL,1)>=1
                CurrAlias=AllAlias(RegLimits1{SubL}(RegL,1):RegLimits1{SubL}(RegL,2));
                CommonAlias=intersect(CurrAlias,upper(GoGeneList{RegL1}));
                fprintf(fid,'PPI %u %u %u\n',RegLimits1{SubL}(RegL,1),RegLimits1{SubL}(RegL,2),round(length(CommonAlias)*100/length(CurrAlias)));
            end
        end
    end
    for NbL=2:PpiNb
        for RegL=1:7
            if RegLimits2{NbL-1}(RegL,2)-RegLimits2{NbL-1}(RegL,1)>=1
                CurrAlias=AllAlias(RegLimits2{NbL-1}(RegL,1):RegLimits2{NbL-1}(RegL,2));
                CommonAlias=intersect(CurrAlias,upper(GoGeneList{RegL1}));
                fprintf(fid,'PPI %u %u %u\n',RegLimits2{NbL-1}(RegL,1),RegLimits2{NbL-1}(RegL,2),round(length(CommonAlias)*100/length(CurrAlias)));
            end
        end
    end    
    fclose(fid)
end

%write highligths for histograms
for RegL=[1,3:6]
    fid=fopen(sprintf('cacn_ppi_highlight_histo_%u.txt',RegL),'w');
    fprintf(fid,'PPI 0 %u fill_color=bq2-%u\n',AliasNb,RegL+7);
    fclose(fid)
end
fid=fopen('cacn_ppi_highlight_histo_2.txt','w');
fprintf(fid,'PPI 0 %u fill_color=bq2-14\n',AliasNb);
fclose(fid)


Colors={'red','green','blue','yellow','purple'};
ColorOrder=[1,2,3,4,5];

%% write ordered alias with info
%construct region info
Region=zeros(AliasNb,1);
for SubL=1:SubUnitNb
    for RegL=1:6
        if RegLimits1{SubL}(RegL,2)-RegLimits1{SubL}(RegL,1)>=1
            Region(RegLimits1{SubL}(RegL,1):RegLimits1{SubL}(RegL,2))=RegL;
        end
    end
end
for NbL=1:length(RegLimits2)
    for RegL=1:6
        if RegLimits2{NbL}(RegL,2)-RegLimits2{NbL}(RegL,1)>=1
            Region(RegLimits2{NbL}(RegL,1):RegLimits2{NbL}(RegL,2))=RegL;
        end
    end
end

%construct GoMatrix
GoMat=uint8(zeros(AliasNb,length(Go{1})+length(Go{2})+length(Go{3})));
GoPos=0;
for BranchL=1:3
    for GoL=1:length(Go{BranchL})
        GoPos=GoPos+1;
        CurrGeneList={};
        CurrGenes=GeneList{BranchL}{GoL};
        CommaPos=findstr(',',CurrGenes);
        if ~isempty(CommaPos)
            CurrGeneList{1}=CurrGenes(1:CommaPos(1)-1);
            for j=1:length(CommaPos)-1
                CurrGeneList{end+1}=CurrGenes(CommaPos(j)+2:CommaPos(j+1)-1);
            end
            CurrGeneList{end+1}=CurrGenes(CommaPos(end)+2:end);
        else
            CurrGeneList{1}=CurrGenes;
        end
        [temp,GenePos,temp]=intersect(AllAlias,upper(CurrGeneList));
        GoMat(GenePos,GoPos)=1;
    end
end

cd('/home/mbellis/sosma/data/cacn/circos/ppi')
fid=fopen('gene_info.txt','w');
%write first header line
fprintf(fid,'Position\tGeneID\tEnsemblID\tRegion\tSubNb')
for SubL=1:SubUnitNb
    fprintf(fid,'\t%s',SubUnits{SubL});
end
fprintf(fid,'\tGO TERMS\n')

%write snd header line
for TabL=1:30
    fprintf(fid,'\t');
end
fprintf(fid,'\tBP TERMS');
for TabL=1:length(Go{1})-1
    fprintf(fid,'\t')
end
fprintf(fid,'MF TERMS');
for TabL=1:length(Go{2})-1
    fprintf(fid,'\t')
end
fprintf(fid,'CC TERMS\n');

%write third header line
for TabL=1:30
    fprintf(fid,'\t');
end
for BranchL=1:3
    for GoL=1:length(Go{BranchL})
        fprintf(fid,'\t%s',Def{BranchL}{GoL});
    end
end
fprintf(fid,'\n');

%write GO information
%Pv1
fprintf(fid,'Pv')
for TabL=1:30
    fprintf(fid,'\t');
end
for BranchL=1:3
    for GoL=1:length(Go{BranchL})
        fprintf(fid,'\t%s',num2str(str2num(Pv1{BranchL}{GoL})));
    end
end
fprintf(fid,'\n');
%Found
fprintf(fid,'Found')
for TabL=1:30
    fprintf(fid,'\t');
end
for BranchL=1:3
    for GoL=1:length(Go{BranchL})
        fprintf(fid,'\t%u',Found{BranchL}(GoL));
    end
end
fprintf(fid,'\n');
%ListSize
fprintf(fid,'ListSize')
for TabL=1:30
    fprintf(fid,'\t');
end
for BranchL=1:3
    for GoL=1:length(Go{BranchL})
        fprintf(fid,'\t%u',ListSize{BranchL}(GoL));
    end
end
fprintf(fid,'\n');
%TotalFound
fprintf(fid,'TotalFound')
for TabL=1:30
    fprintf(fid,'\t');
end
for BranchL=1:3
    for GoL=1:length(Go{BranchL})
        fprintf(fid,'\t%u',TotalFound{BranchL}(GoL));
    end
end
fprintf(fid,'\n');
%TotalSize
fprintf(fid,'TotalSize')
for TabL=1:30
    fprintf(fid,'\t');
end
for BranchL=1:3
    for GoL=1:length(Go{BranchL})
        fprintf(fid,'\t%u',TotalSize{BranchL}(GoL));
    end
end
fprintf(fid,'\n');
fprintf(fid,'\n');


UGeneIds=upper(GeneIds);    
for GeneL=1:length(AllAlias)
    %Position
    fprintf(fid,'%u',GeneL);
    %GeneID
    CurrGene=AllAlias{GeneL};
    fprintf(fid,'\t%s',CurrGene);
    GenePos=strmatch(CurrGene,UGeneIds);
    %EnsemblID
    if ~isempty(GenePos)
        fprintf(fid,'\t%s',EnsIds{GenePos(1)});
    else
        fprintf(fid,'\t');
    end
    %Region
    fprintf(fid,'\t%u',Region(GeneL));
    %SubNb
    fprintf(fid,'\t%u',MaxSubUnit(GeneL));
    %Interacting sub units
    for SubL=1:SubUnitNb
        if PsiMat(GeneL,SubL)>0
            fprintf(fid,'\t1');
        else
            fprintf(fid,'\t');
        end
    end
    %GO terms
    for GoL=1:size(GoMat,2)
        if GoMat(GeneL,GoL)==1
            fprintf(fid,'\t1');
        else
            fprintf(fid,'\t');
        end
    end
    fprintf(fid,'\n');
end
    
    
fclose(fid)    


%% write CIRCOS highlight

%write CIRCOS highlight files (Proteins that interact with only one subunit)
cd('/home/mbellis/sosma/data/cacn/circos/ppi')
fid=fopen('cacn_ppi_highlight_1.txt','w');
SumPos=0;
for SubL=1:SubUnitNb
    Pos=find(PsiMat(:,SubL)==1);
    %SumPos=SumPos+length(Pos)
    if ~isempty(Pos)
        %if length(Pos)>10            
            fprintf(fid,'PPI %u %u fill_color=bq1-%u\n',min(Pos),max(Pos),ColorOrder(CacnType(SubL)));
        %end
    end
end
fclose(fid)

%write CIRCOS highlight files (Proteins that interact with more than one subunit)
MaxSubUnit=max(PsiMat,[],2);
cd('/home/mbellis/sosma/data/cacn/circos/ppi')
fid=fopen('cacn_ppi_highlight_2.txt','w');

for NbL=2:PpiNb
    Pos=find(MaxSubUnit==NbL);
    if ~isempty(Pos)
        fprintf(fid,'PPI %u %u fill_color=link-%u\n',min(Pos),max(Pos),NbL);            
    end
end
fclose(fid)

%write CIRCOS highlight files (Regions)
cd('/home/mbellis/sosma/data/cacn/circos/ppi')
fid=fopen('cacn_ppi_highlight_3.txt','w');

for SubL=1:SubUnitNb
    for RegL=1:7
        if RegLimits1{SubL}(RegL,2)-RegLimits1{SubL}(RegL,1)>=1
            fprintf(fid,'PPI %u %u fill_color=reg-%u\n',RegLimits1{SubL}(RegL,1),RegLimits1{SubL}(RegL,2),RegL);            
        end
    end
end
for NbL=1:length(RegLimits2)
    for RegL=1:7
        if RegLimits2{NbL}(RegL,2)-RegLimits2{NbL}(RegL,1)>=1
            fprintf(fid,'PPI %u %u fill_color=reg-%u\n',RegLimits2{NbL}(RegL,1),RegLimits2{NbL}(RegL,2),RegL);            
        end
    end
end
fclose(fid)


%% write CIRCOS 'chromosome' file

%color coding according to Brewer palette defined as bq1 in CIRCOS configuration file

cd('/home/mbellis/sosma/data/cacn/circos/ppi')
fid=fopen('cacn_ppi_chromosome.txt','w');
PsiMatBin=PsiMat;
PsiMatBin(find(PsiMatBin))=1;
PpiNbList=sum(PsiMatBin);

% for CacnL=1:SubUnitNb
%     if PpiNbList(CacnL)>0
%         fprintf(fid,'chr - %s %s-%s 0 %u %s\n',SubUnits{CacnL},SubUnitShortName{CacnType(CacnL)},SubUnits{CacnL}(length(CacnTypes{CacnType(CacnL)})+1:end),PpiNbList(CacnL),Colors{ColorOrder(CacnType(CacnL))});
%     else
%         fprintf(fid,'chr - %s %s-%s 0 10 %s\n',SubUnits{CacnL},SubUnitShortName{CacnType(CacnL)},SubUnits{CacnL}(length(CacnTypes{CacnType(CacnL)})+1:end),Colors{ColorOrder(CacnType(CacnL))});
%     end
% end
% fprintf(fid,'chr - PPI Interactors 0 %u %s\n',length(AllAlias),BrewerColor{end});

% subuntis size proportional to nb of interactants
for CacnL=1:10
        %fprintf(fid,'chr - %s %s-%s 0 %u bq1-%u\n',SubUnits{CacnL},SubUnitShortName{CacnType(CacnL)},SubUnits{CacnL}(length(CacnTypes{CacnType(CacnL)})+1:end),PpiNbList(CacnL),ColorOrder(CacnType(CacnL)));
        fprintf(fid,'chr - %s %s 0 %u bq1-%u\n',SubUnits{CacnL},SubUnits{CacnL},PpiNbList(CacnL)+PpiNb,ColorOrder(CacnType(CacnL)));
end
for CacnL=11:SubUnitNb
        %fprintf(fid,'chr - %s %s-%s 0 %u bq1-%u\n',SubUnits{CacnL},SubUnitShortName{CacnType(CacnL)},SubUnits{CacnL}(length(CacnTypes{CacnType(CacnL)})+1:end),PpiNbList(CacnL),ColorOrder(CacnType(CacnL)));
        fprintf(fid,'chr - %s %s-%s 0 %u bq1-%u\n',SubUnits{CacnL},SubUnitShortName{CacnType(CacnL)},SubUnits{CacnL}(length(CacnTypes{CacnType(CacnL)})+1:end),PpiNbList(CacnL)+PpiNb,ColorOrder(CacnType(CacnL)));
end
fprintf(fid,'chr - PPI Interactors 0 %u bq1-%u\n',length(AllAlias),ColorOrder(end));
fclose(fid)

%same size for all subunits
% for CacnL=1:10
%         %fprintf(fid,'chr - %s %s-%s 0 %u bq1-%u\n',SubUnits{CacnL},SubUnitShortName{CacnType(CacnL)},SubUnits{CacnL}(length(CacnTypes{CacnType(CacnL)})+1:end),PpiNbList(CacnL),ColorOrder(CacnType(CacnL)));
%         fprintf(fid,'chr - %s %s 0 %u bq1-%u\n',SubUnits{CacnL},SubUnits{CacnL},PpiNb+1,ColorOrder(CacnType(CacnL)));
% end
% for CacnL=11:SubUnitNb
%         %fprintf(fid,'chr - %s %s-%s 0 %u bq1-%u\n',SubUnits{CacnL},SubUnitShortName{CacnType(CacnL)},SubUnits{CacnL}(length(CacnTypes{CacnType(CacnL)})+1:end),PpiNbList(CacnL),ColorOrder(CacnType(CacnL)));
%         fprintf(fid,'chr - %s %s-%s 0 %u bq1-%u\n',SubUnits{CacnL},SubUnitShortName{CacnType(CacnL)},SubUnits{CacnL}(length(CacnTypes{CacnType(CacnL)})+1:end),PpiNb+1,ColorOrder(CacnType(CacnL)));
% end
% fprintf(fid,'chr - PPI Interactors 0 %u bq1-%u\n',length(AllAlias),ColorOrder(end));
% fclose(fid)

%% write CIRCOS configuration file
cd('/home/mbellis/sosma/data/cacn/circos/ppi')
fid=fopen('cacn_ppi.conf','w');
%COLORS
fprintf(fid,'<colors>\n');
fprintf(fid,'<<include /usr/local/circos-0.55/etc/colors.conf>>\n');
fprintf(fid,'<<include /usr/local/circos-0.55/etc/brewer.conf>>\n');
fprintf(fid,'bq1-1 = 251,180,174\n');
fprintf(fid,'bq1-2 = 179,205,227\n');
fprintf(fid,'bq1-3 = 204,235,197\n');
fprintf(fid,'bq1-4 = 222,203,228\n');
fprintf(fid,'bq1-5 = 254,217,166\n');


fprintf(fid,'reg-1 = 0,255,0\n');
fprintf(fid,'reg-2 = 255,255,0\n');
fprintf(fid,'reg-3 = 255,0,255\n');
fprintf(fid,'reg-4 = 255,0,0\n');
fprintf(fid,'reg-5 = 0,255,255\n');
fprintf(fid,'reg-6 = 255,150,0\n');
fprintf(fid,'reg-7 = 100,100,100\n');



%brewer qualitative paired
fprintf(fid,'bq2-1 = 51, 160, 44\n');
fprintf(fid,'bq2-2 = 255, 255, 0\n');
fprintf(fid,'bq2-3 = 106, 61, 154\n');
fprintf(fid,'bq2-4 = 227, 26, 28\n');
fprintf(fid,'bq2-5 = 31, 120, 180\n');
fprintf(fid,'bq2-6 = 255, 127, 0\n');
fprintf(fid,'bq2-7 = 120,120,120\n');
fprintf(fid,'bq2-8 = 178, 223, 138\n');
fprintf(fid,'bq2-9 = 255, 255, 153\n');
fprintf(fid,'bq2-10 = 202, 178, 214\n');
fprintf(fid,'bq2-11 = 251, 154, 153\n');
fprintf(fid,'bq2-12 = 166, 206, 227\n');
fprintf(fid,'bq2-13 = 253, 191, 111\n');
fprintf(fid,'bq2-14 = 200,200,200\n');

%linkcolor
for ColorL=1:PpiNb
    fprintf(fid,'link-%u = %u,%u,%u\n',ColorL,round(255*LinkColorMap(ColorL,1)),round(255*LinkColorMap(ColorL,2)),round(255*LinkColorMap(ColorL,3)));
end
fprintf(fid,'</colors>\n');

%FONTS
fprintf(fid,'<fonts>\n');
fprintf(fid,'<<include /usr/local/circos-0.55/etc/fonts.conf>>\n');
fprintf(fid,'</fonts>\n');

%IMAGE
fprintf(fid,'<image>\n');
%fprintf(fid,'file = cacn_ppi.png\n');
fprintf(fid,'<<include /usr/local/circos-0.55/etc/image.conf>>\n');
fprintf(fid,'</image>\n');

%LINKS
 fprintf(fid,'<links>\n'); 
 fprintf(fid,'<link chain>\n'); 
 fprintf(fid,'ribbon = yes\n');
 fprintf(fid,'file = cacn_ppi_bundle.txt\n');
 fprintf(fid,'bezier_radius = 0r\n');
 fprintf(fid,'radius = 0.93r\n');
 fprintf(fid,'thickness = 1p\n');
 fprintf(fid,'</link>\n');
 fprintf(fid,'</links>\n');

%  fprintf(fid,'<links>\n');
%  fprintf(fid,'<link chain>\n'); 
%  fprintf(fid,'file = ppi_links.txt\n');
%  fprintf(fid,'bezier_radius = 0r\n');
%  fprintf(fid,'radius = 0.9r\n');
%  fprintf(fid,'thickness = 1p\n'); 
%  fprintf(fid,'</link>\n'); 
%  fprintf(fid,'</links>\n');

%INCLUDE 
% fprintf(fid,'<<include ideogram.conf>>\n');
% fprintf(fid,'karyotype = cacn_ppi_chromosome.txt\n');
% fprintf(fid,'chromosomes_units = %u\n',round(sum(PpiNbList)/100));
% fprintf(fid,'chromosomes_display_default = yes\n');
% 
Radius=zeros(6,2);
Radius(1,:)=[1,1.05];
Radius(2,:)=[1.05,1.10];
Radius(3,:)=[1.10,1.15];
Radius(4,:)=[1.15,1.20];
Radius(5,:)=[1.20,1.25];
Radius(6,:)=[1.25,1.30];


%HIGHLIGHTS
fprintf(fid,'<highlights>\n');

%histogram highlight
for i=1:6
    fprintf(fid,'<highlight>\n');
    fprintf(fid,'file = cacn_ppi_highlight_histo_%u.txt\n',i);
    fprintf(fid,'r0 = %0.2fr\n',Radius(i,1)+0.02);
    fprintf(fid,'r1 = %0.2fr\n',Radius(i,2)+0.02);
    fprintf(fid,'</highlight>\n');
end

% fprintf(fid,'<highlight>\n');
% fprintf(fid,'file = cacn_ppi_highlight_1.txt\n');
% fprintf(fid,'r0 = 0.88r\n');
% fprintf(fid,'r1 = 0.93r\n');
% fprintf(fid,'</highlight>\n');

fprintf(fid,'<highlight>\n');
fprintf(fid,'file = cacn_ppi_highlight_2.txt\n');
fprintf(fid,'r0 = 0.88r\n');
fprintf(fid,'r1 = 0.93r\n');
fprintf(fid,'</highlight>\n');

fprintf(fid,'<highlight>\n');
fprintf(fid,'file = cacn_ppi_highlight_3.txt\n');
fprintf(fid,'r0 = 0.96r\n');
fprintf(fid,'r1 = 0.98r\n');
fprintf(fid,'</highlight>\n');

fprintf(fid,'</highlights>\n');

%HISTOGRAMS

fprintf(fid,'<plots>\n');
for i=1:6
fprintf(fid,'<plot>\n');
fprintf(fid,'type = histogram\n');
fprintf(fid,'file = cacn_ppi_histogram_%u.txt\n',i);
fprintf(fid,'min = 0\n');
fprintf(fid,'max = 100\n');
fprintf(fid,'r0 = %0.2fr\n',Radius(i,1)+0.02);
fprintf(fid,'r1 = %0.2fr\n',Radius(i,2)+0.02);
fprintf(fid,'color = bq2-%u\n',i);
fprintf(fid,'thickness = 2p\n');
fprintf(fid,'fill_under = yes\n');
fprintf(fid,'fill_color =  bq2-%u\n',i);
fprintf(fid,'</plot>\n');
end
fprintf(fid,'</plots>\n');

%fprintf(fid,'<<include housekeeping.conf>>\n')
fprintf(fid,'units_ok = bupr\n');
fprintf(fid,'units_nounit = n\n');

fprintf(fid,'<<include ideogram.conf>>\n');
fprintf(fid,'karyotype = cacn_ppi_chromosome.txt\n');
fprintf(fid,'chromosomes_units = %u\n',round(sum(PpiNbList)/100));
fprintf(fid,'chromosomes_display_default = yes\n');
fprintf(fid,'chromosomes_scale = PPI:6\n');
fclose(fid)


%% write CIRCOS ideogram file
cd('/home/mbellis/sosma/data/cacn/circos/ppi')
fid=fopen('ideogram.conf','w');
fprintf(fid,'<ideogram>\n\n');
fprintf(fid,'<spacing>\n\n');
fprintf(fid,'default = 0.5u\n');

% fprintf(fid,'<pairwise CACNA1B;CACNA1C>\n');
% fprintf(fid,'spacing = 2u\n');
% fprintf(fid,'</pairwise>\n');
% 
% 
% fprintf(fid,'<pairwise CACNA1C;CACNA1D>\n');
% fprintf(fid,'spacing = 2u\n');
% fprintf(fid,'</pairwise>\n');
% 
fprintf(fid,'<pairwise Cav3.3;CACNA2D1>\n');
fprintf(fid,'spacing = 2u\n');
fprintf(fid,'</pairwise>\n');

fprintf(fid,'<pairwise CACNA2D4;CACNB1>\n');
fprintf(fid,'spacing = 2u\n');
fprintf(fid,'</pairwise>\n');

fprintf(fid,'<pairwise CACNB4;CACNG1>\n');
fprintf(fid,'spacing = 2u\n');
fprintf(fid,'</pairwise>\n');

fprintf(fid,'<pairwise CACNG8;PPI>\n');
fprintf(fid,'spacing = 4u\n');
fprintf(fid,'</pairwise>\n');

fprintf(fid,'<pairwise PPI;Cav1.1>\n');
fprintf(fid,'spacing = 4u\n');
fprintf(fid,'</pairwise>\n');

fprintf(fid,'</spacing>\n\n');
fprintf(fid,'<<include ideogram.position.conf>>\n');
fprintf(fid,'<<include ideogram.label.conf>>\n');
fprintf(fid,'</ideogram>\n\n');
fclose(fid)

%% write CIRCOS ideogram position
cd('/home/mbellis/sosma/data/cacn/circos/ppi')
fid=fopen('ideogram.position.conf','w');
fprintf(fid,'radius = 0.70r\n');
fprintf(fid,'thickness = 50p\n');
fprintf(fid,'fill = yes\n');
fprintf(fid,'fill_color = black\n');
fprintf(fid,'stroke_thickness = 2\n');
fprintf(fid,'stroke_color = black\n');
fclose(fid)

%% write CIRCOS links

%bundles and links
cd('/home/mbellis/sosma/data/cacn/circos/ppi')
fid=fopen('cacn_ppi_bundle.txt','w');
LinkRank=0;
for SubL=1:SubUnitNb
    Pos=find(PsiMat(:,SubL)==1);
    if ~isempty(Pos)
        LinkRank=LinkRank+1;
        fprintf(fid,'link%ubundle %s %u %u z=%u,color=bq1-%u\n',LinkRank,SubUnits{SubL},round((PpiNbList(SubL)+PpiNb)/2),round((PpiNbList(SubL)+PpiNb)/2),PpiNb+1,ColorOrder(CacnType(SubL)));
        fprintf(fid,'link%ubundle PPI %u %u z=%u,color=bq1-%u\n',LinkRank,min(Pos),max(Pos),PpiNb+1,ColorOrder(CacnType(SubL)));
    end
end

for SubL=1:SubUnitNb
    for PpiL=1:size(PsiMat,1)
        CurrPpiNb=PsiMat(PpiL,SubL);
        if CurrPpiNb>1
            LinkRank=LinkRank+1;
            fprintf(fid,'link%u %s %u %u z=%u,color=link-%u\n',LinkRank,SubUnits{SubL},round((PpiNbList(SubL)+PpiNb)/2+(PpiNbList(SubL)+PpiNb)/PpiNb),round((PpiNbList(SubL)+PpiNb)/2+(PpiNbList(SubL)+PpiNb)/PpiNb),CurrPpiNb,CurrPpiNb);
            fprintf(fid,'link%u PPI %u %u z=%u,color=link-%u\n',LinkRank,PpiL,PpiL,CurrPpiNb,CurrPpiNb);
        end
    end
end
fclose(fid)

%Only bundles
% cd('/home/mbellis/sosma/data/cacn/circos/ppi')
% fid=fopen('cacn_ppi_bundle.txt','w');
% for SubL=1:SubUnitNb
%     Pos=find(PsiMat(:,SubL)==1);
%     if ~isempty(Pos)
%         fprintf(fid,'link%ubundle %s %u %u z=1,color=bq1-%u\n',LinkRank,SubUnits{SubL},3,3,ColorOrder(CacnType(SubL)));
%         fprintf(fid,'link%ubundle PPI %u %u z=1,color=bq1-%u\n',LinkRank,min(Pos),max(Pos),ColorOrder(CacnType(SubL)));
%     end
% end
% LinkRank=0;
% for SubL=1:SubUnitNb
%     MemPpiNb=-1;
%     for PpiL=1:size(PsiMat,1)
%         CurrPpiNb=PsiMat(PpiL,SubL);
%         if CurrPpiNb>1
%             PpiL
%             if CurrPpiNb~=MemPpiNb             
%                 if MemPpiNb>0           
%                     LinkRank=LinkRank+1;
%                     fprintf(fid,'link%u %s %u %u z=%u,color=link-%u\n',LinkRank,SubUnits{SubL},MemPpiNb,MemPpiNb,MemPpiNb,MemPpiNb);
%                     fprintf(fid,'link%u PPI %u %u z=%u,color=link-%u\n',LinkRank,StartPos,PpiL,MemPpiNb,CurrPpiNb);
%                     MemPpiNb=CurrPpiNb;
%                     StartPos=PpiL;
%                 else
%                     MemPpiNb=CurrPpiNb;
%                     StartPos=PpiL;
%                 end
%             end
%         end
%     end
%     if StartPos~=PpiL & CurrPpiNb>1
%         z=z+1;
%         LinkRank=LinkRank+1;
%         fprintf(fid,'link%u %s %u %u z=%u,color=link-%u\n',LinkRank,SubUnits{SubL},MemPpiNb,MemPpiNb,MemPpiNb,MemPpiNb);
%         fprintf(fid,'link%u PPI %u %u z=%u,color=link-%u\n',LinkRank,StartPos,PpiL,MemPpiNb,CurrPpiNb);
%     end
% end
% fclose(fid)


% Individuals links
% cd('/home/mbellis/sosma/data/cacn/circos/ppi')
% fid=fopen('ppi_links.txt','w');
% LinkRank=0;
% for SubL=1:SubUnitNb
%     SubPos=-1;
%     for PpiL=1:size(PsiMat,1)
%         CurrPpiNb=PsiMat(PpiL,SubL);
%         if CurrPpiNb>0
%             LinkRank=LinkRank+1;
%             fprintf(fid,'link%u %s %u %u color=link-%u\n',LinkRank,SubUnits{SubL},CurrPpiNb,CurrPpiNb,CurrPpiNb);
%             fprintf(fid,'link%u PPI %u %u color=link-%u\n',LinkRank,PpiL,PpiL,CurrPpiNb);
%         end
%     end
% end
% fclose(fid)

%% Find ps in the nervous region (5) of m27 in different class
for i=1:7
Ps{i}=find(PsMatrix(:,17)==i & Clu(:,1)==5);
end

