%===========================%
% FUNCTION NET_DISPLAYLIST  %
%===========================%

% NET_DISPLAYLIST displays  Netstensor results for a given list in several networks

%INPUT PARAMETERS

%  1      ChipRank:
%  2     DirSuffix: directory suffix (e.g. _2to22_1p is suffix of directory m27_2to22_1p)
%  2     NetRanks1: network ranks used by PsMapping (information on probe sets)
%  3     NetRanks2: network used for recovering CORR and ANTI values
%  4      ListName:
%  5   PathwayName:
%  6       Species:
%  7     LocalFlag: load A,C,F data from local disk (1) or external disk
%  8     ClassFlag: split probe sets between s (SS,MS,CX(SS),CX(MS),HX(SS)HX(MS))
%                   and rest (All% - s)
%  9  DiscreteFlag: significative correlation values (i.e. >0) are written as 1 in netstensor
%                   files
% 10 RedundantFlag: keep all probesets (1) or only one probeset among redundant probesets
% 11      CorrType: indicates how is calculated correlation between probesets by the rank in
%                   the CorrTypeName list
%                   {'C','C-A','rawC','raw(C-A)','scoreC','scoreC-scoreA'}
% 12       ValFlag: keep in memory correlation values of all networsk
% 13    QLimitFlag: indicates if qlimit has been used to construct networks
% 14       TsnFlag: 1: write files necessary to do NetsTensor analysis
%                   0: load, process and save results of NetsTensor analysis

%OUTPUT FILE

%noqlimit:
%C-A
%net_displaylist(27,'_2to22_1p',[2:22],[119:133],'mouse_krebs_proteasome_mapk',{'Krebs','Proteasome','MAPPK'},'mouse',0,1,0,1,2,1,0,1)
%C
% write TSN files
%net_displaylist(27,'_2to22_1p',[2:22],[119:133],'mouse_krebs_proteasome_mapk',{'Krebs','Proteasome','MAPPK'},'mouse',0,1,0,1,1,1,0,1)
% load TSN results
%net_displaylist(27,'_2to22_1p',[2:22],[119:133],'mouse_krebs_proteasome_mapk',{'Krebs','Prot%easome','MAPPK'},'mouse',0,1,0,1,1,1,0,0)
%net_displaylist(5,'_1p',[17:21,35:64],[89:104],'mouse_krebs_proteasome_mapk',{'Krebs','Prote%asome','MAPPK'},'mouse',0,1,0,1,1,1,0,0)

%net_displaylist('8','_1p',[7:21],[119:133],'mouse_krebs_proteasome_mapk',{'Krebs','Proteasome','MAPPK'},'mouse',0,1,0,1,1,1,0,1)


%qlimit
%net_displaylist(27,[2:22],[2:22],'mouse_krebs_proteasome_mapk',{'Krebs','Proteasome','MAPPK'},'mouse')
%net_displaylist(27,[2:22],[2:22],'mouse_mapk',{'MAPK'},'mouse',0,1,0,1,2,1)
%net_displaylist(8,[7:21],[7:21],'mouse_mapk',{'MAPK'},'mouse',0,1,0,1,2,1)
%net_displaylist(5,[17:21,35:64],[17:21,35:64],'mouse_mapk',{'MAPK'},'mouse',0,1,0,1,2,1)
%net_displaylist(2,[12:32],[12:32],'human_mapk',{'MAPK'},'human',0,1,0,1,2,1)
%net_displaylist(3,[6:41],[6:41],'human_mapk',{'MAPK'},'human',0,1,0,1,2,1)
%net_displaylist(6,[6:20],[6:20],'rat_mapk',{'MAPK'},'rat',0,1,0,1,2,1)

%net_displaylist(8,[24:38],'mouse_krebs_proteasome_mapk','mouse')

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


function net_displaylist(ChipRank,DirSuffix,NetRanks1,NetRanks2,ListName,PathwayName,Species,LocalFlag,ClassFlag,DiscreteFlag,RedundantFlag,CorrType,ValFlag,QLimitFlag,TsnFlag)
global K

%% load ensembl id and gene name
% ex:
% EnsId{1} = ENSMUSG00000085299
% Def{1} = predicted gene, 16627 [Source:MGI Symbol;Acc:MGI:4439551]
% Gene{1} = predicted gene, 16627 [Source:MGI Symbol;Acc:MGI:4439551]

cd(K.dir.gene)
[EnsId,Def,Gene]=textread(sprintf('%s_ensembl.txt',Species),'%s%s%s','delimiter','\t');

%% load pathways
% recover Pathway variable (structure)
% ex :
% Pathway(1) =
%    keggGeneRank: [31x1 double]
%            name: 'Citrate cycle (TCA cycle) - Mus musculus (mouse)'
%     ensGeneRank: [31x1 double]

NetNb=length(NetRanks2);
PsNb=K.chipSet.probesetNb(ChipRank);
if CorrType==1|CorrType==3
    %C or raw(C)
    if QLimitFlag
        %for continuous analysis
        Densities{1}=0.3:0.1:0.6;
        %for discrete analysis
        Densities{2}=0.4:0.1:1.0;
    else
        Densities{1}=0.4:0.1:0.6;
        Densities{2}=0.5:0.1:1.0;
    end
else
    %C-A or raw(C)-raw(A)
    if QLimitFlag
        Densities{1}=0:0.1:0.5;
        Densities{2}=0.1:0.1:0.6;
    else
        Densities{1}=0:10.1:0.5;
        Densities{2}=0.2:0.1:0.6;
    end
end
cd(K.dir.list)
load(sprintf('%s_ens',ListName))
ListNb=length(Pathway);
switch Species
    case 'mouse'
        EnsPrefix='ENSMUSG';
    case 'human'
        EnsPrefix='ENSG';
    case 'rat'
        EnsPrefix='ENSRNOG';
end

%% load list of ps by class
%probe sets are in lists corresponding to classes, sub-classes or combination of classes
%ss: single ps targeting a single gene
%md: multiple ps (redundant ps) targeting a single gene
%cx: subset of ss and md combination existing in complex and hypercomplex classes
%s: union of ss, md and cs
%rest: all - s
ClassName={'all','ss','ms','cx','s','rest'};
CorrTypeName={'corr','diff','rawcorr','rawdiff'};


cd(K.dir.chip)
for ClassL=1:6
    fid=fopen(sprintf('m%u_pslist%u.u32',ChipRank,ClassL),'r','ieee-le');
    ClassPsRank{ClassL}=fread(fid,inf,'uint32');
    fclose(fid);
end


cd(sprintf('/home/mbellis/sosma/data/psawn/mldata/%s/m%u%s',Species,ChipRank,DirSuffix))
eval(sprintf('load m%u_n%u_netnb%u_probenb1_newps_stat_netprc100_pvcorr1',ChipRank,NetRanks1(1),length(NetRanks1)));

%find ps ranks corresponding ot ensGeneRank of each pathway (=> may have multiple probeset
%for some genes
PsRank=cell(ListNb,1);
EnsPos=cell(ListNb,1);
PsClass=cell(ListNb,1);
EnsRank=cell(ListNb,1);
for ListL=1:ListNb
    PsRank{ListL}=[];
    EnsPos{ListL}=[];
    PsClass{ListL}=[];
    EnsRank{ListL}=[];
    for GeneL=1:length(Pathway(ListL).ensGeneRank)
        CurrGene=sprintf('%s%011u',EnsPrefix,Pathway(ListL).ensGeneRank(GeneL));
        GenePos=strmatch(CurrGene,Genes.name,'exact');
        if ~isempty(GenePos)
            %find ps assigned to the current ensembl gene
            CurrPsRank=find(PsMatrix(:,1)==GenePos);
            PsRank{ListL}=[PsRank{ListL};CurrPsRank];
            EnsPos{ListL}=[EnsPos{ListL};repmat(GeneL,length(CurrPsRank),1)];
            EnsRank{ListL}=[EnsRank{ListL};repmat(Pathway(ListL).ensGeneRank(GeneL),length(CurrPsRank),1)];
            for PsL=1:length(CurrPsRank)
                FoundPs=0;
                for ClassL=[2,3,4,6]
                    if ~isempty(find(ClassPsRank{ClassL}==CurrPsRank(PsL)))
                        PsClass{ListL}=[PsClass{ListL};ClassL];
                        FoundPs=1;
                        break
                    end
                end
                if FoundPs==0
                    PsClass{ListL}=[PsClass{ListL};0];
                end
            end
        end
    end
end

%reorder according to classes
for ListL=1:ListNb
    [PsClass{ListL} SortOrder]=sort(PsClass{ListL});
    EnsPos{ListL}=EnsPos{ListL}(SortOrder);
    PsRank{ListL}=PsRank{ListL}(SortOrder);
    EnsRank{ListL}=EnsRank{ListL}(SortOrder);
end



%concatenate all pathways

PsRanks=[];
EnsPoss=[];
PsClasses=[];
EnsRanks=[];
for ListL=1:ListNb
    PsRanks=[PsRanks;PsRank{ListL}];
    EnsPoss=[EnsPoss;EnsPos{ListL}];
    PsClasses=[PsClasses;PsClass{ListL}];
    EnsRanks=[EnsRanks;EnsRank{ListL}];
end

%complete gene information
Defs=cell(length(EnsRanks),1);
Genes=cell(length(EnsRanks),1);
for GeneL=1:length(EnsRanks)
    GenePos=strmatch(sprintf('%s%011u',EnsPrefix,EnsRanks(GeneL)),EnsId,'exact');
    if ~isempty(GenePos)
        Defs{GeneL}=Def{GenePos};
        Genes{GeneL}=Gene{GenePos};
    end
end

if TsnFlag

    if RedundantFlag==0
        %construct a list of non redundant probesets (keep the first in list)
        NrPsIndex=zeros(length(PsRanks),1);
        Offset=0;
        TickOffset=0;
        NrTick=1;
        for ListL=1:ListNb
            UniquePsPos=unique(EnsPos{ListL});
            NrPsRank{ListL}=[];
            for PsL=1:length(UniquePsPos);
                Pos=find(EnsPos{ListL}==UniquePsPos(PsL));
                NrPsIndex(Pos(1)+Offset)=1;
                NrPsRank{ListL}(end+1,1)=PsRank{ListL}(Pos(1));
            end
            Offset=Offset+length(EnsPos{ListL});
            NrTick=[NrTick,round(length(UniquePsPos)/2)+TickOffset,length(UniquePsPos)+TickOffset];
            TickOffset=TickOffset+length(UniquePsPos);
        end
        NrPsIndex=find(NrPsIndex);
        PsRanks=PsRanks(NrPsIndex);
        EnsPoss=EnsPoss(NrPsIndex);
        PsClasses=PsClasses(NrPsIndex);
        EnsRanks=EnsRanks(NrPsIndex);
        Defs=Defs(NrPsIndex);
        Genes=Genes(NrPsIndex);
    end
    clear ClassPsRank LinkedPs NewPs PsBy PsMtrix Stat

    %% LOAD C (A,F) DATA

    Offset=0;
    PathTick=1;
    PathLabel={''};
    for ListL=1:ListNb
        PathTick=[PathTick,round(length(PsRank{ListL})/2)+Offset,length(PsRank{ListL})+Offset];
        Offset=Offset+length(PsRank{ListL});
        PathLabel{1,end+1}=PathwayName{ListL};
        PathLabel{1,end+1}='';
    end

    Offset=0;
    ClassTick=1;
    ClassLabel={''};
    for ListL=1:ListNb
        for ClassL=[2,3,4,6]
            ClassPos=find(PsClass{ListL}==ClassL);
            if  ~isempty(ClassPos)
                ClassTick=[ClassTick,round(length(ClassPos)/2)+Offset,length(ClassPos)+Offset];
                Offset=Offset+length(ClassPos);
                ClassLabel{1,end+1}=num2str(ClassL);
                ClassLabel{1,end+1}='';
            end
        end
    end



    for NetL=1:NetNb
        sprintf('loading net %u',NetL)
        if LocalFlag
            cd(sprintf('/home/mbellis/net/clum%u',ChipRank'))
        else
            cd(fullfile(K.dir.net,sprintf('m%u',ChipRank),sprintf('n%u',NetRanks2(NetL))))
        end
        C=load_data(sprintf('c_m%u_n%u.4mat',ChipRank,NetRanks2(NetL)),'./',PsNb,PsNb,'uint8','ieee-le',PsRanks,PsRanks);
        C=double(C);
        if CorrType==2|CorrType==4
            A=load_data(sprintf('a_m%u_n%u.4mat',ChipRank,NetRanks2(NetL)),'./',PsNb,PsNb,'uint8','ieee-le',PsRanks,PsRanks);
            A=double(A);
        end
        if CorrType==3|CorrType==4
            F=load_data(sprintf('f_m%u_n%u.4mat',ChipRank,NetRanks2(NetL)),'./',PsNb,PsNb,'uint8','ieee-le',PsRanks,PsRanks);
            F=double(F);
        end

        switch CorrType
            case 2
                C=C-A;
            case 3
                C=(F/100).*C;
            case 4
                C=(F/100).*(C-A);
        end

        if NetL==1
            %histogram of  CORR values
            h=figure;
            if RedundantFlag
                set(h,'name',sprintf('m%u n%u %s with redundant ps',ChipRank,NetRanks2(NetL),strrep(ListName,'_',' ')))
            else
                set(h,'name',sprintf('m%u n%u %s without redundant ps',ChipRank,NetRanks2(NetL),strrep(ListName,'_',' ')))
            end
            set(gcf,'color',[1,1,1])
            c=C(:);
            c(find(c==0))=[];
            hist(double(c),100)
            xlabel('CORR')
            ylabel('freq')
            if RedundantFlag
                title(sprintf('m%u n%u %s with redundant ps',ChipRank,NetRanks2(NetL),strrep(ListName,'_',' ')))
            else
                title(sprintf('m%u n%u %s without redundant ps',ChipRank,NetRanks2(NetL),strrep(ListName,'_',' ')))
            end
        end

        %WRITE TENSORNET FILES
        cd(sprintf('/home/mbellis/net/clum%u',ChipRank'))
        fid=fopen(sprintf('m%u_n%u_%s_%s_all.network',ChipRank,NetRanks2(NetL),ListName,CorrTypeName{CorrType}),'w');
        for LineL=1:length(C)
            Pos=find(C(LineL,LineL+1:end)>0);
            if ~isempty(Pos)
                for PosL=1:length(Pos)
                    fprintf(fid,'%u\t%u\t%.2f\n',LineL,Pos(PosL)+LineL,single(C(LineL,LineL+Pos(PosL)))/100);
                end
            end
        end
        fclose(fid)

        if DiscreteFlag
            fid=fopen(sprintf('m%u_n%u_%s_%s_alld.network',ChipRank,NetRanks2(NetL),ListName,CorrTypeName{CorrType}),'w');
            for LineL=1:length(C)
                Pos=find(C(LineL,LineL+1:end)>0);
                if ~isempty(Pos)
                    for PosL=1:length(Pos)
                        fprintf(fid,'%u\t%u\t1.00\n',LineL,Pos(PosL)+LineL);
                    end
                end
            end
            fclose(fid)
        end
        if ClassFlag
            %add subset according to probest class (s and rest)
            for ClassL=1:4
                switch ClassL
                    case {1,3}
                        if ClassL==1
                            fid=fopen(sprintf('m%u_n%u_%s_%s_s.network',ChipRank,NetRanks2(NetL),ListName,CorrTypeName{CorrType}),'w');
                        else
                            if DiscreteFlag
                                fid=fopen(sprintf('m%u_n%u_%s_%s_sd.network',ChipRank,NetRanks2(NetL),ListName,CorrTypeName{CorrType}),'w');
                            end
                        end
                        c=C;
                        Pos=find(PsClasses==6);
                        c(Pos,:)=0;
                        c(:,Pos)=0;
                    case {2,4}
                        if ClassL==2
                            fid=fopen(sprintf('m%u_n%u_%s_%s_rest.network',ChipRank,NetRanks2(NetL),ListName,CorrTypeName{CorrType}),'w');
                        else
                            if DiscreteFlag
                                fid=fopen(sprintf('m%u_n%u_%s_%s_restd.network',ChipRank,NetRanks2(NetL),ListName,CorrTypeName{CorrType}),'w');
                            end
                        end
                        c=C;
                        Pos=find(PsClasses<6);
                        c(Pos,:)=0;
                        c(:,Pos)=0;
                end
                for LineL=1:length(c)
                    Pos=find(c(LineL,LineL+1:end)>0);
                    if ~isempty(Pos)
                        if ClassL>2
                            if DiscreteFlag
                                for PosL=1:length(Pos)
                                    fprintf(fid,'%u\t%u\t1.00\n',LineL,Pos(PosL)+LineL);
                                end
                            end
                        else
                            for PosL=1:length(Pos)
                                fprintf(fid,'%u\t%u\t%.2f\n',LineL,Pos(PosL)+LineL,single(c(LineL,LineL+Pos(PosL)))/100);
                            end
                        end
                    end
                end
            end
            fclose(fid)
        end

        %% NETWORK REPRODUCIBILITY
        % computes the number of networks in wich a particular edge between two probesets has a value
        % greater than a given limit
        % calculate also MeanSortOrder (allows to reorder matrix according to the mean of CORR
        % values)

        Limit=[0:10:60,100];
        if NetL==1
            for LimitL=1:7
                NetFreq{LimitL}=zeros(size(C));
                Pos=find(C>Limit(LimitL)&C<100);
                NetFreq{LimitL}(Pos)=1;
            end
            %mean C
            NetFreq{8}=double(C);
        else
            for LimitL=1:7
                Pos=find(C>Limit(LimitL)&C<100);
                NetFreq{LimitL}(Pos)=NetFreq{LimitL}(Pos)+1;
            end
            NetFreq{8}=NetFreq{8}+double(C);
        end
        if ValFlag
            Corr{NetL}=C;
        end
        %mean C
        NetFreq{8}=NetFreq{8}/NetNb;
    end

    %reorder matrix by using mean




    MeanSortOrder=[];
    Offset=0;
    for ListL=1:ListNb
        if RedundantFlag
            CPos=1:length(PsRank{ListL});
        else
            CPos=1:length(NrPsRank{ListL});
        end
        CSum=sum(NetFreq{8}(CPos,CPos));
        [temp CurrSortOrder]=sort(CSum);
        MeanSortOrder=[MeanSortOrder,fliplr(CurrSortOrder)+Offset];
        if RedundantFlag
            Offset=Offset+length(PsRank{ListL});
        else
            Offset=Offset+length(NrPsRank{ListL});
        end
    end




    %MatchIndex=1:length(MeanSortOrder);
    %[temp InvSortOrder]=sort(MeanSortOrder);
    %MatchIndex=MatchIndex(InvSortOrder);


    h=figure;
    set(gcf,'color',[1,1,1])
    set(h,'name','network reproducibility according to CORR limit')
    hold on
    set(gcf,'color',[1,1,1])
    Colors=colors(colormap,7);
    for LimitL=1:7
        CurrL=NetFreq{LimitL}(:);
        NullPos=find(CurrL==0);
        CurrL(NullPos)=[];
        Values=histc(CurrL,1:NetNb);
        plot(1:NetNb,Values/sum(Values),'color',Colors(LimitL,:))
    end
    set(gca,'box','on')
    title('network reproducibility according to CORR limit')
    legend('>0','>10','>20','>30','>40','>50','>60')
    xlabel('network number')
    ylabel('freq')


    h=figure;
    set(gcf,'color',[1,1,1])
    set(h,'name','CORR > Limit  and MEAN')
    for LimitL=1:8
        subplot(4,2,LimitL)
        h=pcolor(NetFreq{LimitL}(MeanSortOrder,MeanSortOrder));
        %h=pcolor(NetFreq{LimitL});
        set(h,'linestyle','none')
        set(gca,'xtick',ClassTick)
        set(gca,'ytick',PathTick)
        set(gca,'yticklabel',PathLabel)
        set(gca,'yticklabel',ClassLabel)
        set(gca,'tickdir','out')
        if LimitL<8
            title(sprintf('network occurence>%u',Limit(LimitL)))
        else
            title('mean(CORR)')
        end
    end

    if ValFlag
        ColNb=ceil(sqrt(NetNb/2.5));
        RowNb=round(NetNb/ColNb);
        h=figure;
        set(gcf,'color',[1,1,1])
        set(h,'name','network correlation in all networks')
        for NetL=1:NetNb
            subplot(RowNb,ColNb,NetL)
            c=Corr{NetL};
            c(find(c>0))=0;
            h=pcolor(double(c(MeanSortOrder,MeanSortOrder)));
            set(h,'linestyle','none')
            title(sprintf('network %u (%u) - %s',NetL,NetRanks2(NetL),CorrTypeName{CorrType}))
        end
    end

    'stop'




    %write files used by tensor net

    cd(sprintf('/home/mbellis/net/clum%u',ChipRank'))


    %write list of networks
    fid=fopen(sprintf('m%u_%s_netlist_%s_all.txt',ChipRank,ListName,CorrTypeName{CorrType}),'w');
    for NetL=1:NetNb
        fprintf(fid,'m%u_n%u_%s_%s_all\n',ChipRank,NetRanks2(NetL),ListName,CorrTypeName{CorrType});
    end
    fclose(fid)
    if ClassFlag
        fid=fopen(sprintf('m%u_%s_netlist_%s_s.txt',ChipRank,ListName,CorrTypeName{CorrType}),'w');
        for NetL=1:NetNb
            fprintf(fid,'m%u_n%u_%s_%s_s\n',ChipRank,NetRanks2(NetL),ListName,CorrTypeName{CorrType});
        end
        fclose(fid)
        fid=fopen(sprintf('m%u_%s_netlist_%s_rest.txt',ChipRank,ListName,CorrTypeName{CorrType}),'w');
        for NetL=1:NetNb
            fprintf(fid,'m%u_n%u_%s_%s_rest\n',ChipRank,NetRanks2(NetL),ListName,CorrTypeName{CorrType});
        end
        fclose(fid)
    end

    if DiscreteFlag
        fid=fopen(sprintf('m%u_%s_netlist_%s_alld.txt',ChipRank,ListName,CorrTypeName{CorrType}),'w');
        for NetL=1:NetNb
            fprintf(fid,'m%u_n%u_%s_%s_alld\n',ChipRank,NetRanks2(NetL),ListName,CorrTypeName{CorrType});
        end
        fclose(fid)


        if ClassFlag
            fid=fopen(sprintf('m%u_%s_netlist_%s_sd.txt',ChipRank,ListName,CorrTypeName{CorrType}),'w');
            for NetL=1:NetNb
                fprintf(fid,'m%u_n%u_%s_%s_sd\n',ChipRank,NetRanks2(NetL),ListName,CorrTypeName{CorrType});
            end
            fclose(fid)
            fid=fopen(sprintf('m%u_%s_netlist_%s_restd.txt',ChipRank,ListName,CorrTypeName{CorrType}),'w');
            for NetL=1:NetNb
                fprintf(fid,'m%u_n%u_%s_%s_restd\n',ChipRank,NetRanks2(NetL),ListName,CorrTypeName{CorrType});
            end
            fclose(fid)
        end
    end

    %gene ids is used only by NetsTensor
    fid=fopen('gene_ids','w');
    for PsL=1:PsNb
        fprintf(fid,'%u\n',PsL);
    end
    fclose(fid)

    ListNames={'all','s','rest','alld','sd','restd'};
    if ClassFlag
        if DiscreteFlag
            ListPos=1:6;
        else
            ListPos=1:3;
        end
    else
        if DiscreteFlag
            ListPos=1;
        else
            ListPos=[1,4];
        end
    end
    if QLimitFlag
        NetNbs=[round(NetNb/2),NetNb-5:NetNb];
    else
        NetNbs=NetNb;
    end
    MinGene=5;


    for ListL=1:length(ListPos)
        CurrListPos=ListPos(ListL);
        CurrListName=ListNames{CurrListPos};

        fid=fopen(sprintf('tsn%u_%s.sh',CurrListPos,CorrTypeName{CorrType}),'w');
        NetworkPath='.';
        if QLimitFlag
            ResultPath=sprintf('./%s',CurrListName);
        else
            ResultPath=sprintf('./noqlimit/%s',CurrListName);
        end
        DatasetsListFile=sprintf('m%u_%s_netlist_%s_%s.txt',ChipRank,ListName,CorrTypeName{CorrType},CurrListName);
        for NetL=1:length(NetNbs)
            for DensL=1:length(Densities{1})
                ResultFile=sprintf('m%u_n%uton%u_%s_%s_%s_g%u_n%u_d%u',ChipRank,NetRanks2(1),NetRanks2(2),ListName,CurrListName,CorrTypeName{CorrType},...
                    MinGene,NetNbs(NetL),round(Densities{1}(DensL)*100));
                fprintf(fid,'/usr/local/netstensor/netstensor --DatasetsListFile="%s" --minGene=%u --minNet=%u --minDensity=%.2f --networksPath="%s" --resultPath="%s" --ResultFile="%s" &\n',...
                    DatasetsListFile,...
                    MinGene,...
                    NetNbs(NetL),...
                    Densities{1}(DensL),...
                    NetworkPath,...
                    ResultPath,...
                    ResultFile);
            end
        end
        fclose(fid)
    end


    'stop'
else

    %% READ TENSORNET FILE (EACH PATHWAYS)
    NetNbs=[round(NetNb/2),NetNb-5:NetNb];

    ListNames={'all','s','rest','alld','sd','restd'};
    CorrTypeName={'corr','diff','rawcorr','rawdiff'};
    ClassName={'all','s','rest'};

    if QLimitFlag
        NetNbs=[round(NetNb/2),NetNb-5:NetNb];
    else
        NetNbs=NetNb;
    end
    MinGene=5;


    if ClassFlag
        ClassNb=3;
    else
        ClassNb=1;
    end
    if DiscreteFlag
        TypeNb=2;
    else
        TypeNb=1;
    end
    if QLimitFlag
        Dir=sprintf('/home/mbellis/net/clum%u',ChipRank');
    else
        Dir=sprintf('/home/mbellis/net/clum%u/noqlimit',ChipRank');
    end

    for TypeL=1:TypeNb
        Ps{TypeL}={};
        PsSize{TypeL}={};
        Density{TypeL}={};
        NbOfNets{TypeL}={};
        Nets{TypeL}={};
        for ClassL=1:ClassNb
            if TypeL==1
                CurrListName=ListNames{ClassL};
                if TypeL==2
                    CurrListName=ListNames{ClassL+3};
                end
                cd(sprintf('%s/%s',Dir,CurrListName'))
                for NetL=1:length(NetNbs)
                    for DensL=1:length(Densities{TypeL})
                        ResultFile=sprintf('m%u_n%uton%u_%s_%s_%s_g%u_n%u_d%u',ChipRank,NetRanks2(1),NetRanks2(2),ListName,CurrListName,CorrTypeName{CorrType},...
                            MinGene,NetNbs(NetL),round(Densities{1}(DensL)*100));
                        try
                            [Ps{TypeL}{ClassL}{NetL,DensL},PsSize{TypeL}{ClassL}{NetL,DensL},Density{TypeL}{ClassL}{NetL,DensL},NbOfNets{TypeL}{ClassL}{NetL,DensL},Nets{TypeL}{ClassL}{NetL,DensL}]=textread(ResultFile,'%s%u%.2f%u%s','delimiter','\t');
                        catch
                            sprintf('result %u m%u_n%uton%u_%s_%s_%s_g%u_n%u_d%u does not exists',ClassL,ChipRank,NetRanks2(1),NetRanks2(2),ListName,CurrListName,CorrTypeName{CorrType},...
                                MinGene,NetNbs(NetL),round(Densities{1}(DensL)*100))
                        end
                    end
                end
            end
        end


        %eval expressions
        for ClassL=1:ClassNb
            for NetL=1:length(NetNbs)
                for DensL=1:length(Densities{TypeL})
                    try
                        if ~isempty(Ps{TypeL}{ClassL}{NetL,DensL})
                            for CluL=1:length(Ps{TypeL}{ClassL}{NetL,DensL});
                                Ps{TypeL}{ClassL}{NetL,DensL}{CluL}=sort(eval(Ps{TypeL}{ClassL}{NetL,DensL}{CluL}));
                            end
                        else
                            Ps{TypeL}{ClassL}{NetL,DensL}={};
                        end
                    catch
                    end
                end
            end
        end

        %keep clusters present exactly in the minNet number of genes
        PsMinNet{TypeL}={};
        DensityMinNet{TypeL}={};
        for ClassL=1:ClassNb
            for NetL=1:length(NetNbs)
                for DensL=1:length(Densities{TypeL})
                    try
                        if ~isempty(Ps{TypeL}{ClassL}{NetL,DensL})
                            Pos=find(NbOfNets{TypeL}{ClassL}{NetL,DensL}==NetNbs(NetL));
                            PsMinNet{TypeL}{ClassL}{NetL,DensL}=Ps{TypeL}{ClassL}{NetL,DensL}(Pos);
                            DensityMinNet{TypeL}{ClassL}{NetL,DensL}=Density{TypeL}{ClassL}{NetL,DensL}(Pos);
                        else
                            PsMinNet{TypeL}{ClassL}{NetL,DensL}={};
                            DensityMinNet{TypeL}{ClassL}{NetL,DensL}={};
                        end
                    catch
                    end
                end
            end
        end
    end

    %number of clusters
    for TypeL=1:TypeNb
        h=figure;
        if TypeL==1
            set(h,'name','number of clusters - continuous corr')
        else
            set(h,'name','number of clusters - discrete corr')
        end
        for ClassL=1:ClassNb
            set(gcf,'color',[1,1,1])
            Colors=colors(colormap,length(NetNbs));
            Legend={};
            CluNb{ClassL}=zeros(length(NetNbs),length(Densities{TypeL}));
            CluNbMinNet{ClassL}=zeros(length(NetNbs),length(Densities{TypeL}));
            CluSize{ClassL}=zeros(length(NetNbs),length(Densities{TypeL}));
            CluSizeMinNet{ClassL}=zeros(length(NetNbs),length(Densities{TypeL}));
            for NetL=1:length(NetNbs)
                CluNb{ClassL}(NetL,:)=zeros(1,length(Densities{TypeL}));
                for DensL=1:length(Densities{TypeL})
                    try
                        CluNb{ClassL}(NetL,DensL)=length(PsSize{TypeL}{ClassL}{NetL,DensL});
                        PsClu=[];
                        for CluL=1:length(PsSize{TypeL}{ClassL}{NetL,DensL});
                            PsClu=union(PsClu,Ps{TypeL}{ClassL}{NetL,DensL}{CluL});
                        end
                        CluSize{ClassL}(NetL,DensL)=length(PsClu);
                    catch
                        CluSize{ClassL}(NetL,DensL)=0;
                    end
                end
                subplot(ClassNb,4,(ClassL-1)*4+1)
                hold on
                plot(Densities{TypeL},CluNb{ClassL}(NetL,:),'color',Colors((NetL),:))
                subplot(ClassNb,4,(ClassL-1)*4+2)
                hold on
                plot(Densities{TypeL},CluSize{ClassL}(NetL,:),'color',Colors((NetL),:))
                Legend{end+1,1}=sprintf('>=%u networks',NetNbs(NetL));


                CluNbMinNet{ClassL}(NetL,:)=zeros(1,length(Densities{TypeL}));
                for DensL=1:length(Densities{TypeL})
                    try
                        CluNbMinNet{ClassL}(NetL,DensL)=length(PsMinNet{TypeL}{ClassL}{NetL,DensL});
                        PsClu=[];
                        for CluL=1:length(PsMinNet{TypeL}{ClassL}{NetL,DensL});
                            PsClu=union(PsClu,PsMinNet{TypeL}{ClassL}{NetL,DensL}{CluL});
                        end
                        CluSizeMinNet{ClassL}(NetL,DensL)=length(PsClu);
                    catch
                        CluSizeMinNet{ClassL}(NetL,DensL)=0;
                    end
                end
                subplot(ClassNb,4,(ClassL-1)*4+3)
                hold on
                plot(Densities{TypeL},CluNbMinNet{ClassL}(NetL,:),'color',Colors((NetL),:))
                subplot(ClassNb,4,(ClassL-1)*4+4)
                hold on
                plot(Densities{TypeL},CluSizeMinNet{ClassL}(NetL,:),'color',Colors((NetL),:))
            end
            subplot(ClassNb,4,(ClassL-1)*4+1)
            title(sprintf('%s - %s',strrep(ListName,'_',' '),ClassName{ClassL}))
            xlabel('density')
            ylabel('number of clusters')
            %set(gca,'ylim',[0,100])
            set(gca,'box','on')
            subplot(ClassNb,4,(ClassL-1)*4+2)
            xlabel('density')
            ylabel('number of probe sets')
            set(gca,'box','on')

            subplot(ClassNb,4,(ClassL-1)*4+3)
            title(sprintf('%s - %s -minNet',strrep(ListName,'_',' '),ClassName{ClassL}))
            xlabel('density')
            ylabel('number of clusters')
            %set(gca,'ylim',[0,100])
            set(gca,'box','on')
            subplot(ClassNb,4,(ClassL-1)*4+4)
            xlabel('density')
            ylabel('number of probe sets')
            set(gca,'box','on')
        end
        legend(Legend)
    end


    % %union of clusters
    % for ClassL=1:ClassNb
    %     for NetL=1:length(NetNbs)
    %         NetL
    %         %for NetL=1
    %         %for DensL=1:length(Densities)
    %         for DensL=1:length(Densities)
    %             ChangedFlag=1;
    %             TestedClu=Ps{ClassL}{NetL,DensL};
    %             SelPs{ClassL}{NetL,DensL}={};
    %             SelClu={};
    %             Sel=zeros(length(TestedClu),1);
    %             while ChangedFlag
    %                 ChangedFlag=0;
    %                 for CluL1=1:length(TestedClu)-1
    %                     for CluL2=CluL1+1:length(TestedClu)
    %                         if Sel(CluL1)==0|Sel(CluL2)==0
    %                             if length(intersect(TestedClu{CluL1},TestedClu{CluL2}))==length(TestedClu{CluL1})|length(intersect(TestedClu{CluL1},TestedClu{CluL2}))==length(TestedClu{CluL2})
    %                                 SelClu{end+1,1}=union(TestedClu{CluL1},TestedClu{CluL2});
    %                                 Sel(CluL1)=1;
    %                                 Sel(CluL2)=1;
    %                                 ChangedFlag=1;
    %                             end
    %                         end
    %                     end
    %                 end
    %             end
    %             Pos=find(Sel==0);
    %             if ~isempty(Pos)
    %                 for PosL=1:length(Pos)
    %                     SelClu{end+1,1}=TestedClu{Pos(PosL)};
    %                 end
    %             end
    %             if ChangedFlag
    %                 TestedClu=SelClu;
    %                 SelClu={};
    %                 Sel=zeros(length(TestedClu),1);
    %             else
    %                 SelPs{ClassL}{NetL,DensL}=TestedClu;
    %             end
    %         end
    %     end
    % end




    %plot found clusters
    for TypeL=1:TypeNb

        FreqPs{TypeL}=cell(1,ClassNb);
        for ClassL=1:ClassNb
            if length(PsMinNet{TypeL})>=ClassL
                h=figure;
                if TypeL==1
                    set(h,'name',sprintf('Tensornet clusters for %s and %s - continuous corr',ListName,ClassName{ClassL}))
                else
                    set(h,'name',sprintf('Tensornet clusters for %s and %s - discrete corr',ListName,ClassName{ClassL}))
                end
                set(gcf,'color',[1,1,1])

                %CurrC=NetFreq{8}(MeanSortOrder,MeanSortOrder);
                %             MeanC=NetFreq{8};
                %             MeanC(find(MeanC<0))=0;
                %             CurrC=NetFreq{LimitPos};
                %             MeanC=[MeanC;CurrC*5];
                %             YLabel={'','mean(CORR) on all networks','',sprintf('Network reproducibility at CORR>%u',Limit(LimitPos)),''};
                %             YTick=[1,round(length(PsRanks)/2),length(PsRanks),round(1.5*length(PsRanks)),2*length(PsRanks)];
                MeanC=[];
                YTick=1;
                YLabel={''};
                FreqPs{TypeL}{ClassL}=zeros(1,length(PsRanks));
                for NetL=1:length(NetNbs)
                    for DensL=1:length(Densities{TypeL})
                        CurrCluNb=0;
                        try
                            CurrCluNb=length(PsMinNet{TypeL}{ClassL}{NetL,DensL});
                        catch
                        end
                        if CurrCluNb>0&CurrCluNb<100
                            if CurrCluNb>1
                                YTick=[YTick,YTick(end)+round(CurrCluNb/2),YTick(end)+CurrCluNb];
                            else
                                YTick=[YTick,YTick(end)+2,YTick(end)+3];
                            end
                            YLabel{1,end+1}=sprintf('netnb %u dens %.2f',NetNbs(NetL),Densities{TypeL}(DensL));
                            YLabel{1,end+1}='';
                            if CurrCluNb>1
                                for CluL=1:CurrCluNb;
                                    CurrLine=zeros(1,length(PsRanks));
                                    %CurrLine(PsMinNet{TypeL}{ClassL}{NetL,DensL}{CluL})=100-ceil(Step*CluL);
                                    CurrLine(PsMinNet{TypeL}{ClassL}{NetL,DensL}{CluL})=50;
                                    if isempty(MeanC)
                                        MeanC=CurrLine;
                                    else
                                        MeanC=[MeanC;CurrLine];
                                    end
                                    FreqPs{TypeL}{ClassL}(find(CurrLine))=FreqPs{TypeL}{ClassL}(find(CurrLine))+1;
                                end
                            else
                                CurrLine=zeros(3,length(PsRanks));
                                %CurrLine(2,PsMinNet{TypeL}{ClassL}{NetL,DensL}{1})=100-ceil(Step);
                                CurrLine(2,PsMinNet{TypeL}{ClassL}{NetL,DensL}{1})=50;
                                if isempty(MeanC)
                                    MeanC=CurrLine;
                                else
                                    MeanC=[MeanC;CurrLine];
                                end
                                FreqPs{TypeL}{ClassL}(find(CurrLine(2,:)))=FreqPs{TypeL}{ClassL}(find(CurrLine(2,:)))+1;
                            end

                        end
                    end
                end
                CurrLine=zeros(3,length(PsRanks));
                if isempty(MeanC)
                    MeanC=CurrLine;
                else
                    MeanC=[MeanC;CurrLine];
                end
                FreqPs{TypeL}{ClassL}=round(FreqPs{TypeL}{ClassL}*100/max(FreqPs{TypeL}{ClassL}));
                MeanC=[MeanC;repmat(FreqPs{TypeL}{ClassL},3,1)];





                h=pcolor([MeanC,zeros(size(MeanC,1),1);zeros(1,size(MeanC,2)+1)]);
                set(h,'linestyle','none')
                set(gca,'ytick',YTick)
                set(gca,'tickdir','out')
                set(gca,'yticklabel',YLabel)
            end
            %    title(PathwayName{ClassL})
        end
    end

    %Probeset occurence in clusters
    MeanC=zeros(1,length(PsRanks));
    YTick=[];
    YLabel={};
    for ClassL=1:3
        for TypeL=1:TypeNb
            if ~isempty(FreqPs{TypeL}{ClassL})
                if isempty(YTick)
                    YTick=2.5;
                else
                    YTick(1,end+1)=YTick(end)+1;
                end
                if TypeL==2
                    YLabel{1,end+1}=sprintf('%s - discrete',ClassName{ClassL});
                else
                    YLabel{1,end+1}=sprintf('%s - continuous',ClassName{ClassL});
                end
                MeanC=[MeanC;FreqPs{TypeL}{ClassL}];
            end
        end
    end
    MeanC=[MeanC;zeros(1,length(PsRanks))];
    h=figure;
    set(h,'name','Probe set occurence in clusters')
    set(gcf,'color',[1,1,1])
    h=pcolor([MeanC,zeros(size(MeanC,1),1);zeros(1,size(MeanC,2)+1)]);
    set(h,'linestyle','none')
    set(gca,'ytick',YTick)
    set(gca,'yticklabel',YLabel)
    title('Probe set occurence in clusters')




    %% SAVE RESULTS
    cd(sprintf('/home/mbellis/net/clum%u',ChipRank'))
    Res.densities=Densities{1};
    Res.clu=PsMinNet;
    Res.freqPs=FreqPs;
    Res.gene=Genes;
    Res.def=Defs;
    Res.ensRank=EnsRanks;
    Res.ensPos=EnsPoss;
    Res.psClass=PsClasses;
    Res.psRank=PsRanks;
    cd(sprintf('/home/mbellis/net/clum%u',ChipRank))
    eval(sprintf('save %s_m%u_n%uton%u_%s_ntsres Res',ListName,ChipRank,NetRanks2(1),NetRanks2(end),CorrTypeName{CorrType}))



end