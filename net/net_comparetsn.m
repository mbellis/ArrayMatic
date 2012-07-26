%==========================%
% FUNCTION NET_COMPARETSN  %
%==========================%

% NET_COMPARETSN compare two results from NetsTensor clustering

%INPUT PARAMETERS

%  1   ChipRank:
%  2   ListName:
%  3    Species:
%  4 CommonFlag: Exists a correspondance list between two chips models used
%  5 QLimitFlag:
%  6    CluFlag: (0)1 : (do not) process cluster information 
%  7 CluNbLimit: maximal number of clusters allowed to be processed if CluFlag=1
%OUTPUT FILE

%net_comparetsn([5,27],{'mouse_krebs_proteasome_mapk_m5_n89ton104_corr_ntsres','mouse_krebs_proteasome_mapk_m27_n119ton133_corr_ntsres'},{'mouse','mouse'},0,0,1)


%net_comparetsn([8,27],{'mouse_mapk_m8_ntsres','mouse_mapk_m27_ntsres'},{'mouse','mouse'},1)
%net_comparetsn([5,8],{'mouse_mapk_m5_ntsres','mouse_mapk_m8_ntsres'},{'mouse','mouse'},0)
%net_comparetsn([5,27],{'mouse_mapk_m5_ntsres','mouse_mapk_m27_ntsres'},{'mouse','mouse'},0)
%net_comparetsn([5,8,27],{'mouse_mapk_m5_ntsres','mouse_mapk_m8_ntsres','mouse_mapk_m27_ntsres'},{'mouse','mouse','mouse'},0)
%net_comparetsn([2,3],{'human_mapk_m2_ntsres','human_mapk_m3_ntsres'},{'human','human'},0,0)
%net_comparetsn(5,{'mouse_mapk_m5_ntsres'},{'mouse'},0,0)
%net_comparetsn(8,{'mouse_mapk_m8_ntsres'},{'mouse'},0,0)
%net_comparetsn(27,{'mouse_mapk_m27_ntsres'},{'mouse'},0,0)
%net_comparetsn(2,{'human_mapk_m2_ntsres'},{'human'},0,0)
%net_comparetsn(3,{'human_mapk_m3_ntsres'},{'human'},0,0)
%net_comparetsn(6,{'rat_mapk_m6_ntsres'},{'rat'},0,0)


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


function net_comparetsn(ChipRank,ListName,Species,CommonFlag,QLimitFlag,CluFlag,CluNbLimit)
global K

ClassName={'all','s','rest'};
if QLimitFlag
    TypeName={'continue','discrete'};
else
    TypeName={'continue'};
end
TypeNb=length(TypeName);
ChipNb=length(ChipRank);

for ChipL=1:ChipNb
    cd(sprintf('/home/mbellis/net/clum%u',ChipRank(ChipL)'))
    eval(sprintf('load %s',ListName{ChipL}))
    Nts{ChipL}=Res;
end

Index=[];
if CommonFlag
    %match data
    cd(K.dir.affyMetadata)
    load(sprintf('m%u_m%u_commonps',min(ChipRank),max(ChipRank)))
    Index={};
    [ComEnsRank,Idx1,Idx2]=intersect(Nts{1}.ensRank,Nts{2}.ensRank);


    Index=[];
    for PsL=1:length(Nts{1}.psRank)
        Pos1=find(ComPsRank(:,1)==Nts{1}.psRank(PsL));
        if ~isempty(Pos1)
            Pos2=find(Nts{2}.psRank==ComPsRank(Pos1,2));
            if ~isempty(Pos2)
                Index(end+1,1)=PsL;
                Index(end,2)=Pos2;
            end
        end
    end
else
    %construct Index
    if  isequal(Species{1},Species{2})
    else
    end
end


if ~isempty(Index)
    %re-plot found clusters in the same order that in m8
    Offset=0;
    LimitPos=4;
    FreqPs=cell(1,ClassNb);
    for ClassL=1:ClassNb
        h=figure;
        set(h,'name',sprintf('C %s',PathwayName{ClassL}))
        set(gcf,'color',[1,1,1])
        set(h,'name',sprintf('Tensornet clusters for %s',PathwayName{ClassL}));

        CurrC=NetFreq{8}(MeanSortOrder,MeanSortOrder);
        MeanC=CurrC(Offset+1:Offset+length(NrPsRank{ClassL}),Offset+1:Offset+length(NrPsRank{ClassL}));

        CurrC=NetFreq{LimitPos}(MeanSortOrder,MeanSortOrder);
        CurrC=CurrC(Offset+1:Offset+length(NrPsRank{ClassL}),Offset+1:Offset+length(NrPsRank{ClassL}));
        MeanC=[MeanC;CurrC*7];
        YNrTick=[0.5,round(length(NrPsRank{ClassL})/2),length(NrPsRank{ClassL}+0.5),...
            length(NrPsRank{ClassL})+round(length(NrPsRank{ClassL})/2)+0.5,2*length(NrPsRank{ClassL})+0.5];
        YLabel={'','mean(CORR) on all networks','',sprintf('Network reproducibility at CORR>%u',Limit(LimitPos)),''};
        FreqPs{ClassL}=zeros(1,size(CurrC,2));
        for NetL=1:length(NetNbs)
            for DensL=1:length(Densities{TypeL})
                CurrCluNb=length(SelPs{ClassL}{NetL,DensL});
                if CurrCluNb>0&CurrCluNb<100
                    Step=90/CurrCluNb;
                    if CurrCluNb>1
                        YNrTick=[YNrTick,YNrTick(end)+round(CurrCluNb/2),YNrTick(end)+CurrCluNb];
                    else
                        YNrTick=[YNrTick,YNrTick(end)+2,YNrTick(end)+3];
                    end
                    YLabeNetFreq{1,end+1}=sprintf('netnb %u dens %.2f',NetNbs(NetL),Densities{TypeL}(DensL));
                    YLabeNetFreq{1,end+1}='';
                    if CurrCluNb>1
                        for CluL=1:CurrCluNb;
                            CurrLine=zeros(1,size(MeanC,2));
                            CurrLine(SelPs{ClassL}{NetL,DensL}{CluL}-Offset)=100-ceil(Step*CluL);
                            MeanC=[MeanC;CurrLine];
                            FreqPs{ClassL}(find(CurrLine))=FreqPs{ClassL}(find(CurrLine))+1;
                        end
                    else
                        CurrLine=zeros(3,size(MeanC,2));
                        CurrLine(2,SelPs{ClassL}{NetL,DensL}{1}-Offset)=100-ceil(Step);
                        MeanC=[MeanC;CurrLine];
                        FreqPs{ClassL}(find(CurrLine(2,:)))=FreqPs{ClassL}(find(CurrLine(2,:)))+1;
                    end

                end
            end
        end
        CurrLine=zeros(3,size(MeanC,2));
        MeanC=[MeanC;CurrLine];
        FreqPs{ClassL}=round(FreqPs{ClassL}*100/max(FreqPs{ClassL}));
        MeanC=[MeanC;repmat(FreqPs{ClassL},3,1)];


        MeanC=MeanC(:,Index{2}{ClassL});

        SubMean=MeanC(1:length(NrPsRank{ClassL}),:);
        SubMean=SubMean(Index{2}{ClassL},:);
        MeanC(1:length(Index{2}{ClassL}),:)=SubMean;

        SubMean=MeanC(1+length(NrPsRank{ClassL}):2*length(NrPsRank{ClassL}),:);
        SubMean=SubMean(Index{2}{ClassL},:);
        MeanC(1+length(NrPsRank{ClassL}):length(NrPsRank{ClassL})+length(Index{2}{ClassL}),:)=SubMean;


        h=pcolor([MeanC,zeros(length(MeanC),1);zeros(1,size(MeanC,2)+1)]);
        set(h,'linestyle','none')
        set(gca,'ytick',YNrTick)
        set(gca,'tickdir','out')
        set(gca,'yticklabel',YLabel)
        title(PathwayName{ClassL})

        Offset=Offset+length(NrPsRank{ClassL});
    end

    %plot freqPs

    LimitInf=10;

    h=figure;
    set(h,'name','Frequency');
    for TypeL=1:TypeNb
        set(gcf,'color',[1,1,1])
        for ClassL=1:2
            subplot(2,2,(TypeL-1)*2+ClassL)
            hold on
            Freq1=Nts{1}.freqPs{TypeL}{ClassL}(Index(:,1));
            Freq2=Nts{2}.freqPs{TypeL}{ClassL}(Index(:,2));
            [Freq2 SortFreq]=sort(Freq2);
            Freq1=Freq1(SortFreq);
            [Freq1 SortFreq]=sort(Freq1);
            Freq2=Freq2(SortFreq);

            %plot(Freq1,Freq2,'+')
            plot(1:length(Freq1),Freq1,'b-')
            plot(1:length(Freq2),Freq2,'r-')
            set(gca,'box','on')
            if TypeL==1&ClassL==1
                legend(strrep(ListName{1},'_',' '),strrep(ListName{2},'_',' '))
            end
            xlabel(sprintf('ordered on %s',strrep(ListName{1},'_',' ')))
            title(sprintf('%s - %s',ClassName{ClassL},TypeName{TypeL}))
        end
    end

    h=figure;
    set(h,'name','Overlap')
    for TypeL=1:TypeNb
        set(gcf,'color',[1,1,1])
        for ClassL=1:2
            subplot(2,2,(TypeL-1)*2+ClassL)
            hold on
            Range=[0:10:100];
            Overlap=zeros(1,length(Range));
            for FreqL=1:length(Range)
                hold on
                Freq1=Nts{1}.freqPs{TypeL}{ClassL}(Index(:,1));
                Freq2=Nts{2}.freqPs{TypeL}{ClassL}(Index(:,2));
                Pos=find(Freq1>=Range(FreqL));
                Overlap(FreqL)=round(length(find(Freq2(Pos)>LimitInf))*100/length(Pos));
                %Overlap(FreqL)=round(length(find(Freq2(Pos)>Range(FreqL)))*100/length(Pos));
                %Overlap(FreqL)=round(length(find(Freq2(Pos)>Range(FreqL))));
                %Overlap(FreqL)=round(length(find(Freq2(Pos)>LimitInf)));
            end
            %plot(Range,Overlap,'+')
            plot(Range,Overlap)
            set(gca,'box','on')
            xlabel(sprintf('PsFreq in %s',strrep(ListName{1},'_',' ')));
            ylabel(sprintf('Overlap in %s',strrep(ListName{2},'_',' ')));
            title(sprintf('%s - %s',ClassName{ClassL},TypeName{TypeL}))
        end
    end


    %Probeset occurence in clusters
    MeanC=zeros(1,size(Index,1));
    YTick=[];
    YLabel={};
    for ClassL=1:3
        for TypeL=1:TypeNb
            if ~isempty(Nts{1}.freqPs{TypeL}{ClassL})
                if isempty(YTick)
                    YTick=2.5;
                else
                    YTick(1,end+1)=YTick(end)+1;
                end
                if TypeL==1
                    YLabel{1,end+1}=sprintf('%s - discrete - m%u',ClassName{ClassL},ChipRank(1));
                else
                    YLabel{1,end+1}=sprintf('%s - continuous - m%u',ClassName{ClassL},ChipRank(1));
                end
                MeanC=[MeanC;Nts{1}.freqPs{TypeL}{ClassL}(Index(:,1))];
            end
        end
    end
    for ClassL=1:3
        for TypeL=1:TypeNb
            if ~isempty(Nts{2}.freqPs{TypeL}{ClassL})
                if isempty(YTick)
                    YTick=2.5;
                else
                    YTick(1,end+1)=YTick(end)+1;
                end
                if TypeL==1
                    YLabel{1,end+1}=sprintf('%s - discrete - m%u',ClassName{ClassL},ChipRank(2));
                else
                    YLabel{1,end+1}=sprintf('%s - continuous - m%u',ClassName{ClassL},ChipRank(2));
                end
                MeanC=[MeanC;Nts{2}.freqPs{TypeL}{ClassL}(Index(:,2))];
            end
        end
    end
    MeanC=[MeanC;zeros(1,size(Index,1))];
    h=figure;
    set(h,'name','Probe set occurence in clusters')
    set(gcf,'color',[1,1,1])
    h=pcolor([MeanC,zeros(size(MeanC,1),1);zeros(1,size(MeanC,2)+1)]);
    set(h,'linestyle','none')
    set(gca,'ytick',YTick)
    set(gca,'yticklabel',YLabel)
    title('Probe set occurence in clusters')


    LimitInf1=10;
    LimitInf2=10;
    NtsPsRank={};
    NtsGene={};
    NtsDef={};
    NtsPsPos={};
    AllGene={};
    for TypeL=1:TypeNb
        for ClassL=1:3
            if ~isempty(Nts{1}.freqPs{TypeL}{ClassL}) & ~isempty(Nts{2}.freqPs{TypeL}{ClassL})
                ComPos=Index(find(Nts{1}.freqPs{TypeL}{ClassL}(Index(:,1))>LimitInf1&Nts{2}.freqPs{TypeL}{ClassL}(Index(:,2))>LimitInf2),1);
                NtsPsPos{TypeL}{ClassL}=ComPos;
                NtsPsRank{TypeL}{ClassL}=Nts{1}.psRank(ComPos);
                NtsGene{TypeL}{ClassL}=Nts{1}.gene(ComPos);
                NtsDef{TypeL}{ClassL}=Nts{1}.def(ComPos);
                [NtsGene{TypeL}{ClassL},GeneIndex,temp]=unique(NtsGene{TypeL}{ClassL});
                NtsDef{TypeL}{ClassL}=NtsDef{TypeL}{ClassL}(GeneIndex);
                AllGene=union(AllGene,NtsGene{TypeL}{ClassL});
            end
        end
    end
end

if isempty(Index) & ChipNb==2
    LimitInf1=10;
    LimitInf2=10;
    NtsGene={};
    NtsDef={};
    PsNb{1}=zeros(2,3);
    PsNb{2}=zeros(2,3);
    GeneNb{1}=zeros(2,3);
    GeneNb{2}=zeros(2,3);
    for TypeL=1:TypeNb
        for ClassL=1:3
            if ~isempty(Nts{1}.freqPs{TypeL}{ClassL}) & ~isempty(Nts{2}.freqPs{TypeL}{ClassL})
                ComPos1=find(Nts{1}.freqPs{TypeL}{ClassL}>LimitInf1);
                ComPos2=find(Nts{2}.freqPs{TypeL}{ClassL}>LimitInf2);
                PsNb{1}(TypeL,ClassL)=length(ComPos1);
                PsNb{2}(TypeL,ClassL)=length(ComPos2);
                GeneNb{1}(TypeL,ClassL)=length(unique(Nts{1}.gene(ComPos1)));
                GeneNb{2}(TypeL,ClassL)=length(unique(Nts{2}.gene(ComPos2)));
                NtsDef{TypeL}{ClassL}=[Nts{1}.def(ComPos1);Nts{2}.def(ComPos2)];
                NtsGene{TypeL}{ClassL}=[Nts{1}.gene(ComPos1);Nts{2}.gene(ComPos2)];
                InterGene=intersect(Nts{1}.gene(ComPos1),Nts{2}.gene(ComPos2));
                [temp,GenePos,temp]=intersect(NtsGene{TypeL}{ClassL},InterGene);
                NtsGene{TypeL}{ClassL}=NtsGene{TypeL}{ClassL}(GenePos);
                NtsDef{TypeL}{ClassL}=NtsDef{TypeL}{ClassL}(GenePos);
                %AllGene=union(AllGene,intersect(Nts{1}.gene(ComPos1),Nts{2}.gene(ComPos2)));
            end
        end
    end

end

if ChipNb==2
    AllGenes=union(Nts{1}.gene,Nts{2}.gene);
    Title={};
    Title{1}=sprintf('\t\t\t\tContinuous analysis');
    if TypeNb==2
        Title{1}=sprintf('%s\tDiscrete analysis',Title{1});
    end

    Title{2}=sprintf('GeneId\tGeneDef\tExist in m%u\tExist in m%u\tClustered in m%u\tClustered in m%u ',ChipRank(1),ChipRank(2),ChipRank(1),ChipRank(2));
    for i=3:5        
        Title{i}=sprintf('\t\t\t\t\t');
    end
    if TypeNb==2
        Title{2}=sprintf('%s\tClustered in m%u\tClustered in m%u ',Title{2},ChipRank(1),ChipRank(2));
        for i=3:5        
            Title{i}=sprintf('%s\t\t\t',Title{i});
        end
    end
    
    if CluFlag
        if TypeNb==1
            MatPos1=4;
        else
            MatPos1=6;
        end
        ZeroNb=2+TypeNb*2;
        CluNb1=0;
        for TypeL=1:TypeNb
            if TypeL==1
                Title{1}=sprintf('%s\t\tTSN clusters in continuous analysis',Title{1});
            else
                Title{1}=sprintf('%s\t\tTSN clusters in discrete analysis',Title{1});
            end
            for ChipL=1:2
                Title{2}=sprintf('%s\tm%u',Title{2},ChipRank(ChipL));
                for ClassL=1:3
                    Title{3}=sprintf('%s\tClass %s',Title{3},ClassName{ClassL});
                    for DensL=1:length(Nts{ChipL}.clu{TypeL}{ClassL})                        
                        CluNb=length(Nts{ChipL}.clu{TypeL}{ClassL}{1,DensL});                        
                        if CluNb<CluNbLimit
                            CluNb1=CluNb1+CluNb;
                            Title{4}=sprintf('%s\tDensity %u',Title{4},round(Nts{ChipL}.densities(DensL)*100));
                            ZeroNb=ZeroNb+CluNb;
                            for i=1:4
                                Title{i}=sprintf('%s%s',Title{i},repmat('\t',1,CluNb-1));
                            end
                            Title{2}=sprintf('%s\t',Title{2});
                            %Title{3}=sprintf('%s\t',Title{3});
                            for CluL=1:CluNb
                                Title{5}=sprintf('%s\tclu%u',Title{5},CluL);
                            end
                        end
                    end
                end
                if ChipL==1
                    MatPos2(1)=MatPos1;
                    MatPos2(2)=MatPos2(1)+ZeroNb;
                end
            end
        end
        GeneMat=zeros(length(AllGenes),ZeroNb);
    else
        GeneMat=zeros(length(AllGenes),2+TypeNb*2);
    end
    AllDef=cell(length(AllGenes),1);

    for GeneL=1:length(AllGenes)
        for TypeL=1:TypeNb
            for ChipL=1:2
                GenePos=strmatch(AllGenes{GeneL},Nts{ChipL}.gene,'exact');
                if  ~isempty(GenePos)
                    GeneMat(GeneL,ChipL)=1;
                    if isempty(AllDef{GeneL})
                        AllDef{GeneL}=Nts{ChipL}.def{GenePos(1)};
                    end
                end

                %             GenePos=strmatch(AllGenes{GeneL},Nts{2}.gene,'exact');
                %             if  ~isempty(GenePos)
                %                 GeneMat(GeneL,2)=1;
                %                 if isempty(AllDef{GeneL})
                %                     AllDef{GeneL}=Nts{2}.def{GenePos(1)};
                %                 end
                %             end
                %             for ChipL=1:2
                try
                    %CurrGenePos=strmatch(AllGenes{GeneL},NtsGene{TypeL}{ChipL},'exact');
                    if ~isempty(GenePos)
                        GeneMat(GeneL,(TypeL-1)*2+2+ChipL)=1;
                        MatPos=MatPos2(ChipL);
                        if CluFlag
                            for ClassL=1:3
                                for DensL=1:length(Nts{ChipL}.clu{TypeL}{ClassL})
                                    CluNb=length(Nts{ChipL}.clu{TypeL}{ClassL}{1,DensL});
                                    if CluNb<CluNbLimit                                        
                                        for CluL=1:CluNb
                                            MatPos=MatPos+1;
                                            
                                            for i=1:length(GenePos)
                                            if ~isempty(find(Nts{ChipL}.clu{TypeL}{ClassL}{1,DensL}{CluL}==GenePos(i)))
                                                GeneMat(GeneL,MatPos)=1;
                                            end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                catch
                end
            end
        end
    end

    cd(sprintf('/home/mbellis/net/clum%u',ChipRank(1)'))
    fid=fopen(sprintf('gene_%s_%s.txt',ListName{1},ListName{2}),'w');    
    for i=1:5
        fprintf(fid,[Title{i},'\n']);
    end
    for GeneL=1:length(AllGenes)
        fprintf(fid,['%s\t%s',repmat('\t%u',1,size(GeneMat,2)),'\n'],AllGenes{GeneL},AllDef{GeneL},GeneMat(GeneL,:));
    end
    fclose(fid)

else
    LimitInf=10;
    
    NtsGene={};
    PsNb=zeros(2,3);
    GeneNb=zeros(2,3);
    for TypeL=1:TypeNb
        for ClassL=1:3
            if ~isempty(Nts{1}.freqPs{TypeL}{ClassL})
              ComPos1=find(Nts{1}.freqPs{TypeL}{ClassL}>LimitInf);               
                PsNb(TypeL,ClassL)=length(ComPos1);
                GeneNb(TypeL,ClassL)=length(unique(Nts{1}.gene(ComPos1)));
                NtsGene{TypeL}{ClassL}=Nts{1}.gene(ComPos1);
                [NtsGene{TypeL}{ClassL},temp,temp]=unique(NtsGene{TypeL}{ClassL});
            end
        end
    end

    
    [AllGenes,GenePos,temp]=unique(Nts{1}.gene);
    AllDef=Nts{1}.def(GenePos);
    GeneMat=zeros(length(AllGenes),4);             

    for GeneL=1:length(AllGenes)
        for TypeL=1:TypeNb
            try
                if ~isempty(strmatch(AllGenes{GeneL},NtsGene{TypeL}{1},'exact'))
                    GeneMat(GeneL,(TypeL-1)*2+1)=1;
                end
            catch
            end
            try
                if ~isempty(strmatch(AllGenes{GeneL},NtsGene{TypeL}{2},'exact'))
                    GeneMat(GeneL,(TypeL-1)*2+2)=1;
                end
            catch
            end
        end
    end
          
    
    cd(sprintf('/home/mbellis/net/clum%u',ChipRank(1)'))    
    fid=fopen(sprintf('gene_%s.txt',ListName{1}),'w');


    for GeneL=1:length(AllGenes)
        fprintf(fid,'%s\t%s\t%u\t%u\t%u\t%u\n',AllGenes{GeneL},AllDef{GeneL},GeneMat(GeneL,1),...
            GeneMat(GeneL,2),GeneMat(GeneL,3),GeneMat(GeneL,4));
    end
    fclose(fid)

end






