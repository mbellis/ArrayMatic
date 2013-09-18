%TRS_WRITE RESULTS allows to write comparisons result


%INPUT PARAMETERS
%1 GeneList: list of genes
%2 RankFlag: display rank (if =1) or log2(signal) (if =0)

%EXTERNAL FILES

%OUTPUT PARAMETERS


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


%trs_writeresults([],'neuro_gene2psrank')
%trs_writeresults([1:14],'neuro_gene2psrank')

function trs_writeresults(varargin)
global K M P

%load comparisons descriptions
cd(P.dir.data)
load Comp
if nargin>=1
    CompRank=varargin{1};
else
    [CompRank,Ok]=listdlg('liststring',M{1}.compName,'listsize',[400,1200]);
end

%load eventually a gene list
cd(K.dir.list)
ListFlag=0;
if nargin>=2
    GeneList=varargin{2};        
    load(GeneList)
    ListFlag=1;
else
    ListFlag=questdlg('do you want ot restrict results to a gene list ?','','yes','no','yes');
    if isequal(ListFlag,'yes')
        ListFlag=1;
        [GeneList,PathName]=uigetfile('*.mat');
        cd(PathName)
        load(GeneList)
    end
end
if ListFlag
    if ~isempty(GenePos{P.chip.chipRank})
        PsRank=[];
        CurrGenePos=GenePos{P.chip.chipRank};
        %recover gene or ensembl id of interrogated probe sets
        for GeneL=1:length(Gene)
            try
                if ~isempty(CurrGenePos{GeneL})
                    for PsL=1:length(CurrGenePos{GeneL})
                        PsRank(end+1,1)=CurrGenePos{GeneL}(PsL);
                    end
                end
            catch
            end
        end
    else
        PsRank=[];
    end
else
    PsRank=1:P.chip.currProbeSetNb;
end

if ~isempty(PsRank)
    %load probe set and gene names
    cd(K.dir.chip)
    [Ps,Ens,ChipGene]=textread(sprintf('m%u_gene.txt',P.chip.chipRank),'%s%s%s','delimiter','\t');

    Fdr=load_data('Fdr_02.float32le',P.dir.data,P.chip.currProbeSetNb,P.chip.currProbeSetNb,'single','ieee-le',PsRank,CompRank);
    ZVar=load_data('ZVar_02.float32le',P.dir.data,P.chip.currProbeSetNb,P.chip.currProbeSetNb,'single','ieee-le',PsRank,CompRank);
    Sensitivity=load_data('Sensitivity_02.float32le',P.dir.data,P.chip.currProbeSetNb,P.chip.currProbeSetNb,'single','ieee-le',PsRank,CompRank);
    FC=load_data('Fc_02.float32le',P.dir.data,P.chip.currProbeSetNb,P.chip.currProbeSetNb,'single','ieee-le',PsRank,CompRank);
    load TotalVar_02
    TotalVar.inc=TotalVar.inc(CompRank);
    TotalVar.dec=TotalVar.dec(CompRank);
    cd(P.dir.results)
    %fid=fopen(sprintf('aorteHP_%s_results_%s.txt',lower(TYPE),date),'w');
    fid=fopen(sprintf('DRG_crush_results_%s.txt',date),'w');

    CompNames=M{1}.compName(CompRank);

    Header='PsId\tGeneId\tGene';
    CompNb=length(CompNames);
    for i=1:CompNb
        Header=sprintf('%s\t%s\t\t\t',Header,CompNames{i});
    end
    fprintf(fid,'%s\n',Header);

    Header='total inc:';
    for CompL=1:CompNb
        Header=sprintf('%s\t%u\t\t\t',Header,round(TotalVar.inc(CompL)));
    end
    fprintf(fid,'%s\n',Header);

    Header='total dec:';
    for CompL=1:CompNb
        Header=sprintf('%s\t%u\t\t\t',Header,round(TotalVar.dec(CompL)));
    end
    fprintf(fid,'%s\n',Header);

    Header='probe set';
    for i=1:CompNb
        Header=sprintf('%s\tZVar\tFDR\tSens\tFC',Header);
    end
    fprintf(fid,'%s\n',Header);


    for PsL=1:length(PsRank)
        CurrPsRank=PsRank(PsL);
        fprintf(fid,'%s\t%s\t%s',Ps{CurrPsRank},Ens{CurrPsRank},ChipGene{CurrPsRank});
        for CompL=1:CompNb
            fprintf(fid,'\t%.3f\t%.3f',ZVar(PsL,CompL),Fdr(PsL,CompL),Sensitivity(PsL,CompL),FC(PsL,CompL));
        end
        fprintf(fid,'\n');
    end
    fclose(fid)
    
    if ListFlag
        PsRanks=PsRank(1:16);
        Range=[1:7];
        XRange=[1:7];
        Legend={};
        ColNb=floor(sqrt(length(PsRanks)));
        RowNb=ceil(length(PsRanks)/ColNb);
        if ColNb*RowNb<length(PsRanks)
            RowNb=RowNb+1;
        end              
        
        h=figure
        set(gcf,'color',[1,1,1])
        Colors=colors(colormap,length(PsRanks));
        Legend={};
        ColNb=floor(sqrt(length(PsRanks)));
        RowNb=ceil(length(PsRanks)/ColNb);
        if ColNb*RowNb<length(PsRanks)
            RowNb=RowNb+1;
        end
        
        for PsL=1:length(PsRanks)
            subplot(RowNb,ColNb,PsL)
            CurrPsRank=PsRanks(PsL);            
            plot(XRange,ZVar(PsL,Range),':','color',[0,0,0])
            hold on
            plot(XRange,ZVar(PsL,Range),'.','color',[0,0,0])
            Pos=find(abs(Fdr(PsL,Range))<=0.10);
            if~isempty(Pos)
                plot(XRange(Pos),ZVar(PsL,Range(Pos)),'o','color',[0,0,0])
            end
            
            plot(XRange,ZVar(PsL,Range),':','color',Colors(PsL,:),'linewidth',2)
            plot(XRange,ZVar(PsL,Range),'.','color',Colors(PsL,:))
            Pos=find(abs(Fdr(PsL,Range))<=0.10);            
            if~isempty(Pos)
                plot(XRange(Pos),ZVar(PsL,Range(Pos)),'o','color',Colors(PsL,:))
            end
            title([ChipGene{CurrPsRank},' ',strrep(Ps{CurrPsRank},'_',' ')]);
            set(gca,'xtick',[1,4,7])
            set(gca,'xticklabel',{'1h','12h','24h'})
            set(gca,'ylim',[-15,25])
        end
    end
end



        