%TRS_DISPLAYRANK allows to dispklay rank or signal of several probe sets


%INPUT PARAMETERS
%1 GeneList: list of genes
%2 RankFlag: display rank (if =1) or log2(signal) (if =0)
%3 OverlapFlag: display all ps on a single figure

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


%trs_displayrank('neuro_gene2psrank',1,0)
%trs_displayrank('neuro_gene2psrank',1,0)
%trs_displayrank('th_1_gene2psrank',1,1)
%trs_displayrank('crel_1_gene2psrank',1,1)
%trs_displayrank('crel_evi_gene2psrank',1,1)
%trs_displayrank('crel_areb6_gene2psrank',1,1)
%trs_displayrank('crel_nfk_gene2psrank',1,1)
%trs_displayrank('crel_2_gene2psrank',1,1)

function trs_displayrank(GeneList,RankFlag,OverlapFlag)
global P K

if P.flag.loadData==1
    DataRanks=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le');
end

cd(K.dir.chip)
[Ps,Ens,ChipGene]=textread(sprintf('m%u_gene.txt',P.chip.chipRank),'%s%s%s','delimiter','\t');

cd(K.dir.list)
load(GeneList)

Colors=colors(colormap,P.biol.nb);
if ~isempty(GenePos{P.chip.chipRank})
    CurrGenePos=GenePos{P.chip.chipRank};
    %recover gene or ensembl id of interrogated probe sets
    Legend={};
    PsRank=[];
    for GeneL=1:length(Gene)
        try
        if ~isempty(CurrGenePos{GeneL})
            for PsL=1:length(CurrGenePos{GeneL})
                CurrPsRank=CurrGenePos{GeneL}(PsL);
                PsRank(end+1,1)=CurrPsRank;
                Legend{end+1}=Gene{GeneL};
            end                
        end
        catch
        end
    end
    %reorder point according to biolrank and  exprank
    [temp BiolOrder]=sort(P.point.biolRank);
    [temp ExpOrder]=sort(P.point.expRank(BiolOrder));
    PointOrder=BiolOrder(ExpOrder);
    if OverlapFlag
        h=figure;
        set(gcf,'color',[1,1,1])
        set(h,'name',sprintf('all Ps of %s',GeneList))
        hold on
        PsColors=colors(colormap,length(PsRank));      
        hs=[];
    end
    for PsL=1:length(PsRank)
        MinVal=100;
        MaxVal=0;
        CurrPs=PsRank(PsL);
        if OverlapFlag==0
            h=figure;
            set(gcf,'color',[1,1,1])
            set(h,'name',sprintf('Ps %u (%s): %s %s',CurrPs,strrep(Ps{CurrPs},'_','\_'),Ens{CurrPs},ChipGene{CurrPs}))
            hold on
        end
        YTick=1;
        BiolRank=P.point.biolRank(PointOrder(1));
        ExpRank=P.point.expRank(PointOrder(1));
        YTickLabel={};
        YTickLabel{1}=[P.exp.name{ExpRank},'-',P.biol.name{BiolRank}];
        h=line([0,100],[1,1]);
        set(h,'color',Colors(1,:))
        set(h,'linestyle',':')
        CurrPs=PsRank(PsL);
        Start=1;
        Val=[];
        MemVal=0;
        for PointL=1:P.point.nb
            CurrExpRank=P.point.expRank(PointOrder(PointL));
            CurrBiolRank=P.point.biolRank(PointOrder(PointL));                          
            if CurrBiolRank~=BiolRank
                BiolRank=CurrBiolRank;
                YTick(end+1)=PointL;
                if CurrExpRank==ExpRank
                    YTickLabel{end+1}=P.biol.name{BiolRank};
                    h=line([0,100],[PointL,PointL]);
                set(h,'color',Colors(P.point.biolRank(PointOrder(PointL)),:))
                set(h,'linestyle',':')
                else
                    ExpRank=CurrExpRank;
                    YTickLabel{end+1}=[P.exp.name{ExpRank},'-',P.biol.name{BiolRank}];
                    h=line([0,100],[PointL,PointL]);
                    set(h,'color',Colors(P.point.biolRank(PointOrder(PointL)),:))
                    set(h,'linestyle','-')
                end
                MeanVal=mean(Val);                
                Val=[];
                h=line([MemVal,MeanVal],[Start,Start]);
                if OverlapFlag
                    set(h,'color',PsColors(PsL,:))                   
                else
                    set(h,'color',[0,0,0])
                end
                set(h,'linestyle','-')                
                if OverlapFlag                                   
                    if BiolRank==P.point.biolRank(PointOrder(P.point.nb))
                        hs(end+1)=h;
                    end
                end
                h=line([MeanVal,MeanVal],[Start,PointL]);
                if OverlapFlag
                    set(h,'color',PsColors(PsL,:))
                    if P.point.biolRank(PointOrder(PointL))==P.point.biolRank(PointOrder(2))
                        hs(end+1)=h;
                    end
                else
                    set(h,'color',[0,0,0])
                end                
                set(h,'linestyle','-')
                MemVal=MeanVal;
                Start=PointL;
            end
            if RankFlag
                CurrVal=DataRanks(CurrPs,PointOrder(PointL));
            else
                CurrVal=log2(interp1(P.chip.refRank,P.chip.refSignal,DataRanks(CurrPs,PointOrder(PointL))));            
                if CurrVal>MaxVal
                    MaxVal=CurrVal;
                end
                if CurrVal<MinVal
                    MinVal=CurrVal;
                end
            end
            if OverlapFlag==0
                plot(CurrVal,PointL,'.','color',Colors(P.point.biolRank(PointOrder(PointL)),:));
            end
            Val(end+1)=CurrVal;
        end
        MeanVal=mean(Val);
        h=line([MemVal,MeanVal],[Start,Start]);
        if OverlapFlag
            set(h,'color',PsColors(PsL,:))        
        else
            set(h,'color',[0,0,0])
        end
        set(h,'linestyle','-')
        h=line([MeanVal,MeanVal],[Start,PointL]);
        if OverlapFlag
            set(h,'color',PsColors(PsL,:))
        else
            set(h,'color',[0,0,0])
        end        
        set(h,'linestyle','-')
        if OverlapFlag==0|PsL==length(PsRank)
            set(gca,'ytick',YTick)
            set(gca,'yticklabel',YTickLabel)
            set(gca,'tickdir','out')
            set(gca,'box','on')
            set(gca,'ticklength',[0.002,0.025])
            set(gcf,'position',[319,40,600,1770])
            set(gca,'position',[0.42,0.0189,0.485,0.966])

            if RankFlag
                xlabel('rank(signal)')
            else
                xlabel('Log2(Signal)')
                set(gca,'xlim',[floor(MinVal),ceil(MaxVal)])
                set(gca,'xtick',[floor(MinVal):ceil(MaxVal)])
            end
        end
        if OverlapFlag==0
            title(sprintf('Ps%u (%s): %s - %s',CurrPs,strrep(Ps{CurrPs},'_','-'),Ens{CurrPs},strrep(ChipGene{CurrPs},'_','-')))
        end
    end
    if OverlapFlag
        legend(hs,Legend)
    end
end