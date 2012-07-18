cd('/home/mbellis/array2/mbellis/sosma/arraymatic/amcdata/metadata/affy')
load m11_psassignation
cd('/home/mbellis/sosma/raydataraw/2010_JULIEN')
load fdr_zvar_mu_vs_wt_julien_2010_commonaffy
[TestCliGsea,GsClasses,CliGsNames]=gsea_findgeneset(11,OperonZVar,OperonFdr,OperonMeanRanks,0,AbsentIndex,'m11_cliques',K.dir.affyMetadata,'m11_cliques_type1_r5489_c2000',K.dir.affyMetadata,TairId,0,0,'cliques');
[TestPathGsea,GsClasses,PathGsNames]=gsea_findgeneset(11,OperonZVar,OperonFdr,OperonMeanRanks,0,AbsentIndex,'m11_athgsdb11',K.dir.affyMetadata,'athgsdb11_type1_r5489_c2000',K.dir.affyMetadata,TairId,0,0,'athgsdb11');
GeneName='WOX8';
GsName='Cliques';

Gsea=CliGsea;
TestGsea=TestCliGsea;
RoundNb=1;
BiolNb=30;
CompNb=RoundNb*BiolNb+1;
ClassNb=length(GsClasses);

%for each class find the maximal gene set rank
MaxGeneSetRank=zeros(ClassNb,1);
for ClassL=1:ClassNb
    for RoundL=1:RoundNb
        for CompL=1:BiolNb
            if ~isempty(Gsea{RoundL}{CompL}{ClassL})
                MaxGeneSetRank(ClassL)=max(MaxGeneSetRank(ClassL),max(Gsea{RoundL}{CompL}{ClassL}));
            end
            if ~isempty(TestGsea{ClassL})
                MaxGeneSetRank(ClassL)=max(MaxGeneSetRank(ClassL),max(TestGsea{ClassL}));
            end
        end
    end
end
%construct and fill matrix
Gs=cell(ClassNb,1);
GsRank=cell(ClassNb,1);
if RoundNb==1
    Colors=[1,0,0];
else
    Colors=colors(colormap,RoundNb);
end
PosClassNb=length(find(MaxGeneSetRank>0));
LineNb=floor(sqrt(PosClassNb));
ColNb=round(PosClassNb/LineNb);
if ColNb*LineNb<PosClassNb
    LineNb=LineNb+1;
end

h=figure;
set(gcf,'color',[1,1,1])
set(h,'name',sprintf('%s %s',GeneName,GsName))
PlotPos=0;
GsStat=cell(ClassNb,1);
for ClassL=1:ClassNb
    if MaxGeneSetRank(ClassL)>0
        PlotPos=PlotPos+1;
        subplot(LineNb,ColNb,PlotPos)
        Gs{ClassL}=zeros(CompNb,MaxGeneSetRank(ClassL));
        Pos=0;
        for RoundL=1:RoundNb
            for CompL=1:BiolNb
                Pos=Pos+1;
                if ~isempty(Gsea{RoundL}{CompL}{ClassL})
                    Gs{ClassL}(Pos,Gsea{RoundL}{CompL}{ClassL})=1;
                end
            end
            if ~isempty(TestGsea{ClassL})
                Gs{ClassL}(CompNb,TestGsea{ClassL})=1;
            end
        end
        
        %find not used gene sets
        NullPos=find(sum(Gs{ClassL})==0);
        GsRank{ClassL}=find(sum(Gs{ClassL})>0);
        GsStat{ClassL}.gsRank=GsRank{ClassL};        
        Gs{ClassL}(:,NullPos)=[];
        %order results
        Repeated=sum(Gs{ClassL}(1:CompNb-1,:));
        [Temp,SortIndex]=sort(Repeated);
        [Temp,ReverseIndex]=sort(SortIndex);
        Gs{ClassL}=Gs{ClassL}(:,SortIndex);
        GsRank{ClassL}=GsRank{ClassL}(SortIndex);
        NullNb=length(find(Repeated==0));
        PosNb=length(find(Repeated>0));
        
        hold on
        Repeated=[];
        for RoundL=1:RoundNb
            Repeated=[Repeated;sum(Gs{ClassL}((RoundL-1)*BiolNb+1:RoundL*BiolNb,:))];            
        end
        Frequency=sum(Repeated,1)*100/(RoundNb*BiolNb);
        plot(Frequency,'.','color',Colors(RoundL,:))
        GsStat{ClassL}.gsFrequency=Frequency(ReverseIndex);
        Pos=Gs{ClassL}(CompNb,:)>0;
        GsStat{ClassL}.inTest=Pos(ReverseIndex);
        %plot(cumsum(ones(1,length(Repeated)))./length(Repeated)*100,'k')
        Pos=find(Gs{ClassL}(CompNb,:));
        title(sprintf('%s (%u(%u%%)+%u(%u%%))',GsClasses{ClassL},NullNb,max(0,round(NullNb*100/length(Pos))),length(Pos)-NullNb,max(0,round((length(Pos)-NullNb)*100/(size(Repeated,2)-NullNb)))))
        set(gca,'box','on')
        plot(Pos,(cumsum(ones(1,length(Pos)))./length(Pos))*100,'k')
        plot(Pos,(cumsum(ones(1,length(Pos)))./length(Pos))*100,'g+')
        set(gca,'ylim',[0,100])
        set(gca,'xlim',[0,length(Repeated)])        
    end
end


%print results
cd('/home/mbellis/sosma/raydataraw/2010_JULIEN')
load gsea_comp
GeneNames={'WRKY2','WOX8'};
for GeneL=1:length(GeneNames)
    for TypeL=1:2
        if TypeL==1
            eval(sprintf('Comp=%s_cliques',GeneNames{GeneL}));
            GsNames=CliGsNames;
            fid=fopen(sprintf('%s_gsea_cliques_%s.csv',GeneNames{GeneL},date),'w');
        else
            eval(sprintf('Comp=%s_pathways',GeneNames{GeneL}));
            GsNames=PathGsNames;
            fid=fopen(sprintf('%s_gsea_pathways_%s.csv',GeneNames{GeneL},date),'w');
        end        
        for ClassL=1:ClassNb
            if ~isempty(Comp{ClassL})
                fprintf(fid,'%s\tFreq\tInTest\n',GsClasses{ClassL})
                for GsL=1:length(Comp{ClassL}.gsRank)
                    fprintf(fid,'%s\t%u\t%u\n',GsNames{Comp{ClassL}.gsRank(GsL)},ceil(Comp{ClassL}.gsFrequency(GsL)),Comp{ClassL}.inTest(GsL));
                end
            end
        end
        fclose(fid)
    end
end
