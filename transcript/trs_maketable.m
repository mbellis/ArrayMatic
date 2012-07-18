%TRS_MAKETABLE is used in large analysis to read bloc of data and construct matlab table containing the fdr values of comparisons.
%Two types of tables are constructed :
%FIRST TYPE: table with systematic comparisons between selected biological conditions
%            the columnar position of a given comparison between two biological conditions is determined
%            in order to select subsets of comparisons between groups of
%            biologicla conditions
%SECOND TYPE: table with sytematic comparisons between non selected intra-experiments biological conditions
%             the columnar position of a given comparison is undetermined,
%             because the table is used as a whole (no selection of
%             particular columns)

%INPUT PARAMETERS
%ModelRank: chip model rank
%BiolIndexRank: allows to recover index of used biological conditions
%               position in P.net.biolIndex{BiolIndexRank}
%TableStart: the first table to be constructed 
%TableEnd: the last table to be constructed
%          TableStart and TableEnd allow to run several instances of the script in parallel 
%DisplayFlag (0,1): indicates if figures must be displayed 
%                   if DisplayFlag==1, only the figure are displayed
%                   (tables are not created)

%EXTERNAL FILE
%blocX.float32le files containing the results of comparisons
%blocX.mat files containing the description of comparaisons and their location in blocX.float32le

%OUTPUT PARAMETERS
%none

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

function trs_maketable(ModelRank,BiolIndexRank,TableStart,TableEnd,DisplayFlag)
global P
if DisplayFlag==0
    TABLE_SIZE=1000;
    BiolNb=length(find(P.net.biolIndex{BiolIndexRank}));
    TablePos=uint32(zeros(BiolNb));
    Pos=0;
    for BiolL1=1:BiolNb-1
        for BiolL2=BiolL1+1:BiolNb
            Pos=Pos+1;
            TablePos(BiolL1,BiolL2)=Pos;
        end
    end
    CompNb1=Pos;
    TableNb=ceil(P.chip.currProbeSetNb/TABLE_SIZE);
    BlocNb=ceil(P.biol.nb/100);

    %calculate size of first type tables
    cd(P.dir.data)
    BiolList=find(P.net.biolIndex{BiolIndexRank});
    CompNb2=0;
    for BlocL1=1:BlocNb
        for BlocL2=BlocL1:BlocNb
            BlocRank=P.net.blocGrid(BlocL1,BlocL2);
            load(sprintf('bloc%u',BlocRank))
            BiolIndex1=find(BiolList>(BlocL1-1)*100&BiolList<=BlocL1*100);
            BiolList1=BiolList(BiolIndex1);
            BiolIndex2=find(BiolList>(BlocL2-1)*100&BiolList<=BlocL2*100);
            BiolList2=BiolList(BiolIndex2);
            if BlocL1==BlocL2
                CurrCompNb1=length(BiolList1)*(length(BiolList1)-1)/2;
            else
                CurrCompNb1=length(BiolList1)*length(BiolList2);
            end
            CompNb2=CompNb2+size(BlocIndex,1)-CurrCompNb1;
        end
    end
    CurrTableNb=TableEnd-TableStart+1;
    t=0;
    tic
    for TableL=TableStart:TableEnd
        if TableL<TableNb
            Table1=single(zeros(TABLE_SIZE,CompNb1));
            Table2=single(zeros(TABLE_SIZE,CompNb2));
            Range=(TableL-1)*TABLE_SIZE+1:TableL*TABLE_SIZE;
        else
            Table1=single(zeros(mod(P.chip.currProbeSetNb,TABLE_SIZE),CompNb1));
            Table2=single(zeros(mod(P.chip.currProbeSetNb,TABLE_SIZE),CompNb2));
            Range=(TableL-1)*TABLE_SIZE+1:P.chip.currProbeSetNb;
        end
        StartTable2=1;
        for BlocL1=1:BlocNb
            BiolIndex1=find(BiolList>(BlocL1-1)*100&BiolList<=BlocL1*100);
            BiolList1=BiolList(BiolIndex1);
            for BlocL2=BlocL1:BlocNb
                %BlocL2
                BiolIndex2=find(BiolList>(BlocL2-1)*100&BiolList<=BlocL2*100);
                BiolList2=BiolList(BiolIndex2);

                BlocRank=P.net.blocGrid(BlocL1,BlocL2);
                load(sprintf('bloc%u',BlocRank))

                UsedPos=[];
                CurrTablePos1=[];
                CurrCompPos1=[];
                for BiolL1=1:length(BiolList1)
                    Biol1=BiolList1(BiolL1);
                    for BiolL2=1:length(BiolList2)
                        Biol2=BiolList2(BiolL2);
                        if Biol1~=Biol2
                            CompPos=find(BlocIndex(:,1)==Biol1&BlocIndex(:,2)==Biol2);
                            if isempty(CompPos)
                                h=errordlg(sprintf('Bloc1:%u, Bloc2:%u, misses Biol%u-Biol%u',BlocL1,BlocL2,Biol1,Biol2));
                                waitfor(h)
                                error('process canceled')
                            end
                            UsedPos(end+1,1)=CompPos;
                            %read data
                            BiolPos1=find(BiolList==Biol1);
                            BiolPos2=find(BiolList==Biol2);
                            CurrTablePos1=[CurrTablePos1,TablePos(min(BiolPos1,BiolPos2),max(BiolPos1,BiolPos2))];
                            CurrCompPos1=[CurrCompPos1,CompPos];
                        end
                    end
                end
                if ~isempty(CurrCompPos1)
                    Table1(:,CurrTablePos1)=load_data(sprintf('bloc%u.float32le',BlocRank),P.dir.data,P.chip.currProbeSetNb,size(BlocIndex,1),'single','ieee-le',Range,CurrCompPos1);
                end
                CurrCompPos2=1:size(BlocIndex,1);
                %eliminate position already processed
                CurrCompPos2(CurrCompPos1)=[];
                if ~isempty(CurrCompPos2)
                    CurrTablePos2=StartTable2:StartTable2+length(CurrCompPos2)-1;
                    Table2(:,CurrTablePos2)=load_data(sprintf('bloc%u.float32le',BlocRank),P.dir.data,P.chip.currProbeSetNb,size(BlocIndex,1),'single','ieee-le',Range,CurrCompPos2);
                    StartTable2=StartTable2+length(CurrCompPos2);
                end
            end
        end
        eval(sprintf('save m%u_%u_1 Table1',ModelRank,TableL))
        eval(sprintf('save m%u_%u_2 Table2',ModelRank,TableL))
        t=t+toc;
        sprintf('%utables processed upon %u in %uH: estimate that it stay %uH',TableL-TableStart+1,CurrTableNb,round(t/3600),round(((TableEnd-TableL)*t)/((TableL-TableStart+1)*3600)))
        tic
    end

else
    Var{1}=single(zeros(BiolNb));
    Var{2}=single(zeros(BiolNb));
    Var{3}=single(zeros(BiolNb));
    for BlocL1=1:BlocNb
        for BlocL2=BlocL1:BlocNb
            BlocRank=P.net.blocGrid(BlocL1,BlocL2);
            load(sprintf('bloc%u',BlocRank))
            BiolIndex1=find(BiolList>(BlocL1-1)*100&BiolList<=BlocL1*100);
            BiolList1=BiolList(BiolIndex1);
            BiolIndex2=find(BiolList>(BlocL2-1)*100&BiolList<=BlocL2*100);
            BiolList2=BiolList(BiolIndex2);
            for BiolL1=1:length(BiolList1)
                Biol1=BiolList1(BiolL1);
                for BiolL2=1:length(BiolList2)
                    Biol2=BiolList2(BiolL2);
                    if Biol1~=Biol2
                        CompPos=find(BlocIndex(:,1)==Biol1&BlocIndex(:,2)==Biol2);

                        if isempty(CompPos)
                            h=errordlg(sprintf('Bloc1:%u, Bloc2:%u, misses Biol%u-Biol%u',BlocL1,BlocL2,Biol1,Biol2));
                            waitfor(h)
                            error('process canceled')
                        end
                        BiolPos1=find(BiolList==Biol1);
                        BiolPos2=find(BiolList==Biol2);
                        for i=1:3
                            Var{i}(min(BiolPos1,BiolPos2),max(BiolPos1,BiolPos2))=BlocIndex(CompPos,2+(i-1)*2+1);
                            Var{i}(max(BiolPos1,BiolPos2),min(BiolPos1,BiolPos2))=BlocIndex(CompPos,2+(i-1)*2+2);
                        end
                    end
                end
            end
        end
    end
    h=figure;
    set(gcf,'color',[1,1,1])
    image(Var{1}/100);
    title('FDR 0.001')

    h=figure;
    set(gcf,'color',[1,1,1])
    image(Var{2}/100);
    title('FDR 0.01')

    h=figure;
    set(gcf,'color',[1,1,1])
    image(Var{3}/100);
    title('FDR 0.1')

    h=figure;
    set(gcf,'color',[1,1,1])
    subplot(2,3,1)
    Val=triu(Var{1});
    Val=Val(:);
    Val(Val==0)=[];
    hist(Val,100);
    title('INC FDR 0.001')
    set(gcf,'color',[1,1,1])
    subplot(2,3,4)
    Val=tril(Var{1});
    Val=Val(:);
    Val(Val==0)=[];
    hist(Val,100);
    title('DEC FDR 0.001')

    subplot(2,3,2)
    Val=triu(Var{2});
    Val=Val(:);
    Val(Val==0)=[];
    hist(Val,100);
    title('INC FDR 0.01')
    set(gcf,'color',[1,1,1])
    subplot(2,3,5)
    Val=tril(Var{2});
    Val=Val(:);
    Val(Val==0)=[];
    hist(Val,100);
    title('DEC FDR 0.01')

    subplot(2,3,3)
    Val=triu(Var{3});
    Val=Val(:);
    Val(Val==0)=[];
    hist(Val,100);
    title('INC FDR 0.1')
    set(gcf,'color',[1,1,1])
    subplot(2,3,6)
    Val=tril(Var{3});
    Val=Val(:);
    Val(Val==0)=[];
    hist(Val,100);
    title('DEC FDR 0.1')

    h=figure;
    set(gcf,'color',[1,1,1])
    for i=1:3
        subplot(1,3,i)
        plot(triu(Var{i}),triu(Var{i}'),'b.','markersize',3)
        set(gca,'box','on')
        xlabel('nb of increased')
        ylabel('nb of decreased')
        switch i
            case 1
                title('FDR 0.001')
            case 2
                title('FDR 0.01')
            case 3
                title('FDR 0.1')
        end
    end
end



