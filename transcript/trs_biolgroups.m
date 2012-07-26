%TRS_BIOLGROUPS allows to select group biological conditions.
%CVM are constructed by using comparisons between or inside groups

%INPUT PARAMETERS
%       GrpType: (used if NotRedundantFlag=0)
%                'split'
%                'tf'
% RedundantFlag: indicate if not redundant biological conditions mus be searched before
%                constituing the groups
%      varargin: if GrpType = 'tf'
%                Fdr: fdr values (one or two) used for doing network   
%                Sensitivity: sensitivity values (one or two) used for doing network
%                liste of tf name (ex: {'WRKY2','WOX8'})
%                liste of tf id (ex: {'At5g56270','At5g45980'})
%                liste of tf gene name (ex: {'WRKY2','WOX8'})
%                liste of corresponding probe set pos {[3802,3809];[12456]}
%                NonRedundantFlag: if =1 restrict choice to NonRedundant biological condition

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

function trs_biolgroups(GrpType,RedundantFlag,varargin)
global K P
%trs_biolgroups('tf',0,[0.001,0.01],[1,1],{'WRKY2','WOX8'},{'At5g56270','At5g45980'},{'WRKY2','WOX8'},{[3108],[4022]},1)


%% MANUAL SELECTION
if RedundantFlag
    %!!!    NOT REFACTORED CODE !!!
    %distance limit used to select not redundant biol cond (i.e. distance > DistLimit)
    DistLimit=[0.1:0.05:0.8];

    if isequal(Type,'fdr')
        %default value for selecting the significative variation between a biol
        %cond and the median biol cond
        Fdr=0.10;
    end

    Continue='yes';
    while isequal(Continue,'yes')
        AnalType=questdlg(sprintf('Select the list of biological condition\nused to find non redundants'),'','All','List','Complement','List');
        if isequal(AnalType,'All')
            SelBiolRank=P.biol.pointIndex;
        elseif isequal(AnalType,'Complement')
            [SelPos,Ok]=listdlg('promptstring','select the biol cond not to be used','liststring',P.biol.grp.name,'selectionmode','multiple');
            if Ok==1
                ClearBiol=[];
                for SelL=1:length(SelPos)
                    ClearBiol=[ClearBiol;P.biol.grp.scoupleindex{SelPos(SelL)}];
                end
                ClearBiol=unique(ClearBiol);
                SelBiolRank=P.biol.pointIndex;
                for ClearL=1:length(ClearBiol)
                    SelBiolRank(find(SelBiolRank==ClearBiol(ClearL)))=[];
                end
            else
                h=errordlg('process canceled');
                waitfor(h)
                error('process canceled')
            end
        else
            %open a file or select
            [SelFile,SelDir]=uigetfile('*.txt','open the file with biol cond ranks');
            cd(SelDir)
            SelBiolRank=load(SelFile);
        end
        FilterIt=questdlg(sprintf('do you want to filter \nthe biological condition ?'),'','yes','no','yes');
        if isequal(FilterIt,'yes')
            options.Resize='on';
            options.WindowStyle='normal';
            Filter=inputdlg({'logical filter','filter name'},'construct logical filter',1,{'P.biol.type~=''L''&P.biol.status==''H''&P.biol.treated==''N''','{TC}HN'},options);
            eval(sprintf('BiolPos=%s;',Filter{1}))
            BiolPos=find(BiolPos);
            SelNb=length(SelBiolRank);
            for BiolL=length(SelBiolRank):-1:1
                if isempty(find(BiolPos==SelBiolRank(BiolL)))
                    SelBiolRank(BiolL)=[];
                end
            end
            h=warndlg(sprintf('%u out of  %u have been eliminated',SelNb-length(SelBiolRank),SelNb));
            waitfor(h)
        end

        if isequal(AnalType,'All')
            if isequal(FilterIt,'yes')
                GrpName=inputdlg({'give the name for this set of not redundant conditions'},'',1,{sprintf('all on %s (filtered => %s)',Type,Filter{2})});
                GrpName=GrpName{1};
            else
                GrpName=sprintf('all on %s',Type);
            end
        elseif isequal(AnalType,'Complement')
            if isequal(FilterIt,'yes')
                GrpName=inputdlg({'give the name for this set of not redundant conditions'},'',1,{sprintf('all wo nr nb ? on %s (filtered => %s)',Type,Filter{2})});
                GrpName=GrpName{1};
            else
                GrpName=inputdlg({'give the name for this set of not redundant conditions'},'',1,{sprintf('all wo nr nb ? on %s',Type)});
                GrpName=GrpName{1};
            end
        else
            if isequal(FilterIt,'yes')
                GrpName=inputdlg({'give the name for this set of not redundant conditions'},'',1,{sprintf('%s on %s (filtered => %s)',SelFile,Type,Filter{2})});
                GrpName=GrpName{1};
            else
                GrpName=inputdlg({'give the name for this set of not redundant conditions'},'',1,{sprintf('%s on %s',SelFile,Type)});
                GrpName=GrpName{1};
            end
        end

        %construct biol cond list from scoupleindex
        %only chip type 2 (1 = MAS4 not used) is used to detect non redundant points

        %ask for a limit on correlation between probeset rank
        InfCorrLimit=0;
        SupCorrLimit=1;
        if isfield(P.biol,'scouplecorr')
            h=figure;
            subplot(1,3,1)
            hist(P.biol.scouplecorr,1000)
            set(gca,'xlim',[0.8,1])

            subplot(1,3,2)
            hist(P.biol.scouplecorr,1000)
            set(gca,'xlim',[0.88,0.92])
            subplot(1,3,3)
            hist(P.biol.scouplecorr,1000)
            set(gca,'xlim',[0.98,0.99])

            CorrLimit=inputdlg({'inferior corr limit (selected >)';'superior corr limit (selected <)'},'',1,{'0.900','0.985'});
            InfCorrLimit=str2num(CorrLimit{1});
            SupCorrLimit=str2num(CorrLimit{2});
            delete(h)
        end


        UsedBiolList={};
        %UsedSCoupleRank=[];
        UsedBiolRank=[];
        for BiolL=1:length(SelBiolRank)
            BiolRank=SelBiolRank(BiolL);
            BiolPos=find(P.biol.pointIndex==BiolRank);
            %chip analyzed with MAS4 are not used
            if ~isfield(P.biol,'scouplecorr')
                P.biol.scouplecorr=cell(P.chip.nb);
            end
            if P.biol.softgrp(BiolPos)~=1 & ~isequal(P.biol.name{BiolRank},'MedianBiol') &P.biol.scouplecorr(BiolPos)<=SupCorrLimit&P.biol.scouplecorr(BiolPos)>=InfCorrLimit
                %UsedSCoupleRank=[UsedSCoupleRank;BiolL];
                UsedBiolList=[UsedBiolList;{sprintf('C%04u_%s',BiolRank,P.biol.name{BiolRank})}];
                UsedBiolRank=[UsedBiolRank;BiolRank];
            end

        end

        %calculate a distance  (1 - nb of common varying probeset / min of varying
        %probeset nb in each of the two biol cond that are compared that is 1 - similarity)
        UsedBiolNb=length(UsedBiolList);
        if isequal(Type,'rank')|isequal(Type,'log2(signal)')
            cd(K.dir.data)
            RankFile=sprintf('rank_%s_l%04u',lower(strrep(P.chip.name,'K_','K')),P.list.ser);
            load(RankFile);
            RTable=Table{1};
            if size(RTable,2)~=UsedBiolNb
                h=errordlg(sprintf('UsedBiolNb=%u & Size(Table,2)=%u !',UsedBiolNb,size(RTable,2)));
                waitfor(h)
            end
            if isequal(Type,'log2(signal)')
                RTable=interp1(K.chip.ref.rank,K.chip.ref.signal,RTable);
                MinSignal=K.chip.ref.signal(2);
                RTable(find(RTable<=0))=MinSignal;
                RTable=log2(RTable);
            end

            PsMean = mean(RTable');
            PsStd = std(RTable');
            SupLim = PsMean + 1*PsStd;
            InfLim = PsMean - 1*PsStd;
            figure
            plot(PsMean,PsStd,'b.')
            xlabel('Ps Mean')
            ylabel('Ps Std')

            Sup = zeros(size(RTable));
            % set values > SupLim to 1 and values <= SupLim to 0
            % => Bindex of values > SupLim
            for i = 1:size(RTable,1)
                Sup(i,find(RTable(i,:) > SupLim(i)))=1;
            end

            Inf = zeros(size(RTable));
            % set values < InfLim to 1 and values >= InfLim to 0
            % => Bindex of values < InfLim
            for i = 1:size(RTable,1)
                Inf(i,find(RTable(i,:) < InfLim(i)))=1;
            end


            Dist = [];

            for i=1:UsedBiolNb
                BiolSup1=Sup(:,i);
                SupNb1 = length(find(BiolSup1 > 0));
                BiolInf1=Inf(:,i);
                InfNb1 = length(find(BiolInf1 > 0));
                for j=i+1:UsedBiolNb
                    BiolSup2=Sup(:,j);
                    SupNb2 = length(find(BiolSup2 > 0));
                    BiolInf2=Inf(:,j);
                    InfNb2 = length(find(BiolInf2 > 0));
                    Dist=[Dist,1-(length(find(BiolSup1&BiolSup2))+length(find(BiolInf1&BiolInf2)))/min(SupNb1+InfNb1,SupNb2+InfNb2)];
                end
            end
        elseif isequal(Type,'fdr')
            ResFile=sprintf('redond_u%04.f_l%04.f.mat',P.user.ser,P.list.ser);
            cd(K.dir.redond)
            load(ResFile)
            %     if length(O.name)~=UsedBiolNb
            %         h=errordlg(sprintf('UsedBiolNb=%u & Length(O.name)=%u !',UsedBiolNb,length(O.name)));
            %         waitfor(h)
            %     end
            %construct the list of used comparison in O
            OPos=zeros(UsedBiolNb,1);
            RankVal=[];
            for ResL=1:length(O.name)
                CurrBiolRank=str2num(O.r{ResL}.name(5:8));
                if ~isempty(find(UsedBiolRank==CurrBiolRank))
                    OPos(find(UsedBiolRank==CurrBiolRank))=ResL;
                    RankVal=[RankVal;CurrBiolRank];
                end
            end


            PsNb=P.chip.limit{1}(2);
            NullBindex=zeros(PsNb,1);
            Dist = [];
            for BiolL1=1:UsedBiolNb-1
                IncBindex1=NullBindex;
                IncBindex1(O.r{OPos(BiolL1)}.incindex{1})=1;
                FdrBindex1=O.r{OPos(BiolL1)}.fdr{1}<=Fdr;
                IncNb1 = length(find(FdrBindex1&IncBindex1));
                DecBindex1=~NullBindex;
                DecBindex1(O.r{OPos(BiolL1)}.incindex{1})=0;
                DecBindex1(O.r{OPos(BiolL1)}.invindex{1})=0;
                DecNb1 = length(find(FdrBindex1&DecBindex1));
                for BiolL2=BiolL1+1:UsedBiolNb
                    IncBindex2=NullBindex;
                    IncBindex2(O.r{OPos(BiolL2)}.incindex{1})=1;
                    FdrBindex2=O.r{OPos(BiolL2)}.fdr{1}<=Fdr;
                    IncNb2 = length(find(FdrBindex2&IncBindex2));
                    DecBindex2=~NullBindex;
                    DecBindex2(O.r{OPos(BiolL2)}.incindex{1})=0;
                    DecBindex2(O.r{OPos(BiolL2)}.invindex{1})=0;
                    DecNb2 = length(find(FdrBindex2&DecBindex2));
                    Dist=[Dist,1-(length(find(IncBindex1&IncBindex2&FdrBindex1&FdrBindex2))+length(find(DecBindex1&DecBindex2&FdrBindex1&FdrBindex2)))/min(IncNb1+DecNb1,IncNb2+DecNb2)];
                end
            end
        else
            h=errordlg(sprintf('unknown type %s',Type));
            waitfor(h)
        end
        Link=linkage(Dist);


        E=squareform(Dist);
        % indexation !
        [Temp ColSort]=sort(E(:,2));
        E1 = E(ColSort,ColSort);
        ME=mean(E);
        [Temp MeanSort]=sort(ME');
        MeanSort=fliplr(MeanSort);
        E2 = E(MeanSort,MeanSort);

        h0=figure;
        subplot(1,3,1)
        h=pcolor(E);
        set(h,'linestyle','none')
        title('distance (D)')
        subplot(1,3,2)
        h=pcolor(E1);
        set(h,'linestyle','none')
        title('sorted on snd column of D')
        subplot(1,3,3)
        h=pcolor(E2);
        set(h,'linestyle','none')
        title('sorted on mean of each column of D')
        cd(K.dir.sosresult)
        K.dir.restree=fullfile(K.dir.sosresult,sprintf('m%03u',ModelRank));
        try
            cd(sprintf('m%03u',ModelRank))
        catch
            mkdir(sprintf('m%03u',ModelRank))
            cd(sprintf('m%03u',ModelRank))
        end
        if isequal(P.flag.station,'windows')
            saveas(h0,sprintf('%s_values_m%03u_%s',strrep(GrpName,' ','_'),ModelRank,date),'bmp')
        else
            saveas(h0,sprintf('%s_values_m%03u_%s',strrep(GrpName,' ','_'),ModelRank,date),'png')
        end



        %prepare list of biol cond with different level of similarity
        for DistL=1:length(DistLimit)
            %starting list
            List{DistL} = [1];
            for i = 2:UsedBiolNb
                if isempty(find(E(i,List{DistL})< DistLimit(DistL)));
                    %increment list if identity between tested condition and all other
                    %conditions already in the list is less thant threshold
                    List{DistL} = [List{DistL} i];
                end
            end
        end

        NewSelection='yes';
        while isequal(NewSelection,'yes')
            %select a distance limit
            LimitRank=1;
            if length(DistLimit)>1
                NRNb=[];
                Msg='';
                for DistL=1:length(DistLimit)
                    NRNb=[NRNb;length(List{DistL})];
                    Msg=sprintf('%s %u (%.2f : %u NR) ',Msg,DistL,DistLimit(DistL),NRNb(DistL));
                end
                [DistLimit',NRNb]
                LimitRank=inputdlg({'give LimitRank of NR selection'},Msg,1,{'1'},'on');
                LimitRank=str2num(LimitRank{1});
            end

            %visualization of trees

            SList=List{LimitRank};
            CurrUsedBiolList=UsedBiolList;
            for BiolL=1:UsedBiolNb
                if ~isempty(find(SList==BiolL))
                    CurrUsedBiolList{BiolL}=sprintf('%s <==>',CurrUsedBiolList{BiolL});
                else
                    CurrUsedBiolList{BiolL}=sprintf('%s --------',CurrUsedBiolList{BiolL});
                end
            end


            hfig=figure;
            set(hfig,'color',[1,1,1])
            [h,Tree,Order]=dendrogram(Link,0);
            CurrUsedBiolList=CurrUsedBiolList(Order);
            title(Type)
            TreeSize=length(Order);
            set(h,'linewidth',2)
            set(h,'color','k')
            set(gca,'ytick',[1:TreeSize]);
            set(gca,'xlim',[0,1]);
            set(gca,'ylim',[0,TreeSize+1]);
            %set(gca,'yticklabel',CurrUsedBiolList);
            set(gca,'position',[0.50 0.02 0.50 0.92])
            Dend=gca;
            set(Dend,'yticklabel',CurrUsedBiolList);
            LeafNb=length(CurrUsedBiolList);
            TreeSize=inputdlg({sprintf('There are %u leaves. Indicate the number of leaves by image',LeafNb)},'',1,{'50'});
            TreeSize=str2num(TreeSize{1});
            LeafName=CurrUsedBiolList;

            XLimit=get(gca,'xlim');
            if LeafNb>=TreeSize
                set(gca,'ytick',[1:TreeSize]);
                set(gca,'xlim',XLimit);
                set(gca,'ylim',[0,TreeSize+1]);
                set(gca,'yticklabel',LeafName(Order((1:TreeSize))));
                set(gca,'position',[0.50 0.11 0.450 0.815])
                title(sprintf('%s (%s) tree nb 1 - Not redundant selected at fdr<=%.2f - Chip A',strrep(P.chip.name,'_',' '),date,DistLimit(LimitRank)))
            else
                set(gca,'ytick',[1:LeafNb]);
                set(gca,'ylim',XLimit);
                set(gca,'ylim',[0,LeafNb+1]);
                set(gca,'yticklabel',LeafName(Order));
                set(gca,'position',[0.50 0.11 0.450 0.815])
                title(sprintf('%s (%s) - Not redundant selected fdr<=%.2f - Chip A',strrep(P.chip.name,'_',' '),date,DistLimit(LimitRank)))
            end
            Units=get(0,'units');
            set(0,'units','normalized')
            scrsz = get(0,'ScreenSize');
            eval(['set(0,''units'',''',Units,''')'])
            h1=get(gca,'parent');
            set(h1,'units','normalized')
            set(h1,'Position',[0 0 scrsz(3)*0.5  scrsz(4)*0.93])
            cd(K.dir.restree)
            set(gca,'xlim',XLimit);
            SaveIt=questdlg('Do you want to save figures ?','','yes','no','no');
            if isequal(SaveIt,'yes')
                if isequal(P.flag.station,'windows')
                    saveas(h1,sprintf('%s_dist%u_m%03u_%s',strrep(GrpName,' ','_'),round(100*DistLimit(LimitRank)),ModelRank,date),'bmp')
                else
                    saveas(h1,sprintf('%s_dist%u_m%03u_%s',strrep(GrpName,' ','_'),round(100*DistLimit(LimitRank)),ModelRank,date),'png')
                end
            end

            DisplayNb=ceil(LeafNb/TreeSize);
            if DisplayNb>=2
                %for DisplayLoop=[2:DisplayNb-1]
                Display=1;
                Rank=1;
                while Display==1
                    Rank=inputdlg(sprintf('give the rank of sub tree to be displayed\n (>0 & <%u) or cancel (<0 or >=%u)',DisplayNb+1,DisplayNb+1),'',1,{sprintf('%u',Rank+1)});
                    Rank=str2num(Rank{1});
                    if Rank<1|Rank>DisplayNb
                        Answer=questdlg('do you want to cancel','','Yes','No','No');
                        if isequal(Answer,'Yes')
                            Display=0;
                        end
                    else
                        if Rank<DisplayNb
                            DisplayLoop=Rank;
                            set(gca,'ytick',[(DisplayLoop-1)*TreeSize+1:(DisplayLoop*TreeSize)]);
                        else
                            set(gca,'ytick',[(DisplayNb-1)*TreeSize+1:LeafNb]);
                        end
                        set(gca,'xlim',XLimit);
                        if Rank<DisplayNb
                            set(gca,'ylim',[(DisplayLoop-1)*TreeSize,(DisplayLoop*TreeSize)+1]);
                            set(gca,'yticklabel',LeafName(Order((DisplayLoop-1)*TreeSize+1:DisplayLoop*TreeSize)));
                        else
                            set(gca,'ylim',[(DisplayNb-1)*TreeSize,(DisplayNb*TreeSize)+1]);
                            set(gca,'yticklabel',LeafName(Order((DisplayNb-1)*TreeSize+1:LeafNb)));
                        end
                        set(gca,'position',[0.50 0.11 0.450 0.815])
                        title(sprintf('%s (%s) tree nb %u - Not redundant selected at fdr<=0.10 - Chip A',strrep(P.chip.name,'_',' '),date,Rank))
                        h1=get(gca,'parent');
                        set(h1,'units','normalized')
                        set(h1,'Position',[0 0 scrsz(3)*0.50  scrsz(4)*0.93])
                        if isequal(SaveIt,'yes')
                            if isequal(P.flag.station,'windows')
                                if isequal(FilterIt,'yes')
                                    saveas(h1,sprintf('%s_dist%u_%s_%s_%s_nb%u',strrep(GrpName,' ','_'),round(100*DistLimit(LimitRank)),Filter{2},lower(P.chip.name),date,Rank),'bmp')
                                else
                                    saveas(h1,sprintf('%s_dist%u_%s_%s_nb%u',strrep(GrpName,' ','_'),round(100*DistLimit(LimitRank)),lower(P.chip.name),date,Rank),'bmp')
                                end
                            else
                                if isequal(FilterIt,'no')
                                    saveas(h1,sprintf('%s_dist%u_%s_%s_nb%u',strrep(GrpName,' ','_'),round(100*DistLimit(LimitRank)),lower(P.chip.name),date,Rank),'png')
                                else
                                    saveas(h1,sprintf('%s_dist%u_%s_%s_%s_nb%u',strrep(GrpName,' ','_'),round(100*DistLimit(LimitRank)),Filter{2},lower(P.chip.name),date,Rank),'png')
                                end
                            end
                        end
                    end
                end
            end
            NewSelection=questdlg(sprintf('yes = try another distance limit\n(%.02f); no = save results',DistLimit(LimitRank)),'','yes','no','no');
        end
    end
else
    %% SPLIT
    switch GrpType
        case 'split'
            %Construct groups of biological conditions
            if ~isfield(P.biol,'grp')
                P.biol.grp.name={};
                P.biol.grp.usedBiolRanks={};
                P.biol.grp.biolRanks={};
            end


            UsedBiolRanks=find(P.net.biolIndex);
            
            Seed=0;

            BiolNb=length(UsedBiolRanks)+1;
            while BiolNb>length(UsedBiolRanks)
                BiolNb=inputdlg(sprintf('(total=%u) How many biol cond in each group?',length(UsedBiolRanks)),'',1,{'30'});
                BiolNb=str2num(BiolNb{1});
                if BiolNb>UsedBiolRanks
                    h=warndlg(sprintf('BiolNb must be <=%u',length(UsedBiolRanks)));
                    waitfor(h)
                end
            end
            %select the seed for the random process
            while Seed==0
                Seed=inputdlg('Seed of algo (twister != 0)','',1,{'5489'});
                if isempty(Seed)
                    h=errordlg('Process canceled');
                    waitfor(h)
                    error('process canceled')
                else
                    Seed=str2num(Seed{1});
                    if Seed==0
                        h=warndlg('seed must not be equal to 0');
                        waitfor(h)
                    end
                end
            end
            rand('twister',Seed);
            CurrBiolRanks=UsedBiolRanks;
            if BiolNb<=length(CurrBiolRanks)
                if floor(length(CurrBiolRanks)/BiolNb)>1
                    SeveralGrp=questdlg(sprintf('It is possible to make %u groups\nof\n %u biol cond. Do you want that ?',floor(length(CurrBiolRanks)/BiolNb),BiolNb),'','yes','no','yes');
                    if isequal(SeveralGrp,'no')
                        RandPos=randperm(length(CurrBiolRanks));
                        BiolRanks{1}=CurrBiolRanks(RandPos(1:BiolNb));
                    else
                        if mod(length(CurrBiolRanks),BiolNb)>0
                            LastGrp=questdlg(sprintf('Last group has %u condtions.\nDo you want to use it ?',mod(length(CurrBiolRanks),BiolNb)),'','yes','no','no');
                        else
                            LastGrp='no';
                        end
                        for GrpL=1:floor(length(CurrBiolRanks)/BiolNb)
                            RandPos=randperm(length(CurrBiolRanks));
                            BiolRanks{GrpL}=CurrBiolRanks(RandPos(1:BiolNb));
                            CurrBiolRanks(RandPos(1:BiolNb))=[];
                        end
                        if isequal(LastGrp,'yes')
                            BiolRanks{GrpL+1}=CurrBiolRanks;
                        end
                    end
                else
                    RandPos=randperm(length(CurrBiolRanks));
                    BiolRanks{1}=CurrBiolRanks(RandPos(1:BiolNb));
                end
            end


            GrpName=inputdlg({'give the name for this set of not redundant conditions'},'',1,{''});
            GrpName=GrpName{1};
            GrpPos=length(P.biol.grp.name)+1;
            GrpNb=length(BiolRanks);
            ReplaceFlag=0;
            GrpPositions=[];
            for GrpL=1:GrpNb
                %same name could exist but with a different BiolNb
                CurrGrpName=sprintf('%s %ubiolcond grp%u twister_seed%u',GrpName,length(BiolRanks{GrpL}),GrpL,Seed);
                NamePos=strmatch(CurrGrpName,P.biol.grp.name,'exact');
                Correct=0;
                while Correct==0
                    if ~isempty(NamePos)
                        if ReplaceFlag
                            GrpPos=NamePos;
                        else
                            ReplaceIt=questdlg(sprintf('%s exist. Do you want to replace it ?',GrpName),'','yes','no','no');
                            if isequal(ReplaceIt,'yes')
                                GrpPos=NamePos;
                                ReplaceFlag=1;
                                Correct=1;
                            else
                                GrpName=inputdlg({'give the name for this set of not redundant conditions'},'',1,{''});
                                GrpName=GrpName{1};
                                CurrGrpName=sprintf('%s %ubiolcond grp%u twister_seed%u',GrpName,length(BiolRanks{GrpL}),GrpL,Seed);
                                NamePos=strmatch(CurrGrpName,P.biol.grp.name,'exact');
                            end
                        end
                    else
                        Correct=1;
                    end
                end
                GrpPositions(end+1,1)=GrpPos;
                P.biol.grp.usedBiolRanks{GrpPos,1}=UsedBiolRanks;
                P.biol.grp.name{GrpPos,1}=sprintf('%s %ubiolcond grp%u twister_seed%u',GrpName,length(BiolRanks{GrpL}),GrpL,Seed);
                P.biol.grp.biolRanks{GrpPos,1}=sort(BiolRanks{GrpL});
                GrpPos=length(P.biol.grp.name)+1;
            end
            cd(P.dir.project)
            eval(sprintf('save %s P',P.project.name));

            %fill K.net
            BLOC_SIZE=100;
            NET_NB=20;
            %first split created biol grp in two equal part
            %sort GrpPosition to always have GrpPos2>GrpPos1
            GrpPosition=sort(GrpPosition);
            GrpPos1=GrpPosition(1:floor(GrpNb/2));
            GrpPos2=GrpPosition(floor(GrpNb/2)+1:floor(GrpNb/2)+length(GrpPos1));
            NetGrps=zeros(NET_NB,2);
            NetGrp(1:length(GrpPos1),1)=GrpPos1';
            NetGrp(1:length(GrpPos1),2)=GrpPos2';
            NetPos=length(GrpPos1);
            %add random pairs of biological conditions if necessary
            if length(GrpPos1)<NET_NB
                for NetL=1:NET_NB-length(GrpPos1)
                    Pos1=NetGpr(1,1);
                    Pos2=NetGpr(1,2);
                    while ~isempty(find(NetGrp(:,1)==Pos1))&~isempty(find(NetGrp(:,2)==Pos2))
                        RandPos=randperm(GrpNb);
                        Pos1=GrpPositions(RandPos(1));
                        Pos2=GrpPositions(RandPos(2));
                        if Pos1>Pos2
                            temp=Pos1;
                            Pos1=Pos2;
                            Pos2=temp;
                        end
                    end
                    NetPos=NetPos+1;
                    NetGrp(NetPos:1)=Pos1;
                    NetGrp(NetPos:2)=Pos2;
                end
            end

            for NetL=1:NET_NB
                Pos1=NetGrp(NetL,1);
                Pos2=NetGrp(NetL,2);
                BiolRanks1=P.biol.grp.biolRanks{Pos1};
                BiolRanks2=P.biol.grp.biolRanks{Pos2};
                NetName=sprintf('g%us% vs g%us% tws%u',Pos1,length(BiolRanks1),Pos2,length(BiolRanks2),Seed);
                NetPos=strmatch(NetName,K.net{ModelRank}.name,'exact');
                if isempty(NamePos)
                    NetPos=length(K.net{ModelRank})+1;
                    K.net{ModelRank}.name{NetPos,1}=NetName;
                    K.net{ModelRank}.rank(NetPos,1)=max(K.net{ModelRank}.rank)+1;
                    K.net{ModelRank}.biolRank{NetPos,1}{1}=BiolRanks1;
                    K.net{ModelRank}.biolRank{NetPos,1}{2}=BiolRanks2;
                    K.net{ModelRank}.compNb(NetPos,1)=length(K.net{ModelRank}.biolRank{NetPos,1}{1})*length(K.net{ModelRank}.biolRank{NetPos,1}{2});
                    K.net{ModelRank}.fdr(NetPos,:)=[0,001 0.010];
                    K.net{ModelRank}.s(NetPos,:)=[0,1];
                    K.net{ModelRank}.blocNb(NetPos,1)=ceil(PsNb/BLOC_SIZE);
                    K.net{ModelRank}.blocSize(NetPos,1)=BLOC_SIZE;
                    K.net{ModelRank}.comment{NetPos,1}='';
                else
                    K.net{ModelRank}.biolRank{NetPos,1}{1}=BiolRanks1;
                    K.net{ModelRank}.biolRank{NetPos,1}{2}=BiolRanks2;
                    K.net{ModelRank}.compNb(NetPos,1)=length(K.net{ModelRank}.biolRank{NetPos,1}{1})*length(K.net{ModelRank}.biolRank{NetPos,1}{2});
                    K.net{ModelRank}.netMade(NetPos,1)=0;
                end
                K.net{ModelRank}.netMade(NetPos,1)=0;
            end
            cd(K.dir.common)
            Tempo=K.net;
            save netlist Tempo
            clear Tempo
            %% TF
        case 'tf'

            if nargin>2
                Fdr=varargin{1};
                Sensitivity=varargin{2};
                Tf=varargin{3};
                TfId=varargin{4};
                TfGeneSymbol=varargin{5};
                TfPs=varargin{6};
                StatFlag=0;
                TargetFlag=0;
                TfNb=length(Tf);
                %nb of probe set that target the TF
                TfPsNb=zeros(TfNb,1);
                for TfL=1:TfNb
                    TfPsNb(TfL)=length(TfPs{TfL});
                end
                NonRedundantFlag=0;
                if nargin>=9
                    NonRedundantFlag=varargin{7};
                end                
            else
                TfPs={};
                StatFlag=1;
                TargetFlag=1;
                NonRedundantFlag=1;
            end
            %calculate mean ranks for biological data
            if isfield(P.biol,'scouplenb')
                BiolNb=P.biol.scouplenb;
                BiolRanks=P.biol.scoupleindex{1};
                PointRanks=P.biol.scouple{1};
                PsNb=P.chip.limit{1}(2);
                DataRanks=single(zeros(PsNb,BiolNb));

                %construct psxmean(rank) matrix
                cd(K.dir.work)
                for BiolL=1:BiolNb
                    CurrPoints=PointRanks{BiolRanks(BiolL)};
                    eval(sprintf('load Data%04u',CurrPoints(1)))
                    Data=CData{1}.rank;
                    eval(sprintf('load Data%04u',CurrPoints(2)))
                    Data=[Data,CData{1}.rank];
                    DataRanks(:,BiolL)=mean(Data,2);
                end
            else
                BiolRanks=find(P.biol.used);
                BiolNb=length(BiolRanks);
                PointRanks=P.biol.pairs(BiolRanks);
                PsNb=P.chip.currProbeSetNb;
                DataRanks=uint8(zeros(PsNb,BiolNb));
                %construct psxmean(rank) matrix
                cd(P.dir.data)
                AllData=load_data('DataRanks.float32le','./',PsNb,P.point.nb,'single','ieee-le');
                for BiolL=1:BiolNb
                    CurrPoints=PointRanks{BiolL};
                    DataRanks(:,BiolL)=uint8(mean(AllData(:,CurrPoints),2));
                end
            end

            if StatFlag
                %load stat on probe sets in order to use a sub set of ps classes
                [FileName,FileDir]=uigetfile('*.mat','open probe set stat file or cancel if not used');
                if isequal(class(FileName),'double')
                    StatFlag=0;
                else
                    StatFlag=1;
                    cd(FileDir)
                    load(FileName)
                    clear NewPs
                end
            end
            if TargetFlag
                %load tf gene interractions
                [FileName,FileDir]=uigetfile('*.txt','open tf gene interactions or cancel if not known');
                if isequal(class(FileName),'double')
                    TargetFlag=0;
                else
                    cd(FileDir)
                    [Tfs,Targets,temp]=textread(FileName,'%s%s%s','delimiter','\t');
                    TargetFlag=1;
                    %eliminate nonmus genes
                    MusPos=strmatch('ENSMUSG',Targets);
                    Tfs=Tfs(MusPos);
                    Targets=Targets(MusPos);

                    %load tf name
                    [FileName,FileDir]=uigetfile('*.txt','open tf names');
                    cd(FileDir)
                    [Tf,TfId,TfGeneSymbol]=textread(FileName,'%s%s%s','delimiter','\t');
                    %load probe set order
                    [FileName,FileDir]=uigetfile('*.mat','open chip probeset');
                    cd(FileDir)
                    load(FileName)
                end
            end

            TfNb=length(Tf);
            if  isempty(TfPs)
                %probe sets that target the TF
                TfPs=cell(TfNb,1);
                %nb of probe set that target the TF
                TfPsNb=zeros(TfNb,1);
            end
            if StatFlag
                %list of nb of probes targeting the TF
                TfPNb=cell(TfNb,1);
            end
            if TargetFlag
                %nb of genes assigned to at least one probe set that are targeted by the TF
                TargetNb=zeros(TfNb,1);
                %Ensembl id of targeted genes
                TargetIds=cell(TfNb,1);
                %Probe set ranks of the targeted genes
                TargetPs=cell(TfNb,1);
                for TfL=1:TfNb
                    CurrTfId=TfId{TfL};
                    CurrTfId=str2num(CurrTfId(8:end));
                    PsRanks=find(PsInfo(:,1)==CurrTfId);
                    if ~isempty(PsRanks)
                        TfPs{TfL}=PsRanks;
                        TfPsNb(TfL)=length(PsRanks);
                        TfPNb{TfL}=PsInfo(PsRanks,9);
                    end
                    TargetPos=strmatch(Tf{TfL},Tfs,'exact');
                    TargetNb(TfL)=length(TargetPos);
                    TargetIds{TfL}=[];
                    TargetPs{TfL}=[];
                    for TargetL=1:TargetNb(TfL)
                        CurrTarget=Targets{TargetPos(TargetL)};
                        CurrTarget=str2num(CurrTarget(8:end));
                        PsRanks=find(PsInfo(:,1)==CurrTarget);
                        if ~isempty(PsRanks)
                            TargetIds{TfL}(end+1,1)=CurrTarget;
                            TargetPs{TfL}=[TargetPs{TfL};PsRanks];
                        end
                    end
                end
            end
            %select TF for which exists a probe set
            AbsentPos=find(TfPsNb==0);
            %list absent TF
            TfGeneSymbol(AbsentPos)

            %keep present ones
            Tf(AbsentPos)=[];
            TfGeneSymbol(AbsentPos)=[];
            TfPs(AbsentPos)=[];
            TfPsNb(AbsentPos)=[];
            if StatFlag
                TfPNb(AbsentPos)=[];
            end
            if TargetFlag
                TargetNb(AbsentPos)=[];
                TargetIds(AbsentPos)=[];
                TargetPs(AbsentPos)=[];
                %sort by descendent nb of target
                [TargetNb,SortIndex]=sort(TargetNb,'descend');
                Tf=Tf(SortIndex);
                TfGeneSymbol=TfGeneSymbol(SortIndex);
                TfPs=TfPs(SortIndex);
                TfPsNb=TfPsNb(SortIndex);
                TfPNb=TfPNb(SortIndex);
                TargetIds=TargetIds(SortIndex);
                TargetPs=TargetPs(SortIndex);
            end
            TfNb=length(Tf);
            %transform cell into vector
            MTf=Tf;
            MTfGeneSymbol=TfGeneSymbol;
            MTfPs=TfPs;
            MTfPsNb=TfPsNb;
            if StatFlag
                MTfPNb=TfPNb;
            end
            if TargetFlag
                MTargetNb=TargetNb;
                MTargetIds=TargetIds;
                MTargetPs=TargetPs;
            end
            Tf={};
            TfGeneSymbol={};
            TfPs=[];
            TfPsNb=[];
            if StatFlag
                TfPNb=[];
            end
            if TargetFlag
                TargetNb=[];
                TargetId={};
                TargetPs={};
                TargetPsNb=[];
            end

            for TfL=1:TfNb
                for PsL=1:MTfPsNb(TfL)
                    Tf{end+1,1}=MTf{TfL};
                    TfGeneSymbol{end+1,1}=MTfGeneSymbol{TfL};
                    TfPs(end+1,1)=MTfPs{TfL}(PsL);
                    TfPsNb(end+1,1)=MTfPsNb(TfL);
                    if StatFlag
                        TfPNb(end+1,1)=MTfPNb{TfL}(PsL);
                    end
                    if TargetFlag
                        TargetNb(end+1,1)=MTargetNb(TfL);
                        TargetId{end+1,1}=MTargetIds{TfL};
                        TargetPs{end+1,1}=MTargetPs{TfL};
                        TargetPsNb(end+1,1)=length(MTargetPs{TfL});
                    end
                end
            end
            %sort by decreasing probe nb and increasing ps nb
            [TfPsNb,SortIndex]=sort(TfPsNb);
            Tf=Tf(SortIndex);
            TfGeneSymbol=TfGeneSymbol(SortIndex);
            TfPs=TfPs(SortIndex);
            if StatFlag
                TfPNb=TfPNb(SortIndex);
            end
            if TargetFlag
                TargetNb=TargetNb(SortIndex);
                TargetId=TargetId(SortIndex);
                TargetPs=TargetPs(SortIndex);
                TargetPsNb=TargetPsNb(SortIndex);
            end

            if StatFlag
                [TfPNb,SortIndex]=sort(TfPNb,'descend');
                Tf=Tf(SortIndex);
                TfGeneSymbol=TfGeneSymbol(SortIndex);
                TfPs=TfPs(SortIndex);
                TfPsNb=TfPsNb(SortIndex);
                if TargetFlag
                    TargetNb=TargetNb(SortIndex);
                    TargetId=TargetId(SortIndex);
                    TargetPs=TargetPs(SortIndex);
                    TargetPsNb=TargetPsNb(SortIndex);
                end
            end

            if NonRedundantFlag
                if isfield(P.biol,'scouplenb')
                    %select up and down among the non redundant set of biological conditions
                    BiolRanks=P.biol.nr{1}.startbiol{4};
                    BiolPos=zeros(length(BiolRanks),1);
                    %find the position in the Psxmean(Ranks) matrix
                    for BiolL=1:length(BiolRanks)
                        BiolPos(BiolL)=find(P.biol.scoupleindex{1}==BiolRanks(BiolL));
                    end
                    %extract the relevant part of the matrix
                    DataRanks=DataRanks(:,BiolPos);
                else
                    BiolRanks=find(P.net.biolIndex);
                    UsedBiolRanks=zeros(P.biol.nb,1);
                    UsedBiolRanks(BiolRanks)=1;
                    UsedBiolRanks=UsedBiolRanks(find(P.biol.used));
                    DataRanks=DataRanks(:,find(UsedBiolRanks));
                end
            end

            %do statistics on each TF
            TfNb=length(Tf);
            %the 20th,25th and 30th rank percentile
            RankLimit=cell(TfNb,1);
            %limit for percentiles
            Limit=[[20,80];[25,75];[30,70]];
            for TfL=1:TfNb
                Ranks=DataRanks(TfPs(TfL),:);
                %eliminate zero values
                NullPos=find(Ranks<=0);
                Ranks(NullPos)=[];
                RankLimit{TfL}=zeros(3,2);
                for LimitL=1:3
                    RankLimit{TfL}(LimitL,1)=length(find(Ranks<=Limit(LimitL,1)));
                    RankLimit{TfL}(LimitL,2)=length(find(Ranks<=Limit(LimitL,2)));
                end
            end

            Colors=colors(colormap,TfNb);
            BiolNb=length(BiolRanks);
            h=figure;
            set(gcf,'color',[1,1,1])
            set(h,'name','TF PROFILES')
            hold on
            for TfL=1:TfNb
                Ranks=DataRanks(TfPs(TfL),:);
                plot(1:BiolNb,sort(Ranks),'color',Colors(TfL,:))
            end
            set(gca,'box','on')


            %             %select a limit and construct P.biol.nr informations
            %             PercName={'20th';'25th';'30th'};
            %             GrpSize=30;
            %             FirstNrPos=length(P.biol.nr{1})+1;
            %             NetNb=0;
            %             UsedPs=[];
            %             for TfL=1:TfNb
            %                 Pos1=find(RankLimit{TfL}(:,1)>=GrpSize);
            %                 if ~isempty(Pos1)
            %                     Pos1=Pos1(1);
            %                     Pos2=find(RankLimit{TfL}(:,2)>=GrpSize);
            %                     if ~isempty(Pos2)
            %                         Pos2=Pos2(1);
            %                         NetNb=NetNb+1;
            %                         UsedPs(end+1,1)=TfL;
            %                         Ranks=DataRanks(TfPs(TfL),:);
            %                         NrPos=length(P.biol.nr{1})+1;
            %                         P.biol.nr{1}.startbiol{NrPos}=P.biol.nr{1}.startbiol{4};
            %                         P.biol.nr{1}.name{NrPos}=sprintf('%s DOWN ps%u %sperc',TfGeneSymbol{TfL},TfPs(TfL),PercName{Pos1});
            %                         P.biol.nr{1}.scoupleindex{NrPos}=BiolRanks(find(Ranks<=RankLimit{TfL}(Pos1,1)));
            %                         P.biol.nr{1}.type{NrPos}='up/down';
            %                         P.biol.nr{1}.corrlimit{NrPos}=P.biol.nr{1}.corrlimit{4};
            %                         NrPos=length(P.biol.nr{1})+1;
            %                         P.biol.nr{1}.startbiol{NrPos}=P.biol.nr{1}.startbiol{4};
            %                         P.biol.nr{1}.name{NrPos}=sprintf('%s UP ps%u %sperc',TfGeneSymbol{TfL},TfPs(TfL),PercName{Pos1});
            %                         P.biol.nr{1}.scoupleindex{NrPos}=BiolRanks(find(Ranks>=RankLimit{TfL}(Pos1,2)));
            %                         P.biol.nr{1}.type{NrPos}='up/down';
            %                         P.biol.nr{1}.corrlimit{NrPos}=P.biol.nr{1}.corrlimit{4};
            %                     end
            %                 end
            %             end


            GrpSize=30;
            if isfield(P.biol,'scouplenb')
                FirstNrPos=length(P.biol.nr{1}.name)+1;
            else
                FirstNrPos=length(P.biol.grp.name)+1;
            end
            NetNb=0;
            UsedPs=[];
            for TfL=1:TfNb
                Ranks=DataRanks(TfPs(TfL),:);
                [temp,SortIndex]=sort(Ranks);
                NullNb=length(find(Ranks==0));
                NetNb=NetNb+1;
                UsedPs(end+1,1)=TfL;
                if isfield(P.biol,'scouplenb')
                    NrPos=length(P.biol.nr{1}.name)+1;
                    P.biol.nr{1}.startbiol{NrPos}=P.biol.nr{1}.startbiol{4};
                    P.biol.nr{1}.name{NrPos}=sprintf('%s DOWN ps%u',TfGeneSymbol{TfL},TfPs(TfL));
                    P.biol.nr{1}.scoupleindex{NrPos}=sort(BiolRanks(SortIndex(1+NullNb:GrpSize+NullNb)));
                    P.biol.nr{1}.type{NrPos}='up/down';
                    P.biol.nr{1}.corrlimit{NrPos}=P.biol.nr{1}.corrlimit{4};
                    NrPos=length(P.biol.nr{1}.name)+1;
                    P.biol.nr{1}.startbiol{NrPos}=P.biol.nr{1}.startbiol{4};
                    P.biol.nr{1}.name{NrPos}=sprintf('%s UP ps%u',TfGeneSymbol{TfL},TfPs(TfL));
                    P.biol.nr{1}.scoupleindex{NrPos}=sort(BiolRanks(SortIndex(end-GrpSize+1:end)));
                    P.biol.nr{1}.type{NrPos}='up/down';
                    P.biol.nr{1}.corrlimit{NrPos}=P.biol.nr{1}.corrlimit{4};
                else
                    NrPos=length(P.biol.grp.name)+1;
                    P.biol.grp.name{NrPos}=sprintf('%s DOWN ps%u',TfGeneSymbol{TfL},TfPs(TfL));
                    P.biol.grp.biolRanks{NrPos}=sort(BiolRanks(SortIndex(1+NullNb:GrpSize+NullNb)));
                    NrPos=length(P.biol.grp.name)+1;
                    P.biol.grp.name{NrPos}=sprintf('%s UP ps%u',TfGeneSymbol{TfL},TfPs(TfL));
                    P.biol.grp.biolRanks{NrPos}=sort(BiolRanks(SortIndex(end-GrpSize+1:end)));
                end
            end
            if isfield(P.biol,'scouplenb')
                cd(K.dir.point)
                eval(sprintf('save %s P',P.file.point))
            else
                cd(P.dir.project)
                eval(sprintf('save %s P',P.project.name))
            end

            %fill K.net
            BLOC_SIZE=1000;
            if isfield(P.biol,'scouplenb')
                for NetL=1:NetNb
                    NetPos=length(K.net{ModelRank}.name)+1;
                    K.net{ModelRank}.name{NetPos,1}=sprintf('%s ps%u',TfGeneSymbol{UsedPs(NetL)},TfPs(UsedPs(NetL)));
                    K.net{ModelRank}.rank(NetPos,1)=max(K.net{ModelRank}.rank)+1;
                    K.net{ModelRank}.biolRank{NetPos,1}{1}=P.biol.nr{1}.scoupleindex{FirstNrPos};
                    K.net{ModelRank}.biolRank{NetPos,1}{2}=P.biol.nr{1}.scoupleindex{FirstNrPos+1};
                    FirstNrPos=FirstNrPos+2;
                    K.net{ModelRank}.compNb(NetPos,1)=length(K.net{ModelRank}.biolRank{NetPos,1}{1})*length(K.net{ModelRank}.biolRank{NetPos,1}{2});
                    K.net{ModelRank}.fdr(NetPos,:)=Fdr;
                    K.net{ModelRank}.s(NetPos,:)=Sensitivity;
                    K.net{ModelRank}.blocNb(NetPos,1)=ceil(PsNb/BLOC_SIZE);
                    K.net{ModelRank}.blocSize(NetPos,1)=BLOC_SIZE;
                    K.net{ModelRank}.comment{NetPos,1}='';
                    K.net{ModelRank}.netMade(NetPos,1)=0;
                end
            else
                ModelRank=P.chip.chipRank;
                for NetL=1:NetNb
                    NetPos=length(K.net{ModelRank}.name)+1;
                    K.net{ModelRank}.name{NetPos,1}=sprintf('%s ps%u',TfGeneSymbol{UsedPs(NetL)},TfPs(UsedPs(NetL)));
                    K.net{ModelRank}.rank(NetPos,1)=max(K.net{ModelRank}.rank)+1;
                    K.net{ModelRank}.biolRank{NetPos,1}{1}=P.biol.grp.biolRanks{FirstNrPos};
                    K.net{ModelRank}.biolRank{NetPos,1}{2}=P.biol.grp.biolRanks{FirstNrPos+1};
                    FirstNrPos=FirstNrPos+2;
                    K.net{ModelRank}.compNb(NetPos,1)=length(K.net{ModelRank}.biolRank{NetPos,1}{1})*length(K.net{ModelRank}.biolRank{NetPos,1}{2});
                    K.net{ModelRank}.fdr(NetPos,:)=Fdr;
                    K.net{ModelRank}.s(NetPos,:)=Sensitivity;
                    K.net{ModelRank}.blocNb(NetPos,1)=ceil(PsNb/BLOC_SIZE);
                    K.net{ModelRank}.blocSize(NetPos,1)=BLOC_SIZE;
                    K.net{ModelRank}.comment{NetPos,1}='';
                    K.net{ModelRank}.netMade(NetPos,1)=1;
                end
            end
            cd(K.dir.common)
            Tempo=K.net;
            save netlist Tempo
            clear Tempo
    end
end
