%TRS_NOTREDUNDANT allows to select non redundant biological conditions.
% A distance between biological conditions is calculated, which allows to 
% select conditions that ere not too similar

%INPUT PARAMETERS
%DataType: either point distances or biological conditions distances
%VARARGIN
%TRank: Several types of distances are stored in Tree
%By default point distances and biological conditions distances are stored 
%in T{1} and T[2} respectively. But another T position can be used for both types/

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

function [BiolIndex]=trs_notredundant(DataType,varargin)
global K P T



switch DataType
    case 'point distances' 
        if nargin==1
            TRank=1;
        else
            TRank=varargin{1};
        end
        BiolIndex=T{TRank}.biolIndex;
        P.net.biolIndexes{TRank}=zeros(P.biol.nb,1);
        cd(P.dir.data)
        load Tree
        PointDistances=squareform(T{1}.distances);
        BiolDistances=squareform(T{TRank}.distances);
        %construct distribution of intra biological conditions point distances        
        PairsDistances=[];
        for BiolL=1:length(BiolIndex)
            CurrBiol=BiolIndex(BiolL);
            PointRanks=P.biol.pairs{CurrBiol};            
            PairsDistances=[PairsDistances;PointDistances(PointRanks(1),PointRanks(2))];            
        end
        
        %construct distribution of intra-experiment and inter biological conditions point distance
        ExpDistances=[];
        for ExpL=1:P.exp.nb
            CurrBiols=unique(P.point.biolRank(find(P.point.expRank==ExpL)));
            if length(CurrBiols)>1
                for BiolL1=1:length(CurrBiols)-1
                    BiolPos1=find(BiolIndex==CurrBiols(BiolL1));
                    if ~isempty(BiolPos1)
                        PointRanks1=P.biol.pairs{CurrBiols(BiolL1)};
                        for BiolL2=BiolL1+1:length(CurrBiols)
                            BiolPos2=find(BiolIndex==CurrBiols(BiolL2));
                            if ~isempty(BiolPos2)
                                PointRanks2=P.biol.pairs{CurrBiols(BiolL2)};
                                ExpDistances=[ExpDistances;mean([PointDistances(PointRanks1(1),PointRanks2(1)),...
                                    PointDistances(PointRanks1(1),PointRanks2(2)),...
                                    PointDistances(PointRanks1(2),PointRanks2(1)),...
                                    PointDistances(PointRanks1(2),PointRanks2(2))])];
                            end
                        end
                    end
                end
            end
        end
                
        PairsDistances=round(PairsDistances(:));
        ExpDistances=round(ExpDistances(:));
        BiolDistances=round(BiolDistances(:));

        h=figure;
        set(gcf,'color',[1,1,1])
        set(h,'name',sprintf('biological conditions distances_%s_%s',P.project.name,date))
        hold on
        BinVal=histc(ExpDistances,[1:1:max(ExpDistances)]);
        plot([1:1:max(ExpDistances)],BinVal*100/sum(BinVal),'b')
        plot([1:1:max(ExpDistances)],BinVal*100/sum(BinVal),'r.')
        BinVal=histc(PairsDistances,[1:1:max(PairsDistances)]);
        plot([1:1:max(PairsDistances)],BinVal*100/sum(BinVal),'c')
        plot([1:1:max(PairsDistances)],BinVal*100/sum(BinVal),'m.','markersize',3)
        BinVal=histc(BiolDistances,[1:1:max(BiolDistances)]);
        plot([1:1:max(BiolDistances)],BinVal*100/sum(BinVal),'g')
        plot([1:1:max(BiolDistances)],BinVal*100/sum(BinVal),'k.','markersize',3)

        set(gca,'box','on')
        set(gca,'yticklabel','')
        xlabel('distance between two points')
        ylabel('frequency for intra-biological(cyan-magenta), intra-experiment (blue-red) and inter-experiment (black-green)')
        title(sprintf('biological conditions distances %s - %s',strrep(T{TRank}.selType,' ','-'),P.project.name))
        cd(P.dir.resTree)
        set_figsize('1024px')
        saveas(h,sprintf('biological_condition_distances_%s_plot_%s_%s.png',strrep(T{TRank}.selType,' ','-'),P.project.name,date),'png')

        %test several limits
        BiolDistances=squareform(T{TRank}.distances);
        Limits=[60:10:150];
        BiolPositions=cell(length(Limits),1);
        PointIndex=cell(length(Limits),1);
        for LimitL=1:length(Limits)
            CurrLimit=Limits(LimitL);
            for ExpL=1:length(P.exp.name)
                PointRanks=P.exp.pointIndex{ExpL};
                BiolRanks=unique(P.point.biolRank(PointRanks));
                if length(BiolRanks)==1
                    BiolPos=find(BiolIndex==BiolRanks);
                    if ~isempty(BiolPos)
                        BiolPositions{LimitL}=[BiolPositions{LimitL};BiolPos];
                        PointIndex{LimitL}=[PointIndex{LimitL};P.biol.pairs{BiolRanks}'];
                    end
                else
                    FindFlag=0;
                    Linked=zeros(length(BiolRanks));
                    LinkedNb=0;
                    for BiolL1=1:length(BiolRanks)-1
                        if ~isempty(find(BiolIndex==BiolRanks(BiolL1)))
                            PointRanks1=P.biol.pairs{BiolRanks(BiolL1)};
                            for BiolL2=BiolL1+1:length(BiolRanks)
                                if ~isempty(find(BiolIndex==BiolRanks(BiolL2)))
                                    PointRanks2=P.biol.pairs{BiolRanks(BiolL2)};
                                    CurrDistance=mean([PointDistances(PointRanks1(1),PointRanks2(1)),...
                                        PointDistances(PointRanks1(1),PointRanks2(2)),...
                                        PointDistances(PointRanks1(2),PointRanks2(1)),...
                                        PointDistances(PointRanks1(2),PointRanks2(2))]);
                                    if CurrDistance>CurrLimit
                                        Linked(BiolL1,BiolL2)=1;
                                        LinkedNb=LinkedNb+1;
                                    end
                                end
                            end
                        end


                        if LinkedNb>0
                            if length(Linked)<=8
                                Selected=maximalCliques(Linked);
                                if ~isempty(Selected)
                                    Selected=Selected{1};
                                else
                                    Selected=[];
                                end
                            else
                                cd(P.dir.resTree)
                                fid=fopen('cliquer.txt','w');
                                fprintf(fid,'c rank 1\n');
                                fprintf(fid,'p edge %u %u\n',length(Linked),sum(sum(Linked)));
                                for PsL=1:length(Linked)
                                    Pos=find(Linked(PsL,:));
                                    if ~isempty(Pos)
                                        for PosL=1:length(Pos)
                                            fprintf(fid,'e %u %u\n',PsL,Pos(PosL));
                                        end
                                    end
                                end
                                fclose(fid);
                                eval(sprintf('! /usr/local/cliquer/cl -su %s/cliquer.txt>%s/cliquer_res.txt;',P.dir.resTree,P.dir.resTree))
                                Selected=load('cliquer_res.txt');

                            end
                            if ~isempty(Selected)
                                FindFlag=1;
                                for BiolL=1:length(Selected)
                                    BiolPos=find(BiolIndex==BiolRanks(Selected(BiolL)));
                                    BiolPositions{LimitL}=[BiolPositions{LimitL};BiolPos];
                                    PointIndex{LimitL}=[PointIndex{LimitL};P.biol.pairs{BiolRanks(Selected(BiolL))}'];
                                end
                            end
                        end
                    end
                    if FindFlag==0
                        %take only one biologicalconditon in this experiment
                        for BiolL=1:length(BiolRanks)
                            BiolPos=find(BiolIndex==BiolRanks(BiolL));
                            if ~isempty(BiolPos)
                                BiolPositions{LimitL}=[BiolPositions{LimitL};BiolPos];
                                PointIndex{LimitL}=[PointIndex{LimitL};P.biol.pairs{BiolRanks(BiolL)}'];
                                break
                            end
                        end
                    end
                end
            end
            BiolPositions{LimitL}=unique(BiolPositions{LimitL});
            PointIndex{LimitL}=unique(PointIndex{LimitL});           
        end

        %select pairs used for constructing network
        List={'DIST   NB OF BIOL COND'};
        for LimitL=1:length(Limits)
            List{end+1,1}=sprintf('%u    %u',Limits(LimitL),length(BiolPositions{LimitL}));
        end
        Select=0;
        while Select<=1
            Select=listdlg('liststring',List,'selectionmode','single','initialvalue',2,'promptstring','select the minimal distance between biol cond used for network construction');
        end

        CurrBiolDistances=BiolDistances(BiolPositions{Select-1},BiolPositions{Select-1});
        if max(T{TRank}.distances)<=1
            Factor=1/60;
        else
            Factor=round(median(CurrBiolDistances(:))/50);
        end
        h=figure;
        set(h,'color',[1,1,1])
        set(h,'name',sprintf('%u DIFFERENT BIOL COND AT LIMIT>%u- %s',length(BiolPositions{Select-1}),Limits(Select-1),P.project.name))
        image(CurrBiolDistances./Factor)
        title(sprintf('%u DIFFERENT BIOL COND AT LIMIT>%u - %s',length(BiolPositions{Select-1}),Limits(Select-1),P.project.name))
        set_figsize('960px')
        cd(P.dir.resTree)
        saveas(h,sprintf('selectedbiol_pointsel_%s_%u_%s_%s.png',T{TRank}.selType,Limits(Select-1),P.project.name,date),'png')

        P.net.biolIndexes{TRank}(BiolIndex(BiolPositions{Select-1}))=1;
        cd(P.dir.project)
        eval(sprintf('save %s P',P.project.name))
        %points

    case {'biological conditions distances'}
        if nargin==1
            TRank=2;
        else
            TRank=varargin{1};
        end

        BiolIndex=T{TRank}.biolIndex;
        BiolDistances=squareform(T{TRank}.distances);
        P.net.biolIndexes{TRank}=zeros(P.biol.nb,1);
        cd(P.dir.data)
        load Tree
        PointDistances=squareform(T{1}.distances);
        %construct distribution of intra biological conditions point distances
        PairsDistances=[];
        for BiolL=1:length(BiolIndex)
            CurrBiol=BiolIndex(BiolL);
            PointRanks=P.biol.pairs{CurrBiol};
            PairsDistances=[PairsDistances;PointDistances(PointRanks(1),PointRanks(2))];
        end

        %construct distribution of intra-experiment and inter biological conditions point distance
        ExpDistances=[];
        for ExpL=1:P.exp.nb
            CurrBiols=unique(P.point.biolRank(find(P.point.expRank==ExpL)));
            if length(CurrBiols)>1
                for BiolL1=1:length(CurrBiols)-1
                    BiolPos1=find(BiolIndex==CurrBiols(BiolL1));
                    if ~isempty(BiolPos1)
                        for BiolL2=BiolL1+1:length(CurrBiols)
                            BiolPos2=find(BiolIndex==CurrBiols(BiolL2));
                            if ~isempty(BiolPos2)
                                ExpDistances=[ExpDistances;BiolDistances(BiolPos1,BiolPos2)];
                            end
                        end
                    end
                end
            end
        end

        PairsDistances=round(PairsDistances(:));
        ExpDistances=round(ExpDistances(:)*100);
        BiolDistances=round(BiolDistances(:)*100);

        h=figure;
        set(gcf,'color',[1,1,1])
        set(h,'name',sprintf('biological conditions distances_%s_%s',P.project.name,date))
        hold on
        BinVal=histc(ExpDistances,[1:1:max(ExpDistances)]);
        plot([1:1:max(ExpDistances)],BinVal*100/sum(BinVal),'b')
        plot([1:1:max(ExpDistances)],BinVal*100/sum(BinVal),'r.')
        BinVal=histc(PairsDistances,[1:1:max(PairsDistances)]);
        plot([1:1:max(PairsDistances)],BinVal*100/sum(BinVal),'c')
        plot([1:1:max(PairsDistances)],BinVal*100/sum(BinVal),'m.','markersize',3)
        BinVal=histc(BiolDistances,[1:1:max(BiolDistances)]);
        plot([1:1:max(BiolDistances)],BinVal*100/sum(BinVal),'g')
        plot([1:1:max(BiolDistances)],BinVal*100/sum(BinVal),'k.','markersize',3)

        set(gca,'box','on')
        set(gca,'yticklabel','')
        xlabel('distance between two biological conditions')
        ylabel('frequency for intra-biological(cyan-magenta), intra-experiment (blue-red) and inter-experiment (black-green)')
        title(sprintf('point and biological conditions distances - %s',P.project.name))
        cd(P.dir.resTree)
        set_figsize('1024px')
        saveas(h,sprintf('biological_condition_distances_%s_plot_%s_%s.png',strrep(T{TRank}.selType,' ','-'),P.project.name,date),'png')

        %test several limits
        BiolDistances=squareform(T{TRank}.distances*100);
        MedianDist=Median(round(T{TRank}.distances*10))*10;
        Limits=[MedianDist:10:90];
        BiolPositions=cell(length(Limits),1);
        PointIndex=cell(length(Limits),1);
        for LimitL=1:length(Limits)
            CurrLimit=Limits(LimitL);
            for ExpL=1:length(P.exp.name)
                PointRanks=P.exp.pointIndex{ExpL};
                BiolRanks=unique(P.point.biolRank(PointRanks));
                if length(BiolRanks)==1
                    BiolPos=find(BiolIndex==BiolRanks);
                    if ~isempty(BiolPos)
                        BiolPositions{LimitL}=[BiolPositions{LimitL};BiolPos];
                        PointIndex{LimitL}=[PointIndex{LimitL};P.biol.pairs{BiolRanks}'];
                    end
                else
                    FindFlag=0;
                    Linked=zeros(length(BiolRanks));
                    LinkedNb=0;
                    for BiolL1=1:length(BiolRanks)-1
                        BiolPos1=find(BiolIndex==BiolRanks(BiolL1));
                        if ~isempty(BiolPos1)
                            for BiolL2=BiolL1+1:length(BiolRanks)
                                BiolPos2=find(BiolIndex==BiolRanks(BiolL2));
                                if ~isempty(BiolPos2)
                                    CurrDistance=BiolDistances(BiolPos1,BiolPos2);
                                    if CurrDistance>CurrLimit
                                        Linked(BiolL1,BiolL2)=1;
                                        LinkedNb=LinkedNb+1;
                                    end
                                end
                            end
                        end
                        if LinkedNb>0
                            if length(Linked)<=8
                                Selected=maximalCliques(Linked);
                                if ~isempty(Selected)
                                    Selected=Selected{1};
                                else
                                    Selected=[];
                                end
                            else
                                cd(P.dir.resTree)
                                fid=fopen('cliquer.txt','w');
                                fprintf(fid,'c rank 1\n');
                                fprintf(fid,'p edge %u %u\n',length(Linked),sum(sum(Linked)));
                                for PsL=1:length(Linked)
                                    Pos=find(Linked(PsL,:));
                                    if ~isempty(Pos)
                                        for PosL=1:length(Pos)
                                            fprintf(fid,'e %u %u\n',PsL,Pos(PosL));
                                        end
                                    end
                                end
                                fclose(fid);
                                eval(sprintf('! /usr/local/cliquer/cl -su %s/cliquer.txt>%s/cliquer_res.txt;',P.dir.resTree,P.dir.resTree))
                                Selected=load('cliquer_res.txt');

                            end
                            if ~isempty(Selected)
                                FindFlag=1;
                                for BiolL=1:length(Selected)
                                    BiolPos=find(BiolIndex==BiolRanks(Selected(BiolL)));
                                    BiolPositions{LimitL}=[BiolPositions{LimitL};BiolPos];
                                    PointIndex{LimitL}=[PointIndex{LimitL};P.biol.pairs{BiolRanks(Selected(BiolL))}'];
                                end
                            end
                        end
                    end
                    if FindFlag==0
                        %take only one biologicalconditon in this experiment
                        for BiolL=1:length(BiolRanks)
                            BiolPos=find(BiolIndex==BiolRanks(BiolL));
                            if ~isempty(BiolPos)
                                BiolPositions{LimitL}=[BiolPositions{LimitL};BiolPos];
                                PointIndex{LimitL}=[PointIndex{LimitL};P.biol.pairs{BiolRanks(BiolL)}'];
                                break
                            end
                        end
                    end
                end
            end
            BiolPositions{LimitL}=unique(BiolPositions{LimitL});
            PointIndex{LimitL}=unique(PointIndex{LimitL});
        end

        %select pairs used for constructing network
        List={'DIST   NB OF BIOL COND'};
        for LimitL=1:length(Limits)
            List{end+1,1}=sprintf('%u    %u',Limits(LimitL),length(BiolPositions{LimitL}));
        end
        Select=0;
        while Select<=1
            Select=listdlg('liststring',List,'selectionmode','single','initialvalue',2,'promptstring','select the minimal distance between biol cond used for network construction');
        end

        CurrBiolDistances=BiolDistances(BiolPositions{Select-1},BiolPositions{Select-1});
        if max(T{TRank}.distances<=1)
            Factor=1/60;
        else
            Factor=round(median(CurrBiolDistances(:))/50);
        end
        h=figure;
        set(h,'color',[1,1,1])
        set(h,'name',sprintf('%u DIFFERENT BIOL COND AT LIMIT>%u- %s',length(BiolPositions{Select-1}),Limits(Select-1),P.project.name))
        image(CurrBiolDistances./Factor)
        title(sprintf('%u DIFFERENT BIOL COND AT LIMIT>%u - %s',length(BiolPositions{Select-1}),Limits(Select-1),P.project.name))
        set_figsize('960px')
        cd(P.dir.resTree)
        saveas(h,sprintf('selectedbiol_biolsel_%s,%u_%s_%s.png',T{TRank}.selType,Limits(Select-1),P.project.name,date),'png')

        P.net.biolIndexes{TRank}(BiolIndex(BiolPositions{Select-1}))=1;
        cd(P.dir.project)
        eval(sprintf('save %s P',P.project.name))
        %points
       
    case 'variation distances'

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
                [SelPos,Ok]=listdlg('promptstring','select the biol cond not to be used','liststring',P.biol.nr{1}.name,'selectionmode','multiple');
                if Ok==1
                    ClearBiol=[];
                    for SelL=1:length(SelPos)
                        ClearBiol=[ClearBiol;P.biol.nr{1}.scoupleindex{SelPos(SelL)}];
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
                    NrName=inputdlg({'give the name for this set of not redundant conditions'},'',1,{sprintf('all on %s (filtered => %s)',Type,Filter{2})});
                    NrName=NrName{1};
                else
                    NrName=sprintf('all on %s',Type);
                end
            elseif isequal(AnalType,'Complement')
                if isequal(FilterIt,'yes')
                    NrName=inputdlg({'give the name for this set of not redundant conditions'},'',1,{sprintf('all wo nr nb ? on %s (filtered => %s)',Type,Filter{2})});
                    NrName=NrName{1};
                else
                    NrName=inputdlg({'give the name for this set of not redundant conditions'},'',1,{sprintf('all wo nr nb ? on %s',Type)});
                    NrName=NrName{1};
                end
            else
                if isequal(FilterIt,'yes')
                    NrName=inputdlg({'give the name for this set of not redundant conditions'},'',1,{sprintf('%s on %s (filtered => %s)',SelFile,Type,Filter{2})});
                    NrName=NrName{1};
                else
                    NrName=inputdlg({'give the name for this set of not redundant conditions'},'',1,{sprintf('%s on %s',SelFile,Type)});
                    NrName=NrName{1};
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
                hist(P.biol.scouplecorr{1},1000)
                set(gca,'xlim',[0.8,1])

                subplot(1,3,2)
                hist(P.biol.scouplecorr{1},1000)
                set(gca,'xlim',[0.88,0.92])
                subplot(1,3,3)
                hist(P.biol.scouplecorr{1},1000)
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
                if P.biol.softgrp{1}(BiolPos)~=1 & ~isequal(P.biol.name{BiolRank},'MedianBiol') &P.biol.scouplecorr{1}(BiolPos)<=SupCorrLimit&P.biol.scouplecorr{1}(BiolPos)>=InfCorrLimit
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
            K.dir.restree=fullfile(K.dir.sosresult,sprintf('m%03u',ModelRank))
            try
                cd(sprintf('m%03u',ModelRank))
            catch
                mkdir(sprintf('m%03u',ModelRank))
                cd(sprintf('m%03u',ModelRank))
            end
            if isequal(P.flag.station,'windows')
                saveas(h0,sprintf('%s_values_m%03u_%s',strrep(NrName,' ','_'),ModelRank,date),'bmp')
            else
                saveas(h0,sprintf('%s_values_m%03u_%s',strrep(NrName,' ','_'),ModelRank,date),'png')
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
                        saveas(h1,sprintf('%s_dist%u_m%03u_%s',strrep(NrName,' ','_'),round(100*DistLimit(LimitRank)),ModelRank,date),'bmp')
                    else
                        saveas(h1,sprintf('%s_dist%u_m%03u_%s',strrep(NrName,' ','_'),round(100*DistLimit(LimitRank)),ModelRank,date),'png')
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
                                        saveas(h1,sprintf('%s_dist%u_%s_%s_%s_nb%u',strrep(NrName,' ','_'),round(100*DistLimit(LimitRank)),Filter{2},lower(P.chip.name),date,Rank),'bmp')
                                    else
                                        saveas(h1,sprintf('%s_dist%u_%s_%s_nb%u',strrep(NrName,' ','_'),round(100*DistLimit(LimitRank)),lower(P.chip.name),date,Rank),'bmp')
                                    end
                                else
                                    if isequal(FilterIt,'no')
                                        saveas(h1,sprintf('%s_dist%u_%s_%s_nb%u',strrep(NrName,' ','_'),round(100*DistLimit(LimitRank)),lower(P.chip.name),date,Rank),'png')
                                    else
                                        saveas(h1,sprintf('%s_dist%u_%s_%s_%s_nb%u',strrep(NrName,' ','_'),round(100*DistLimit(LimitRank)),Filter{2},lower(P.chip.name),date,Rank),'png')
                                    end
                                end
                            end
                        end
                    end
                end
                NewSelection=questdlg(sprintf('yes = try another distance limit\n(%.02f); no = save results',DistLimit(LimitRank)),'','yes','no','no');
            end



            %BiolRank=P.biol.pointIndex(UsedSCoupleRank(SList));
            Continue=0;
            ReplaceIt='no';


            if ~isfield(P.biol,'nr')
                for ChipL=1:P.chip.nb
                    P.biol.nr{ChipL}.name={};
                    P.biol.nr{ChipL}.startbiol={};
                    P.biol.nr{ChipL}.name={};
                    P.biol.nr{ChipL}.type={};
                    P.biol.nr{ChipL}.corrlimit={};
                    P.biol.nr{ChipL}.scoupleindex={};
                end
            end
            NrPos=length(P.biol.nr{1}.name)+1;
            if isequal(AnalType,'All')
                NrName=sprintf('all on %s at dist > %.2f',Type,DistLimit(LimitRank));
            elseif isequal(AnalType,'Complement')
                NrName=inputdlg({'give the name for this set of not redundant conditions'},'',1,{sprintf('all wo nr nb 1 on %s at dist > %.2f',Type,DistLimit(LimitRank))});
                NrName=NrName{1};
            else
                NrName=inputdlg({'give the name for this set of not redundant conditions'},'',1,{sprintf('%s on %s at dist > %.2f',SelFile,Type,DistLimit(LimitRank))});
                NrName=NrName{1};
            end
            NamePos=strmatch(NrName,P.biol.nr{1}.name,'exact');
            Correct=0;
            while Correct==0
                if ~isempty(NamePos)
                    ReplaceIt=questdlg(sprintf('%s exist. Do you want to replace it ?',NrName),'','yes','no','no');
                    if isequal(ReplaceIt,'yes')
                        NrPos=NamePos;
                        Correct=1;
                    else
                        NrName=inputdlg({'give the name for this set of not redundant conditions'},'',1,{''});
                        NrName=NrName{1};
                        NamePos=strmatch(NrName,P.biol.nr{1}.name,'exact');
                    end
                else
                    Correct=1;
                end
            end

            UsedBiolRankMem=UsedBiolRank;
            UsedBiolRank=UsedBiolRank(SList);
            Seed=0;
            BiolNb=length(UsedBiolRank);
            if isequal(ReplaceIt,'no')
                LessBiol=questdlg(sprintf('there exist %u selected biol cond.\n Do you want less ?',length(UsedBiolRank)),'','yes','no','no');
                if isequal(LessBiol,'yes')
                    BiolNb=inputdlg('How many biol cond ?','',1,{'32'});
                    BiolNb=str2num(BiolNb{1});
                    Seed=inputdlg('Seed of algo (twister != 0)','',1,{'5489'});
                    Seed=str2num(Seed{1});
                    if Seed==0
                        h=errordlg('seed must not = 0. Process canceled!')
                        error('process canceled')
                    end
                    rand('twister',Seed);
                    if BiolNb<length(UsedBiolRank)
                        if floor(length(UsedBiolRank)/BiolNb)>1
                            SeveralGrp=questdlg(sprintf('It is possible to make %u groups\nof\n %u biol cond. Do you want that ?',floor(length(UsedBiolRank)/BiolNb),BiolNb),'','yes','no','yes');
                            if isequal(SeveralGrp,'no')
                                RandPos=randperm(length(UsedBiolRank));
                                BiolRanks{1}=UsedBiolRank(RandPos(1:BiolNb));
                            else
                                if mod(length(UsedBiolRank),BiolNb)>0
                                    LastGrp=questdlg(sprintf('Last group has %u condtions.\nDo you want to use it ?',mod(length(UsedBiolRank),BiolNb)),'','yes','no','no');
                                else
                                    LastGrp='no';
                                end
                                for i=1:floor(length(UsedBiolRank)/BiolNb)
                                    RandPos=randperm(length(UsedBiolRank));
                                    BiolRanks{i}=UsedBiolRank(RandPos(1:BiolNb));
                                    UsedBiolRank(RandPos(1:BiolNb))=[];
                                end
                                if isequal(LastGrp,'yes')
                                    BiolRanks{i+1}=UsedBiolRank;
                                end
                            end
                        else
                            RandPos=randperm(length(UsedBiolRank));
                            BiolRanks{1}=UsedBiolRank(RandPos(1:BiolNb));
                        end
                    end
                else
                    BiolRanks{1}=UsedBiolRank;
                end
            else
                BiolRanks{1}=UsedBiolRank;
            end
            %P.biol.nr{1}.startbiol{1}=P.biol.pointIndex(UsedSCoupleRank);
            for i=1:length(BiolRanks)
                P.biol.nr{1}.startbiol{NrPos,1}=UsedBiolRankMem;
                if Seed==0
                    P.biol.nr{1}.name{NrPos,1}=sprintf('%s %u %01u',NrName,length(BiolRanks{i}),i);
                else
                    P.biol.nr{1}.name{NrPos,1}=sprintf('%s %u %01u twister_seed%u',NrName,length(BiolRanks{i}),i,Seed);
                end
                P.biol.nr{1}.type{NrPos,1}=Type;
                P.biol.nr{1}.corrlimit{NrPos,1}=[InfCorrLimit,SupCorrLimit];
                P.biol.nr{1}.scoupleindex{NrPos,1}=sort(BiolRanks{i});
                if isequal(FilterIt,'yes')
                    P.biol.nr{1}.filter{NrPos,1}=Filter{1};
                end
                NrPos=NrPos+1;
            end
            cd(K.dir.point)
            eval(sprintf('save %s P -v7',P.file.point))
            Continue=questdlg('Do you want to do another selection ?','','yes','no','yes');
        end
end