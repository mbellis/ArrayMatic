%=====================
% FUNCTION GEO_METADB
%=====================

% GEO_METADB allows to display structure of Geometadb.sqlite database
% loaded from Metzer lab at http://gbnci.abcc.ncifcrf.gov/geo/ and to display particular GPL extracted 
% from this database 
% 
%INPUT PARAMETERS
% 1 Action:
%       open geo metadb
%           open GEOmetadb.sqlite located in K.dir.geoMetadata
%       close
%           close GEO metadb
%       view tables and fields
%           view tables and fields existing in GEOmetadb.sqlite
%       write exemplar records
%       view technology
%       view species
%       display a GPL
%       print a GPL
%       import GSE
%       modif version
%       find biological conditions
%       import biological conditions
%
% FUNCTIONS
%
% VERIF
%
% PRINT

% EXTERNAL TOOL
% uses mksqlite.m developped by by Martin Kortmann <mail@kortmann.de>
% (sources downloaded at developer.berlios.de/projects/mksqlite/)

%geo_metadb('display a GPL',[3533,11095,11097,201,3921,4910,4911,4912,4913,4914,4915,4916,5082,5175,5188,570,571,6244,80,8300,91,96,10740,11096,11098,32,339,5811,1261,6096,6193,75,81,8321,1355,85,341,6247,88,9199,6543,6194],1)
%geo_metadb('display a GPL',[96,10740,11096,11098,32,339,5811,1261,6096,6193,75,8321,1355,85,341,6247,88,9199,6543,6194],1)
%geo_metadb('display a GPL',[82,86,87],1)

%geo_metadb('find biological conditions','GPL198',{'ATH1','(A|a)rabidopsis','(A.)?thaliana'})
%geo_metadb('GSE dictionary',[11095,11097,201,3921,4910,4911,4912,4913,4914,4915,4916,5082,5175,5188,570,571,6244,80,8300,91,96,10740,11096,11098,32,339,5811,1261,6096,6193,75,81,8321,1355,85,341,6247,88,9199,6543,6194])
%geo_metadb('GSE dictionary',[11095,11097,201,3921,4910,4911,4912,4913,4914,4915,4916,5082,5175,5188,570,571,6244,80,8300,91,96,10740,11096,11098,32,339,5811,1261,6096,6193,75,81,8321,1355,85,341,6247,88,9199,6543,6194],'([Bb]rains? ?)|([Nn]eur\w* ?)|([Nn]erv\w* ?)|[Cc]hannel\w* ?)|([Pp]urkinje ?)|([Cc]erebell\w+ ?)|([Cc]erebr\w+ ?)|([Aa]tax\w+ ?)|(polyglu\w+ ?)','neuro')
%geo_metadb('GSE dictionary',[11095,11097,201,3921,4910,4911,4912,4913,4914,4915,4916,5082,5175,5188,570,571,6244,80,8300,91,96,10740,11096,11098,32,339,5811,1261,6096,6193,75,81,8321,1355,85,341,6247,88,9199,6543,6194],'[Pp]urkinje ?','purkinje')
%geo_metadb('GSE dictionary',[11095,11097,201,3921,4910,4911,4912,4913,4914,4915,4916,5082,5175,5188,570,571,6244,80,8300,91,96,10740,11096,11098,32,339,5811,1261,6096,6193,75,81,8321,1355,85,341,6247,88,9199,6543,6194],'([Tt]h[12])|([Tt]h17)|([Tt]hreg)|([Tt][Cc][Rr])','th')
%geo_metadb('GSE dictionary',[3533,11095,11097,201,3921,4910,4911,4912,4913,4914,4915,4916,5082,5175,5188,570,571,6244,80,8300,91,96,10740,11096,11098,32,339,5811,1261,6096,6193,75,81,8321,1355,85,341,6247,88,9199,6543,6194],'[Aa]xotomy|DRG\w*|(([Aa]xon\w* ?|[Nn]eur\w* ?|[Nn]erv\w* ?)([Oo]utgrowth ?|[Rr]egeneration ?|[Rr]egrowth ?))','nerve_regeneration')
%geo_metadb('GSE dictionary',[3533,11095,11097,201,3921,4910,4911,4912,4913,4914,4915,4916,5082,5175,5188,570,571,6244,80,8300,91,96,10740,11096,11098,32,339,5811,1261,6096,6193,75,81,8321,1355,85,341,6247,88,9199,6543,6194],'([Tt]h[12])|([Tt]h17)|([Tt]hreg)|([Tt][Cc][Rr])','th')
%geo_metadb('GSE dictionary',[11095,11097,201,3921,4910,4911,4912,4913,4914,4915,4916,5082,5175,5188,570,571,6244,80,8300,91,96,10740,11096,11098,32,339,5811,1261,6096,6193,75,81,8321,1355,85,341,6247,88,9199,6543,6194],'DRG','drg')
%geo_metadb('GSE special symbols',[11095,11097,201,3921,4910,4911,4912,4913,4914,4915,4916,5082,5175,5188,570,571,6244,80,8300,91,96,10740,11096,11098,32,339,5811,1261,6096,6193,75,81,8321,1355,85,341,6247,88,9199,6543,6194],'[\|\$\[\]\)(}{_''"#`~&^%*?<>?@=???]',1,1,30,30)
%geo_metadb('GSE special symbols',[11095,11097,201,3921,4910,4911,4912,4913,4914,4915,4916,5082,5175,5188,570,571,6244,80,8300,91,96,10740,11096,11098,32,339,5811,1261,6096,6193,75,81,8321,1355,85,341,6247,88,9199,6543,6194],'?',1,1,30,30)
%geo_metadb('GSE special symbols',[339],'?',1,1,30,30)
%geo_metadb('GSE special symbols',[339],sprintf('%c',char(65533)),1,1,30,30)
%geo_metadb('GSE special symbols',[11095,11097,201,3921,4910,4911,4912,4913,4914,4915,4916,5082,5175,5188,570,571,6244,80,8300,91,96,10740,11096,11098,32,339,5811,1261,6096,6193,75,81,8321,1355,85,341,6247,88,9199,6543,6194],'\S\.\S|\S;\S|\S:\S|\S,\S|\S!\S|\S\?\S|',3,2,30,30)
%geo_metadb('print GSEs',[11095,11097,201,3921,4910,4911,4912,4913,4914,4915,4916,5082,5175,5188,570,571,6244,80,8300,91,96,10740,11096,11098,32,339,5811,1261,6096,6193,75,81,8321,1355,85,341,6247,88,9199,6543,6194],'([Tt]h[12])|([Tt]h17)|([Tt]hreg)|([Tt][Cc][Rr])','Th_cells')
%geo_metadb('print GSEs',[81,85,88,339,570,1261,1355],'DRG\w*|(([Nn]eur\w* ?|[Nn]erv\w* ?)([Rr]egeneration ?|[Rr]egrowth ?))','nerve_regeneration')
%geo_metadb('print GSEs',[81,85,88,96,339,341,570,1261,1355,6193,6247],'[Aa]xotomy|DRG\w*|(([Aa]xon\w* ?|[Nn]eur\w* ?|[Nn]erv\w* ?)([Oo]utgrowth ?|[Rr]egeneration ?|[Rr]egrowth ?))','nerve_regeneration')
%geo_metadb('print GSEs',[81,85,91,96,201,339,570,571,1261,1355,3533,3921,5811,6244,8300,8321,11096,11098],'([Tt]h[12])|([Tt]h17)|([Tt]hreg)|([Tt][Cc][Rr])','th')

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


function varargout=geo_metadb(Action,varargin)
global K

%first version used '30-Jun-2010'
%snd version use VERSION='17-Jul-2010';
%third VERSION='10-Dec-2011';
VERSION='14-Jul-2012';
switch Action

    case 'open geo metadb'
%% OPEN GEO METADB
        cd(K.dir.geoMetadata)
        K.tmp.geoDbId=mksqlite(0,'open','GEOmetadb.sqlite');

    case 'close geo metadb'
%% CLOSE GEO METADB
        Ok=VERIF();
        if Ok
            cd(K.dir.geoMetadata)
            mksqlite(K.tmp.geoDbId,'close');
            K.tmp.geoDbId=0;
        end


    case 'view tables and fields'
%% FIELDS
        Ok=VERIF();
        if Ok
            List = mksqlite( K.tmp.geoDbId,'SELECT name FROM sqlite_master WHERE type=''table'' ORDER BY ''name''');
            Table={};
            for ListL=1:length(List)
                Table{ListL,1}=List(ListL).name;
            end
            Sel=listdlg('liststring',Table,'selectionmode','single');
            while Sel
                eval(sprintf('Fields = mksqlite(K.tmp.geoDbId,''PRAGMA TABLE_INFO(%s)'');',Table{Sel}))
                if ~isempty(Fields)
                    FieldList={};
                    FieldList{end+1}=sprintf('Fields of %s',Table{Sel});
                    FieldList{end+1}=' No.  | Name                      | Type';
                    FieldList{end+1}='------+---------------------------+----------------------';

                    for FieldL = 1:size( Fields, 1 )
                        FieldList{end+1}=sprintf( '  %2u  | %-25s | %s \n', FieldL, Fields(FieldL).name, Fields(FieldL).type);
                    end
                    listdlg('liststring',FieldList,'selectionmode','single','listsize',[300,500]);
                end
                Sel=listdlg('liststring',Table,'selectionmode','single');
            end
        end

    case 'write exemplar records'
%% EXAMPLE
        Gpl=varargin{1};
        Gse=varargin{2};
        cd(K.dir.geoMetadata)
        GseFid=fopen(sprintf('tables4%s.txt',Gse),'w');
        fprintf(GseFid,'\n*******************************************\t')
        %TABLE GPL
        fprintf(GseFid,'TABLE gpl\n');
        eval(sprintf('Selection=mksqlite(''SELECT * FROM gpl WHERE gpl=''''%s'''''')',Gpl))
        PRINT(GseFid,Selection)

        %TABLE GSE
        fprintf(GseFid,'TABLE gse\n');
        eval(sprintf('Selection=mksqlite(''SELECT * FROM gse WHERE gse=''''%s'''''')',Gse))
        PRINT(GseFid,Selection)

        %TABLE GSE_GSM
        fprintf(GseFid,'TABLE gse_gsm\n');
        eval(sprintf('Selection=mksqlite(''SELECT * FROM gse_gsm WHERE gse=''''%s'''''')',Gse))
        PRINT(GseFid,Selection)

        %TABLE GSM with first GSM of the previous selection
        fprintf(GseFid,'TABLE gsm\n');
        eval(sprintf('Selection=mksqlite(''SELECT * FROM gsm WHERE gsm=''''%s'''''')',Selection(1).gsm))
        PRINT(GseFid,Selection)

        %TABLE GDS of THE CURRENT GSE if posssible otherwise the last one
        eval(sprintf('Selection=mksqlite(''SELECT * FROM gds WHERE gpl=''''%s'''''')',Gpl))
        for SelL=1:length(Selection)
            if isequal(Selection(SelL).gse,Gse)
                break
            end
        end
        fprintf(GseFid,'TABLE gds\n');
        eval(sprintf('Selection=mksqlite(''SELECT * FROM gds WHERE gpl=''''%s'''' AND gse=''''%s'''''')',Gpl,Selection(SelL).gse))
        PRINT(GseFid,Selection)

        %TABLE GSE_GPL
        fprintf(GseFid,'TABLE gse_gpl\n');
        eval(sprintf('Selection=mksqlite(''SELECT * FROM gse_gpl WHERE gpl=''''%s'''' LIMIT 1'')',Gpl))
        PRINT(GseFid,Selection)

        %TABLE GDS_SUBSET (first record)
        fprintf(GseFid,'TABLE gds_subset\n');
        Selection=mksqlite('SELECT * FROM gds_subset LIMIT 1');
        PRINT(GseFid,Selection)

        %TABLE SMATRIX (first record)
        fprintf(GseFid,'TABLE sMatrix\n');
        eval(sprintf('Selection=mksqlite(''SELECT * FROM sMatrix WHERE gpl=''''%s'''''')',Gpl))
        PRINT(GseFid,Selection)

        %TABLE GEODB_COLUMN_DESC
        fprintf(GseFid,'TABLE geodb_column_desc\n');
        Selection=mksqlite('SELECT * FROM geodb_column_desc LIMIT 1');
        PRINT(GseFid,Selection)

        %TABLE GEOCONVERT
        fprintf(GseFid,'TABLE geoConvert\n');
        Selection=mksqlite('SELECT * FROM geoConvert LIMIT 1');
        PRINT(GseFid,Selection)

        %TABLE METAINFO
        fprintf(GseFid,'TABLE metaInfo\n');
        Selection=mksqlite('SELECT * FROM metaInfo LIMIT 1');
        PRINT(GseFid,Selection)
        fclose(GseFid)

    case 'view technologies'
%% TECHNOLOGY
        DisplayFlag=1;
        if nargin==2
            DisplayFlag=varargin{1};
        end
        VERIF()
        TechSel=mksqlite('SELECT technology FROM GPL');
        Technologies={};
        for SelL=1:length(TechSel)
            if ~isempty(TechSel(SelL).technology)
                Technologies{end+1,1}=TechSel(SelL).technology;
            end
        end
        Technologies=unique(Technologies);
        if DisplayFlag
            Temp=listdlg('liststring',Technologies,'PromptString','Available Technologies','ListSize',[500,600]);
        end
        if nargout==1
            varargout{1}=Technologies;
        end
        geo_metadb('close geo metadb');

    case 'view species'
%% SPECIES
        Ok=VERIF();
        if Ok
            %FIND EXISTING SPECIES IN GPL
            Species={};
            Sel=mksqlite('SELECT organism FROM gpl');
            Rank=0;
            for SpecL=1:length(Sel)
                CurrSpecies=Sel(SpecL).organism;
                if ~isempty(CurrSpecies)
                    Rank=Rank+1;
                    Species{Rank,1}=CurrSpecies;
                end
            end
            Species=unique(Species);

            %SELECT ONE OR SEVERAL SPECIES
            SpeciesSel=listdlg('liststring',Species,'PromptString','select one or several species','ListSize',[500,600]);
            if ~isempty(SpeciesSel)
                %PROCESS EACH SELECTED SPECIES

                Manufacturer={};
                Technology={};
                Model={};
                GseNb=[];
                for SelL=1:length(SpeciesSel)

                    %RECOVER THE GPL TARGETING THIS SPECIES
                    eval(sprintf('Gpls=mksqlite(''SELECT gpl,title,manufacturer,technology FROM gpl WHERE organism=''''%s'''''');',Species{SpeciesSel(SelL)}))            %process each Gpl of the current species

                    %RECOVER INFORMATION FOR EACH GPL

                    for GplL=1:length(Gpls)
                        GseDates{GplL,1}={};
                        GseNb{GplL,1}=[];
                        GsmNb{GplL,1}=[];
                        Model{GplL,1}=Gpls(GplL).title;
                        Manufacturer{GplL,1}=Gpls(GplL).manufacturer;
                        Technology{GplL,1}=Gpls(GplL).technology;

                        %RECOVER THE GSE MADE WITH THE CURRENT GPL
                        eval(sprintf('Gses=mksqlite(''SELECT gse FROM gse_gpl WHERE gpl=''''%s'''''');',Gpls(GplL).gpl))

                        %proceed each Gse of the current Gpl
                        for GseL=1:length(Gses)

                            %RECOVER SUBMISSION DATE
                            eval(sprintf('Date=mksqlite(''SELECT submission_date FROM gse WHERE gse=''''%s'''''');',Gses(GseL).gse))

                            %RECOVER THE NUMBER OF GSM FOR THE CURRENT GSE
                            eval(sprintf('Gsms=mksqlite(''SELECT gsm FROM gse_gsm WHERE gse=''''%s'''''');',Gses(GseL).gse))

                            DatePos=strmatch(Date.submission_date,GseDates{GplL},'exact');
                            if isempty(DatePos)
                                GseDates{GplL,1}{end+1,1}=Date.submission_date;
                                GseNb{GplL,1}(end+1,1)=1;
                                GsmNb{GplL,1}(end+1,1)=length(Gsms);
                            else
                                GseNb{GplL,1}(DatePos)=GseNb{GplL}(DatePos)+1;
                                GsmNb{GplL,1}(DatePos)=GsmNb{GplL}(DatePos)+length(Gsms);
                            end
                        end
                    end

                    %put a string in  missing Manufacturer
                    for ManL=1:length(Manufacturer)
                        if isempty(Manufacturer{ManL})
                            Manufacturer{ManL}='-';
                        end
                    end

                    %PROCESS MANUFACTURERS

                    Manufacturers=unique(Manufacturer);
                    %recover existing technologies in GEO
                    Technologies=geo_metadb('view technologies',0);

                    %recover the number of Gse and Gsm for each Manufacturer
                    ManGseNb=zeros(length(Manufacturers),1);
                    ManGsmNb=zeros(length(Manufacturers),1);
                    for ManL=1:length(Manufacturers)

                        %RECOVER THE GPL OF THE CURRENT MANUFACTURER
                        if isequal(Manufacturers{ManL},'-')
                            eval(sprintf('Gpls=mksqlite(''SELECT gpl FROM gpl WHERE manufacturer is NULL AND organism="%s"'');',strrep(strrep(Manufacturers{ManL},'''',' '),'"',''),Species{SpeciesSel(SelL)}))
                        else                     
                            eval(sprintf('Gpls=mksqlite(''SELECT gpl FROM gpl WHERE manufacturer="%s" AND organism="%s"'');',strrep(strrep(Manufacturers{ManL},'''',' '),'"',''),Species{SpeciesSel(SelL)}))                            
                        end

                        for GplL=1:length(Gpls)

                            %RECOVER THE GSE OF THE CURRENT GPL
                            eval(sprintf('Gses=mksqlite(''SELECT gse FROM gse_gpl WHERE gpl=''''%s'''''');',Gpls(GplL).gpl))
                            ManGseNb(ManL)=ManGseNb(ManL)+ length(Gses);

                            if length(Gses)>0
                                for GseL=1:length(Gses)

                                    %RECOVER THE GSM OF THE CURRENT GSE
                                    eval(sprintf('Gsms=mksqlite(''SELECT gsm FROM gse_gsm WHERE gse=''''%s'''''');',Gses(GseL).gse))

                                    ManGsmNb(ManL)=ManGsmNb(ManL)+ length(Gsms);
                                end
                            end
                        end
                    end

                    %sort over the number of GSM
                    [ManGsmNb SortOrder]=sort(ManGsmNb,'descend');
                    ManGseNb=ManGseNb(SortOrder);
                    Manufacturers=Manufacturers(SortOrder);
                    if sum(ManGseNb)>0 & sum(ManGsmNb)>0

                        %construct the displayed list of manufacturers
                        ManDisp=cell(length(Manufacturers),1);
                        for ManL=1:length(Manufacturers)
                            ManDisp{ManL}=sprintf('%s (%u GSE - %u GSM)',Manufacturers{ManL},ManGseNb(ManL),ManGsmNb(ManL));
                        end

                        %print Manufacturer information
                        fid=fopen(sprintf('%s_gsenb_%s.txt',strrep(Species{SpeciesSel(SelL)},' ','_'),VERSION),'w');
                        fprintf(fid,'Manufacturer\tGSE_NB\tGMS_NB\n')
                        for ManL=1:length(Manufacturers)
                            fprintf(fid,'%s\t%u\t%u\n',Manufacturers{ManL},ManGseNb(ManL),ManGsmNb(ManL));
                        end
                        fclose(fid)

                        %display information for choosing which keywords use to select
                        %manufactorers that will be processed individually (all others
                        %will be gathered and names 'other')
                        ManSel=listdlg('liststring',ManDisp,'PromptString','Write down the Keywords you want to use to select manufacturers','ListSize',[500,600]);
                        ManKeywords=inputdlg('enter list of the Keywords /nyou wrote down with this format/n  {'','',''}','',1);
                        ManKeywords=ManKeywords{1};
                        eval(sprintf('ManKeywords=%s;',ManKeywords))
                        %find the manufacturers corresponding to the selected keywords
                        ManPos={};
                        FoundPos=[];
                        for KeyL=1:length(ManKeywords)
                            ManPos{KeyL}=[];
                            for ManL=1:length(Manufacturers)
                                if ~isempty(findstr(upper(ManKeywords{KeyL}),upper(Manufacturers{ManL})))
                                    ManPos{KeyL}=[ManPos{KeyL};ManL];
                                    FoundPos=[FoundPos;ManL];
                                end
                            end
                        end
                        %add a category with gathering not selected manufacturers
                        NotFoundPos=setdiff(1:length(Manufacturers),FoundPos);
                        ManKeywords{end+1}='others';
                        ManPos{end+1}=NotFoundPos;

                        %RECOVER AND DISPLAY INFORMATIONS FOR THE SELECTED MANUFACTURERS

                        TechNb=length(Technologies);
                        %existing technologies
                        ExisTech=zeros(length(Technologies),1);
                        for TechL=1:TechNb
                            for ManL=1:length(ManKeywords)
                                for PosL=1:length(ManPos{ManL})
                                    ManIndex=strmatch(Manufacturers{ManPos{ManL}(PosL)},Manufacturer,'exact');
                                    for IndexL=1:length(ManIndex)
                                        if isequal(Technology{ManIndex(IndexL)},Technologies{TechL})
                                            if sum(GseNb{ManIndex(IndexL)})>0
                                                ExistTech(TechL)=1;
                                                break
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        ExistTech=find(ExistTech);
                        TechNb=length(ExistTech);
                        XPlot=ceil(sqrt(TechNb));
                        YPlot=floor(TechNb/XPlot);
                        if XPlot*YPlot==TechNb
                            XPlot=ceil(sqrt(TechNb+1));
                            YPlot=floor((TechNb+1)/XPlot);
                        end
                        if XPlot*YPlot<TechNb
                            YPlot=YPlot+1;
                        end

                        h=figure;
                        set(h,'name',Species{SpeciesSel(SelL)})
                        ManGseNbs=zeros(length(ManKeywords),1);
                        ManGsmNbs=zeros(length(ManKeywords),1);
                        %                     CurrColor=colormap;
                        %                     ColorNb=size(CurrColor,1);
                        %                     Colors=zeros(TechNb,3);
                        %                     ColorStep=floor(ColorNb/TechNb);
                        %                     for TechL=1:TechNb
                        %                         if mod((TechL-1)*ColorStep+1,ColorNb)
                        %                             Colors(TechL,:)=CurrColor(mod((TechL-1)*ColorStep+1,ColorNb),:);
                        %                         else
                        %                             Colors(TechL,:)=CurrColor(end,:);
                        %                         end
                        %                     end
                        Colors=colors(colormap,length(ManKeywords));
                        for TechL=1:TechNb
                            subplot(YPlot,XPlot,TechL)
                            hold on
                            %                Colors='mrygcbkmrygcbk';
                            MinDate=now;
                            MaxDate=0;
                            for ManL=1:length(ManKeywords)
                                CurrDates=[];
                                CurrGseNb=[];
                                CurrGsmNb=[];
                                %find gse and gsm fot the current manufacturer of the
                                %current technology (if exists)
                                for PosL=1:length(ManPos{ManL})
                                    ManIndex=strmatch(Manufacturers{ManPos{ManL}(PosL)},Manufacturer,'exact');
                                    for IndexL=1:length(ManIndex)
                                        if isequal(Technology{ManIndex(IndexL)},Technologies{ExistTech(TechL)})
                                            CurrGseNb=[CurrGseNb;GseNb{ManIndex(IndexL)}];
                                            CurrGsmNb=[CurrGsmNb;GsmNb{ManIndex(IndexL)}];
                                            CurrDates=[CurrDates;datenum(datevec(GseDates{ManIndex(IndexL)},'yyyy-mm-dd'))];
                                        end
                                    end
                                end
                                %plot current manufacturer of the current technology
                                if ~isempty(CurrGseNb)
                                    MinDate=min(MinDate,min(CurrDates));
                                    MaxDate=max(MaxDate,max(CurrDates));
                                    ManGseNbs(ManL)=ManGseNbs(ManL)+sum(CurrGseNb);
                                    ManGsmNbs(ManL)=ManGsmNbs(ManL)+sum(CurrGsmNb);
                                    [CurrDate,DateOrder]=sort(CurrDates);
                                    CurrGseNb=CurrGseNb(DateOrder);
                                    CurrGsmNb=CurrGsmNb(DateOrder);
                                    CurrGseNb=cumsum(CurrGseNb);
                                    CurrGsmNb=cumsum(CurrGsmNb);
                                    plot(CurrDate,CurrGseNb,'color',Colors(ManL,:),'marker','.')
                                    plot(CurrDate,CurrGsmNb,'color',Colors(ManL,:),'marker','+','markersize',3)
                                    plot(CurrDate,CurrGseNb,'color',Colors(ManL,:))
                                    plot(CurrDate,CurrGsmNb,'color',Colors(ManL,:))
                                end

                            end
                            %                datetick('x')
                            MinDate=datevec(MinDate);
                            MinDate=datenum([MinDate(1) 01 01]);
                            MaxDate=datevec(MaxDate);
                            MaxDate=datenum([MaxDate(1)+1 01 01]);
                            set(gca,'xlim',[MinDate,MaxDate])
                            set(gca,'xtick',[datenum(MinDate(1)),datenum(MaxDate(1))])
                            datetick('x','keeplimits','keepticks')
                            title(sprintf('%s',Technologies{ExistTech(TechL)}))
                            set(gca,'box','on')
                            set(gca,'yscale','log')
                        end

                        %put legend
                        subplot(YPlot,XPlot,TechNb+1)
                        hold on
                        Legend=cell(1,length(ManKeywords));
                        for ManL=1:length(ManKeywords)
                            Legend{ManL}=sprintf('%s (%u gse - %u gsm)',ManKeywords{ManL},ManGseNbs(ManL),ManGsmNbs(ManL));
                        end
                        for ManL=1:length(ManKeywords)
                            plot(0,0,'color',Colors(ManL,:))
                        end
                        legend(Legend,'location','BestOutside')
                        set(gca,'yscale','log')
                        set(gca,'visible','off')
                        set(gcf,'color',[1,1,1])
                        set_figsize('1280px')
                        saveas(h,sprintf('%s_GEO_%s',strrep(Species{SpeciesSel(SelL)},' ','_'),VERSION),'png')
                    end
                end %current species
            end
            geo_metadb('close geo metadb');
        end
        
    case 'read a GPL'
%% READ A GPL
        %process GPLs
        if nargin<2
            h=errordlg('misses GPL');
            waitfor(h)
            error('process canceled')
        end
        GPLs=varargin{1};
        if nargin==3
            KeepFlag=1;
            KeepGpl=varargin{2};
        else
            KeepFlag=0;
        end
        Ok=VERIF();
        if Ok
            try
            feature('DefaultCharacterSet','UTF-8')
            GseFields={'title';...
                'status';...
                'submissionDate';...
                'lastUpdateDate';...
                'pubmedId';...
                'summary';...
                'type';...
                'contributor';...
                'webLink';...
                'overallDesign';...
                'repeats';...
                'repeatsSampleList';...
                'variable';...
                'variableDescription';...
                'contact';...
                'supplementaryFile'};
            CurrGseFields={'title';...
                'status';...
                'submission_date';...
                'last_update_date';...
                'pubmed_id';...
                'summary';...
                'type';...
                'contributor';...
                'web_link';...
                'overall_design';...
                'repeats';...
                'repeats_sample_list';...
                'variable';...
                'variable_description';...
                'contact';...
                'supplementary_file'};

            GsmFields={'gsm';...
                'Id';...
                'title';...
                'seriesId';...
                'gpl';...
                'status';...
                'submissionDate';...
                'lastUpdateDate';...
                'type';...
                'sourceNameCh1';...
                'organismCh1';...
                'characteristicsCh1';...
                'moleculeCh1';...
                'labelCh1';...
                'treatmentProtocolCh1';...
                'extractProtocolCh1';...
                'labelProtocolCh1';...
                'sourceNameCh2';...
                'organismCh2';...
                'characteristicsCh2';...
                'moleculeCh2';...
                'labelCh2';...
                'treatmentProtocolCh2';...
                'extractProtocolCh2';...
                'labelProtocolCh2';...
                'hybProtocol';...
                'description';...
                'dataProcessing';...
                'contact';...
                'supplementaryFile';...
                'dataRowCount';...
                'channelCount'};
            CurrGsmFields={'gsm';...
                'ID';...              
                'title';...
                'series_id';...
                'gpl';...
                'status';...
                'submission_date';...
                'last_update_date';...
                'type';...
                'source_name_ch1';...
                'organism_ch1';...
                'characteristics_ch1';...
                'molecule_ch1';...
                'label_ch1';...
                'treatment_protocol_ch1';...
                'extract_protocol_ch1';...
                'label_protocol_ch1';...
                'source_name_ch2';...
                'organism_ch2';...
                'characteristics_ch2';...
                'molecule_ch2';...
                'label_ch2';...
                'treatment_protocol_ch2';...
                'extract_protocol_ch2';...
                'label_protocol_ch2';...
                'hyb_protocol';...
                'description';...
                'data_processing';...
                'contact';...
                'supplementary_file';...
                'data_row_count';...
                'channel_count'};
            FactorFields={'gdsIds';...
                'gdsRanks';...
                'factorNames';...
                'factorValues'};


            cd(K.dir.geoMetadata)
            LogFid=fopen(sprintf('display_gpl_%s_%s.log',VERSION,Date),'a');
            GdsFid=fopen(sprintf('gds_%s.log',VERSION),'a');
            fprintf(GdsFid,'Gpl\tAction\tGds\tGsm\n');
            GseFid=fopen(sprintf('gse_%s.log',VERSION),'a');
            fprintf(GseFid,'Gpl\tGse\tField\tOldValue\tNewValue\n');
            GsmFid=fopen(sprintf('gsm_%s.log',VERSION),'a');
            fprintf(GsmFid,'Gpl\tGse\tGsm\tField\tOldValue\tNewValue\n');
            %process each GPL
            for GplL=1:length(GPLs)
                try
                    GPL=GPLs{GplL};
                catch
                    GPL=GPLs(GplL);
                end
                %modify eventuelly to use the right format
                if isnumeric(GPL)
                    GPL=sprintf('GPL%u',GPL);
                end
                CurrGplRank=str2num(GPL(4:end));
                cd(K.dir.geoMetadata)
                if exist(sprintf('%s.mat',GPL),'file')
                    if KeepFlag==0
                        CurrKeepGpl=questdlg(sprintf('%s GPL already exists. Do you want to update it (otherwise cleared)',GPL),'','yes','no','yes');
                        waitfor(CurrKeepGpl)
                        if isequal(CurrKeepGpl,'yes')
                            eval(sprintf('load %s.mat;',GPL))
                            CurrKeepGpl=1;
                        else
                            CurrKeepGpl=0;
                        end
                    else
                        if KeepGpl
                            try
                                eval(sprintf('load %s.mat;',GPL))
                                CurrKeepGpl=1;
                            catch
                                CurrKeepGpl=0;
                            end
                        end
                    end
                else
                    CurrKeepGpl=0;
                end
                if CurrKeepGpl
                    NewGds=0;
                    NewGse=0;
                    NewGsm=0;
                else
                    %initiate empty structures
                    Gse=[];
                    Gsm=[];
                    Gds.gsmRank=[];
                    Gds.gdsIds=[];
                    Gds.metadb=[];         
                end


                %SCAN GDS TO RECOVER EVENTUAL FACTORS
                Gdss=mksqlite('SELECT * FROM gds_subset');

                %load or construct Gpl which contains ,for each Gpl, the list of Gdss
                % positions
                try
                    eval(sprintf('load GDS_%s.mat',VERSION))
                catch
                    Gpl=[];
                    Gpl.gplRank=[];
                    Gpl.gdsPos={};
                    for GdsL=1:length(Gdss)
                        %verify GPL
                        eval(sprintf('CurrGpl=mksqlite(''SELECT  gpl FROM gds WHERE ID=%u'');',Gdss(GdsL).ID));
                        if ~isempty(CurrGpl)
                            CurrGpl=str2num(CurrGpl.gpl(4:end));
                            GplPos=find(Gpl.gplRank==CurrGpl);
                            if isempty(GplPos)
                                Gpl.gplRank(end+1,1)=CurrGpl;
                                Gpl.gdsPos{end+1,1}=GdsL;
                            else
                                Gpl.gdsPos{GplPos,1}(end+1)=GdsL;
                            end
                        end
                    end
                    eval(sprintf('save GDS_%s.mat Gpl',VERSION))
                end

                if CurrKeepGpl
                    %construct list of GSM for each GDS
                    RegisteredGds=[];
                    RegisteredGds.gdsRank=[];
                    RegisteredGds.gsm={};
                    for GsmL=1:length(Gds.gsmRank)
                        for GdsL=1:length(Gds.gdsIds{GsmL})
                            GdsPos=find(RegisteredGds.gdsRank==Gds.gdsIds{GsmL}(GdsL));
                            if isempty(GdsPos)
                                RegisteredGds.gdsRank(end+1)=Gds.gdsIds{GsmL}(GdsL);
                                RegisteredGds.gsm{end+1,1}=Gds.gsmRank(GsmL);
                            else
                                RegisteredGds.gsm{GdsPos}(end+1)=Gds.gsmRank(GsmL);
                            end
                        end
                    end
                end

                %process each gds and write the list of gds ids in each found
                %gsm
                GplPos=find(Gpl.gplRank==CurrGplRank);
                if ~isempty(GplPos)
                    for GdsL=1:length(Gpl.gdsPos{GplPos})
                        GdsPos=Gpl.gdsPos{GplPos}(GdsL);
                        %recover GSMs
                        Gsms=regexp(Gdss(GdsPos).sample_id,'(?<=GSM)\d*(?=,?)','match');
                        CurrGsms=[];
                        for GsmL=1:length(Gsms)
                            CurrGsms(GsmL)=str2num(Gsms{GsmL});
                        end
                        CurrGds=Gdss(GdsPos).ID;
                        if CurrKeepGpl
                            %recover registered GSMs
                            RegisteredPos=find(RegisteredGds.gdsRank==CurrGds);
                            if ~isempty(RegisteredPos)
                                %test if some GSM have been removed
                                RemovedGsms=setdiff(RegisteredGds.gsm{RegisteredPos},CurrGsms);
                                if ~isempty(RemovedGsms)
                                    for GsmL=1:length(RemovedGsms)
                                        CurrGsmRank=RemovedGsms{GsmL};
                                        fprintf(GdsFid,'%u\tremoved\t%u\t%u\n',CurrGplRank,CurrGds,CurrGsmRank);
                                        GsmPos=find(Gds.gsmRank==CurrGsmRank);
                                        if ~isempty(GsmPos)
                                            CurrGdss=setdiff(Gds.gdsIds{GsmPos},CurrGds);
                                            if isempty(CurrGds)
                                                Gds.gsmRank(GsmPos)=[];
                                                Gds.gdsIds{GsmPos}=[];
                                                Gds.metadb{GsmPos}=[];
                                            else
                                                Gds.gdsIds{GsmPos}=CurrGdss;
                                                Gds.metadb{GsmPos}=VERSION;
                                            end
                                        end
                                    end
                                end
                                %test if some GSM have been addeed
                                AddedGsms=setdiff(CurrGsms,RegisteredGds.gsm{RegisteredPos});
                                if ~isempty(AddedGsms)
                                    for GsmL=1:length(AddedGsms)
                                        CurrGsmRank=AddedGsms{GsmL};
                                        fprintf(GdsFid,'%u\tadded\t%u\t%u\n',CurrGplRank,CurrGds,CurrGsmRank);
                                        GsmPos=find(Gds.gsmRank==CurrGsmRank);
                                        if isempty(GsmPos)
                                            GsmPos=length(Gds.gsmRank)+1;
                                            Gds.gsmRank(GsmPos,1)=CurrGsmRank;
                                            Gds.gdsIds{GsmPos}=CurrGds;
                                            Gds.metadb{GsmPos}=VERSION;
                                        else
                                            Gds.gdsIds{GsmPos}=sort([Gds.gdsIds{GsmPos},CurrGds]);
                                        end
                                    end
                                end
                            else
                                for GsmL=1:length(CurrGsms)
                                    CurrGsmRank=CurrGsms(GsmL);
                                    GsmPos=find(Gds.gsmRank==CurrGsmRank);
                                    if isempty(GsmPos)
                                        GsmPos=length(Gds.gsmRank)+1;
                                        Gds.gsmRank(GsmPos,1)=CurrGsmRank;
                                        Gds.gdsIds{GsmPos}=CurrGds;
                                        Gds.metadb{GsmPos}=VERSION;
                                    else
                                        Gds.gdsIds{GsmPos}=sort([Gds.gdsIds{GsmPos},CurrGds]);
                                    end
                                end
                            end
                        else
                            for GsmL=1:length(CurrGsms)
                                CurrGsmRank=CurrGsms(GsmL);
                                GsmPos=find(Gds.gsmRank==CurrGsmRank);
                                if isempty(GsmPos)
                                    GsmPos=length(Gds.gsmRank)+1;
                                    Gds.gsmRank(GsmPos,1)=CurrGsmRank;
                                    Gds.gdsIds{GsmPos}=CurrGds;
                                    Gds.metadb{GsmPos}=VERSION;
                                else
                                    Gds.gdsIds{GsmPos}=sort([Gds.gdsIds{GsmPos},CurrGds]);
                                end
                            end
                        end
                    end
                end
          

                %recover all Gses
                eval(sprintf('Gses=mksqlite(''SELECT * FROM sMatrix WHERE gpl=''''%s'''' ORDER BY ID'');',GPL));
                GseNb=length(Gses);
                %order GSE
                GseRanks=zeros(length(Gses),1);
                for GseL=1:length(Gses)
                    GseRanks(GseL)=str2num(Gses(GseL).gse(4:end));
                end
                [GseRanks,GseOrder]=sort(GseRanks);
                %process each Gse
                if CurrKeepGpl
                    %need MemGse to calculate new values before comparing them to old ones
                    MemGse=[];
                    MemGse.outGplGsmNb=[];
                    MemGse.gseRank=[];
                end

                for GseL=1:GseNb
                    GseP=GseOrder(GseL);
                    %gse name ('GSEXXX')
                    %RECOVER INFORMATION OF THE CURRENT GSE
                    eval(sprintf('CurrGse=mksqlite(''SELECT * FROM gse WHERE gse=''''%s'''''');',Gses(GseP).gse));

                    if ~isempty(CurrGse)
                        %if isempty, it means that the Gse is private and not
                        %yet released
                        CurrGseRank=str2num(Gses(GseP).gse(4:end));
                        if ~CurrKeepGpl
                            RegisterIt=1;
                        else
                            if length(find(Gse.gseRank==CurrGseRank))>1
                                fprintf(GseFid,'%u\t%u\tGsePos\t%u\t%u\n',CurrGplRank,CurrGseRank,length(find(Gse.gseRank==CurrGseRank)),Gses(GseP).ID);
                            end
                            GsePos=find(Gse.gseRank==CurrGseRank&Gse.id==Gses(GseP).ID);
                            if isempty(GsePos)
                                RegisterIt=1;
                                NewGse=NewGse+1;
                            else
                                RegisterIt=0;
                                ChangedFlag=0;
                                if Gse.id(GsePos)~=Gses(GseP).ID;
                                    fprintf(GseFid,'%u\t%u\tid\t%u\t%u\n',CurrGplRank,CurrGseRank,Gse.id(GsePos),Gses(GseP).ID)
                                    ChangedFlag=1;
                                end
                                %verify the number of GSM
                                eval(sprintf('Gsms=mksqlite(''SELECT * FROM gse_gsm WHERE gse=''''%s'''''');',Gses(GseP).gse));
                                %keep this information to determine the number GsmNb                        
                                if Gse.gsmNb(GsePos)~=length(Gsms);
                                    fprintf(GseFid,'%u\t%u\tgsmNb\t%u\t%u\n',CurrGplRank,CurrGseRank,Gse.gsmNb(GsePos),length(Gsms))
                                    Gse.gsmNb(GsePos)=length(Gsms);
                                    ChangedFlag=1;
                                end
                                if Gse.gsmCount(GsePos)~=Gses(GseP).GSM_Count
                                    fprintf(GseFid,'%u\t%u\tgsmCount\t%u\t%u\n',CurrGplRank,CurrGseRank,Gse.gsmCount(GsePos),Gses(GseP).GSM_Count);
                                    Gse.gsmCount(GsePos)=Gses(GseP).GSM_Count;
                                    ChangedFlag=1;
                                end
                                %Gse.outGplGsmNb(GsePos)=0;


                                for FieldL=1:length(GseFields)
                                    if eval(sprintf('~isequal(Gse.%s{GsePos},CurrGse.%s)',GseFields{FieldL},CurrGseFields{FieldL}))                                        
                                        try
                                            fprintf(GseFid,'%u\t%u\t%s\t%s\t%s\n',CurrGplRank,CurrGseRank,GseFields{FieldL},num2str(eval(sprintf('Gse.%s{GsePos}',GseFields{FieldL}))),num2str(eval(sprintf('CurrGse.%s',CurrGseFields{FieldL}))));
                                        catch
                                            fprintf(GseFid,'%u\t%u\t%s\t%s\t%s\n',CurrGplRank,CurrGseRank,GseFields{FieldL},eval(sprintf('Gse.%s{GsePos}',GseFields{FieldL})),eval(sprintf('CurrGse.%s',CurrGseFields{FieldL})));
                                        end
                                        eval(sprintf('Gse.%s{GsePos}=CurrGse.%s;',GseFields{FieldL},CurrGseFields{FieldL}));
                                        ChangedFlag=1;
                                    end
                                end
                                if ~isequal(Gse.lastUpdateDate{GsePos},Gses(GseP).Last_Update_Date)
                                    fprintf(GseFid,'%u\t%u\tlastUpdateDate\t%s\t%s\n',CurrGplRank,CurrGseRank,Gse.lastUpdateDate{GsePos},Gses(GseP).Last_Update_Date);   
                                    try
                                        Gse.lastUpdateDate{GsePos}=Gses(GseP).Last_Update_Date;
                                    catch
                                        Gse.lastUpdateDate(GsePos)='';
                                        Gse.lastUpdateDate{GsePos}='';
                                        Gse.lastUpdateDate{GsePos}=Gses(GseP).Last_Update_Date;
                                    end
                                    ChangedFlag=1;
                                end
                                if ChangedFlag
                                    Gse.metadb{GsePos}=VERSION;
                                end
                            end
                        end
                        if RegisterIt
                            if isempty(Gse)
                                GsePos=1;
                            else
                                GsePos=length(Gse.gse)+1;
                            end
                            Gse.gse{GsePos}=Gses(GseP).gse;
                            Gse.gseRank(GsePos)=CurrGseRank;
                            Gse.id(GsePos)=Gses(GseP).ID;
                            %verify the number of GSM
                            eval(sprintf('Gsms=mksqlite(''SELECT * FROM gse_gsm WHERE gse=''''%s'''''');',Gses(GseP).gse));
                            %keep this information to determine the number GsmNb
                            Gse.gsmNb(GsePos)=length(Gsms);
                            Gse.gsmCount(GsePos)=Gses(GseP).GSM_Count;
                            Gse.outGplGsmNb(GsePos)=0;

                            for FieldL=1:length(GseFields)
                                eval(sprintf('Gse.%s{GsePos}=CurrGse.%s;',GseFields{FieldL},CurrGseFields{FieldL}));
                            end
                            Gse.lastUpdateDate{GsePos}=Gses(GseP).Last_Update_Date;
                            Gse.metadb{GsePos}=VERSION;
                            Gse.isBiol(GsePos)=0;
                        end

                        %process each GSM
                        if CurrKeepGpl
                            %need MemGsm to calculate new values before comparing them to old ones
                            MemGsm=[];
                            MemGsm.gsmRank=[];
                            MemGsm.gseRanks={};
                            MemGsm.gseNb=[];
                        end
                        for GsmL=1:length(Gsms)
                            %RECOVER INFO OF EACH GSM
                            eval(sprintf('CurrGsm=mksqlite(''SELECT * FROM gsm WHERE gsm=''''%s'''''');',Gsms(GsmL).gsm));
                            if ~isempty(CurrGsm)
                                CurrGsmRank=str2num(CurrGsm.gsm(4:end));

                                % try
                                if isequal(CurrGsm.gpl,GPL)
                                    %some GSE may have GSM belonging to several GPL
                                    GsmPos=0;
                                    if ~isempty(Gsm)
                                        %uses index 1, because sometimes the same
                                        %record is repeated in CurrGsm
                                        if length(CurrGsm)>1
                                            for i=1:length(CurrGsm)
                                                if ~isequal(CurrGsm(i).gsm,CurrGsm(1).gsm)
                                                    h=warndlg(sprintf('different gsm in %s',Gsms(GsmL).gsm));
                                                    waitfor(h)
                                                end
                                            end
                                            CurrGsm=CurrGsm(1);
                                        end
                                        GsmPos=strmatch(CurrGsm.gsm,Gsm.gsm,'exact');
                                    end
                                    if ~CurrKeepGpl
                                        RegisterIt=1;
                                    else
                                        if isempty(GsmPos)
                                            RegisterIt=1;
                                            NewGsm=NewGsm+1;
                                        else
                                            RegisterIt=0;
                                            ChangedFlag=0;
                                            try
                                                MemGsm.gseRanks{GsmPos}=sort([MemGsm.gseRanks{GsmPos},Gse.gseRank(GsePos)]);
                                                MemGsm.gseNb(GsmPos)=MemGsm.gseNb(GsmPos)+1;
                                            catch
                                                MemGsm.gseRanks{GsmPos}=Gse.gseRank(GsePos);
                                                MemGsm.gseNb(GsmPos)=1;
                                            end
                                            MemGsm.gsmRank(GsmPos)=CurrGsmRank;
                                            for FieldL=1:length(GsmFields)
                                                try
                                                    if eval(sprintf('~isequal(Gsm.%s{GsmPos},CurrGsm.%s)',GsmFields{FieldL},CurrGsmFields{FieldL}))
                                                        try
                                                            fprintf(GsmFid,'%u\t%u\t%u\t%s\t%s\t%s\n',CurrGplRank,CurrGseRank,CurrGsmRank,GsmFields{FieldL},num2str(eval(sprintf('Gsm.%s{GsmPos}',GsmFields{FieldL}))),num2str(eval(sprintf('CurrGsm.%s',CurrGsmFields{FieldL}))));
                                                        catch
                                                            fprintf(GsmFid,'%u\t%u\t%u\t%s\t%s\t%s\n',CurrGplRank,CurrGseRank,CurrGsmRank,GsmFields{FieldL},eval(sprintf('Gsm.%s{GsmPos}',GsmFields{FieldL})),eval(sprintf('CurrGsm.%s',CurrGsmFields{FieldL})));
                                                        end
                                                        eval(sprintf('Gsm.%s{GsmPos}=CurrGsm.%s;',GsmFields{FieldL},CurrGsmFields{FieldL}));
                                                        ChangedFlag=1;
                                                    end
                                                catch
                                                    if eval(sprintf('~isequal(Gsm.%s(GsmPos),CurrGsm.%s)',GsmFields{FieldL},CurrGsmFields{FieldL}))
                                                        fprintf(GsmFid,'%u\t%u\t%u\t%s\t%u\t%u\n',CurrGplRank,CurrGseRank,CurrGsmRank,GsmFields{FieldL},eval(sprintf('Gsm.%s(GsmPos)',GsmFields{FieldL})),eval(sprintf('CurrGsm.%s',CurrGsmFields{FieldL})));
                                                        eval(sprintf('Gsm.%s(GsmPos)=CurrGsm.%s;',GsmFields{FieldL},CurrGsmFields{FieldL}));
                                                        ChangedFlag=1;
                                                    end
                                                end
                                            end




                                            %RECOVER EVENTUAL FACTORS
                                            GdsPos=find(Gds.gsmRank==Gsm.gsmRank(GsmPos));

                                            if ~isempty(GdsPos)
                                                GdsIds=Gds.gdsIds{GdsPos};
                                                MemGsm.gdsIds{GsmPos}=GdsIds;
                                                for GdsL=1:length(GdsIds)
                                                    GdsId=GdsIds(GdsL);
                                                    %RECOVER THE CURRENT GDS
                                                    eval(sprintf('CurrGds=mksqlite(''SELECT * FROM gds_subset WHERE ID=%u'');',GdsId))
                                                    GdsRank=str2num(CurrGds.gds(4:end));
                                                    NewFlag=0;
                                                    if isfield(MemGsm,'gdsRanks')
                                                        if length(MemGsm.gdsRanks)>=GsmPos
                                                            if isempty(MemGsm.gdsRanks{GsmPos})
                                                                NewFlag=1;
                                                            else
                                                                %in general exist several factors
                                                                MemGsm.gdsIds{GsmPos}=unique([MemGsm.gdsIds{GsmPos},GdsId]);
                                                                %may also exist several Gds
                                                                MemGsm.gdsRanks{GsmPos}=unique([MemGsm.gdsRanks{GsmPos},GdsRank]);
                                                                %find if the current factor already
                                                                %exists

                                                                FactorPos=strmatch(CurrGds.type,MemGsm.factorNames{GsmPos},'exact');
                                                                if isempty(FactorPos)
                                                                    MemGsm.factorNames{GsmPos}{end+1}=CurrGds.type;
                                                                    MemGsm.factorValues{GsmPos}{end+1}=CurrGds.description;
                                                                else
                                                                    %verify that the value is the same
                                                                    if iscell(MemGsm.factorValues{GsmPos}{FactorPos})
                                                                        %exist already
                                                                        %several values
                                                                        if isempty(strmatch(CurrGds.description,MemGsm.factorValues{GsmPos}{FactorPos},'exact'))
                                                                            MemGsm.factorValues{GsmPos}{FactorPos}{end+1,1}=CurrGds.description;
                                                                        end
                                                                    else
                                                                        if ~isequal(MemGsm.factorValues{GsmPos}{FactorPos},CurrGds.description)
                                                                            %transform the content of MemGsm.factorValues{GsmPos} into a cell
                                                                            MemVal=MemGsm.factorValues{GsmPos}{FactorPos};
                                                                            MemGsm.factorValues{GsmPos}{FactorPos}={};
                                                                            MemGsm.factorValues{GsmPos}{FactorPos}{1}=MemVal;
                                                                            MemGsm.factorValues{GsmPos}{FactorPos}{end+1,1}=CurrGds.description;
                                                                            %h=warndlg(sprintf('There exist several values for factor %s in %s',CurrGds.type,MemGsm.gsm{GsmPos}));
                                                                            %waitfor(h)
                                                                        end
                                                                    end
                                                                end

                                                            end
                                                        else
                                                            NewFlag=1;
                                                        end
                                                    else
                                                        NewFlag=1;
                                                    end
                                                    if NewFlag==1
                                                        MemGsm.gdsRanks{GsmPos}=GdsRank;
                                                        MemGsm.gdsIds{GsmPos}=GdsId;
                                                        MemGsm.factorNames{GsmPos}{1}=CurrGds.type;
                                                        MemGsm.factorValues{GsmPos}{1}=CurrGds.description;
                                                    end
                                                end
                                                MemGsm.gdsNb(GsmPos)=length(MemGsm.gdsRanks{GsmPos});


                                                for FieldL=1:length(FactorFields)
                                                    if eval(sprintf('~isequal(MemGsm.%s{GsmPos},Gsm.%s{GsmPos})',FactorFields{FieldL},FactorFields{FieldL}))
                                                        try
                                                            fprintf(GsmFid,'%u\t%u\t%u\t%s\t%s\t%s\n',CurrGplRank,CurrGseRank,CurrGsmRank,FactorFields{FieldL},num2str(eval(sprintf('MemGsm.%s{GsmPos}',FactorFields{FieldL}))),num2str(eval(sprintf('Gsm.%s{GsmPos}',FactorFields{FieldL}))));
                                                        catch
                                                            fprintf(GsmFid,'%u\t%u\t%u\t%s\t%s\t%s\n',CurrGplRank,CurrGseRank,CurrGsmRank,FactorFields{FieldL},eval(sprintf('MemGsm.%s{GsmPos}',FactorFields{FieldL})),eval(sprintf('Gsm.%s{GsmPos}',FactorFields{FieldL})));
                                                        end
                                                        eval(sprintf('Gsm.%s{GsmPos}=MemrGsm.%s{GsmPos};',FactorFields{FieldL},FactorFields{FieldL}));
                                                        ChangedFlag=1;
                                                    end
                                                end
                                                if ChangedFlag
                                                    Gsm.metadb{GsmPos}=VERSION;
                                                end
                                            end

                                        end
                                    end
                                    if RegisterIt
                                        if GsmPos>0
                                            %increment Gse nb (several Gse reference the same GSM)
                                            try
                                                Gsm.gseRanks{GsmPos}=sort([Gsm.gseRanks{GsmPos},Gse.gseRank(GsePos)]);
                                                Gsm.gseNb(GsmPos)=Gsm.gseNb(GsmPos)+1;
                                            catch
                                                Gsm.gseRanks{GsmPos}=Gse.gseRank(GsePos);
                                                Gsm.gseNb(GsmPos)=1;
                                            end
                                        end
                                        if isempty(Gsm)
                                            GsmPos=1;
                                        else
                                            GsmPos=length(Gsm.gsm)+1;
                                        end

                                        %Gsm.gsm{GsmPos}=CurrGsm.gsm;
                                        Gsm.gsmRank(GsmPos)=CurrGsmRank;
                                        for FieldL=1:length(GsmFields)
                                            try
                                                eval(sprintf('Gsm.%s{GsmPos}=CurrGsm.%s;',GsmFields{FieldL},CurrGsmFields{FieldL}));
                                            catch
                                                eval(sprintf('Gsm.%s(GsmPos)=CurrGsm.%s;',GsmFields{FieldL},CurrGsmFields{FieldL}));
                                            end
                                        end


                                        Gsm.gseRanks{GsmPos}=Gse.gseRank(GsePos);
                                        Gsm.gseNb(GsmPos)=1;
                                        Gsm.metadb{GsmPos}=VERSION;
                                        Gsm.isBiol(GsmPos)=0;

                                        %RECOVER EVENTUAL FACTORS
                                        GdsPos=find(Gds.gsmRank==Gsm.gsmRank(GsmPos));

                                        if ~isempty(GdsPos)
                                            GdsIds=Gds.gdsIds{GdsPos};
                                            Gsm.gdsIds{GsmPos}=GdsIds;
                                            for GdsL=1:length(GdsIds)
                                                GdsId=GdsIds(GdsL);
                                                %RECOVER THE CURRENT GDS
                                                eval(sprintf('CurrGds=mksqlite(''SELECT * FROM gds_subset WHERE ID=%u'');',GdsId))
                                                GdsRank=str2num(CurrGds.gds(4:end));
                                                NewFlag=0;
                                                if isfield(Gsm,'gdsRanks')
                                                    if length(Gsm.gdsRanks)>=GsmPos
                                                        if isempty(Gsm.gdsRanks{GsmPos})
                                                            NewFlag=1;
                                                        else
                                                            %in general exist several factors
                                                            Gsm.gdsIds{GsmPos}=unique([Gsm.gdsIds{GsmPos},GdsId]);
                                                            %may also exist several Gds
                                                            Gsm.gdsRanks{GsmPos}=unique([Gsm.gdsRanks{GsmPos},GdsRank]);
                                                            %find if the current factor already
                                                            %exists

                                                            FactorPos=strmatch(CurrGds.type,Gsm.factorNames{GsmPos},'exact');
                                                            if isempty(FactorPos)
                                                                Gsm.factorNames{GsmPos}{end+1}=CurrGds.type;
                                                                Gsm.factorValues{GsmPos}{end+1}=CurrGds.description;
                                                            else
                                                                %verify that the value is the same
                                                                if iscell(Gsm.factorValues{GsmPos}{FactorPos})
                                                                    %exist already
                                                                    %several values
                                                                    if isempty(strmatch(CurrGds.description,Gsm.factorValues{GsmPos}{FactorPos},'exact'))
                                                                        Gsm.factorValues{GsmPos}{FactorPos}{end+1,1}=CurrGds.description;
                                                                    end
                                                                else
                                                                    if ~isequal(Gsm.factorValues{GsmPos}{FactorPos},CurrGds.description)
                                                                        %transform the content of Gsm.factorValues{GsmPos} into a cell
                                                                        MemVal=Gsm.factorValues{GsmPos}{FactorPos};
                                                                        Gsm.factorValues{GsmPos}{FactorPos}={};
                                                                        Gsm.factorValues{GsmPos}{FactorPos}{1}=MemVal;
                                                                        Gsm.factorValues{GsmPos}{FactorPos}{end+1,1}=CurrGds.description;
                                                                        %h=warndlg(sprintf('There exist several values for factor %s in %s',CurrGds.type,Gsm.gsm{GsmPos}));
                                                                        %waitfor(h)
                                                                    end
                                                                end
                                                            end

                                                        end
                                                    else
                                                        NewFlag=1;
                                                    end
                                                else
                                                    NewFlag=1;
                                                end
                                                if NewFlag==1
                                                    Gsm.gdsRanks{GsmPos}=GdsRank;
                                                    Gsm.gdsIds{GsmPos}=GdsId;
                                                    Gsm.factorNames{GsmPos}{1}=CurrGds.type;
                                                    Gsm.factorValues{GsmPos}{1}=CurrGds.description;
                                                end
                                            end
                                            Gsm.gdsNb(GsmPos)=length(Gsm.gdsRanks{GsmPos});
                                        end
                                    end
                                else
                                    %indicate that exist other gpl
                                    if CurrKeepGpl
                                        try
                                            MemGse.outGplGsmNb(GsePos)=Gse.outGplGsmNb(GsePos)+1;
                                        catch
                                            MemGse.outGplGsmNb(GsePos)=1;
                                        end
                                        MemGse.gseRank(GsePos)=CurrGseRank;
                                    else
                                        Gse.outGplGsmNb(GsePos)=Gse.outGplGsmNb(GsePos)+1;
                                    end
                                end
                            end
                        end %of GsmL
                    end %if ~isempty

                    %process MemGsm
                    if CurrKeepGpl
                        GsmPos=find(MemGsm.gsmRank);
                        for GsmL=1:length(GsmPos)
                            CurrGsmPos=GsmPos(GsmL);
                            if length(Gsm.gseRanks{CurrGsmPos})==length(MemGsm.gseRanks{CurrGsmPos})
                                if Gsm.gseRanks{CurrGsmPos}~=MemGsm.gseRanks{CurrGsmPos}
                                    fprintf(GsmFid,'%u\t%u\t%u\tgseRanks\t%s\t%s\n',CurrGplRank,CurrGseRank,Gsm.gsmRank(CurrGsmPos),num2str(MemGsm.gseRanks{CurrGsmPos}),num2str(Gsm.gseRanks{CurrGsmPos}));
                                    Gsm.gseRanks{CurrGsmPos}=MemGsm.gseRanks{CurrGsmPos};
                                    Gsm.metadb{CurrGsmPos}=VERSION;
                                end
                            else
                                fprintf(GsmFid,'%u\t%u\t%u\tgseRanks\t%s\t%s\n',CurrGplRank,CurrGseRank,Gsm.gsmRank(CurrGsmPos),num2str(MemGsm.gseRanks{CurrGsmPos}),num2str(Gsm.gseRanks{CurrGsmPos}));
                                Gsm.gseRanks{CurrGsmPos}=MemGsm.gseRanks{CurrGsmPos};
                                Gsm.metadb{CurrGsmPos}=VERSION;
                            end
                            try
                                if Gsm.gseNb(CurrGsmPos)~=MemGsm.gseNb(CurrGsmPos)
                                    fprintf(GsmFid,'%u\t%u\t%u\tgseNb\t%u\t%u\n',CurrGplRank,CurrGseRank,Gsm.gsmRank(CurrGsmPos),MemGsm.gseNb(CurrGsmPos),Gsm.gseNb(CurrGsmPos));
                                    Gsm.gseNb(CurrGsmPos)=MemGsm.gseNb(CurrGsmPos);
                                    Gsm.metadb{CurrGsmPos}=VERSION;
                                end
                            catch
                                %patch to recover a previous error (Gsm.gseNb=1 instead of
                                %Gsm.gseNb(GsmPos)=1
                                fprintf(GsmFid,'%u\t%u\t%u\tgseNb\t%u\t-\n',CurrGplRank,CurrGseRank,Gsm.gsmRank(CurrGsmPos),MemGsm.gseNb(CurrGsmPos));
                                Gsm.gseNb(CurrGsmPos)=MemGsm.gseNb(CurrGsmPos);
                                Gsm.metadb{CurrGsmPos}=VERSION;
                            end
                        end
                    end
                end %of GseL
                %process MemGse and MemGsm
                if CurrKeepGpl
                    GsePos=find(MemGse.gseRank);
                    for GseL=1:length(GsePos)
                        CurrGsePos=GsePos(GseL);
                        if Gse.outGplGsmNb(CurrGsePos)~=MemGse.outGplGsmNb(CurrGsePos)
                            fprintf(GseFid,'%u\t%u\toutGplGsmNb\t%u\t%u\n',CurrGplRank,Gse.gseRank(CurrGsePos),MemGse.outGplGsmNb(CurrGsePos),Gse.outGplGsmNb(CurrGsePos));
                            Gse.outGplGsmNb(CurrGsePos)=MemGse.outGplGsmNb(CurrGsePos);
                            Gse.metadb{CurrGsePos}=VERSION;
                        end
                    end
                end



                %FIND CONTRIBUTORS
                Contributor.lastName={};
                Contributor.firstNames={};
                Contributor.gses={};
                Contributor.gseRanks={};
                Contributor.gseIndexes={};
                Contributor.gseNb=[];
                Contributor.metadb=[];
                for GseL=1:length(Gse.gse)
                    CurrNames=sprintf(';%s;',Gse.contributor{GseL});
                    SepPos=findstr(';',CurrNames);
                    for NameL=1:length(SepPos)-1
                        CurrName=CurrNames(SepPos(NameL)+1:SepPos(NameL+1)-1);
                        ColonPos=findstr(',',CurrName);
                        if isempty(ColonPos)
                            ColonPos=findstr(' ',CurrName);
                        end
                        if isempty(ColonPos)
                            LastName=CurrName;
                        else
                            FirstNames=CurrName(1:ColonPos(end)-1);
                            LastName=upper(CurrName(ColonPos(end)+1:end));
                        end
                        if ~isempty(LastName)
                            NamePos=strmatch(LastName,Contributor.lastName,'exact');
                            if isempty(NamePos)
                                Contributor.lastName{end+1}=LastName;
                                Contributor.firstNames{end+1}=FirstNames;
                                Contributor.gses{end+1}{1}=Gse.gse{GseL};
                                Contributor.gseRanks{end+1}=Gse.gseRank(GseL);
                                Contributor.gseIndexes{end+1}=GseL;
                                Contributor.gseNb(end+1)=1;
                                Contributor.metadb{end+1}=VERSION;
                            else
                                if isempty(find(Contributor.gseRanks{NamePos}==Gse.gseRank(GseL)))
                                    [Contributor.gseRanks{NamePos},SortOrder]=sort([Contributor.gseRanks{NamePos},Gse.gseRank(GseL)]);
                                    Contributor.gses{NamePos}{end+1}=Gse.gse{GseL};
                                    Contributor.gses{NamePos}=Contributor.gses{NamePos}(SortOrder);
                                    Contributor.gseIndexes{NamePos}=[Contributor.gseIndexes{NamePos},GseL];
                                    Contributor.gseIndexes{NamePos}=Contributor.gseIndexes{NamePos}(SortOrder);
                                    Contributor.gseNb(NamePos)=Contributor.gseNb(NamePos)+1;
                                end
                            end
                        end
                    end
                end
                cd(K.dir.geoMetadata)
                fid=fopen(sprintf('contributors_%s_%s.txt',GPL,VERSION),'w');
                fprintf(fid,'first name\tlast name\tgse\tgse nb\n')
                for ContL=1:length(Contributor.gses)
                    fprintf(fid,'%s\t%s',strrep(strrep(Contributor.firstNames{ContL},sprintf('\t'),''),',',' '),Contributor.lastName{ContL});
                    fprintf(fid,'\t%s',Contributor.gses{ContL}{1});
                    if length(Contributor.gses{ContL})>1
                        for GseL=2:length(Contributor.gses{ContL})
                            fprintf(fid,'/%s',Contributor.gses{ContL}{GseL});
                        end
                    end
                    fprintf(fid,'\t%u\n',Contributor.gseNb(ContL));
                end
                fclose(fid)


                Pb=find(Gse.gsmCount~=Gse.gsmNb);
                Pb=find(Gse.gsmNb(Pb)-Gse.gsmCount(Pb)~=Gse.outGplGsmNb(Pb));
                if ~isempty(Pb)
                    fid=fopen(sprintf('missing_gsm_%s_%s.txt',GPL,VERSION),'w');
                    fprintf(fid,'gse\tgsm count\tgsm nb\tin other gpl\n');
                    for PbL=1:length(Pb)
                        fprintf(fid,'%u\t%u\t%u\t%u\n',Gse.gseRank(Pb(PbL)),Gse.gsmCount(Pb(PbL)),Gse.gsmNb(Pb(PbL)),Gse.outGplGsmNb(Pb(PbL)));
                    end
                    fclose(fid)
                end


                Gsm.gseRank=zeros(length(Gsm.gsm),1);
                for GsmL=1:length(Gsm.gsm)
                    if length(Gsm.gseRanks{GsmL})==1
                        Gsm.gseRank(GsmL)=Gsm.gseRanks{GsmL};
                    else
                        Gsm.gseRank(GsmL)=Gsm.gseRanks{GsmL}(1);
                    end
                end


                ModifiedGse=strmatch(VERSION,Gse.metadb,'exact');
                ModifiedGsm=strmatch(VERSION,Gsm.metadb,'exact');
                if CurrKeepGpl==0
                    NewGse=length(ModifiedGse);
                    NewGsm=length(ModifiedGsm);
                end
                %h=warndlg(sprintf('%u new, %u updated',NewGse,length(ModifiedGse)-NewGse),'GSE');
                %waitfor(h)
                %h=warndlg(sprintf('%u new, %u updated',NewGsm,length(ModifiedGsm)-NewGsm),'GSM');
                %waitfor(h)                
                %cd(K.dir.geoMetadata)
                %sprintf('GPL%u: %u new GSE, %u updated',CurrGplRank,NewGse,length(ModifiedGse)-NewGse)
                %sprintf('GPL%u: %u new GSM, %u updated',CurrGplRank,NewGsm,length(ModifiedGsm)-NewGsm)
                fprintf(LogFid,'%u\t%u\t%u\n',CurrGplRank,NewGse,length(ModifiedGse)-NewGse)
                fprintf(LogFid,'%u\t%u\t%u\n',CurrGplRank,NewGsm,length(ModifiedGsm)-NewGsm)
                eval(sprintf('save %s Gsm Gse Gds Contributor',GPL))
            end
            catch
                fclose(LogFid);
                fclose(GdsFid);
                fclose(GseFid);
                fclose(GsmFid);

                h=warndlg(sprintf('GPL%u has failed',CurrGplRank));
                waitfor(h)
                'stop'
                feature('DefaultCharacterSet','windows-1252')
                geo_metadb('close geo metadb');
            end
        end


%% GSE special symbols
        %find neighbourhood of special charater to set up recover_words function used in GSE
        %dictionnary
        %writes two list by reordering list on words that are just ahead and behind of the symbol
        %geo_metadb('GSE special symbols',[339],sprintf('%c',char(65533)),1,1,30,30)
    case 'GSE special symbols'
        if nargin<3
            h=errordlg('GSE special symbols needs 6 parameters: a list of GPL in numerical format, a regular expression pattern, length of pattern, position of symbol and length of left and right neighbhourood');
            waitfor(h)
            error('process canceled')
        end
        feature('DefaultCharacterSet','UTF-8')
        GPLs=varargin{1};
        Exp=varargin{2};
        MatchSize=varargin{3};
        SymbolPos=varargin{4};
        LeftSize=varargin{5};
        RightSize=varargin{6};
        StrSize=LeftSize+MatchSize+RightSize;
        SymbolPos=LeftSize+SymbolPos;
        %recover all occurences of the symbol with text on both side
        WordList='';
        for GplL=1:length(GPLs)
            %search if current GPL exists in K.chip
            ChipRank=strmatch(sprintf('GPL%u',GPLs(GplL)),K.chip.geoName,'exact');
            if isempty(ChipRank)
                h=warndlg(sprintf('GPL%u does not exists',GPLs(GPL)));
                waitfor(h)
            else
                GPL=sprintf('GPL%u.mat',GPLs(GplL));
                cd(K.dir.geoMetadata)
                if ~exist(GPL,'file')
                    h=warndlg(sprintf('%s does not exists',GPL));
                    waitfor(h)
                else
                    eval(sprintf('load %s',GPL))
                    for GseL=1:length(Gse.summary)
                        %for GseL=1:20
                        String=Gse.summary{GseL};
                        if ~isempty(String)
                            %[Match,Index]=regexp(String,'[\|\$\[\]\)(}{_''"#`~&^%*?<>?@=???]','match','start');
                            %[Match,Index]=regexp(String,'\S\.\S|\S;\S|\S:\S|\S,\S|\S!\S|\S\?\S|','match','start');                            
                            [Match,Index]=regexp(String,Exp,'match','start');                            
                            if ~isempty(Index)
                                for IndexL=1:length(Index)
                                    CurrIndex=Index(IndexL);                                    
                                    CurrStr=String(CurrIndex-min(LeftSize,CurrIndex-1):CurrIndex+MatchSize-1+min(RightSize,length(String)-CurrIndex-MatchSize+1));
                                    CurrStr=regexprep(CurrStr,'\s',' ');
                                    if length(CurrStr)~=StrSize
                                        %fill missing position with _ character
                                        CurrIndex=min(LeftSize,CurrIndex-1)+1;
                                        CurrStr=[repmat('_',1,LeftSize-CurrIndex+1),CurrStr,repmat('_',1,RightSize+MatchSize-length(CurrStr)+CurrIndex-1)];                                        
                                    end
                                    WordList(end+1,:)=CurrStr;
                                end
                            end
                        end
                    end
                end
            end
        end
        'stop'
        %split WordList according to the symbols recovered
        %list of unique symbols
        Symbol=unique(WordList(:,SymbolPos));
        WordLists=cell(length(Symbol),1);
        for SymbolL=1:length(Symbol)
            CurrSymbolPos=find(WordList(:,SymbolPos)==Symbol(SymbolL));            
            WordLists{SymbolL}=WordList(CurrSymbolPos,:);                       
        end

        %clear 'words' which are not significative (e.g. , before a tabulation, or , in a number)
        %filtering for , symbol
        if ~isempty(find(Symbol==','))
            ItemPos=find(Symbol==',');
            %tabulation are replaced by ;
            % , before a tabulation is not significative
            ClearPos=find(WordLists{ItemPos}(:,SymbolPos+1)==';');
            WordLists{ItemPos}(ClearPos,:)=[];                        
            %numbers
            Idx1=regexp(WordLists{ItemPos}(:,SymbolPos-1)','\d','start');
            Idx2=regexp(WordLists{ItemPos}(:,SymbolPos+1)','\d','start');
            Idx=intersect(Idx1,Idx2);
            ClearPos=[];
            for IdL=1:length(Idx)
                if regexp(WordLists{ItemPos}(Idx(IdL),:),'[\(\s><]\d+,\d+[\s\)]','start')<=SymbolPos&regexp(WordLists{ItemPos}(Idx(IdL),:),'[\(\s><]\d+,\d+[\s\)]','end')>=SymbolPos
                    ClearPos(end+1)=IdL;
                end
            end
            WordLists{ItemPos}(Idx(ClearPos),:)=[];            
            %reference page
            Idx1=regexp(WordLists{ItemPos}(:,SymbolPos-1)','\d','start');
            Idx2=regexp(WordLists{ItemPos}(:,SymbolPos+1)','\d','start');
            Idx=intersect(Idx1,Idx2);
            ClearPos=[];
            for IdL=1:length(Idx)
                if regexp(WordLists{ItemPos}(Idx(IdL),:),' \d+(,\d+)+ (?!tetra)','start')<=SymbolPos&regexp(WordLists{ItemPos}(Idx(IdL),:),' \d+(,\d+)+ (?!tetra)','end')>=SymbolPos
                    ClearPos(end+1)=IdL;
                end
            end
            WordLists{ItemPos}(Idx(ClearPos),:)=[];
        end
            

        %filtering for ; symbol
        if ~isempty(find(Symbol==';'))
            ItemPos=find(Symbol==';');
            %tabulation
            ClearPos=find(WordLists{ItemPos}(:,SymbolPos+1)==';');
            WordLists{ItemPos}(ClearPos,:)=[];            
            %numbers
            Idx1=regexp(WordLists{ItemPos}(:,SymbolPos-1)','\d','start');
            Idx2=regexp(WordLists{ItemPos}(:,SymbolPos+1)','\d','start');
            Idx=intersect(Idx1,Idx2);
            WordLists{ItemPos}(Idx,:)=[];
            %chromosomal bands:
            Idx1=regexp(WordLists{ItemPos}(:,SymbolPos-1)','\d','start');
            Idx2=regexp(WordLists{ItemPos}(:,SymbolPos+1)','[p,q]','start');
            Idx=intersect(Idx1,Idx2);
            WordLists{ItemPos}(Idx,:)=[];
            %translocation X
            Idx1=regexp(WordLists{ItemPos}(:,SymbolPos-1)','X','start');
            Idx2=regexp(WordLists{ItemPos}(:,SymbolPos+1)','\d','start');
            Idx=intersect(Idx1,Idx2);
            WordLists{ItemPos}(Idx,:)=[];
        end

            
        %filtering for : symbol
        if ~isempty(find(Symbol==':'))
            ItemPos=find(Symbol==':');
            %tabulation
            ClearPos=find(WordLists{ItemPos}(:,SymbolPos+1)==';');
            WordLists{ItemPos}(ClearPos,:)=[];            
            %nummbers
            Idx1=regexp(WordLists{ItemPos}(:,SymbolPos-1)','\d','start');
            Idx2=regexp(WordLists{ItemPos}(:,SymbolPos+1)','\d','start');
            Idx=intersect(Idx1,Idx2);
            WordLists{ItemPos}(Idx,:)=[];
            %ftp:
            Idx1=regexp(WordLists{ItemPos}(:,SymbolPos-3)','f','start');
            Idx2=regexp(WordLists{ItemPos}(:,SymbolPos-2)','t','start');
            Idx3=regexp(WordLists{ItemPos}(:,SymbolPos-1)','p','start');
            Idx=intersect(Idx1,intersect(Idx2,Idx3));
            WordLists{ItemPos}(Idx,:)=[];
            %http:
            Idx1=regexp(WordLists{ItemPos}(:,SymbolPos-4)','h','start');
            Idx2=regexp(WordLists{ItemPos}(:,SymbolPos-3)','t','start');
            Idx3=regexp(WordLists{ItemPos}(:,SymbolPos-2)','t','start');
            Idx4=regexp(WordLists{ItemPos}(:,SymbolPos-1)','p','start');
            Idx=intersect(Idx1,intersect(Idx2,intersect(Idx3,Idx4)));   
            WordLists{ItemPos}(Idx,:)=[]; 
        end
        %filtering for . symbol
        if ~isempty(find(Symbol=='.'))
            ItemPos=find(Symbol=='.');
            %tabulation
            ClearPos=find(WordLists{ItemPos}(:,SymbolPos+1)==';');
            WordLists{ItemPos}(ClearPos,:)=[];
            %comma
            ClearPos=find(WordLists{ItemPos}(:,SymbolPos+1)==',');
            WordLists{ItemPos}(ClearPos,:)=[];
            %closing parenthesis
            ClearPos=find(WordLists{ItemPos}(:,SymbolPos+1)==')');
            WordLists{ItemPos}(ClearPos,:)=[];
            %nummbers
            Idx1=regexp(WordLists{ItemPos}(:,SymbolPos-1)','\d','start');
            Idx2=regexp(WordLists{ItemPos}(:,SymbolPos+1)','\d','start');
            Idx=intersect(Idx1,Idx2);
            WordLists{ItemPos}(Idx,:)=[];
            % .txt
            Idx1=regexp(WordLists{ItemPos}(:,SymbolPos+1)','t','start');
            Idx2=regexp(WordLists{ItemPos}(:,SymbolPos+2)','x','start');
            Idx3=regexp(WordLists{ItemPos}(:,SymbolPos+3)','t','start');
            Idx=intersect(Idx1,intersect(Idx2,Idx3));
            WordLists{ItemPos}(Idx,:)=[];
            %www.
            Idx1=regexp(WordLists{ItemPos}(:,SymbolPos-3)','w','start');
            Idx2=regexp(WordLists{ItemPos}(:,SymbolPos-2)','w','start');
            Idx3=regexp(WordLists{ItemPos}(:,SymbolPos-1)','w','start');
            Idx=intersect(Idx1,intersect(Idx2,Idx3));
            WordLists{ItemPos}(Idx,:)=[];
            %.com
            Idx1=regexp(WordLists{ItemPos}(:,SymbolPos+1)','c','start');
            Idx2=regexp(WordLists{ItemPos}(:,SymbolPos+2)','o','start');
            Idx3=regexp(WordLists{ItemPos}(:,SymbolPos+3)','m','start');
            Idx=intersect(Idx1,intersect(Idx2,Idx3));
            WordLists{ItemPos}(Idx,:)=[];
            %.gov
            Idx1=regexp(WordLists{ItemPos}(:,SymbolPos+1)','g','start');
            Idx2=regexp(WordLists{ItemPos}(:,SymbolPos+2)','o','start');
            Idx3=regexp(WordLists{ItemPos}(:,SymbolPos+3)','v','start');
            Idx=intersect(Idx1,intersect(Idx2,Idx3));
            WordLists{ItemPos}(Idx,:)=[];
             %.edu
            Idx1=regexp(WordLists{ItemPos}(:,SymbolPos+1)','e','start');
            Idx2=regexp(WordLists{ItemPos}(:,SymbolPos+2)','d','start');
            Idx3=regexp(WordLists{ItemPos}(:,SymbolPos+3)','u','start');
            Idx=intersect(Idx1,intersect(Idx2,Idx3));
            WordLists{ItemPos}(Idx,:)=[];
             %.org
            Idx1=regexp(WordLists{ItemPos}(:,SymbolPos+1)','o','start');
            Idx2=regexp(WordLists{ItemPos}(:,SymbolPos+2)','r','start');
            Idx3=regexp(WordLists{ItemPos}(:,SymbolPos+3)','g','start');
            Idx=intersect(Idx1,intersect(Idx2,Idx3));
            WordLists{ItemPos}(Idx,:)=[];
            %.nlm
            Idx1=regexp(WordLists{ItemPos}(:,SymbolPos+1)','n','start');
            Idx2=regexp(WordLists{ItemPos}(:,SymbolPos+2)','l','start');
            Idx3=regexp(WordLists{ItemPos}(:,SymbolPos+3)','m','start');
            Idx=intersect(Idx1,intersect(Idx2,Idx3));
            WordLists{ItemPos}(Idx,:)=[];
            %.nih
            Idx1=regexp(WordLists{ItemPos}(:,SymbolPos+1)','n','start');
            Idx2=regexp(WordLists{ItemPos}(:,SymbolPos+2)','i','start');
            Idx3=regexp(WordLists{ItemPos}(:,SymbolPos+3)','h','start');
            Idx=intersect(Idx1,intersect(Idx2,Idx3));
            WordLists{ItemPos}(Idx,:)=[];
            %.pdf
            Idx1=regexp(WordLists{ItemPos}(:,SymbolPos+1)','p','start');
            Idx2=regexp(WordLists{ItemPos}(:,SymbolPos+2)','d','start');
            Idx3=regexp(WordLists{ItemPos}(:,SymbolPos+3)','f','start');
            Idx=intersect(Idx1,intersect(Idx2,Idx3));
            WordLists{ItemPos}(Idx,:)=[];
            %.html
            Idx1=regexp(WordLists{ItemPos}(:,SymbolPos+1)','h','start');
            Idx2=regexp(WordLists{ItemPos}(:,SymbolPos+2)','t','start');
            Idx3=regexp(WordLists{ItemPos}(:,SymbolPos+3)','m','start');
            Idx4=regexp(WordLists{ItemPos}(:,SymbolPos+4)','l','start');
            Idx=intersect(Idx1,intersect(Idx2,intersect(Idx3,Idx4)));
            WordLists{ItemPos}(Idx,:)=[];
            %e.g.
            Idx1=regexp(WordLists{ItemPos}(:,SymbolPos-1)','e','start');
            Idx2=regexp(WordLists{ItemPos}(:,SymbolPos+1)','g','start');
            Idx3=regexp(WordLists{ItemPos}(:,SymbolPos+2)','.','start');
            Idx=intersect(Idx1,intersect(Idx2,Idx3));
            WordLists{ItemPos}(Idx,:)=[];
            %i.e.
            Idx1=regexp(WordLists{ItemPos}(:,SymbolPos-1)','i','start');
            Idx2=regexp(WordLists{ItemPos}(:,SymbolPos+1)','e','start');
            Idx3=regexp(WordLists{ItemPos}(:,SymbolPos+2)','.','start');
            Idx=intersect(Idx1,intersect(Idx2,Idx3));
            WordLists{ItemPos}(Idx,:)=[];            
        end

        for SymbolL=1:length(Symbol)
            %replace tabulation
            TabPos=find(WordLists{SymbolL}==char(9));
            WordLists{SymbolL}(TabPos)=' ';            
            %reorder list on words that are just ahead and behind of the symbol
            %construct list of words
            FirstWord=cell(size(WordLists{SymbolL},1),1);
            SndWord=cell(size(WordLists{SymbolL},1),1);
            for WordL=1:size(WordLists{SymbolL},1)

                CurrWord=regexp(WordLists{SymbolL}(WordL,1:SymbolPos-1),'(?<=\s)(\S+$)|(\s$)','match');
                if ~isempty(CurrWord)
                    FirstWord{WordL}=strtrim(CurrWord{1});
                else
                    FirstWord{WordL}='';
                end
                CurrWord=regexp(WordLists{SymbolL}(WordL,SymbolPos+1:end),'(^\S+)(?=\s)|(^\s)','match');
                if ~isempty(CurrWord)
                    SndWord{WordL}=strtrim(CurrWord{1});
                else
                    SndWord{WordL}='';
                end

            end
            [temp SortOrder1]=sort(FirstWord);
            [temp SortOrder2]=sort(SndWord);
            
            WordLists{SymbolL}=[WordLists{SymbolL}(SortOrder1,:),repmat(char(9),length(SortOrder1),1),WordLists{SymbolL}(SortOrder2,:)];
        end
        
        cd(K.dir.geoMetadata)
        for SymbolL=1:length(Symbol)            
            f=fopen(sprintf('gsesymbols_char%03u_%s.txt',double(Symbol(SymbolL)),Date),'w');
            for StrL=1:size(WordLists{SymbolL},1)
                fprintf(f,'%s\n',WordLists{SymbolL}(StrL,:));
            end
            fclose(f)
        end
        feature('DefaultCharacterSet','windows-1252')
        

%% GSE dictionary
%find words in GSE summary field
% varargin:
% 1 List of GPL
% dispensable:
% 2 Regular expression to be searched
% 3 Output name
    case 'GSE dictionary'
        if nargin<2
            h=errordlg('GSE dictionary needs a list of GPL in numerical format');
            waitfor(h)
            error('process canceled')
        end
        ExpFlag=0;
        if nargin==4
            ExpFlag=1;
            SearchedName=varargin{3};
            %search for a regular expression
            Exp=varargin{2};
        end
        %load existing vocabulatory if exist
        Words={};
        ClearedWords={};
        Symbols={};

        try
            cd(K.dir.geoMetadata)
            [Words,WordsNb]=textread('words.txt','%s%u');
            ClearedWords=textread('cleared.txt','%s');
            [Symbols,SymbolsNb]=textread('symbols.txt','%s%u');
            WordFlag=1;
        catch
            WordFlag=0;
        end

        feature('DefaultCharacterSet','UTF-8')
        GPLs=varargin{1};
        WordList={};
        WordNb=[];
        WordSize=[];
        GseNb=0;
        FoundGseNb=0;
        if ExpFlag
            WordByGse=sparse([]);
        end
        
        %process each GPL
        for GplL=1:length(GPLs)
            %search if current GPL exists in K.chip
            ChipRank=strmatch(sprintf('GPL%u',GPLs(GplL)),K.chip.geoName,'exact');
            if isempty(ChipRank)
                h=warndlg(sprintf('GPL%u does not exists',GPLs(GPL)));
                waitfor(h)
            else
                GPL=sprintf('GPL%u.mat',GPLs(GplL));
                cd(K.dir.geoMetadata)
                if ~exist(GPL,'file')
                    h=warndlg(sprintf('%s does not exists',GPL));
                    waitfor(h)
                else
                    eval(sprintf('load %s',GPL))
                    GseNb=GseNb+length(Gse.summary);
                    %process each GSE of the current GPL
                    for GseL=1:length(Gse.summary)
                        %for GseL=1:20
                        String=[Gse.title{GseL},' ',Gse.summary{GseL},' ',Gse.overallDesign{GseL}];
                        if ~isempty(String)
                            if ExpFlag
                                if isempty(regexp(String,Exp,'match'))
                                    String='';
                                else
                                    FoundGseNb=FoundGseNb+1;
                                end
                            end
                            if ~isempty(String)                               
                                [CurrWordList,CurrWordNb,CurrWordSize]=recover_words(String,' ',1,0);
                                if isempty(WordList)
                                    WordList=CurrWordList;
                                    WordNb=CurrWordNb;
                                    WordSize=CurrWordSize;
                                    if ExpFlag
                                        GplList=GPLs(GplL);
                                        GseList=Gse.gseRank(GseL);
                                        WordByGse=ones(length(WordList),1);
                                    end
                                else
                                    %merge with previous GSE analysis
                                    %Intersection
                                    [KnownWord,Index1,Index2]=intersect(WordList,CurrWordList);
                                    if ExpFlag
                                        GplList(end+1)=GPLs(GplL);
                                        GseList(end+1)=Gse.gseRank(GseL);
                                        GsePos=size(WordByGse,2)+1;
                                    end
                                    if ~isempty(KnownWord)
                                        WordNb(Index1)=WordNb(Index1)+CurrWordNb(Index2);
                                        if ExpFlag
                                            WordByGse(Index1,GsePos)=1;
                                        end
                                    end
                                    %union
                                    [NewWord,Index]=setdiff(CurrWordList,WordList);
                                    if ~isempty(NewWord)
                                        WordList=[WordList;CurrWordList(Index,:)];
                                        WordNb=[WordNb;CurrWordNb(Index)];
                                        WordSize=[WordSize;CurrWordSize(Index)];
                                        if ExpFlag
                                            WordByGse(end+1:end+length(Index),GsePos)=1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        [WordList,Index]=sort(WordList);
        WordNb=WordNb(Index);
        WordSize=WordSize(Index);
        if ExpFlag
            WordByGse=WordByGse(Index,:);
        end
        if ExpFlag
            [WordList Index]=setdiff(WordList,ClearedWords);
            WordNb=WordNb(Index);
            WordSize=WordSize(Index);
            WordByGse=WordByGse(Index,:);
            [temp,Index1]=intersect(WordList,Words);
            [temp,Index2]=intersect(WordList,Symbols);
            Index1=setdiff(Index1,Index2);
            Index3=setdiff([1:length(WordList)],union(Index1,Index2));
            WordList=[WordList(Index1);WordList(Index2);WordList(Index3)];
            WordNb=[WordNb(Index1);WordNb(Index2);WordNb(Index3)];
            WordSize=[WordSize(Index1);WordSize(Index2);WordSize(Index3)];
            WordByGse=[WordByGse(Index1,:);WordByGse(Index2,:);WordByGse(Index3,:)];
            
            %reorder columns
            ItemNb=sum(WordByGse);
            [ItemNb SortIndex]=sort(ItemNb);
            ItemNb=fliplr(ItemNb);
            SortIndex=fliplr(SortIndex);
            WordByGse=WordByGse(:,SortIndex);
            GplList=GplList(SortIndex);
            GseList=GseList(SortIndex);            
            

            cd(K.dir.geoMetadata)
            FileNb=ceil((size(WordByGse,2)+4)/1000);
            for FileL=1:FileNb
                Range=[(FileL-1)*1000+1:min(size(WordByGse,2),FileL*1000)];
                fid=fopen(sprintf('gsewords_%s_%s_%u.txt',SearchedName,Date,FileL),'w');
                fprintf(fid,['\t\t\t\t',repmat('\t%u',1,length(Range)),'\n'],Range);
                fprintf(fid,['\t\t\t\tGPL',repmat('\t%u',1,length(Range)),'\n'],GplList(Range));
                fprintf(fid,['\t\t\t\tGSE',repmat('\t%u',1,length(Range)),'\n'],GseList(Range));
                fprintf(fid,['\t\t\t\tsum',repmat('\t%u',1,length(Range)),'\n'],ItemNb(Range));
                fprintf(fid,'word/symbol\t#found\t#total\t%%\n')

                ColNb=length(Range);
                Format=['%s',repmat('\t%u',1,ColNb+4),'\n'];
                for WordL=1:length(Index1)
                    CurrWord=WordList{WordL};
                    WordPos=strmatch(CurrWord,Words,'exact');
                    fprintf(fid,Format,CurrWord,WordNb(WordL),WordsNb(WordPos),round(WordNb(WordL)*100/WordsNb(WordPos)),WordSize(WordL),WordByGse(WordL,Range));
                end
                for WordL=length(Index1)+1:length(Index1)+length(Index2)
                    CurrWord=WordList{WordL};
                    WordPos=strmatch(CurrWord,Symbols,'exact');
                    fprintf(fid,Format,CurrWord,WordNb(WordL),SymbolsNb(WordPos),round(WordNb(WordL)*100/SymbolsNb(WordPos)),WordSize(WordL),WordByGse(WordL,Range));
                end
                for WordL=length(Index1)+length(Index2)+1:length(Index1)+length(Index2)+length(Index3)
                    fprintf(fid,Format,WordList{WordL},WordNb(WordL),0,0,WordSize(WordL),WordByGse(WordL,Range));
                end
                fclose(fid)
            end
        else
            Status=zeros(length(WordList),1);
            if WordFlag

                [temp Index]=intersect(WordList,Words);
                Status(Index)=1;
                [temp Index]=intersect(WordList,Symbols);
                Status(Index)=2;
                [temp Index]=intersect(WordList,ClearedWords);
                Status(Index)=3;
            end
            cd(K.dir.geoMetadata)
            eval(sprintf('save gsewords_%s GPLs WordList WordNb WordSize Status;',Date))
            %fid=fopen(sprintf('gsewords_%s.txt',Date),'w');
            fid=fopen(sprintf('gse_overalldesign_%s.txt',Date),'w');
            for WordL=1:length(WordList)
                fprintf(fid,'%s\t%u\t%u\t%u\n',WordList{WordL},WordNb(WordL),WordSize(WordL),Status(WordL));
            end
            fclose(fid)
        end
        
        feature('DefaultCharacterSet','windows-1252')
        
%% PRINT GSEs
    case 'print GSEs'
        if nargin<4
            h=errordlg('print GSEs needs  a list of GPL in numerical format, a list of Keywords in cell format {'','','',...} and a name for output file');
            waitfor(h)
            error('process canceled')                
        end
        GPLs=varargin{1};
        Exp=varargin{2};
        FileName=varargin{3};
        %scan each GPL
        cd(K.dir.geoMetadata)
        fid=fopen(sprintf('%s_%s.txt',FileName,Date),'w');
        %print date
        fprintf(fid,'%s\n',date);
        %print regular expression
        fprintf(fid,'%s\n',Exp);
        fprintf(fid,'\n\n\n');
        for GplL=1:length(GPLs)
            ChipRank=strmatch(sprintf('GPL%u',GPLs(GplL)),K.chip.geoName,'exact');
            if isempty(ChipRank)
                h=warndlg(sprintf('GPL%u does not exists',GPLs(GplL)));
                waitfor(h)
            else
                %information on current GPL
                fprintf(fid,'GPL%u: %s (m%u)\n\n',GPLs(GplL),K.chip.name{ChipRank},ChipRank);
                GPL=sprintf('GPL%u.mat',GPLs(GplL));
                cd(K.dir.geoMetadata)
                if ~exist(GPL,'file')
                    h=warndlg(sprintf('%s does not exists',GPL));
                    fprintf(fid,'no GPL.mat file for this chip\n\n');
                    waitfor(h)
                else
                    eval(sprintf('load %s',GPL))
                end
                FoundGse=0;
                FieldsFlag=0;
                for GseL=1:length(Gse.summary)
                    String=[Gse.title{GseL},' ',Gse.summary{GseL},' ',Gse.overallDesign{GseL}];
                    if ~isempty(String)
                        if ~isempty(regexp(String,Exp,'match'))
                            FoundGse=1;
                            if FieldsFlag==0
                                FieldsFlag=1;
                                %print current GSEs information
                                fprintf(fid,'gse\t');
                                fprintf(fid,'title\t');
                                fprintf(fid,'status\t');
                                fprintf(fid,'submited\t');
                                fprintf(fid,'updated\t');
                                fprintf(fid,'pubmed\t');
                                fprintf(fid,'summary\t');
                                fprintf(fid,'type\t');
                                fprintf(fid,'url\t');
                                fprintf(fid,'design\t');
                                fprintf(fid,'repeats\t');
                                fprintf(fid,'list\t');
                                fprintf(fid,'variable\t');
                                fprintf(fid,'description\t');
                                fprintf(fid,'contact\t');
                                fprintf(fid,'supplement\t');
                                fprintf(fid,'gsmnb\t');
                                fprintf(fid,'gsmcount\t');
                                fprintf(fid,'outgpl\t\n');
                            end
                            %PRINT GSE INFORMATION
                            fprintf(fid,'%u\t',Gse.gseRank(GseL));
                            fprintf(fid,'%s\t',strrep(Gse.title{GseL},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gse.status{GseL},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gse.submissionDate{GseL},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gse.lastUpdateDate{GseL},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gse.pubmedId{GseL},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gse.summary{GseL},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gse.type{GseL},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gse.webLink{GseL},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gse.overallDesign{GseL},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gse.repeats{GseL},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gse.repeatsSampleList{GseL},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gse.variable{GseL},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gse.variableDescription{GseL},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gse.contact{GseL},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gse.supplementaryFile{GseL},sprintf('\t'),' '));
                            fprintf(fid,'%u\t',Gse.gsmNb(GseL));
                            fprintf(fid,'%u\t',Gse.gsmCount(GseL));
                            fprintf(fid,'%u\t',Gse.outGplGsmNb(GseL));
                            fprintf(fid,'\n')
                        end %of isempty(regexp(String,Exp,'match'))
                    end
                end %of GseL
                if FoundGse==0
                    fprintf(fid,'no keywords in this chip\n\n');
                end
            end %isempty ChipRank
        end %of GplL
        fclose(fid)



%% PRINT A GPL
    case 'print a GPL'
        Ok=VERIF();
        if nargin<2
            h=errordlg('misses GPL');
            waitfor(h)
            error('process canceled')
        end
        GPL=varargin{1};
        %modify eventuelly to use the right format
        if isnumeric(GPL)
            GPL=sprintf('GPL%u',GPL);
        end
        if nargin==3
            ControlFlag=varargin{2};
        end
        if Ok
            cd(K.dir.geoMetadata)
            eval(sprintf('load %s',GPL))

            %RECOVER SPECIES
            eval(sprintf('Species=mksqlite(''SELECT organism FROM gpl WHERE gpl=''''%s'''''');',GPL))

            %RECOVER FACTOR NAMES
            FactorNames={};
            FactorValues={};
            for GsmL=1:length(Gsm.gsm)
                if length(Gsm.factorNames)>=GsmL
                    if ~isempty(Gsm.factorNames{GsmL})
                        for FactorL=1:length(Gsm.factorNames{GsmL})
                            FactorNames{end+1,1}=Gsm.factorNames{GsmL}{FactorL};
                            FactorValues{end+1,1}=Gsm.factorValues{GsmL}{FactorL};
                        end
                    end
                else
                    break
                end
            end
            FactorNames=unique(FactorNames);
            FValues={};
            for FactorL=1:length(FactorValues)
                if iscell(FactorValues{FactorL})
                    for i=1:length(FactorValues{FactorL})
                        FValues{end+1,1}=FactorValues{FactorL}{i};
                    end
                else
                    FValues{end+1,1}=FactorValues{FactorL};
                end
            end
            FactorValues={};
            FactorValues=unique(FValues);

            %PRINT
            fid=fopen(sprintf('%s_%s_%s.txt',strrep(Species.organism,' ','_'),GPL,VERSION),'w');

            %write header
            %species name
            fprintf(fid,'%s\t%s\n',Species.organism,VERSION);
            %field names
            fprintf(fid,'GSE')
            for i=1:19
                fprintf(fid,'\t')
            end
            fprintf(fid,'GSM')
            for i=1:12
                fprintf(fid,'\t')
            end
            fprintf(fid,'FACTORS\n')

            %GSE
            fprintf(fid,'gse\t');
            fprintf(fid,'title\t');
            fprintf(fid,'status\t');
            fprintf(fid,'submited\t');
            fprintf(fid,'updated\t');
            fprintf(fid,'pubmed\t');
            fprintf(fid,'summary\t');
            fprintf(fid,'type\t');
            fprintf(fid,'url\t');
            fprintf(fid,'design\t');
            fprintf(fid,'repeats\t');
            fprintf(fid,'list\t');
            fprintf(fid,'variable\t');
            fprintf(fid,'description\t');
            fprintf(fid,'contact\t');
            fprintf(fid,'supplement\t');
            fprintf(fid,'gsmnb\t');
            fprintf(fid,'gsmcount\t');
            fprintf(fid,'outgpl\t');
            %GSM
            fprintf(fid,'gsm\t');
            fprintf(fid,'title\t');
            fprintf(fid,'source\t');
            fprintf(fid,'characteristics\t');
            fprintf(fid,'molecule\t');
            fprintf(fid,'label\t');
            fprintf(fid,'treatment\t');
            fprintf(fid,'extraction\t');
            fprintf(fid,'labelling\t');
            fprintf(fid,'hybridation\t');
            fprintf(fid,'description\t');
            fprintf(fid,'algorithm');
            %FACTORS
            for FactorL=1:length(FactorNames)
                fprintf(fid,'\t%s',FactorNames{FactorL});
            end
            fprintf(fid,'\n');


            [GseRanks,GseOrder]=sort(Gse.gseRank);
            %GseRanks=Gse.gseRank;
            %GseOrder=1:length(Gse.gse);
            %for ExpL=1:1
            for ExpL=1:length(Gse.gse)
                GseL=GseOrder(ExpL);
                if GseRanks(ExpL)>0
                    %PRINT GSE INFORMATION
                    fprintf(fid,'%u\t',Gse.gseRank(GseL));
                    fprintf(fid,'%s\t',strrep(Gse.title{GseL},sprintf('\t'),' '));
                    fprintf(fid,'%s\t',strrep(Gse.status{GseL},sprintf('\t'),' '));
                    fprintf(fid,'%s\t',strrep(Gse.submissionDate{GseL},sprintf('\t'),' '));
                    fprintf(fid,'%s\t',strrep(Gse.lastUpdateDate{GseL},sprintf('\t'),' '));
                    fprintf(fid,'%s\t',strrep(Gse.pubmedId{GseL},sprintf('\t'),' '));
                    fprintf(fid,'%s\t',strrep(Gse.summary{GseL},sprintf('\t'),' '));
                    fprintf(fid,'%s\t',strrep(Gse.type{GseL},sprintf('\t'),' '));
                    fprintf(fid,'%s\t',strrep(Gse.webLink{GseL},sprintf('\t'),' '));
                    fprintf(fid,'%s\t',strrep(Gse.overallDesign{GseL},sprintf('\t'),' '));
                    fprintf(fid,'%s\t',strrep(Gse.repeats{GseL},sprintf('\t'),' '));
                    fprintf(fid,'%s\t',strrep(Gse.repeatsSampleList{GseL},sprintf('\t'),' '));
                    fprintf(fid,'%s\t',strrep(Gse.variable{GseL},sprintf('\t'),' '));
                    fprintf(fid,'%s\t',strrep(Gse.variableDescription{GseL},sprintf('\t'),' '));
                    fprintf(fid,'%s\t',strrep(Gse.contact{GseL},sprintf('\t'),' '));
                    fprintf(fid,'%s\t',strrep(Gse.supplementaryFile{GseL},sprintf('\t'),' '));
                    fprintf(fid,'%u\t',Gse.gsmNb(GseL));
                    fprintf(fid,'%u\t',Gse.gsmCount(GseL));
                    fprintf(fid,'%u\t',Gse.outGplGsmNb(GseL));
                    fprintf(fid,'\n')
                    %PRINT GSM INFORMATION
                    GsmPos=find(Gsm.gseRank==Gse.gseRank(GseL));
                    if ~isempty(GsmPos)
                        for GsmL=1:length(GsmPos)
                            for i=1:19
                                fprintf(fid,'\t');
                            end
                            fprintf(fid,'%u\t',Gsm.gsmRank(GsmPos(GsmL)));
                            fprintf(fid,'%s\t',strrep(Gsm.title{GsmPos(GsmL)},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gsm.sourceNameCh1{GsmPos(GsmL)},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gsm.characteristicsCh1{GsmPos(GsmL)},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gsm.moleculeCh1{GsmPos(GsmL)},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gsm.labelCh1{GsmPos(GsmL)},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gsm.treatmentProtocolCh1{GsmPos(GsmL)},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gsm.extractProtocolCh1{GsmPos(GsmL)},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gsm.labelProtocolCh1{GsmPos(GsmL)},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gsm.hybProtocol{GsmPos(GsmL)},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gsm.description{GsmPos(GsmL)},sprintf('\t'),' '));
                            fprintf(fid,'%s\t',strrep(Gsm.dataProcessing{GsmPos(GsmL)},sprintf('\t'),' '));

                            %PRINT FACTORS IF EXIST
                            if length(Gsm.factorNames)>=GsmPos(GsmL)
                                if ~isempty(Gsm.factorNames{GsmPos(GsmL)})
                                    [Factors,Order]=sort(Gsm.factorNames{GsmPos(GsmL)});
                                    Values=Gsm.factorValues{GsmPos(GsmL)}(Order);
                                    ListPos=1;
                                    for FactorL=1:length(Factors)
                                        for i=ListPos:length(FactorNames)
                                            if isequal(Factors{FactorL},FactorNames{i})
                                                if iscell(Values{FactorL})
                                                    fprintf(fid,'%s',strrep(Values{FactorL}{1},sprintf('\t'),' '));
                                                    for j=2:size(Values{FactorL},1)
                                                        fprintf(fid,' & %s',strrep(Values{FactorL}{j},sprintf('\t'),' '));
                                                    end
                                                else
                                                    fprintf(fid,'%s',strrep(Values{FactorL},sprintf('\t'),' '));
                                                end
                                                if FactorL==length(Factors)
                                                    fprintf(fid,'\n');
                                                else
                                                    fprintf(fid,'\t');
                                                end
                                                ListPos=i+1;
                                                break
                                            else
                                                fprintf(fid,'\t');
                                            end
                                        end
                                    end
                                else
                                    for i=1:length(FactorNames)-1
                                        fprintf(fid,'\t');
                                    end
                                    fprintf(fid,'\n');
                                end
                            else
                                for i=1:length(FactorNames)-1
                                    fprintf(fid,'\t');
                                end
                                fprintf(fid,'\n');
                            end
                        end % of GsmL
                    end %if ~isempty
                end %if GseRank>0
            end %of GseL
            fclose(fid)
        end %of Ok
        
%% PRINT A SCRIPT TO IMPORT GSE 
    case 'import gse'

        VERIF()
        if nargin~=2
            h=errordlg('print import script needs GPL parameter')
            waitfor(h)
            error('process canceled')
        end

        GPL=varargin{1}
        if isnumeric(GPL)
            GPL=sprintf('GPL%u',GPL);
        end
        cd(K.dir.geoMetadata')
        eval(sprintf('load %s',GPL))

        fid=fopen(sprintf('load_%s.h',GPL),'w');

        fprintf(fid,'#!/bin/bash\n');
        fprintf(fid,sprintf('cd /usr/data/experiments/geo/GPL%s\n',GPL));
        SupNb=0;
        for GseL=1:length(Gse.gse)
            if ~isempty(Gse.supplementaryFile)
                if ischar(Gse.supplementaryFile{GseL})
                    Gse.gse{GseL}
                    SupNb=SupNb+1;
                    fprintf(fid,sprintf('ncftpget ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/series/%s/%s_RAW.tar\n',Gse.gse{GseL},Gse.gse{GseL}));
                end
            end
        end
        fclose(fid)
        h=warndlg(sprintf('exist %u GSE with CEL files among %u',SupNb,length(Gse.gse)));
        waitfor(h)
%% MODIF VERSION
    case 'modif version'
        VERIF()
        if nargin~=2
            h=errordlg('print import script needs GPL parameter')
            waitfor(h)
            error('process canceled')
        end

        GPL=varargin{1};
        if isnumeric(GPL)
            GPL=sprintf('GPL%u',GPL);
        end
        cd(K.dir.geoMetadata')
        eval(sprintf('load %s',GPL))

        for GdsL=1:length(Gds.gsmRank)
            Gds.metadb{GdsL}=VERSION;
        end
        for GseL=1:length(Gse.gse)
            Gse.metadb{GseL}=VERSION;
            Gse.isBiol(GseL,1)=0;
        end
        for GsmL=1:length(Gsm.gsm)
            Gsm.metadb{GsmL}=VERSION;
            Gsm.isBiol(GsmL,1)=0;
        end
        for ContL=1:length(Contributor.firstNames)
            Contributor.metadb{GsmL}=VERSION;
        end

        eval(sprintf('save %s Gsm Gse Gds Contributor',GPL))



%% FIND BIOLOGICAL CONDITIONS
    case 'find biological conditions'
        %VERIF()
        if nargin<2
            h=errordlg('print import script needs GPL parameter')
            waitfor(h)
            error('process canceled')
        end
        ReplaceFlag=0;
        if nargin==3
            ReplaceFlag=1;
            for VarL=1:length(varargin{2})
                Replace{VarL}=sprintf('(^|_|,\\s|,|\\s)+%s(_|\\s|,\\s|,|$)+',varargin{2}{VarL});
            end
        end

        GPL=varargin{1};
        if isnumeric(GPL)
            GPL=sprintf('GPL%u',GPL);
        end
        cd(K.dir.geoMetadata')
        eval(sprintf('load %s',GPL))
        
        
        
        Update=questdlg('Do you want to overwrite existing bioliogical conditions informations?','','no','yes','no');
        if isequal(Update,'yes')
           Gse.isBiol=zeros(length(Gse.gse),1);
           Gse.isFactor=zeros(length(Gse.gse),1);
           Gsm.isBiol=zeros(length(Gsm.gsm),1);
           Gsm.isFactor=zeros(length(Gsm.gsm),1);
           Gsm.replicate=zeros(length(Gsm.gsm),1);
           for GsmL=1:length(Gsm)
               Gsm.biolName{GsmL}='';        
           end               
           
        end

        %process GSM with factors
        for GseL=1:length(Gse.gse)
            if Gse.isBiol(GseL)==0
                GsmPos=find(Gsm.gseRank==Gse.gseRank(GseL));
                if ~isempty(GsmPos)
                    %verify if exist some factor
                    if length(Gsm.factorNames)>=GsmPos(end)
                        if length(Gsm.factorNames{GsmPos(1)})>0
                            %construct compound biological name                            
                            BiolNames={};
                            GsmBiolNames=cell(length(GsmPos),1);
                            for GsmL=1:length(GsmPos)
                                BiolName='*';
                                for FactorL=1:length(Gsm.factorNames{GsmPos(GsmL)})
                                    BiolName=[BiolName,strrep(Gsm.factorValues{GsmPos(GsmL)}{FactorL},' ','-'),'*'];
                                end
                                GsmBiolNames{GsmL}=BiolName;
                                if isempty(strmatch(BiolName,BiolNames,'exact'))
                                    BiolNames{end+1,1}=BiolName;
                                end
                            end
                            %find unique biol names
                            BiolNames=unique(BiolNames);
                            Replicate=zeros(length(BiolNames),1);
                            AllFound=1;
                            for GsmL=1:length(GsmPos)
                                if isequal(GsmBiolNames{GsmL},'*')
                                    Gsm.biolName{GsmPos(GsmL),1}='';
                                    Gsm.isBiol(GsmPos(GsmL))=0;
                                    Gsm.isFactor(GsmPos,1)=0;                                
                                    AllFound=0;
                                else
                                    Gsm.biolName{GsmPos(GsmL),1}=GsmBiolNames{GsmL};
                                    BiolPos=strmatch(GsmBiolNames{GsmL},BiolNames,'exact');
                                    Replicate(BiolPos)=Replicate(BiolPos)+1;
                                    Gsm.replicate(GsmPos(GsmL),1)=Replicate(BiolPos);
                                    Gsm.isFactor(GsmPos,1)=1;                                   
                                    Gsm.isBiol(GsmPos(GsmL))=1;
                                end
                            end

                            %fprintf(fid,'%s\n',Gse.gse{GseL});
                            if AllFound
                                Gse.isBiol(GseL)=1;                                
                                Gse.isFactor(GseL)=1;
%                                 for GsmL=1:length(GsmPos)
%                                     fprintf(fid,'\t%s\t0\t%s\t%s\t%u\n',Gsm.gsm{GsmPos(GsmL)},Gsm.title{GsmPos(GsmL)},Gsm.biolName{GsmPos(GsmL)},Gsm.replicate(GsmPos(GsmL)));
%                                 end
                            else
                                Gse.isBiol(GseL)=-1;
                                Gse.isFactor(GseL)=0;
                                %fprintf(fid,'some GSM have no factor value and are processed below\n')
                            end
                        else
                            Gsm.isFactor(GsmPos,1)=0;
                            Gse.isFactor(GseL)=0;
                        end
                    else
                        Gsm.isFactor(GsmPos,1)=0;                        
                        Gse.isFactor(GseL)=0;
                    end
                end
            end
        end

        %uses regular expressions for GSM without factors
        RegExprs={'(_\d+,\s|,\s|,|\s/|/|\s-|\s_|-|_|\s)+(biological)?(\s|_|-)*rep\w*(\s|No.\s*)*\d+(,\s|,|\s|-|_|$)+';
            '(1st|2nd|3rd)(biological)?(\s|_)*rep\w*(\s|_|$)+';
            '(#\s|#|-|_|,\s|,|\s)*\d+$';
            '(,\s|,|\s)sample(\s)*\d+(,\s|,|\s|-|_|$)+';
            '(,\s|,|\s)No.(\s)*\d+(,|\s|-|_|$)+'};       
        RegNb=zeros(length(RegExprs),1);
        UsedNb=zeros(length(RegExprs),1);
        NotFound=[];
        Verif=[];
        
        
        for RegL=1:length(RegExprs)
            %a Gse is marked isFactor=0 if some Gsm do not have factor and some other have factors
            %GSM with factors have foundBiol=1
            GsePos=find(Gse.isBiol==0&Gse.isFactor==0);
            if isempty(GsePos)
                break
            else
                for GseL=1:length(GsePos)
                    GsmPos=find(Gsm.gseRank==Gse.gseRank(GsePos(GseL)));
                    if ~isempty(GsmPos)
                        GsmBiolNames=cell(length(GsmPos),1);
                        Result=cell(length(GsmPos),1);
                        FoundReg=0;
                        NotFound=0;
                        for GsmL=1:length(GsmPos)
                            GsmBiolNames{GsmL}=Gsm.title{GsmPos(GsmL)};
                            if Gsm.isBiol(GsmPos(GsmL))~=1
                                if ReplaceFlag
                                    for RepL=1:length(Replace)
                                        GsmBiolNames{GsmL}=regexprep(GsmBiolNames{GsmL},Replace{RepL},'');
                                    end
                                end
                                %replace ExpName_ExpNb-PointNb_ current in some experiment names by ''%
                                %GsmBiolNames{GsmL}=regexprep(GsmBiolNames{GsmL},'^\w+_(A|A-)?\d+(-|_)(\d+(\s|_|-)?)?','');
                                GsmBiolNames{GsmL}=regexprep(GsmBiolNames{GsmL},'^\w+\s?\w+_\d+(-|_)\d+(\s|_|-)?','');
                                %replace - by _ to catch some common error
                                %(JA-1 et JA_1 indicating the same biological condition)
                                GsmBiolNames{GsmL}=regexprep(GsmBiolNames{GsmL},'-','_');
                                [start,stop,match]=regexpi(GsmBiolNames{GsmL},RegExprs{RegL},'start','end','match');
                                if ~isempty(start)
                                    FoundReg=1;
                                    RegNb(RegL)=RegNb(RegL)+1;
                                    %replace match by empty string and ignore
                                    %trailing spaces                                   
                                    GsmBiolNames{GsmL}=regexprep(regexprep(GsmBiolNames{GsmL},RegExprs{RegL},'','ignorecase'),'\s+$','');
                                    Result{GsmL}=sprintf('%s *** %s',GsmBiolNames{GsmL},match{1});
                                else
                                    NotFound=1;
                                    %assume that it is a singleton (no
                                    %replicate)
                                    GsmBiolNames{GsmL}=regexprep(GsmBiolNames{GsmL},'\s+$','');
                                    Result{GsmL}=sprintf('NOT FOUND in %s',GsmBiolNames{GsmL});
                                end
                            end
                        end
                        if FoundReg
                            %display results
                            %Sel=listdlg('liststring',Result,'promptstring',sprintf('%s reg %s - select one if list is correct',Gse.gse{GsePos(GseL)},RegExprs{RegL}),'listsize',[500,500]);
                            %if ~isempty(Sel)&NotFound==0
                                %record biological conditions
                                BiolNames=unique(GsmBiolNames);
                                Replicate=zeros(length(BiolNames),1);
                                Result=cell(length(GsmPos),1);
                                %fprintf(fid,'%s\n',Gse.gse{GsePos(GseL)});
                                for GsmL=1:length(GsmPos)
                                    if Gsm.isBiol(GsmPos(GsmL))~=1
                                        Gsm.biolName{GsmPos(GsmL),1}=GsmBiolNames{GsmL};
                                        BiolPos=strmatch(GsmBiolNames{GsmL},BiolNames,'exact');
                                        Replicate(BiolPos)=Replicate(BiolPos)+1;
                                        Gsm.replicate(GsmPos(GsmL),1)=Replicate(BiolPos);
                                        Gsm.isBiol(GsmPos(GsmL),1)=1;
                                        Result{GsmL}=sprintf('%s _ rep%u',GsmBiolNames{GsmL},Replicate(BiolPos));
                                        %                                     fprintf(fid,'\t%s\t%u\t%s\t%s\t%u\n',Gsm.gsm{GsmPos(GsmL)},NotFound,Gsm.title{GsmPos(GsmL)},GsmBiolNames{GsmL},Replicate(BiolPos));
                                    end
                                end
                                Gse.isBiol(GsePos(GseL))=1;
                                '========================='
                                Gse.gse{GsePos(GseL)}
                                Result
                                '========================='
                                if NotFound
                                    Verif=[Verif;GsePos(GseL)];
                                end
                                UsedNb(RegL)=UsedNb(RegL)+length(GsmPos);
                                NotFound=setdiff(NotFound,GsePos(GseL));
                            %end
                        else
                            NotFound=[NotFound;GsePos(GseL)];
                        end
                    end
                end
            end
        end
        %add NotFound
%         if ~isempty(NotFound)
%             for NotL=1:length(NotFound)
%                 GsmPos=find(Gsm.gseRank==NotFound(NotL));
%                 if ~isempty(GsmPos)
%                     fprintf(fid,'%s\n',Gse.gse{NotFound(NotL)});
%                     for GsmL=1:length(GsmPos)
%                         fprintf(fid,'\t\t%s\n',Gsm.title{GsmPos(GsmL)})
%                     end
%                 else
%                     fprintf(fid,'%s has no GSM\n',Gse.gse{NotFound(NotL)});
%                 end
%             end
%         end
        
        fid=fopen(sprintf('BiolCond_%s_%s.txt',GPL,date),'w');        
        GsePos=find(Gse.isFactor);
        if~isempty(GsePos)
            fprintf(fid,'GSM WITH FACTORS\n\n');
            fprintf(fid,'GSE\tGSM\tVerif\tSampleName\tBiolName\tReplicate\n\n');
            for GseL=1:length(GsePos)
                fprintf(fid,'%s\n',Gse.gse{GsePos(GseL)});
                GsmPos=find(Gsm.gseRank==Gse.gseRank(GsePos(GseL)));
                if ~isempty(GsmPos)
                    for GsmL=1:length(GsmPos)
                        if Gsm.isFactor(GsmPos(GsmL))
                            fprintf(fid,'\t%s\t0\t%s\t%s\t%u\n',Gsm.gsm{GsmPos(GsmL)},Gsm.title{GsmPos(GsmL)},Gsm.biolName{GsmPos(GsmL)},Gsm.replicate(GsmPos(GsmL)));
                        else
                            fprintf(fid,'\t%s\t1\t%s\t%s\t%u\n',Gsm.gsm{GsmPos(GsmL)},Gsm.title{GsmPos(GsmL)},Gsm.biolName{GsmPos(GsmL)},Gsm.replicate(GsmPos(GsmL)));
                        end
                    end
                else
                    fprintf(fid,'\tNO GSM\n');
                end
            end
        end

        GsePos=find(Gse.isFactor==0);
        if~isempty(GsePos)

            fprintf(fid,'\n\n');
            fprintf(fid,'GSM WITHOUT FACTORS\n\n');
            fprintf(fid,'Used regular expressions\n');
            for RegL=1:length(RegExprs)
                fprintf(fid,'%s\n',RegExprs{RegL});
            end
            fprintf(fid,'\n\n');

            fprintf(fid,'GSE\tGSM\tVerif\tSampleName\tBiolName\tReplicate\n\n');
            for GseL=1:length(GsePos)
                fprintf(fid,'%s\n',Gse.gse{GsePos(GseL)});
                GsmPos=find(Gsm.gseRank==Gse.gseRank(GsePos(GseL)));
                if ~isempty(GsmPos)
                    for GsmL=1:length(GsmPos)
                        if Gsm.isBiol(GsmPos(GsmL))
                            fprintf(fid,'\t%s\t0\t%s\t%s\t%u\n',Gsm.gsm{GsmPos(GsmL)},Gsm.title{GsmPos(GsmL)},Gsm.biolName{GsmPos(GsmL)},Gsm.replicate(GsmPos(GsmL)));
                        else
                            fprintf(fid,'\t%s\t1\t%s\t%s\t%u\n',Gsm.gsm{GsmPos(GsmL)},Gsm.title{GsmPos(GsmL)},Gsm.biolName{GsmPos(GsmL)},Gsm.replicate(GsmPos(GsmL)));
                        end
                    end
                else
                    fprintf(fid,'\tNO GSM\n');
                end
            end
        end

        fclose(fid)
        eval(sprintf('save %s Gsm Gse Gds Contributor',GPL))

%% IMPORT BIOLOGICAL CONDITIONS        
    case 'import biological conditions'
        
        cd(K.dir.geoMetadata')
        [GplFile,GplDir]=uigetfile('*.mat','select a GPL file');
        GPL=regexp(GplFile,'^.*(?=.mat)','match');
        GPL=GPL{1};
        cd(GplDir)
        load(GplFile)
        cd(K.dir.geoMetadata)
        [FileName,FileDir]=uigetfile('BiolCond*.txt','select the biological name list');
        cd(FileDir)
        [GsmName,BiolName,Rep]=textread(FileName,'%s%s%u','delimiter','\t'); 
        for BiolL=1:length(GsmName)
            GsmPos=strmatch(GsmName{BiolL},Gsm.gsm,'exact');
            if ~isempty(GsmPos)
                Gsm.biolName{GsmPos}=BiolName{BiolL};
                Gsm.replicate(GsmPos)=Rep(BiolL);
                Gsm.isBiol(GsmPos)=1;
                GsePos=find(Gse.gseRank==Gsm.gseRank(GsmPos));
                if ~isempty(GsePos)
                Gse.isBiol(GsePos)=1;
                else
                    h=errordlg(sprintf('GSE%u does not exist',Gsm.gseRank(GsmPos)));
                    waitfor(h)
                    error('process canceled')
                end
            else
                h=warndlg(sprintf('GSM%u does not exist',GsmName{BiolL}));
                waitfor(h)                
            end
        end
        cd(GplDir)
        eval(sprintf('save %s Gsm Gse Gds Contributor',GPL))



end %switch Action



%% function VERIF()
function varargout=VERIF()
global K
Ok=0;
if isfield(K,'tmp')
    if isfield(K.tmp,'geoDbId')
        if K.tmp.geoDbId>0
            Ok=1;
        end
    end
end

if Ok==0
    geo_metadb('open geo metadb')
    Ok=1; 
end
if nargout==1
    varargout{1}=Ok;
end

%% function PRINT()
function PRINT(Fid,Selection)
Fields=fields(Selection);
for FieldL=1:length(Fields)
    eval(sprintf('Value=Selection.%s;',Fields{FieldL}))
    if ~isempty(Value)
        if isequal(class(Value),'char')
            fprintf(Fid,'%s = %s\n',Fields{FieldL},Value);
        else
            fprintf(Fid,'%s = %f\n',Fields{FieldL},Value);
        end
    else
        fprintf(Fid,'%s = []\n',Fields{FieldL});
    end
end
fprintf(Fid,'\n*******************************************\t')

