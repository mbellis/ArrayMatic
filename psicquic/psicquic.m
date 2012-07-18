%====================%
% FUNCTION PSICQUIC  %
%====================%

% PSICQUIC recover Psicquic information

%INPUT PARAMETERS

%  1     Action: subprogram used
%                psimat: interrogate PSICQUIC databases individually and recover results
%                        in CurrInfo variable
%                analysis: eventually merge several CurrInfo in a single Info file which 
%                          is processed to recover restults into matrix and structures 
%                score: display scores
%                output: write results
%  2   ListType: type of list:
%                unique : a one column text file with gene name (select all interactors
%                         for each gene)       
%                pair : a two column text fiel with gene name (select only paired
%                       interactors)   
%                direct : a list of gene passed as parameter (replace File) 
%  3        Dir: directory to load and save files
%  4       File: 
%  5     Output:
%  6 varargin[1}: for  action psimat: First gene of the the range of genes processed (Start)
%                 for action analysis: First genes (Starts) if exists several files 
%                 for action score: name of Authors whose interactions must be deleted
%  7 varargin[2}: for  action psimat: Last gene of the the range of genes processed (Stop)
%                 for action analysis: Last genes (Stops) if exists several files 

% MITAB FORMAT
%  1 Unique identifier for interactor A, represented as databaseName:ac, where databaseName is
%    the name of the corresponding database as defined in the PSI-MI controlled vocabulary,
%    and ac is the unique primary identifier of the molecule in the database. Identifiers from
%    multiple databases can be separated by "|". It is recommended that proteins be identified
%    by stable identifiers such as their UniProtKB or RefSeq accession number.
%  2 Unique identifier for interactor B.
%  3 Alternative identifier for interactor A, for example the official gene symbol as defined
%    by a recognised nomenclature committee. Representation as databaseName:identifier.
%    Multiple identifiers separated by "|".
%  4 Alternative identifier for interactor B.
%  5 Aliases for A, separated by "|". Representation as databaseName:identifier.
%    Multiple identifiers separated by "|".
%  6 Aliases for B.
%  7 Interaction detection methods, taken from the corresponding PSI-MI controlled
%  Vocabulary,
%    and represented as darabaseName:identifier(methodName), separated by "|".
%  8 First author surname(s) of the publication(s) in which this interaction has been shown,
%    optionally followed by additional indicators, e.g. "Doe-2005-a". Separated by "|".
%  9 Identifier of the publication in which this interaction has been shown. Database name
%    taken from the PSI-MI controlled vocabulary, represented as databaseName:identifier.
%    Multiple identifiers separated by "|".
% 10 NCBI Taxonomy identifier for interactor A. Database name for NCBI taxid taken from the
%    PSI-MI controlled vocabulary, represented as databaseName:identifier (typicaly
%    databaseName is set to 'taxid'). Multiple identifiers separated by "|".
%    Note: In this column, the databaseName:identifier(speciesName) notation is only there
%    for consistency. Currently no taxonomy identifiers other than NCBI taxid are anticipated,
%    apart from the use of -1 to indicate "in vitro", -2 to indicate "chemical synthesis",
%    -3 indicates "unknown", -4 indicates "in vivo" and -5 indicates "in silico".
% 11 NCBI Taxonomy identifier for interactor B.
% 12 Interaction types, taken from the corresponding PSI-MI controlled vocabulary,
%    and represented as dataBaseName:identifier(interactionType), separated by "|".
% 13 Source databases and identifiers, taken from the corresponding PSI-MI controlled
%    vocabulary, and represented as databaseName:identifier(sourceName).
%    Multiple source databases can be separated by "|".
% 14 Interaction identifier(s) in the corresponding source database, represented by
%    databaseName:identifier
% 15 Confidence score. Denoted as scoreType:value. There are many different types of
%    confidence score, but so far no controlled vocabulary. Thus the only current
%    recommendation is to use score types consistently within one source.
%    Multiple scores separated by "|".

% LINE_ANALYSIS

% Interaction matrix
%  1 gene or gene pair rank in imput list
%  2 Psicquiq service rank
%  3 first gene not found
%  4 second gene (if gene pairs used) not found
%  5 main first gene name (if gene pairs used) or current interactor position in MainName list
%  6 main second gene name position in MainName list (if gene pairs used)
%  7 source code 
%  8 species of first gene (if gene pairs used) or of current interactor
%  9 species of second gene (if gene pairs used)
% 10 interaction score1
% 11 interaction score2
% 12 interaction score3
% 13 first author
% 14 pubmed position in Dict.pubmed
% 15 omim position in Dict.omim
% 16 imex position in Dict.omim
% 17 interaction type1
% 18 interaction type2
% 19 interaction type3
% 20 interaction type4
% 21 detection method type1
% 22 detection method type2


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


%psicquic('output','direct','/home/mbellis/sosma/docu/canal_calcique/psicquic',...
% {'Cav1.1';'Cav1.2';'Cav1.3';'Cav1.4';'Cav2.1';'Cav2.2';...
%    'Cav2.3';'Cav3.1';'Cav3.2';'Cav3.3';...
%     'Alpha2-D1';'Alpha2-D2';'Alpha2-D3';'Alpha2-D4';...
%     'Beta-1';'Beta-2';'Beta-3';'Beta-4';...
%     'Gamma-1';'Gamma-2';'Gamma-3';'Gamma-4';'Gamma-5';'Gamma-6';'Gamma-7';'Gamma-8'},'cav');


%psicquic('analysis','direct','/home/mbellis/sosma/docu/canal_calcique/psicquic',...
% {'CACNA1S';'CACNA1C';'CACNA1D';'CACNA1F';'CACNA1A';'CACNA1B';...
%    'CACNA1E';'CACNA1G';'CACNA1H';'CACNA1I';...
%    'CACNA2D1';'CACNA2D2';'CACNA2D3';'CACNA2D4';...
%    'CACNB1';'CACNB2';'CACNB3';'CACNB4';'CACNG1';'CACNG2';...
%    'CACNG3';'CACNG4';'CACNG5';'CACNG6';'CACNG7';'CACNG8.txt'},'cav');




%psicquic('mitab','direct','/home/mbellis/sosma/docu/canal_calcique/psicquic',...
% {'CACNA1S';'CACNA1C';'CACNA1D';'CACNA1F';'CACNA1A';'CACNA1B';...
%    'CACNA1E';'CACNA1G';'CACNA1H';'CACNA1I';...
%    'CACNA2D1';'CACNA2D2';'CACNA2D3';'CACNA2D4';...
%    'CACNB1';'CACNB2';'CACNB3';'CACNB4';'CACNG1';'CACNG2';...
%    'CACNG3';'CACNG4';'CACNG5';'CACNG6';'CACNG7';'CACNG8.txt'},'cav');


%File=fullfile('E:','sosma','docu','canal_calcique','psicquic','abul_litterature_list.txt');
%File=fullfile('home/mbellis','sosma','docu','canal_calcique','psicquic','abul_litterature_list.txt');

% psicquic('score','pair','E:\sosma\docu\canal_calcique\psicquic','abul_litterature_list.txt','abul','Abul')
% psicquic('analysis','pair','E:\sosma\docu\canal_calcique\psicquic','abul_litterature_list.txt','abul',[1,58,116,173],[57,115,172,231])
% psicquic('mitab','pair','E:\sosma\docu\canal_calcique\psicquic','abul_litterature_list.txt')
% psicquic('mitab','direct','E:\sosma\docu\canal_calcique\psicquic',{'AAK1'})
% psicquic('mitab','pair','/home/mbellis/sosma/docu/canal_calcique/psicquic','abul_litterature_list.txt','abul',1,57)
% psicquic('mitab','pair','/home/mbellis/sosma/docu/canal_calcique/psicquic','abul_litterature_list.txt','abul',58,115)
% psicquic('mitab','pair','/home/mbellis/sosma/docu/canal_calcique/psicquic','abul_litterature_list.txt','abul',116,172)
% psicquic('mitab','pair','/home/mbellis/sosma/docu/canal_calcique/psicquic','abul_litterature_list.txt','abul',173,231)


function psicquic(Action,ListType,Dir,File,Output,varargin)



switch ListType
    case 'unique'
        [Gene1]=textread(fullfile(Dir,File),'%s','delimiter','\t');
    case 'pair'
        [Gene1,Gene2]=textread(fullfile(Dir,File),'%s%s','delimiter','\t');
    case 'direct'
        Gene1=File;
        ListType='unique';
end
GeneNb=length(Gene1);

%% MITAB
switch Action

    case 'mitab'

          if nargin==7
            Start=varargin{1};
            Stop=varargin{2};
        else
            if isequal(ListType,'direct')
                Start=1;
                Stop=length(File);
            else
                Start=1;
                Stop=length(Gene1);
            end
        end

        %% GET LIST OF SERVICES


        Query='http://www.ebi.ac.uk/Tools/webservices/psicquic/registry/registry?action=STATUS&format=xml';
        Registry=xmlread(Query);
        ServiceNb=Registry.getElementsByTagName('service').getLength;
        ServiceList=Registry.getElementsByTagName('service');
        Services.name=cell(ServiceNb,1);
        Services.soapUrl=cell(ServiceNb,1);
        Services.restUrl=cell(ServiceNb,1);
        Services.restExample=cell(ServiceNb,1);
        Services.active=cell(ServiceNb,1);
        Services.count=zeros(ServiceNb,1);
        Services.version=cell(ServiceNb,1);
        Services.organizationUrl=cell(ServiceNb,1);
        Services.restricted=cell(ServiceNb,1);
        % Note that the item list index is zero-based.
        Fields={'name','soapUrl','restUrl','restExample','active','count','version','organizationUrl','restricted'};
        %Fields={'name','soapUrl','restUrl','restExample','active','version','organizationUrl','restricted'};
        for ItemL = 0:ServiceNb-1
            CurrService = ServiceList.item(ItemL);
            for FieldL=1:length(Fields)
                try
                    CurrField=Fields{FieldL};
                    Val=CurrService.getElementsByTagName(CurrField);
                    if isequal(CurrField,'count')
                        eval(sprintf('Services.%s(ItemL+1)=str2num(char(Val.item(0).getFirstChild.getData));',CurrField));
                    else
                        eval(sprintf('Services.%s{ItemL+1}=char(Val.item(0).getFirstChild.getData);',CurrField));
                    end
                catch
                    sprintf('failed at ItemL=%u and CurrField=%s',ItemL,CurrField)
                end
            end
        end  % End FOR


        %% RECOVER ADDITIONAL INFORMATION ON EACH SERVICES
        % ActivePos=strmatch('true',Services.active,'exact');
        % NotRestrictedPos=strmatch('false',Services.restricted,'exact');
        % UsedPos=intersect(ActivePos,NotRestrictedPos);
        % ServiceNb=length(Services.name);
        % Services.properties=cell(ServiceNb,1);
        % Services.formats=cell(ServiceNb,1);
        % Services.implementation=cell(ServiceNb,1);
        % for ServL=1:length(UsedPos)
        %     CurrServ=UsedPos(ServL);
        %     Query=sprintf('%sproperties',Services.restUrl{CurrServ});
        %     try
        %         Services.properties{CurrServ}=urlread(Query);
        %     catch
        %     end
        %     Query=sprintf('%sformats',Services.restUrl{CurrServ});
        %     try
        %         Services.formats{CurrServ}=urlread(Query);
        %     catch
        %     end
        %     Query=sprintf('%sversion',Services.restUrl{CurrServ});
        %     try
        %         Services.implementation{CurrServ}=urlread(Query);
        %     catch
        %     end
        % end


        %% RECOVER MITAB FORMAT

        %matrix with coded information on all interactants





        ActivePos=strmatch('true',Services.active,'exact');
        NotRestrictedPos=strmatch('false',Services.restricted,'exact');
        UsedPos=intersect(ActivePos,NotRestrictedPos);
        ServiceNb=length(Services.name);

        CurrInfo=cell(GeneNb,1);
        for GeneL=Start:Stop
            CurrGene1=Gene1{GeneL};
            if isequal(ListType,'pair')
                CurrGene2=Gene2{GeneL};
            end
            CurrInfo{GeneL}.count=zeros(ServiceNb,1);
            CurrInfo{GeneL}.txt=cell(ServiceNb,1);
            CurrInfo{GeneL}.time=zeros(ServiceNb,1);
            if GeneL==Start
                CurrInfo{GeneL}.serviceName=Services.name;
                CurrInfo{GeneL}.serviceActive=Services.active;
                CurrInfo{GeneL}.serviceUsed=UsedPos;
            end
            for ServL=1:length(UsedPos)
                tic
                [GeneL ServL]
                CurrServ=UsedPos(ServL);
                AskIt=0;
                if isequal(ListType,'pair')
                    Query=sprintf('%squery/alias:%s?format=count',Services.restUrl{CurrServ},CurrGene1);
                    Count1=str2num(READ_URL(Query));
                    Query=sprintf('%squery/alias:%s?format=count',Services.restUrl{CurrServ},CurrGene2);
                    Count2=str2num(READ_URL(Query));
                    if Count1>0 & Count2>0
                        AskIt=1;
                    end
                else
                    Query=sprintf('%squery/alias:%s?format=count',Services.restUrl{CurrServ},CurrGene1);
                    Count1=str2num(READ_URL(Query));
                    CurrInfo{GeneL}.count(CurrServ)=Count1;
                    if Count1>0
                        AskIt=1;
                    end
                end

                if AskIt
                    if isequal(ListType,'pair')
                        Query=sprintf('%squery/alias:(%s%%20AND%%20%s)?format=count',Services.restUrl{CurrServ},CurrGene1,CurrGene2);
                        CurrInfo{GeneL}.count(CurrServ)=str2num(READ_URL(Query));
                        Query=sprintf('%squery/alias:(%s%%20AND%%20%s)',Services.restUrl{CurrServ},CurrGene1,CurrGene2);
                    else
                        Query=sprintf('%squery/alias:(%s)',Services.restUrl{CurrServ},CurrGene1);
                    end
                    if CurrInfo{GeneL}.count(CurrServ)>0
                        CurrInfo{GeneL}.txt{CurrServ}=READ_URL(Query);
                    end
                end
                CurrInfo{GeneL}.time(CurrServ)=toc;
            end
        end

        cd(Dir)
        eval(sprintf('save %s-%u-%u CurrInfo',Output,Start,Stop))
        
%% ANALYSIS
    case 'analysis'
         if nargin==7
            Starts=varargin{1};
            Stops=varargin{2};
        else
            if isequal(ListType,'direct')
                Starts=1;
                Stops=length(File);
            else
                Starts=1;
                Stops=length(Gene1);
            end
        end
        %merge result files
        cd(Dir)
        if exist(sprintf('%s.mat',Output),'file')
            eval(sprintf('load %s',Output))
        else
            
             for FileL=1:length(Starts)
                eval(sprintf('load %s-%u-%u',Output,Starts(FileL),Stops(FileL)))
                if FileL==1
                    Info=CurrInfo;
                else
                    for GeneL=Starts(FileL):Stops(FileL)
                        Info{GeneL}.count=CurrInfo{GeneL}.count;
                        Info{GeneL}.txt=CurrInfo{GeneL}.txt;
                        Info{GeneL}.time=CurrInfo{GeneL}.time;
                        if GeneL==Starts(FileL)
                            Info{GeneL}.serviceName=CurrInfo{GeneL}.serviceName;
                            Info{GeneL}.serviceActive=CurrInfo{GeneL}.serviceActive;
                            Info{GeneL}.serviceUsed=CurrInfo{GeneL}.serviceUsed;
                        end
                    end
                end
             end

             if length(Starts)==1
                 ServicesName{1}=Info{Starts(1)}.serviceName;
                 ServicesActive{1}=Info{Starts(1)}.serviceActive;
                 ServicesUsed{1}=Info{Starts(1)}.serviceUsed;
             else
                 ServicesName=cell(length(Starts),1);
                 ServicesActive=cell(length(Starts),1);
                 ServicesUsed=cell(length(Starts),1);
                 %construct common reference services
                 DiffNameFlag=0;
                 DiffUsedFlag=0;
                 for FileL1=1:length(Starts)-1
                     for FileL2=FileL1+1:length(Starts)
                         if ~isempty(setdiff(Info{Starts(FileL1)}.serviceName,Info{Starts(FileL2)}.serviceName))|...
                                 ~isempty(setdiff(Info{Starts(FileL2)}.serviceName,Info{Starts(FileL1)}.serviceName))
                             DiffNameFlag=1;
                         end
                         if ~isempty(setdiff(Info{Starts(FileL1)}.serviceUsed,Info{Starts(FileL2)}.serviceUsed))|...
                                 ~isempty(setdiff(Info{Starts(FileL2)}.serviceUsed,Info{Starts(FileL1)}.serviceUsed))
                             DiffUsedFlag=1;
                         end
                     end
                 end
                 if DiffNameFlag==0 & DiffUsedFlag==1
                     %reorder fields
                     for FileL=1:length(Starts)
                         [temp SortOrder]=sort(Info{Starts(FileL)}.serviceName);
                         Info{Starts(FileL)}.serviceName=Info{Starts(FileL)}.serviceName(SortOrder);
                         Used=zeros(length(Info{Starts(FileL)}.serviceName),1);
                         Used(Info{Starts(FileL)}.serviceUsed)=1;
                         Info{Starts(FileL)}.serviceUsed=Used(SortOrder);
                         Info{Starts(FileL)}.serviceActive=Info{Starts(FileL)}.serviceActive(SortOrder);
                         for GeneL=Starts(FileL):Stops(FileL)
                             Info{GeneL}.count=Info{GeneL}.count(SortOrder);
                             Info{GeneL}.txt=Info{GeneL}.txt(SortOrder);
                             Info{GeneL}.time=Info{GeneL}.time(SortOrder);
                         end
                     end
                 end
             end
            eval(sprintf('save %s Info',Output))
        end

        %count field number
        ServiceNb=length(Info{1}.count);
        FieldNb=zeros(ServiceNb,1);
        for GeneL=1:GeneNb
            ServPos=find(Info{GeneL}.count>0);
            for ServL=1:length(ServPos)
                CurrServ=ServPos(ServL);
                LineEnd=regexp(Info{GeneL}.txt{CurrServ},'\n');
                CurrFieldNb=length(find(regexp(Info{GeneL}.txt{CurrServ},'\t')<LineEnd(1)))+1;
                if CurrFieldNb>FieldNb(CurrServ)
                    if FieldNb(CurrServ)>0
                        [GeneL CurrServ FieldNb(CurrServ)]
                    end                        
                    FieldNb(CurrServ)=CurrFieldNb;
                end
            end
        end
        
        %count number of results
        ResNb=zeros(ServiceNb,1);
        NullResNb{1}=zeros(ServiceNb,1);
        NullResNb{2}=zeros(ServiceNb,1);
        NullResNb{3}=zeros(ServiceNb,1);
        for GeneL=1:GeneNb
            ServPos=find(Info{GeneL}.count>0);
            ResNb(ServPos)=ResNb(ServPos)+Info{GeneL}.count(ServPos);
            for ResL=1:3
            ServPos=find(Info{GeneL}.count==-ResL);
            NullResNb{ResL}(ServPos)=NullResNb{ResL}(ServPos)+1;
            end
        end
        
        
        %analysis of txt content
        cd(Dir)
        LogFile=fopen(sprintf('import_log_%s.txt',date),'w');
        fprintf(LogFile,'Service\tLineRank\tGeneNameA\tGeneNameB\tAliasA\tAliasB\n');

        %alias id for interrogating genes
        Gene.ids={};
        Gene.idTypes={};
        %dictionaries of identificator
        Dict.uniprotkbId={};
        Dict.uniprotId={};
        Dict.complexId={};
        Dict.rogId={};
        Dict.irogId={};
        Dict.crogId={};
        Dict.icrogId={};
        Dict.stringId={};
        Dict.refseqId={};
        Dict.ensemblId={};
        Dict.ensemblGeneId={};        
        Dict.entrezId={};
        Dict.locuslinkId={};
        Dict.hgncId={};
        Dict.bindId={};
        Dict.omimId={};
        Dict.pdbId={};
        Dict.intactId={};
        Dict.irefindexId={};
        Dict.respdbId={};
        Dict.unknownId={};
        %dictionary for InteractionType
        Dict.interaction={};
        Dict.interactionCode=[];
        %dictionary for Species
        Dict.species={};
        Dict.speciesCode=[];
        %dictionary for DetectionMethod
        Dict.method={};
        Dict.methodCode=[];
        %dictionary for Source
        Dict.source={};
        Dict.sourceCode=[];
        %dictionary for PUBMED
        Dict.pubmedCode=[];
        %dictionary for Authors
        Dict.author={};
        %PUBMED matrix
        Pubmed=[];        
        %dictionary for OMIM
        Dict.omimCode=[];        
        %OMIM matrix
        Omim=[];
        %dictionary for IMEX
        Dict.imexCode=[];
        %IMEX matrix
        Imex=[];
        %ID matrix
        if isequal(ListType,'pair')
            ID=cell(2,1);        
            ID{1}=[];
        else
            ID=cell(1);
            ID{2}=[];
        end
        %Interac,tion matrix
        Interaction=[];
        MainAlias={};
        for GeneL=1:GeneNb
            Gene.ids{GeneL}={};
            Gene.idTypes{GeneL}=[];
            Pos=find(Info{GeneL}.count>0);
            for ServL=1:length(Pos)
                CurrServ=Pos(ServL);
                ServName=Info{1}.serviceName{CurrServ};
                %recover information              
                LinePos=[0,regexp(Info{GeneL}.txt{CurrServ},'\n')];
                if LinePos(end)<length(Info{GeneL}.txt{CurrServ})
                    LinePos(end+1)=length(Info{GeneL}.txt{CurrServ})+1;
                end
                LineNb=length(LinePos)-1;
                for LineL=1:LineNb
                    CurrLine=Info{GeneL}.txt{CurrServ}(LinePos(LineL)+1:LinePos(LineL+1)-1);
                    TabPos=regexp(CurrLine,'\t');
                    if length(TabPos)>14
                        %in Intact and UniProt, there exist an extended format with 31 fields
                        %the first 15 fields are the standard PSIMAT fields
                        CurrLine=CurrLine(1:TabPos(15)-1);
                    end
                    if isequal(ListType,'pair')
                        [Dict,Gene,Interaction,Pubmed,Omim,Imex,ID,MainAlias]=LINE_ANALYSIS(GeneL,ListType,CurrServ,ServName,CurrLine,Gene1{GeneL},LogFile,Dict,Gene,Interaction,Pubmed,Omim,Imex,ID,MainAlias,Gene2{GeneL});
                    else
                        [Dict,Gene,Interaction,Pubmed,Omim,Imex,ID,MainAlias]=LINE_ANALYSIS(GeneL,ListType,CurrServ,ServName,CurrLine,Gene1{GeneL},LogFile,Dict,Gene,Interaction,Pubmed,Omim,Imex,ID,MainAlias);
                    end
                end
            end
        end
        fclose(LogFile)
        cd(Dir)
%         IdType={};
%         for i=1:length(Dict.idType)
%             Val=regexp(Dict.idType{i},'.*(?=:)','match');
%             if ~isempty(Val)
%                 if isempty(strmatch(Val{1},IdType,'exact'))
%                     IdType{end+1,1}=Val{1};
%                 end
%             else
%                 if isempty(strmatch(Dict.idType{i},IdType,'exact'))
%                     IdType{end+1,1}=Dict.idType{i};
%                 end
%             end
%         end

        
        eval(sprintf('save %s Info Interaction Dict MainAlias Imex Omim Pubmed ID Gene',Output))

%% SCORE
    case 'score'
        cd(Dir)
        if exist(sprintf('%s.mat',Output),'file')
            eval(sprintf('load %s',Output))
        else
            h=errordlg('run analysis first');
            waitfor(h)
            error('process canceled')
        end
        'stop'
        if nargin==6
            Author=varargin{1};
            AuthorPos=strmatch(Author,Dict.author);
            %eliminate interaction deposited by this author
            if ~isempty(AuthorPos)
                if length(AuthorPos)==1
                    InterPos=find(Interaction(:,13)==AuthorPos);
                    if ~isempty(InterPos)
                        Interaction(InterPos,:)=[];
                        h=warndlg(sprintf('%u line(s) corresponding to %s author deleted from Interaction',length(InterPos),Dict.author{AuthorPos}));
                        waitfor(h)
                    end
                else
                    h=errordlg('exist %u authors corresponding to %s',length(AuthorPos),Author);
                    waitfor(h)
                    error('process canceled')
                    %to be developped
                end
            end
        end
        
        %number of gene found
        GeneNb=length(Gene1);
        FoundGene=unique(Interaction(:,1));
        FoundGeneNb=length(FoundGene);
        
        %score distribution
        h=figure;
        set(gcf,'color',[1,1,1])
        set(h,'name','SCORE DISTRIBUTION')
        UsedServ=unique(Interaction(:,2));
        ServNb=length(UsedServ);
        if ~isempty(strmatch('iRefIndex',Info{1}.serviceName(UsedServ),'exact'))
            PlotNb=ServNb+2;
        end
        RowNb=ceil(sqrt(PlotNb));
        ColNb=floor(PlotNb/RowNb);
        if RowNb*ColNb<PlotNb
            ColNb=ColNb+1;
        end
        Offset=0;

        for ServL=1:ServNb
            subplot(RowNb,ColNb,ServL+Offset)
            CurrServ=UsedServ(ServL);
            ServPos=find(Interaction(:,2)==CurrServ);
            Score1=Interaction(ServPos,10);
            if isequal(Info{1}.serviceName{CurrServ},'iRefIndex')
                Score2=Interaction(ServPos,11);
                Score3=Interaction(ServPos,12);
                Offset=2;
            end
            if length(find(Score1>0))>0
                hist(Score1(find(Score1>0)),round(length(find(Score1>0))/10))
            end
            title(sprintf('%s: %u interactions for %u pairs',...
                Info{1}.serviceName{CurrServ},...
                length(ServPos),...
                length(unique(Interaction(ServPos,5:6),'rows'))))          
            xlabel(sprintf('%u with score>0',length(find(Score1>0))))
            set(gca,'box','on')
            %[length(unique(Score1))*100/length(find(Score1>0)), length(Score1),length(find(Score1>0)),min(Score1(find(Score1>0))),max(Score1)]
            if isequal(Info{1}.serviceName{CurrServ},'iRefIndex')
                subplot(RowNb,ColNb,ServL+1)
                if length(find(Score1>0))>0
                    hist(Score2(find(Score2>0)),round(length(find(Score2>0))/10))
                end
                title(sprintf('%s: %u interactions for %u pairs',...
                    Info{1}.serviceName{CurrServ},...
                    length(ServPos),...
                    length(unique(Interaction(ServPos,5:6),'rows'))))
                xlabel(sprintf('%u with score>0',length(find(Score2>0))))
                set(gca,'box','on')
                %[length(unique(Score2))*100/length(find(Score2>0)), length(Score2),length(find(Score2>0)),min(Score2(find(Score2>0))),max(Score2)]
                subplot(RowNb,ColNb,ServL+2)
                if length(find(Score1>0))>0
                    hist(Score3(find(Score3>0)),round(length(find(Score3>0))/10))
                end
                title(sprintf('%s: %u interactions for %u pairs',...
                    Info{1}.serviceName{CurrServ},...
                    length(ServPos),...
                    length(unique(Interaction(ServPos,5:6),'rows'))))                
                xlabel(sprintf('%u with score>0',length(find(Score3>0))))
                set(gca,'box','on')
                %[length(unique(Score3))*100/length(find(Score3>0)), length(Score3),length(find(Score3>0)),min(Score3(find(Score3>0))),max(Score3)]
            end
        end
        
        %methods with high and low scores
        %recover inferior and superior limits
        UsedServ=unique(Interaction(:,2));
        ServNb=length(UsedServ);
        MinMaxScore=[];
        for ServL=1:ServNb
            CurrServ=UsedServ(ServL);
            ServPos=find(Interaction(:,2)==CurrServ);
            if isequal(Info{1}.serviceName{CurrServ},'iRefIndex')
                ScoreNb=3;
            else
                ScoreNb=1;
            end
            for ScoreL=1:ScoreNb
                Score=Interaction(ServPos,9+ScoreL);
                Score=Score(find(Score>0));
                if length(Score>50)
                    Score=sort(Score);
                    MinScore=Score(round(length(Score)*0.25));
                    MaxScore=Score(round(length(Score)*0.75));
                    MinMaxScore(end+1,:)=[CurrServ,ScoreL,MinScore,round(length(find(Score<=MinScore))*100/length(Score)),...
                        MaxScore,round(length(find(Score>=MaxScore))*100/length(Score))];
                end
            end
        end
        
        %methods frequency        
        MethodNb=zeros(length(Dict.method),size(MinMaxScore,1));
        MethodMin=MethodNb;
        MethodMax=MethodNb;
        for ServL=1:size(MinMaxScore,1)
            CurrServ=MinMaxScore(ServL,1);
            ServPos=find(Interaction(:,2)==CurrServ);
            Score=Interaction(ServPos,9+MinMaxScore(ServL,2));
            ScorePos=find(Score>0);
            Score=Score(ScorePos);
            ServPos=ServPos(ScorePos);
            for MethodL=1:length(Dict.method)
%                 MethodNb(MethodL,ServL)=length(find(Interaction(ServPos,21)==Dict.methodCode(MethodL)));
%                 if MethodNb(MethodL,ServL)>0
%                     MinPos=find(Score<=MinMaxScore(ServL,3));
%                     MethodMin(MethodL,ServL)=round(length(find(Interaction(ServPos(MinPos),21)==Dict.methodCode(MethodL)))*100/MethodNb(MethodL,ServL));
%                     MaxPos=find(Score>=MinMaxScore(ServL,5));
%                     MethodMax(MethodL,ServL)=round(length(find(Interaction(ServPos(MaxPos),21)==Dict.methodCode(MethodL)))*100/MethodNb(MethodL,ServL));
%                 end
                MethodNb(MethodL,ServL)=length(find(Interaction(ServPos,21)==MethodL));
                if MethodNb(MethodL,ServL)>0
                    MinPos=find(Score<=MinMaxScore(ServL,3));
                    MethodMin(MethodL,ServL)=round(length(find(Interaction(ServPos(MinPos),21)==MethodL))*100/MethodNb(MethodL,ServL));
                    MaxPos=find(Score>=MinMaxScore(ServL,5));
                    MethodMax(MethodL,ServL)=round(length(find(Interaction(ServPos(MaxPos),21)==MethodL))*100/MethodNb(MethodL,ServL));
                end                
            end            
        end  
        %eliminate not used methods
        ClearPos=find(sum(MethodNb,2)==0);
        KeepPos=find(sum(MethodNb,2)>0);
        MethodNb(ClearPos,:)=[];
        MethodMin(ClearPos,:)=[];
        MethodMax(ClearPos,:)=[];
        Dict.method(KeepPos)
        
        %interaction frequency
        InteractionNb=zeros(length(Dict.interaction),size(MinMaxScore,1));
        InteractionMin=InteractionNb;
        InteractionMax=InteractionNb;
        for ServL=1:size(MinMaxScore,1)
            CurrServ=MinMaxScore(ServL,1);
            ServPos=find(Interaction(:,2)==CurrServ);
            Score=Interaction(ServPos,9+MinMaxScore(ServL,2));
            ScorePos=find(Score>0);
            Score=Score(ScorePos);
            ServPos=ServPos(ScorePos);
            for InteractionL=1:length(Dict.interaction)
                InteractionNb(InteractionL,ServL)=length(find(Interaction(ServPos,17)==InteractionL));
                if InteractionNb(InteractionL,ServL)>0
                MinPos=find(Score<=MinMaxScore(ServL,3));
                InteractionMin(InteractionL,ServL)=round(length(find(Interaction(ServPos(MinPos),17)==InteractionL))*100/InteractionNb(InteractionL,ServL));
                MaxPos=find(Score>=MinMaxScore(ServL,5));
                InteractionMax(InteractionL,ServL)=round(length(find(Interaction(ServPos(MaxPos),17)==InteractionL))*100/InteractionNb(InteractionL,ServL));                
                end
            end            
        end  
        %eliminate not used methods
        ClearPos=find(sum(InteractionNb,2)==0);
        KeepPos=find(sum(InteractionNb,2)>0);
        InteractionNb(ClearPos,:)=[];
        InteractionMin(ClearPos,:)=[];
        InteractionMax(ClearPos,:)=[];
        Dict.interaction(KeepPos)
        
        %% OUTPUT
    case 'output'
        cd(Dir)
        if exist(sprintf('%s.mat',Output),'file')
            eval(sprintf('load %s',Output))
        else
            h=errordlg('run analysis first');
            waitfor(h)
            error('process canceled')
        end
        Muller={};
        [Muller{1}.id,Muller{1}.name,Muller{1}.abondance,Muller{1}.cav21,Muller{1}.cav22,Muller{1}.cav23]=textread('muller1.txt','%s%s%s%s%s%s','delimiter','\t');
        [Muller{2}.id,Muller{2}.name,Muller{2}.abondance,Muller{2}.cav21,Muller{2}.cav22,Muller{2}.cav23]=textread('muller2.txt','%s%s%s%s%s%s','delimiter','\t');

        %calculate repetitivity
        IdTypeNb=size(ID{1},2);
        for InterL=1:size(Interaction,1)
            InterPos=find(Interaction(:,1)==Interaction(InterL,1));
            IdemPos=[];
            for IdL=1:IdTypeNb
                if ID{1}(InterL,IdL)>0
                    IdemPos=union(IdemPos,find(ID{1}(InterPos,IdL)==ID{1}(InterL,IdL)));
                end
            end
            Interaction(InterL,4)=length(IdemPos);
        end

        IdDict={'uniprotkbId';...
            'uniprotId';...
            'entrezId';...
            'locuslinkId';...
            'stringId';...
            'refseqId';...
            'ensemblId';...
            'ensemblGeneId';...
            'hgncId';...
            'pdbId';...
            'respdbId';...
            'intactId';...
            'irefindexId';...
            'bindId';...
            'complexId';...
            'rogId';...
            'irogId';...
            'crogId';...
            'icrogId';...
            'omimId';...
            'unknownId'};

        %don't process genemania
        GeneMania=strmatch('genemania',lower(Info{1}.serviceName),'exact');

        % %write data
        % fid=fopen(sprintf('%s_info.txt',Output),'w');
        % fprintf(fid,[repmat('\t',1,8),'MULLER\t\t\t\tMULLERnonb\t\t\t\','\tINTERACTOR\n']);
        % fprintf(fid,'RANK\tSUBUNIT\tEMPTY\tREPEAT\tSOURCE1\tSOURCE2\tREFBIB\tSPECIES\tAbondance\tCav2.1\tCav2.2\tCav2.3\tAbondance\tCav2.1\tCav2.2\tCav2.3\tUNIPROTKB\tUNIPROTKB\tENTREZ\tENTREZ\tSTRING\tREFSEQ\tENSEMBL\tENSEMBL\tHGNC\tPDB\tPDB\tINTACT\tIREFINDEX\tBIND\tCOMPLEX\tROGID\tIROGID\tCROGID\tICROGID\tOMIM\tUNKNOWN\tINTERACTION TYPE\t\t\tDETECTION METHOD\n');
        % Format=[repmat('\t%u',1,size(Interaction,2)-16),'\n'];
        % FoundMuller{1}=[];
        % FoundMuller{2}=[];
        % for InterL=1:size(Interaction,1)
        %     if Interaction(InterL,2)~=GeneMania
        %         %find Muller intersection
        %         if ID{1}(InterL,1)>0
        %             MulPos=[0,0];
        %             for MulL=1:2
        %                 CurrMulPos=strmatch(Dict.uniprotkbId{ID{1}(InterL,1)},Muller{MulL}.id,'exact');
        %                 if ~isempty(CurrMulPos)
        %                     MulPos(MulL)=CurrMulPos;
        %                     if isempty(find(FoundMuller{MulL}==MulPos(MulL)))
        %                         FoundMuller{MulL}(end+1,1)=MulPos(MulL);
        %                     end
        %                 end
        %             end
        %         end
        %         PubmedPos=find(Pubmed(:,InterL));
        %         CurrPubmed='';
        %         if ~isempty(PubmedPos)
        %             if length(PubmedPos)==1
        %                 CurrPubmed=num2str(PubmedPos);
        %             else
        %                 for PubL=1:length(PubmedPos)-1
        %                     CurrPubmed=[CurrPubmed,num2str(PubmedPos(PubL)),','];
        %                 end
        %                 CurrPubmed=[CurrPubmed,num2str(PubmedPos(end))];
        %             end
        %         end
        %         fprintf(fid,'%u\t%s\t%u\t%u\t%u\t%u\t%s\t%u',InterL,Gene1{Interaction(InterL,1)},Interaction(InterL,[3,4,2,7]),CurrPubmed,Interaction(InterL,8));
        %         for MulL=1:2
        %             if MulPos(MulL)>0
        %                 fprintf(fid,'\t%s\t%s\t%s\t%s',Muller{MulL}.abondance{MulPos(MulL)},...
        %                     Muller{MulL}.cav21{MulPos(MulL)},Muller{MulL}.cav22{MulPos(MulL)},...
        %                     Muller{MulL}.cav23{MulPos(MulL)});
        %             else
        %                 fprintf(fid,'\t\t\t\t');
        %             end
        %         end
        %         for IdL=1:length(IdDict)
        %             CurrId=ID{1}(InterL,IdL);
        %             if CurrId==0
        %                 fprintf(fid,'\t');
        %             else
        %                 eval(sprintf('Id=Dict.%s{CurrId};',IdDict{IdL}))
        %                 fprintf(fid,'\t%s',Id);
        %             end
        %         end
        %         fprintf(fid,Format,Interaction(InterL,17:end));
        %     end
        % end
        % fclose(fid)


        %write data without repetition
        fid=fopen(sprintf('%s_info.txt',Output),'w');
        fprintf(fid,[repmat('\t',1,8),'MULLER\t\t\t\tMULLERnonb\t\t\t\','\tINTERACTOR\n']);
        fprintf(fid,'RANK\tSUBUNIT\tEMPTY\tREPEAT\tSOURCE1\tSOURCE2\tREFBIB\tSPECIES\tAbondance\tCav2.1\tCav2.2\tCav2.3\tAbondance\tCav2.1\tCav2.2\tCav2.3\tUNIPROTKB\tUNIPROTKB\tENTREZ\tENTREZ\tSTRING\tREFSEQ\tENSEMBL\tENSEMBL\tHGNC\tPDB\tPDB\tINTACT\tIREFINDEX\tBIND\tCOMPLEX\tROGID\tIROGID\tCROGID\tICROGID\tOMIM\tUNKNOWN\tINTERACTION TYPE\tDETECTION METHOD\tSCORE1\tSCORE2\tSCORE3\n');
        FoundMuller{1}=[];
        FoundMuller{2}=[];
        Processed=zeros(size(Interaction,1),1);
        for InterL=1:size(Interaction,1)
            if Processed(InterL)==0
                RepeatNb=Interaction(InterL,4);
                if RepeatNb>1
                    InterPos=find(Interaction(:,1)==Interaction(InterL,1));
                    IdemPos=[];
                    for IdL=1:IdTypeNb
                        if ID{1}(InterL,IdL)>0
                            IdemPos=union(IdemPos,find(ID{1}(InterPos,IdL)==ID{1}(InterL,IdL)));
                        end
                    end
                else
                    IdemPos=InterL;
                end
                Processed(IdemPos)=1;

                if Interaction(InterL,2)~=GeneMania
                    %find Muller intersection
                    if ID{1}(InterL,1)>0
                        MulPos=[0,0];
                        for MulL=1:2
                            CurrMulPos=strmatch(Dict.uniprotkbId{ID{1}(InterL,1)},Muller{MulL}.id,'exact');
                            if ~isempty(CurrMulPos)
                                MulPos(MulL)=CurrMulPos;
                                if isempty(find(FoundMuller{MulL}==MulPos(MulL)))
                                    FoundMuller{MulL}(end+1,1)=MulPos(MulL);
                                end
                            end
                        end
                    end
                    %process pubmed references
                    PubmedPos=find(Pubmed(:,InterL));
                    CurrPubmed='';
                    if ~isempty(PubmedPos)
                        CurrPubmed='/';
                        for PubL=1:length(PubmedPos)
                            CurrPubmed=[CurrPubmed,num2str(PubmedPos(PubL)),'/'];
                        end
                    end
                    %process source 2
                    SourceCode=[];
                    for IdemL=1:length(IdemPos)
                        SourceCode=union(SourceCode,Interaction(IdemPos(IdemL),7));
                    end
                    SourceCode=setdiff(SourceCode,0);
                    Sources='/';
                    for i=1:length(SourceCode)
                        Sources=[Sources,num2str(SourceCode(i)),'/'];
                    end
                    %process interaction types and method types
                    InterCode=[];
                    MethodCode=[];
                    for IdemL=1:length(IdemPos)
                        for ColL=17:19
                            InterCode=union(InterCode,Interaction(IdemPos(IdemL),ColL));
                        end
                        for ColL=20:size(Interaction,2)
                            MethodCode=union(MethodCode,Interaction(IdemPos(IdemL),ColL));
                        end
                    end
                    InterCode=setdiff(InterCode,0);
                    MethodCode=setdiff(MethodCode,0);
                    if ~isempty(InterCode)
                        Interactions='/';
                        for i=1:length(InterCode)
                            Interactions=[Interactions,num2str(InterCode(i)),'/'];
                        end
                    else
                        Interactions='';
                    end
                    if ~isempty(MethodCode)
                        Methods='/';
                        for i=1:length(MethodCode)
                            Methods=[Methods,num2str(MethodCode(i)),'/'];
                        end
                    else
                        Methods='';
                    end
                    %process scores
                    for ScoreL=1:3
                        ScoreCode=[];
                        for IdemL=1:length(IdemPos)
                            ScoreCode=union(ScoreCode,Interaction(IdemPos(IdemL),9+ScoreL));
                        end
                        ScoreCode=setdiff(ScoreCode,-1);
                        if ~isempty(ScoreCode)
                            Scores{ScoreL}='/';
                            for i=1:length(ScoreCode)
                                Scores{ScoreL}=[Scores{ScoreL},num2str(ScoreCode(i)),'/'];
                            end
                        else
                            Scores{ScoreL}='';
                        end
                    end

                    fprintf(fid,'%u\t%s\t%u\t%u\t%u\t%s\t%s\t%u',InterL,Gene1{Interaction(InterL,1)},Interaction(InterL,[3,4,2]),Sources,CurrPubmed,Interaction(InterL,8));
                    for MulL=1:2
                        if MulPos(MulL)>0
                            fprintf(fid,'\t%s\t%s\t%s\t%s',Muller{MulL}.abondance{MulPos(MulL)},...
                                Muller{MulL}.cav21{MulPos(MulL)},Muller{MulL}.cav22{MulPos(MulL)},...
                                Muller{MulL}.cav23{MulPos(MulL)});
                        else
                            fprintf(fid,'\t\t\t\t');
                        end
                    end
                    for IdL=1:length(IdDict)
                        CurrId=ID{1}(InterL,IdL);
                        if CurrId==0
                            fprintf(fid,'\t');
                        else
                            eval(sprintf('Id=Dict.%s{CurrId};',IdDict{IdL}))
                            fprintf(fid,'\t%s',Id);
                        end
                    end
                    fprintf(fid,'\t%s\t%s\t%s\t%s\t%s\n',Interactions,Methods,Scores{1},Scores{2},Scores{3});
                end
            end
        end
        fclose(fid)

        for MulL=1:2
            NotFoundMuller{MulL}=1:length(Muller{MulL}.id);
            NotFoundMuller{MulL}=setdiff(NotFoundMuller{MulL},FoundMuller{MulL});
            if ~isempty(NotFoundMuller{MulL})
                fid=fopen(sprintf('nofoundmuller_%u.txt',MulL),'w');
                for IdL=1:length(NotFoundMuller{MulL})
                    Pos=NotFoundMuller{MulL}(IdL);
                    fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\n',Muller{MulL}.id{Pos},Muller{MulL}.name{Pos},Muller{MulL}.abondance{Pos},Muller{MulL}.cav21{Pos},Muller{MulL}.cav22{Pos},Muller{MulL}.cav23{Pos});
                end
                fclose(fid)
            end
        end
        'stop'

%% PUBMED
    case 'pubmed'
url='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=14529713'
% end
      cd(Dir)
        if exist(sprintf('%s.mat',Output),'file')
            eval(sprintf('load %s',Output))
        else
            h=errordlg('run analysis first');
            waitfor(h)
            error('process canceled')
        end

        fid=fopen(sprintf('%s_pubmed1.txt',Output),'w');
        fprintf(fid,'pubrank\tsubunits\tauthors\ttitle\tabstract\tjournal\tyear\tvolume\tissue\tpages\tpubid\n')
        for PubL=1:length(Dict.pubmedCode)
            if Dict.pubmedCode(PubL)>0
                Url=sprintf('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=%u',Dict.pubmedCode(PubL));
                Ref=urlread(Url);                
                
                GenePos=unique(Interaction(find(Interaction(:,14)==PubL),1));                
                Authors=regexp(Ref,'(?<=names ml {\n)[^}]*(?=})','match');
                if isempty(Authors)
                    Authors=regexp(Ref,'(?<=name ml ")[^"]*(?=")','match');
                end
                if ~isempty(Authors)
                    Authors=regexp(Authors{1},'(?<=")[^",]*(?=")','match');
                else
                    Authors{1}='-';
                end
                Title=regexp(Ref,'(?<=title {\n *name ")[^}]*(?=")','match');
                Abstract=regexp(Ref,'(?<=abstract ")[^"]*(?=")','match');
                Journal=regexp(Ref,'(?<=from journal [^"]*")[^"]*(?=")','match');
                Year=regexp(Ref,'(?<= year )\d{4,4}','match');
                Volume=regexp(Ref,'(?<=volume ")[^"]*(?=")','match');
                Issue=regexp(Ref,'(?<=issue ")[^"]*(?=")','match');
                Pages=regexp(Ref,'(?<=pages ")[^"]*(?=")','match');
                
                fprintf(fid,'%u\t',PubL);                
                if ~isempty(GenePos)
                    fprintf(fid,Gene1{GenePos(1)});
                    if length(GenePos)>1
                        for i=2:length(GenePos)
                            fprintf(fid,', %s',Gene1{GenePos(i)});
                        end
                    end
                end
                fprintf(fid,'\t');
                if ~isempty(Authors)
                    fprintf(fid,Authors{1});
                    if length(Authors)>1
                        for i=2:length(Authors)
                            fprintf(fid,', %s',Authors{i});
                        end
                    end
                end
                fprintf(fid,'\t');
                if ~isempty(Title)
                    fprintf(fid,'%s\t',regexprep(Title{1},'\n',''));
                end
                if ~isempty(Abstract)
                    fprintf(fid,'%s\t',regexprep(Abstract{1},'\n',''));
                end
                if ~isempty(Journal)
                    fprintf(fid,'%s\t',Journal{1});
                end
                if ~isempty(Year)
                    fprintf(fid,'%s\t',Year{1});
                end
                if ~isempty(Volume)
                    fprintf(fid,'%s\t',Volume{1});
                end
                if ~isempty(Issue)
                    fprintf(fid,'%s\t',Issue{1});
                end
                if ~isempty(Pages)
                    fprintf(fid,'%s\t',Pages{1});
                end
                fprintf(fid,'%u\n',Dict.pubmedCode(PubL));
            end
        end
        fclose(fid)


%% MAKE MATRIX AND FIGURE OF INTERACTANTS
    case 'matrix'
% PsiMat=zeros(length(AllAlias),length(FileList));
% for FileL=1:length(FileList)
%     CurrAliases=unique(MainAlias{FileL});
%     for AliasL=1:length(CurrAliases)
%         PsiMat(strmatch(CurrAliases{AliasL},AllAlias,'exact'),FileL)=1;
%     end
% end
% %codage according to repetition
% for i=1:length(PsiMat)
%     PsiMat(i,find(PsiMat(i,:)))=sum(PsiMat(i,:));
% end
% %reorder data starting from last subunit
% for FileL=length(FileList):-1:1
%     [Temp SortOrder]=sort(PsiMat(:,FileL));
%     SortOrder=flipud(SortOrder);
%     PsiMat=PsiMat(SortOrder,:);
%     AllAlias=AllAlias(SortOrder);
% end
% 
% SubUnits=FileList;
% for FileL=1:length(FileList)
%     SubUnits{FileL}=FileList{FileL}(1:length(FileList{FileL})-4);
% end
% 
% figure
% set(gcf,'color',[1,1,1])
% image(PsiMat'*floor(255/size(PsiMat,2)))
% set(gca,'ytick',[1:length(SubUnits)]+0.5)
% set(gca,'yticklabel',SubUnits)
% title('INTERACTANTS OF CALCIUM CHANEL SUBUNITS')
% map=colormap;
% map(1,:)=[1,1,1];
% colormap(map)
% set(gca,'xgrid','on')
% 

end


%% LINE_ANALYSIS
function [Dict,Gene,Interaction,Pubmed,Omim,Imex,ID,MainAlias]=LINE_ANALYSIS(GeneL,ListType,Service,ServName,Line,FirstGeneName,LogFile,Dict,Gene,Interaction,Pubmed,Omim,Imex,ID,MainAlias,varargin)



[NameA,NameB,AltIdA,AltIdB,AliasesA,AliasesB,DetectionMethod,FirstAuthor,Biblio,SpeciesA,SpeciesB,InteractionType,Source,InteractionId,ConfidenceScore]=strread(Line,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter','\t');

GeneName{1}=FirstGeneName;
if isequal(ListType,'pair')
    GeneName{2}=varargin{1};
    AliasNb=2;
else
    AliasNb=1;
end

%recover information
InfoPos=size(Interaction,1)+1;

%process all Gene Identificators
IdTypes={'uniprotkb';...
    'uniprotkb';...
    'entrez gene/locuslink';...
    'entrezgene/locuslink';...    
    'string';...
    'refseq';... 
    'ensembl';...
    'ensembl';...    
    'hgnc';...
    'pdb';...
    'rcsb pdb';...
    'intact';...
    'irefindex';...
    'bind complexid';...    
    'complex';...
    'rogid';...
    'irogid';...
    'crogid';...
    'icrogid';...    
    'omim';...
    'unknown'};

IdDict={'uniprotkbId';...
    'uniprotId';...
    'entrezId';...    
    'locuslinkId';...        
    'stringId';...    
    'refseqId';...
    'ensemblId';...
    'ensemblGeneId';...
    'hgncId';...
    'pdbId';...
    'respdbId';...    
    'intactId';...
    'irefindexId';...
    'bindId';...    
    'complexId';...
    'rogId';...
    'irogId';...
    'crogId';...
    'icrogId';...
    'omimId';...
    'unknownId'};

Ids=cell(2,1);
Ids{1}={};
Ids{2}={};
IdTs=cell(2,1);
IdTs{1}=[];
IdTs{2}=[];

for IdL=1:3
    switch IdL
        case 1
            String{1}=NameA{1};
            String{2}=NameB{1};
        case 2
            String{1}=AltIdA{1};
            String{2}=AltIdB{1};
        case 3
            String{1}=AliasesA{1};
            String{2}=AliasesB{1};
    end
    for NameL=1:2
        CurrString=String{NameL};
        if ~isempty(findstr(CurrString,':'))
            SepPos=regexp(CurrString,'\|');
            SepPos=[0,SepPos,length(CurrString)+1];
            for SepL=1:length(SepPos)-1
                CurrName=CurrString(SepPos(SepL)+1:SepPos(SepL+1)-1);
                %recover current gene id
                Id=regexp(CurrName,'(?<=:).*','match');
                Id=Id{1};
                %extract identificator
                Pos=findstr(Id,'(');
                if ~isempty(Pos)
                    Id=Id(1:Pos(end)-1);
                end
                %recover current type of reference
                IdType=regexp(CurrName,'.*(?=:)','match');
                IdType=IdType{1};
                if isequal(lower(IdType),'ensembl')
                    CurrId=regexp(Id,'(?<=\d*\.?)ENS.*','match');
                    if ~isempty(CurrId)
                        Id=CurrId{1};
                    end
                end
                if isequal(lower(IdType),'uniprotkb')
                    if ~isempty(findstr('_',Id))
                        Id=regexp(Id,'.*(?=_)','match');  
                        Id=Id{1};
                    end
                end
%                 if isequal(lower(IdType),'string')
%                     %remove Ensembl ID
%                     CurrId=regexp(Id,'(?<=\d*\.?)ENS.*','match');
%                     if ~isempty(CurrId)
%                         Id='';
%                     end
%                 end                
                if ~isempty(Id)
                    try
                    if isempty(strmatch(lower(Id),lower(Ids{NameL}),'exact'))
                        IdPos=length(Ids{NameL})+1;
                        Ids{NameL}{IdPos,1}=Id;

                        IdTypePos=strmatch(lower(IdType),IdTypes,'exact');
                        if~isempty(IdTypePos)
                            %exist two types of uniprotkb id
                            if isequal(lower(IdType),'uniprotkb')
                                if length(Id)==6 & regexp(Id,'^[A-Z][0-9A-Z]{5,5}')==1
                                    IdTs{NameL}(IdPos,1)=1;
                                else
                                    IdTs{NameL}(IdPos,1)=2;
                                end
                                %exist two types of uniprotkb id
                            elseif isequal(lower(IdType),'ensembl')
                                if ~isempty(regexp(Id,'^ENS'))
                                    IdTs{NameL}(IdPos,1)=7;
                                else
                                    IdTs{NameL}(IdPos,1)=8;
                                end
                            elseif isequal(lower(IdType),'entrez gene/locuslink')|isequal(lower(IdType),'entrezgene/locuslink')
                                if IdL<=2
                                    IdTs{NameL}(IdPos,1)=2;
                                else
                                    IdTs{NameL}(IdPos,1)=4;
                                end
                            else
                                IdTs{NameL}(IdPos,1)=IdTypePos;
                            end

                        else
                            %not referenced identificator type
                            'stop'
                        end
                    end
                    catch
                        'stop'
                    end
                end
            end
        end
    end
end


GeneName1=GeneName{1};
if AliasNb==2
    GeneName2=GeneName{2};
end
EmptyFlag=[];
EmptyFlag(1)=0;
EmptyFlag(2)=0;
%procesess first gene
CurrGenePos1=strmatch(lower(GeneName{1}),lower(Ids{1}),'exact');
if ~isempty(CurrGenePos1)
    GenePos(1)=1;
    CurrIdType1=IdTs{1}(CurrGenePos1);
    if AliasNb==1
        GenePos(2)=2;
        %find if there exists for the snd gene an similar identificator
        CurrGenePos2=find(IdTs{2}==CurrIdType1);
        if ~isempty(CurrGenePos2)          
            CurrGenePos2=CurrGenePos2(1);
        else
            %take for the second gene the lowest identificator rank in the IdTypes list
            [CurrIdType2 CurrGenePos2]=min(IdTs{2});
        end
        try
        GeneName2=Ids{2}{CurrGenePos2};
        catch
            'stop'
        end
    end
else
    CurrGenePos1=strmatch(lower(GeneName{1}),lower(Ids{2}),'exact');
    if ~isempty(CurrGenePos1)
        GenePos(1)=2;
        CurrIdType1=IdTs{2}(CurrGenePos1);
        if AliasNb==1
            GenePos(2)=1;
            CurrGenePos2=find(IdTs{1}==CurrIdType1);
            if ~isempty(CurrGenePos2)
                CurrGenePos2=CurrGenePos2(1);
            else
                [CurrIdType2 CurrGenePos2]=min(IdTs{1});
            end
            GeneName2=Ids{1}{CurrGenePos2};
        end
    else
        %first gene not found
        EmptyFlag(1)=1;
        if AliasNb==1
            GenePos(1)=1;
            GenePos(2)=2;
            [CurrIdType1 CurrGenePos1]=min(IdTs{1});
            %find if there exists for the snd gene an similar identificator
            CurrGenePos2=find(IdTs{2}==CurrIdType1);
            if ~isempty(CurrGenePos2)
                CurrGenePos2=CurrGenePos2(1);
            else
                %take for the second gene the lowest identificator rank in the IdTypes list
                [CurrIdType2 CurrGenePos2]=min(IdTs{2});
            end
            GeneName2=Ids{2}{CurrGenePos2};
        end
    end
end
%process eventually second gene
if AliasNb==2
    if EmptyFlag(1)==0
        if GenePos(1)==1
            GenePos(2)=2;
            CurrGenePos2=strmatch(lower(GeneName{2}),lower(Ids{2}),'exact');
            if isempty(CurrGenePos2)               
                EmptyFlag(2)=1;
                [CurrIdType2 CurrGenePos2]=min(IdTs{2});
                GeneName2=Ids{2}{CurrGenePos2};
            end
        elseif GenePos(1)==2
            GenePos(2)=1;
            CurrGenePos2=strmatch(lower(GeneName{2}),lower(Ids{1}),'exact');
            if isempty(CurrGenePos2)               
                EmptyFlag(2)=1;
                [CurrIdType2 CurrGenePos2]=min(IdTs{1});
                GeneName2=Ids{1}{CurrGenePos2};
            end
        end
    else       
        CurrGenePos2=strmatch(lower(GeneName{2}),lower(Ids{1}),'exact');
        if ~isempty(CurrGenePos2)
            GenePos(1)=2;
            GenePos(2)=1;            
            [CurrIdType1 CurrGenePos1]=min(IdTs{2});
            GeneName1=Ids{2}{CurrGenePos1};
        else
            GenePos(1)=1;
            GenePos(2)=2;
            [CurrIdType1 CurrGenePos1]=min(IdTs{1});
            GeneName1=Ids{1}{CurrGenePos1};
            CurrGenePos2=strmatch(lower(GeneName{2}),lower(Ids{2}),'exact');            
            if isempty(CurrGenePos2)               
                EmptyFlag(2)=1;
                [CurrIdType2 CurrGenePos2]=min(IdTs{2});
                GeneName2=Ids{2}{CurrGenePos2};
            end
        end
    end
end

%fill first gene aliases
if AliasNb==1 & EmptyFlag(1)==0
    for IdL=1:length(Ids{GenePos(1)})
        try
        if isempty(strmatch(Ids{GenePos(1)}{IdL},Gene.ids{GeneL},'exact'))
            Gene.ids{GeneL}{end+1,1}=Ids{GenePos(1)}{IdL};
            Gene.idTypes{GeneL}(end+1,1)=IdTs{GenePos(1)}(IdL);
        end
        catch
            'stop'
        end
    end
end


%fill MainAlias
AliasPos=[];
CurrAliasPos=strmatch(lower(GeneName1),lower(MainAlias),'exact');
if isempty(CurrAliasPos)
    MainAlias{end+1,1}=GeneName1;
    AliasPos(1)=length(MainAlias);
else
    AliasPos(1)=CurrAliasPos;
end

CurrAliasPos=strmatch(lower(GeneName2),lower(MainAlias),'exact');
if isempty(CurrAliasPos)
    MainAlias{end+1,1}=GeneName2;
    AliasPos(2)=length(MainAlias);
    else
    AliasPos(2)=CurrAliasPos;
end



%fill ID and Dict
if AliasNb==1    
    for IdL=1:length(Ids{GenePos(2)})        
        eval(sprintf('IdPos=strmatch(Ids{GenePos(2)}{IdL},Dict.%s,''exact'');',IdDict{IdTs{GenePos(2)}(IdL)}));
        if isempty(IdPos)
            eval(sprintf('IdPos=length(Dict.%s)+1;',IdDict{IdTs{GenePos(2)}(IdL)}));
        end
        eval(sprintf('Dict.%s{IdPos,1}=Ids{GenePos(2)}{IdL};',IdDict{IdTs{GenePos(2)}(IdL)}));
        ID{1}(InfoPos,IdTs{GenePos(2)}(IdL))=IdPos;
    end
end    
if AliasNb==2
    for AliasL=1:2
        for IdL=1:length(Ids{GenePos(AliasL)})
            eval(sprintf('IdPos=strmatch(Ids{GenePos(AliasL)}{IdL},Dict.%s,''exact'');',IdDict{IdTs{GenePos(AliasL)}(IdL)}));
            if isempty(IdPos)
                eval(sprintf('IdPos=length(Dict.%s)+1;',IdDict{IdTs{GenePos(AliasL)}(IdL)}));
            end
            eval(sprintf('Dict.%s{IdPos,1}=Ids{GenePos(AliasL)}{IdL};',IdDict{IdTs{GenePos(AliasL)}(IdL)}));
            ID{AliasL}(InfoPos,IdTs{GenePos(AliasL)}(IdL))=IdPos;
        end
    end
end


%Gene position in Gene1 list
Interaction(InfoPos,1)=GeneL;

%Database (aka Service in Psicquiq)  code
Interaction(InfoPos,2)=Service;

%Is gene name field empty
Interaction(InfoPos,3)=EmptyFlag(1);
if AliasNb==2
    Interaction(InfoPos,4)=EmptyFlag(2);
end

if AliasNb==1
    Interaction(InfoPos,5)=AliasPos(2);
else
    Interaction(InfoPos,5)=AliasPos(1);
    Interaction(InfoPos,6)=AliasPos(2);
end
%examples
%psi-mi:"MI:0468"(HPRD)|psi-mi:"MI:0463"(GRID)|psi-mi:"MI:0467"(reactome)
Val=Source{1};
Pos1=findstr('MI:',Val);
Pos2=findstr('"',Val);
Pos3=findstr('(',Val);
Pos4=findstr(')',Val);
if ~isempty(Pos1)
    for PosL=1:length(Pos1)
        %find code (if MI: exists, then code is terminated by ")
        CodePos=find(Pos2>Pos1(PosL));
        CurrCode=str2num(Val(Pos1(PosL)+3:Pos2(CodePos(1))-1));
        %verify if it exists in the dictionary
        DictPos=find(Dict.sourceCode==CurrCode);
        if isempty(DictPos)
            DictPos=length(Dict.source)+1;
            Dict.sourceCode(DictPos,1)=CurrCode;
            %find name (if exist, in parenthesis, just after ")
            NamePos=find(Pos3==Pos2(CodePos(1))+1);
            if ~isempty(NamePos)
                Dict.source{DictPos,1}=Val(Pos3(NamePos)+1:Pos4(NamePos)-1);
            else
                Dict.source{DictPos,1}='-';
            end
        else
            if isequal(Dict.source{DictPos},'-')
                NamePos=find(Pos3==Pos2(CodePos(1))+1);
                if ~isempty(NamePos)
                    Dict.source{DictPos}=Val(Pos3(NamePos)+1:Pos4(NamePos)-1);
                end
            end
        end
        Interaction(InfoPos,7)=DictPos;
    end
end

%species information
for SpL=1:2
    if SpL==1
        Val=SpeciesA{1};
    else
        Val=SpeciesB{1};
    end
    Pos1=findstr('taxid:',Val);
    Pos2=findstr('(',Val);
    Pos3=findstr(')',Val);
    if isempty(Pos1)
        CurrCode=0;
        CurrName='-';
    else
        if ~isempty(Pos2)
            CurrCode=str2num(Val(Pos1+6:Pos2-1));
            CurrName=Val(Pos2+1:Pos3-1);
        else
            CurrCode=str2num(Val(Pos1+6:end));
            CurrName='-';
        end
    end
    CurrDictPos=find(Dict.speciesCode==CurrCode);
    if isempty(CurrDictPos)
        DictPos(SpL)=length(Dict.speciesCode)+1;
        Dict.speciesCode(DictPos(SpL),1)=CurrCode;
        Dict.species{DictPos(SpL),1}=CurrName;
    else
        DictPos(SpL)=CurrDictPos;
        %replace eventually by scientific name
        if isempty(findstr(' ',Dict.species{DictPos(SpL)}))
            if isempty(findstr(' ',CurrName))
                Dict.species{DictPos(SpL),1}=CurrName;
            end
        end
    end
end
if AliasNb==1
    if GenePos(1)==1
        Interaction(InfoPos,8)=DictPos(2);
    else
        Interaction(InfoPos,8)=DictPos(1);
    end
else
        Interaction(InfoPos,8)=DictPos(1);
        Interaction(InfoPos,9)=DictPos(2);
end

%SCORE
Interaction(InfoPos,10)=-1;
Interaction(InfoPos,11)=-1;
Interaction(InfoPos,12)=-1;
switch lower(ServName)
    case {'bind','biogrid','apid','virhostnet'}
        %no score
    case 'genemania'
        try
            Score=regexp(ConfidenceScore{1},'(?<=function-based confidence:).*','match');
            Interaction(InfoPos,10)=eval(Score{1});
        catch
        end
    case 'intact'
        try
            Score=regexp(ConfidenceScore{1},'(?<=intact-miscore:)\d+\.?\d+','match');
            Interaction(InfoPos,10)=str2num(Score{1});
        catch
        end
    case 'irefindex'
        try
            Score=regexp(ConfidenceScore{1},'(?<=lpr:)\d+','match');
            Interaction(InfoPos,10)=str2num(Score{1});
            Score=regexp(ConfidenceScore{1},'(?<=hpr:)\d+','match');
            Interaction(InfoPos,11)=str2num(Score{1});
            Score=regexp(ConfidenceScore{1},'(?<=np:)\d+','match');
            Interaction(InfoPos,12)=str2num(Score{1});
        catch
        end
    case 'mint'
        try
            Score=regexp(ConfidenceScore{1},'(?<=mint-score:)\d+\.?\d+','match');
            Interaction(InfoPos,10)=str2num(Score{1});
        catch
        end
    case 'spike'
        try
            Score=regexp(ConfidenceScore{1},'(?<=SpikeIntegrity:)\d+','match');
            Interaction(InfoPos,10)=str2num(Score{1});
        catch
        end
    case 'string'
        try
            Score=regexp(ConfidenceScore{1},'(?<=score:)\d+','match');
            Interaction(InfoPos,10)=str2num(Score{1});
        catch
        end
    case 'uniprot'
        try
            Score=regexp(ConfidenceScore{1},'(?<=intact-miscore:)\d+.\d+','match');
            Interaction(InfoPos,10)=str2num(Score{1});
        catch
        end
    otherwise
        h=errordlg(sprintf('%s case not processed in prog psciqui',ServName));
        waitfor(h)
        error('process canceled')        
end

%FIRST AUTHOR
Interaction(InfoPos,13)=0;
if ~isempty(findstr('et al',FirstAuthor{1}))
    Author=regexp(FirstAuthor{1},'.+(?= et al)','match');
else
    Author=regexp(FirstAuthor{1},'.+(?=\()','match');
end
if ~isempty(Author)
    Author=strtrim(Author{1});
    if ~isequal(Author,'-')
        %verify if it exists in the dictionary
        DictPos=strmatch(Author,Dict.author,'exact');
        if isempty(DictPos)
            DictPos=length(Dict.author)+1;
            Dict.author{DictPos,1}=Author;
        end
        Interaction(InfoPos,13)=DictPos;
    end
end

%BIBLIOGRAPHY (PUBMED,OMIM,IMEX) => in Matrix
Val=Biblio{1};
Interaction(InfoPos,14)=0;
Interaction(InfoPos,15)=0;
Interaction(InfoPos,16)=0;
if ~isequal(Val,'-')
    for PubL=1:3
        if PubL==1
            BibName='pubmed:';
            CurrDict='Dict.pubmed';
            BibMat='Pubmed';
        elseif PubL==2
            BibName='omim:';
            CurrDict='Dict.omim';
            BibMat='Omim';
        else
            BibName='imex:IM-';
            CurrDict='Dict.imex';
            BibMat='Imex';
        end
        Pos1=findstr(BibName,Val);
        Pos2=findstr('|',Val);
        for PosL=1:length(Pos1)
            %find code
            CodePos=find(Pos2>Pos1(PosL));
            if ~isempty(CodePos)
                CurrCode=str2num(Val(Pos1(PosL)+length(BibName):Pos2(CodePos(1))-1));
            else
                CurrCode=str2num(Val(Pos1(PosL)+length(BibName):end));
            end
            %process only if it is a number
            if ~isempty(CurrCode)
                %verify if it exists in the dictionary
                eval(sprintf('DictPos=find(%sCode==CurrCode);',CurrDict));
                if isempty(DictPos)
                    eval(sprintf('DictPos=length(%sCode)+1;',CurrDict));
                    eval(sprintf('%sCode(DictPos,1)=CurrCode;',CurrDict));
                end
                %fill matrix
                eval(sprintf('%s(DictPos,InfoPos)=1;',BibMat));
                %fill Interaction
                Interaction(InfoPos,13+PubL)=DictPos;
            end
        end
    end
end


%interaction type  position
%example:
%psi-mi:"MI:0407"(direct interaction)|psi-mi:"MI:0403"(colocalization)
%|psi-mi:"MI:0915"(physical association)
Val=InteractionType{1};
Pos1=findstr('MI:',Val);
Pos2=findstr('"',Val);
Pos3=findstr('(',Val);
Pos4=findstr(')',Val);
if ~isempty(Pos1)
    if length(Pos1)>4
        h=errordlg(sprintf('there exist %u detection methods => change Interaction size',length(Pos1)));
        waitfor(h)
    end
    for PosL=1:length(Pos1)
        %find code (if MI: exists, then code is terminated by ")
        CodePos=find(Pos2>Pos1(PosL));
        CurrCode=str2num(Val(Pos1(PosL)+3:Pos2(CodePos(1))-1));
        %verify if it exists in the dictionary
        DictPos=find(Dict.interactionCode==CurrCode);
        if isempty(DictPos)
            DictPos=length(Dict.interaction)+1;
            Dict.interactionCode(DictPos,1)=CurrCode;
            %find name (if exist, in parenthesis, just after ")
            NamePos=find(Pos3==Pos2(CodePos(1))+1);
            if ~isempty(NamePos)
                Dict.interaction{DictPos,1}=Val(Pos3(NamePos)+1:Pos4(NamePos)-1);
            else
                Dict.interaction{DictPos,1}='-';
            end
        end
        Interaction(InfoPos,16+PosL)=DictPos;
    end
end
%detection method position
%examples
%psi-mi:"MI:0006"(anti bait coimmunoprecipitation)|psi-mi:"MI:0007"(anti
%tag coimmunoprecipitation)
%psi-mi:"MI:10010"
Val=DetectionMethod{1};
Pos1=findstr('MI:',Val);
Pos2=findstr('"',Val);
Pos3=findstr('(',Val);
Pos4=findstr(')',Val);
if ~isempty(Pos1)
    for PosL=1:length(Pos1)
        %find code (if MI: exists, then code is terminated by ")
        CodePos=find(Pos2>Pos1(PosL));
        CurrCode=str2num(Val(Pos1(PosL)+3:Pos2(CodePos(1))-1));
        %verify if it exists in the dictionary
        DictPos=find(Dict.methodCode==CurrCode);
        if isempty(DictPos)
            DictPos=length(Dict.method)+1;
            Dict.methodCode(DictPos,1)=CurrCode;
            %find name (if exist, in parenthesis, just after ")
            NamePos=find(Pos3==Pos2(CodePos(1))+1);
            if ~isempty(NamePos)
                Dict.method{DictPos,1}=Val(Pos3(NamePos)+1:Pos4(NamePos)-1);
            else
                Dict.method{DictPos,1}='-';
            end
        else
            if isequal(Dict.method{DictPos},'-')
                NamePos=find(Pos3==Pos2(CodePos(1))+1);
                if ~isempty(NamePos)
                    Dict.method{DictPos}=Val(Pos3(NamePos)+1:Pos4(NamePos)-1);
                end
            end
        end
        Interaction(InfoPos,19+PosL)=DictPos;
    end
end



%% READ_URL
function Read=READ_URL(Query)
Read='-1';
try
    LoopNb=0;
    Continue=1;
    while Continue
        LoopNb=LoopNb+1;
        try
            Read=urlread(Query);
            Continue=0;
        catch
            if LoopNb==10
                Continue=0;
                s=lasterror;
                if ~isempty(findstr('Error using',s.message))
                    Read='-2';
                else
                    Read='-3';
                end
            end
        end
    end
catch
    sprintf('%s not processed',Query)
end
