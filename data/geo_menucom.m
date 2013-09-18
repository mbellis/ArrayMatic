% processes GEO menu calls
%
%=======================
% FUNCTION GEO_MENUCOM
%=======================

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

function varargout=geo_menucom(Action,varargin)
global K P
switch Action
    case 'read a GPL'
%% READ A GPL        
        if nargin==1
            VERIF()
            [Species,GPL]=SELECT_GPL('multiple');
            if ~isempty(GPL)
                KeepFlag=questdlg('Do you want to keep existing GPL mat files ?','','yes','no','cancel','yes');
                if isequal(KeepFlag,'yes')
                    geo_metadb('read a GPL',GPL,1)
                elseif isequal(KeepFlag,'no')
                    geo_metadb('read a GPL',GPL,0)                    
                end
            end
        end
        
    case 'control a GPL'
%% CONTROL A GPL        
        if nargin==1
            
            GPL=inputdlg('Indicate GPL rank','');
            GPL=str2num(GPL{1});
            if isnumeric(GPL)
                geo_metadb('control a GPL',GPL,1)
            else
                h=warndlg('GPL rank must be numeric');
                waitfor(h)
            end
        end     
        
    case 'print a GPL'
%% PRINT A GPL                %% PRINT A GPL
        VERIF
        [Species,GPL]=SELECT_GPL('single');
        if ~isempty(GPL)
            geo_metadb('print a GPL',GPL)
        end

    case 'find biological conditions'        
%% FIND BIOLOGICAL CONDITIONS
        VERIF
        [Species,GPL]=SELECT_GPL('single');
        if ~isempty(GPL)
            geo_metadb('find biological conditions',GPL)
        end

    case 'import data'
%% IMPORT DATA        
        cd(K.dir.geoMetadata)
        [FileName,FileDir]=uigetfile('GPL*.mat','select a GPL');
        if isnumeric(FileName)
            h=warndlg('use ''display a GPL'' to create the GPL you want to use');
            waitfor(h)
        else
            cd(FileDir)
            load(FileName)
            if ~isfield(Gse,'imported')
                Gse.imported=zeros(length(Gse.gse),1);
            end
            FileName=regexp(FileName,'(?=.)^GPL\d+','match');
            FileName=FileName{1};
            %load data            
            GplDir=regexp(FileName,'GPL\d+','match');
            GplDir=fullfile(K.dir.geoExperiments,GplDir{1});
            if ~exist(GplDir,'dir')
                mkdir(GplDir)
            end
            cd(GplDir)
            %select GSE to be imported
            [GseSel,GseOK] = listdlg('ListString',Gse.gse,'SelectionMode','multiple','ListSize',[600,600],'Name','FTP load of GSE','PromptString','Select GSE(s)');
            SuppFilePos=[];
            for GseL=1:length(Gse.supplementaryFile)
                if ~isempty(Gse.supplementaryFile{GseL})
                    SuppFilePos(end+1,1)=GseL;
                end
            end
            if length(Gse.imported)<length(Gse.gsmNb)
                Gse.imported=[Gse.imported;zeros(length(Gse.gsmNb)-length(Gse.imported),1)];
            end
            GseSel=find(Gse.imported==0&Gse.gsmNb'>0);
            GseSel=intersect(GseSel,SuppFilePos);
            if GseOK & length(GseSel)>0     
                %eliminate doublons
                [temp,KeepIndex,temp]=unique(Gse.gseRank(GseSel));
                GseSel=GseSel(KeepIndex);
                for GseL=1:length(GseSel)                    
                    CurrGsePos=GseSel(GseL);
                    CurrGseRank=Gse.gseRank(CurrGsePos);
                    GseUpdatePos=find(Gse.gseRank==CurrGseRank);
                    GsmNb=length(find(Gsm.gseRank==CurrGseRank));
                    fprintf('\nGSE%u %u GSM (%u/%u)\n',CurrGseRank,GsmNb,GseL,length(GseSel))
                    %test if exist directory
                    %do not try to import it no GSM
                    
                    OutGsmNb=Gse.outGplGsmNb(CurrGsePos);
                    if GsmNb>0
                        if OutGsmNb<50
                            cd(GplDir)
                            if ~exist(sprintf('GSE%u_RAW.tar',CurrGseRank),'file')
                                Gse.imported(GseUpdatePos)=0;
                                try
                                    eval(sprintf('!wget -nv ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/series/GSE%u/GSE%u_RAW.tar',CurrGseRank,CurrGseRank));
                                    if exist(sprintf('GSE%u_RAW.tar',CurrGseRank),'file')
                                        Gse.imported(GseUpdatePos)=1;
                                        if isempty(Gse.supplementaryFile{CurrGsePos})
                                            for i=1:length(GseUpdatePos)
                                                Gse.supplementaryFile{GseUpdatePos(i)}=sprintf('ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/series/GSE%u/GSE%u_RAW.tar',CurrGseRank,CurrGseRank);
                                            end
                                        end
                                    else
                                        fprintf('GSE%u_RAW.tar not imported\n',CurrGseRank)
                                        Gse.imported(GseUpdatePos)=0;
                                        Gse.analyzed(GseUpdatePos)=0;
                                    end
                                catch
                                    fprintf('GSE%u_RAW.tar not imported\n',CurrGseRank)
                                    Gse.imported(GseUpdatePos)=0;
                                    Gse.analyzed(GseUpdatePos)=0;
                                end
                                close(Ftp)
                            else
                                fprintf('GSE%u already exists\n',CurrGseRank)
                                Gse.imported(GseUpdatePos)=1;
                                Gse.analyzed(GseUpdatePos)=0;
                            end
                        else
                            if exist(sprintf('%s_RAW.tar',CurrGseRank),'file')
                                eval(sprintf('delete %s_RAW.tar',CurrGseRank))
                                fprintf('%s_RAW.tar deleted\n',CurrGseRank)
                            end
                            GseDir=fullfile(GplDir,sprintf('GSE%u',CurrGseRank));
                            mkdir(GseDir)
                            cd(GseDir)                            
                            CurrGsmPos=find(Gsm.gseRank==CurrGseRank);
                            ImportedNb=0;
                            for GsmL=1:GsmNb
                                CurrGsmRank=Gsm.gsmRank(CurrGsmPos(GsmL));
                                FtpInfo=Gsm.supplementaryFile{CurrGsmPos(GsmL)};
                                if ~isempty(FtpInfo)
                                    ImportFlag=1;
                                    if ~isempty(findstr('CEL',FtpInfo))
                                        FtpFileName=regexp(FtpInfo,'(?<=suppl/).+(?=CEL)','match');                                                                
                                        FtpFileName=[FtpFileName{1},'CEL.gz'];
                                    elseif ~isempty(findstr('cel',FtpInfo))
                                        FtpFileName=regexp(FtpInfo,'(?<=suppl/).+(?=cel)','match');                                                                
                                        FtpFileName=[FtpFileName{1},'cel.gz'];
                                    else
                                        ImportFlag=0;
                                        fprintf('GSM%u not imported\n',CurrGsmRank);
                                    end
                                    if ImportFlag
                                        eval(sprintf('!wget -nv %s',FtpInfo));
                                        if exist(FtpFileName,'file')
                                            if ~isequal(FtpFileName,sprintf('GSM%u.CEL.gz',CurrGsmRank))
                                                eval(sprintf('!mv %s GSM%u.CEL.gz',FtpFileName,CurrGsmRank))
                                            end
                                            ImportedNb=ImportedNb+1;
                                        else
                                            fprintf('GSM%u not imported\n',CurrGsmRank);
                                        end
                                    end

                                end
                            end
                            if ImportedNb>0
                                Gse.imported(GseUpdatePos)=1;
                                Gse.analyzed(GseUpdatePos)=0;
                                if ImportedNb~=GsmNb
                                    fprintf('GSE%u: %u GSM imported instead of %u\n',CurrGseRank,ImportedNb,GsmNb)
                                end
                            else
                                Gse.imported(GseUpdatePos)=0;
                                Gse.analyzed(GseUpdatePos)=0;
                            end
                        end
                    end
                end
                cd(K.dir.geoMetadata)
                eval(sprintf('save %s Gsm Gse Gds Contributor',FileName))
            end
        end

    case 'do RDN analysis'
%% DO RDN ANALYSIS     
[ChipRank,ChipPos,Type,ProbeSetNb,Gpl,CompName,Chromosomes,Success] =select_chip();
if Success==0
    h=warndlg('RDN analysis canceled');
    waitfor(h)
else
    ChipName=K.chip.shortName{ChipPos};
    MyChipName=sprintf('m%u',ChipRank);
    cd(K.dir.geoMetadata)
    [FileName,FileDir]=uigetfile(sprintf('%s.mat',Gpl),'select a Gpl');
    if isnumeric(FileName)
        h=warndlg('use ''display a GPL'' to create the GPL you want to use');
        waitfor(h)
    else
        %import chip set informationif necessary
        cd(K.dir.affyMetadata)
        if ~exist(sprintf('%s_cdf.mat',MyChipName),'file')|~exist(sprintf('%s_gin.mat',MyChipName),'file')|~exist(sprintf('%s_seq.mat',MyChipName),'file')|~exist(sprintf('%s_pa.mat',MyChipName),'file')
            affy_rdn('recover chip info','','', ChipName, MyChipName,K.dir.affyMetadata,fullfile(K.dir.affyChipData,MyChipName,'libfiles'));
        end

        cd(FileDir)
        load(FileName)
        %load data
        cd(K.dir.geoExperiments)
        if ~exist(Gpl,'dir')
            h=warndlg(sprintf('%s does not exist in %s',Gpl,K.dir.data));
            waitfor(h)
        else
            UpdateFlag=questdlg('Do you want to overwrite existing analysis','','no','yes','no');
            if isequal(UpdateFlag,'yes')
                UpdateFlag=1;
            else
                UpdateFlag=0;
            end
            cd(Gpl);
            GplDir=cd;
            GsePos=find(Gse.imported);
            for GseL=1:length(GsePos)
                cd(GplDir)
                %create a directory and decompress files
                GSE=Gse.gse{GsePos(GseL)};
                GseRank=Gse.gseRank(GsePos(GseL));
                TarFile=sprintf('%s_RAW.tar',GSE);
                if ~exist(TarFile,'file')
                    h=warndlg(sprintf('%s does not exist in %s',TarFile,Gpl));
                    waitfor(h)
                else
                    if ~exist(GSE,'dir')
                        mkdir(GSE)
                    end
                    cd(GSE)
                    GseDir=cd;
                    if ~exist(sprintf('rdn_%s.txt',GSE),'file')|UpdateFlag
                        cd(GplDir)
                        %decompress
                        untar(TarFile,GSE)
                        cd(GseDir)
                        %delete EXP files
                        delete('*.EXP.gz')
                        %eliminate GSM that does not belong to the current
                        %Gpl
                        GsmPos=find(Gsm.gseRank==GseRank);
                        DoIt=1;
                        BiolNames={};
                        AllBiolNames={};
                        Replicates=[];
                        CelNames={};
                        for GsmL=1:length(GsmPos)
                            if ~isequal(Gsm.gpl{GsmPos(GsmL)},Gpl)
                                delete(sprintf('%s.CEL.gz',Gsm.gsm{GsmPos(GsmL)}))
                            else
                                %verify it exists
                                if ~exist(sprintf('%s.CEL.gz',Gsm.gsm{GsmPos(GsmL)}),'file')&~exist(sprintf('%s.cel.gz',Gsm.gsm{GsmPos(GsmL)}),'file')
                                    h=warndlg(sprintf('%s does not exist - process canceled',Gsm.gsm{GsmPos(GsmL)}));
                                    waitfor(h)
                                    DoIt=0;
                                    break
                                else

                                    if isempty(strmatch(Gsm.biolName{GsmPos(GsmL)},BiolNames))
                                        BiolNames{end+1,1}=Gsm.biolName{GsmPos(GsmL)};
                                    end
                                    AllBiolNames{end+1,1}=Gsm.biolName{GsmPos(GsmL)};
                                    Replicates=[Replicates;Gsm.replicate(GsmPos(GsmL))];
                                    CelNames{end+1,1}=sprintf('%s.CEL',Gsm.gsm{GsmPos(GsmL)});
                                    if exist(sprintf('%s.CEL.gz',Gsm.gsm{GsmPos(GsmL)}),'file')
                                        gunzip(sprintf('%s.CEL.gz',Gsm.gsm{GsmPos(GsmL)}));
                                    else
                                        gunzip(sprintf('%s.cel.gz',Gsm.gsm{GsmPos(GsmL)}));
                                    end
                                end
                            end
                        end
                        %delete gz files
                        delete *.gz
                        if DoIt
                            %construct sample name list and order list
                            SampleNames={};
                            PrintOrder=[];
                            [CelNames CelOrder]=sort(CelNames);
                            AllBiolNames=AllBiolNames(CelOrder);
                            Replicates=Replicates(CelOrder);
                            for BiolL=1:length(BiolNames)
                                RepPos=strmatch(BiolNames{BiolL},AllBiolNames,'exact');
                                [temp RepOrder]=sort(Replicates(RepPos));
                                for RepL=1:length(RepOrder)
                                    SampleNames{end+1,1}=sprintf('%s_r%u',BiolNames{BiolL},RepL);
                                    PrintOrder=[PrintOrder,RepPos(RepOrder(RepL))];
                                end
                            end
                            %write a description file used by
                            %simpleaffy
                            fid=fopen(sprintf('%s_description.txt',GSE),'w');
                            fprintf(fid,' CEL_FILE BIOL_COND\n')
                            for CelL=1:length(CelNames)
                                if exist(CelNames{CelL},'file')
                                    fprintf(fid,'%s %s_%u\n',CelNames{CelL},AllBiolNames{CelL},Replicates(CelL));
                                elseif exist(strrep(CelNames{CelL},'CEL','cel'),'file')
                                    fprintf(fid,'%s %s_%u\n',strrep(CelNames{CelL},'CEL','cel'),AllBiolNames{CelL},Replicates(CelL));
                                end
                            end
                            fclose(fid)

                            sprintf('doing %s (%u upon %u)',GSE,GseL,length(GsePos))
                            %analyze RDN
                            if ~exist(sprintf('%s_cel.mat',GSE),'file')
                                affy_rdn('import cel files',GSE,GseDir);
                            end
                            affy_rdn('do rdn analysis',GSE,GseDir,ChipName,MyChipName,K.dir.affyMetadata,fullfile(K.dir.affyChipData,ChipName,'LibFiles'));
                            affy_rdn('print rdn signals',GSE,GseDir,ChipName,MyChipName,K.dir.affyMetadata,fullfile(K.dir.affyChipData,ChipName,'LibFiles'),SampleNames,PrintOrder);
                        end
                    end
                end
            end
        end
    end
end


    case 'do RMA analysis'
%% DO RMA ANALYSIS            
[ChipRank,ChipPos,Type,ProbeSetNb,Gpl,CompName,Chromosomes,Success] =select_chip();
% if ChangeFlag==1 => uses GSM.supplementaryFile toe eventually change the GSM file name
ChangeFlag=1;
if Success==0
    h=warndlg('RMA analysis canceled');
    waitfor(h)
else
    cd(K.dir.geoMetadata)
    [FileName,FileDir]=uigetfile(sprintf('%s.mat',Gpl),'select a GPL');
    if isnumeric(FileName)
        h=warndlg('use ''display a GPL'' to create the GPL you want to use');
        waitfor(h)
    else
        cd(FileDir)
        load(FileName)
        cd(K.dir.geoExperiments)
        if ~exist(Gpl,'dir')
            h=warndlg(sprintf('%s does not exist in %s',Gpl,K.dir.geoExperiments));
            waitfor(h)
        else
            JadeFlag=questdlg('Do you want to make analysis on Jade','','no','yes','yes');
            if isequal(JadeFlag,'yes')
                JadeFlag=1;
                QFid=fopen(sprintf('%s.pbs',Gpl),'w');
                fprintf(QFid,'#PBS -S /bin/bash\n');
                fprintf(QFid,'#PBS -N %s\n',Gpl);
                fprintf(QFid,'#PBS -e %s.mclerr\n',Gpl);
                fprintf(QFid,'#PBS -o %s.mcllog\n',Gpl);
                fprintf(QFid,'#PBS -l walltime=2:00:00\n');
                fprintf(QFid,'#PBS -l select=1:ncpus=8:mpiprocs=9\n');
                fprintf(QFid,'module unload intel\n');
                fprintf(QFid,'module load intel/11.1.072\n');
                fprintf(QFid,'module load numpy/1.5.1\n');
                fprintf(QFid,'module load scipy/0.8.0\n');
                fprintf(QFid,'module load pserie/9.3.26\n');
                fprintf(QFid,'export MPI_DSM_DISTRIBUTE=ON\n');
                fprintf(QFid,'cd /scratch/cinbell/data/experiments/geo/%s\n',Gpl);
                fprintf(QFid,'mpiexec pserie_load_balanced < %s.sh\n',Gpl);
                fclose(QFid);
                ShFid=fopen(sprintf('%s.sh',Gpl),'w');
            else
                JadeFlag=0;
            end
            UpdateFlag=questdlg('Do you want to overwrite existing analysis','','no','yes','no');
            if isequal(UpdateFlag,'yes')
                UpdateFlag=1;
                GsePos=find(Gse.imported);
            else
                UpdateFlag=0;
                GsePos=find(Gse.imported==1 & Gse.analyzed==0);
            end
            %eliminate doublons
            GsePos=sort(GsePos);
            [temp,KeepIndex,temp]=unique(Gse.gseRank(GsePos));
            GsePos=GsePos(KeepIndex);
            %eliminate Gse without supplementary data
            ClearPos=[];
            for GseL=1:length(GsePos)                
                if isempty(Gse.supplementaryFile{GsePos(GseL)})
                    ClearPos(end+1,1)=GseL;
                end
            end
            if ~isempty(ClearPos)
                GsePos(ClearPos)=[];
            end
            cd(Gpl);
            GplDir=cd;
            if JadeFlag
                DataDir=sprintf('/scratch/cinbell/data/experiments/geo/%s',Gpl);
            else
                DataDir=GplDir;
            end
            %create list of Gse to be processed in the current round
            if ~isempty(GsePos)
                MaxGsmNb=inputdlg('max nb of GSM','',1,{'400'});
                MaxGsmNb=str2num(MaxGsmNb{1});
                Processed=0;
                for GseL=1:length(GsePos)
                    %create a directory and decompress files
                    CurrGsePos=GsePos(GseL);                    
                    GSE=Gse.gse{CurrGsePos};                    
                    CurrGse=GSE;
                    GseRank=Gse.gseRank(CurrGsePos);
                    GseUpdatePos=find(Gse.gseRank==GseRank);
                    %                     TarFile=sprintf('%s_RAW.tar',GSE);
                    %                     if ~exist(TarFile,'file')
                    %                         %h=warndlg(sprintf('%s does not exist in %s',TarFile,Gpl));
                    %                         %waitfor(h)
                    %                         fprintf('%s does not exist in %s\n',TarFile,Gpl);
                    %                         Gse.imported(GsePos(GseL))=0;
                    %                     else
                    %create R script                    
                    CurrGseRank=GseRank;
                    DoIt=1;
                    if JadeFlag==0
                        cd(GplDir)
                        if exist(CurrGse,'dir')
                            cd(CurrGse)
                            if exist(sprintf('%s_rma.txt',CurrGse),'file')
                                DoIt=0;
                                Gse.analyzed(GseUpdatePos,1)=1;
                            end
                        end
                    end
                    GsmPos=find(Gsm.gseRank==CurrGseRank);
%                     if isfield(Gse,'trueGsmNb')
%                         GsmNb=Gse.trueGsmNb(CurrGsePos);
%                     else
                        GsmNb=length(GsmPos);
%                     end
                    OutGsmNb=Gse.outGplGsmNb(CurrGsePos);                                                        
                    if GsmNb>0 & GsmNb<MaxGsmNb
                        Processed=Processed+1;
                        if OutGsmNb<50
                            GseFlag=1;
                        else
                            GseFlag=0;
                        end                                        
                        cd(GplDir)
                        fid=fopen(sprintf('%s_%u.r',Gpl,GseRank),'w');
                        fprintf('%s_%u.r\n',Gpl,GseRank);
                        if JadeFlag
                            fprintf(ShFid,'/home/cinbell/mywork/R-2.12.2/bin/R CMD BATCH ./%s_%u.r\n',Gpl,GseRank);
                        end
                        fprintf(fid,'library(affy)\n');

                        %untar GSM files
                        fprintf(fid,'setwd("%s")\n',DataDir);
                        if GseFlag
                            fprintf(fid,'try(system("mkdir %s"),silent=TRUE)\n',CurrGse);
                            fprintf(fid,'system("tar -xf %s_RAW.tar --directory=./%s")\n',CurrGse,CurrGse);
                        end
                        fprintf(fid,'setwd("./%s")\n',CurrGse);
                        %rename files with standardized form GSMxxxx.CEL.gz
                        if JadeFlag
                            fprintf(fid,'try(system("/home/cinbell/prename ''s/gsm/GSM/'' gsm*"),silent=TRUE)\n');                        
                            fprintf(fid,'try(system("/home/cinbell/prename ''s/cel/CEL/'' *cel*"),silent=TRUE)\n');
                            %sometimes GSM have something between the number and .CEL.gz
                            fprintf(fid,'try(system("/home/cinbell/prename ''s/_.+CEL/\\\\.CEL/'' *CEL*"),silent=TRUE)\n');
                        else
                        fprintf(fid,'try(system("rename ''s/gsm/GSM/'' gsm*"),silent=TRUE)\n');                        
                        fprintf(fid,'try(system("rename ''s/cel/CEL/'' *cel*"),silent=TRUE)\n');
                        %sometimes GSM have something between the number and .CEL.gz
                        fprintf(fid,'try(system("rename ''s/_.+CEL/\\\\.CEL/'' *CEL*"),silent=TRUE)\n');
                        end
                        %eventually remove unwanted GSM files
                        for GsmL=1:length(GsmPos)
                            SupplFile=Gsm.supplementaryFile{GsmPos(GsmL)};
                            if ~isempty(findstr('CEL',SupplFile))
                                GsmFileName=regexp(SupplFile,'(?<=n/GSM\d+/).+(?=CEL)','match');
                                GsmFileName=[GsmFileName{1},'CEL.gz'];
                            else                                
                                GsmFileName=regexp(SupplFile,'(?<=n/GSM\d+/).+(?=cel)','match');
                                try
                                    GsmFileName=[GsmFileName{1},'cel.gz'];
                                catch
                                    GsmFileName='';
                                end
                            end                            
                            if ~isempty(GsmFileName)
                                if ChangeFlag
                                    if ~isempty(findstr('suppl/',GsmFileName))
                                        GsmFileName=regexp(GsmFileName,'(?<=suppl/).+','match');
                                        GsmFileName=GsmFileName{1};
                                    end
                                    if ~isequal(Gsm.gpl{GsmPos(GsmL)},Gpl)
                                        fprintf(fid,'unlink("%s")\n',GsmFileName);
                                    else
                                        %gunzip files
                                        fprintf(fid,'system("gunzip %s")\n',GsmFileName);
                                        GsmFileName=regexp(GsmFileName,'.+(?=.gz)','match');
                                        GsmFileName=GsmFileName{1};
                                        if ~isequal(GsmFileName,sprintf('GSM%u.CEL',Gsm.gsmRank(GsmPos(GsmL))))
                                            fprintf(fid,'system("mv %s GSM%u.CEL")\n',GsmFileName,Gsm.gsmRank(GsmPos(GsmL)));
                                        end
                                    end
                                else
                                    if ~isequal(Gsm.gpl{GsmPos(GsmL)},Gpl)
                                        fprintf(fid,'unlink("GSM%u.CEL.gz")\n',Gsm.gsmRank(GsmPos(GsmL)));
                                    else
                                        fprintf(fid,'system("gunzip GSM%u.CEL.gz")\n',Gsm.gsmRank(GsmPos(GsmL)));
                                    end
                                end
                            end
                        end
                        %clear all other gz files
                        fprintf(fid,'unlink("*.gz")\n');

                        %load cel files
                        fprintf(fid,'data=ReadAffy()\n');
                        fprintf(fid,'result=rma(data)\n');
                        fprintf(fid,'write.table(assayData(result)$exprs,file="%s_rma.txt",col.names=sampleNames(phenoData(result)),row.names=featureNames(featureData(result)),sep="\\t",na="NaN",dec=".",quote=FALSE)\n',CurrGse);
                        fprintf(fid,'try(unlink("*.CEL"),silent=TRUE)\n');
                        fprintf(fid,'try(unlink("*.cel"),silent=TRUE)\n');
                        fprintf(fid,'try(remove(data,result),silent=TRUE)\n');
                        %                                 if JadeFlag
                        %                                     %write result file
                        %                                     fprintf(fid,'try(system("mkdir /store/cinbell/data/experiments/geo/%s/%s"),silent=TRUE)\n',Gpl,CurrGse);
                        %                                     fprintf(fid,'try(system("mv * /store/cinbell/data/experiments/geo/%s/%s"),silent=TRUE)\n',Gpl,CurrGse);
                        %                                 end
                        fprintf(fid,'quit("no")\n');
                        fclose(fid);
                    else
                        fprintf('GSE%u %u GsmPos %u TrueGsmNb %u GsmNb \n',GseRank,length(GsmPos),Gse.trueGsmNb(CurrGsePos),Gse.gsmNb(CurrGsePos))
                    end                    
                end %for GseL
            end %if ~isempty(GsePos)
            if JadeFlag
                fclose(ShFid);
            end
        end %if ~exist(Gpl,'dir')
    end %if isnumeric(FileName)
    if JadeFlag==0
        cd(FileDir)
        cd(K.dir.geoMetadata)
        eval(sprintf('save %s Gsm Gse Gds Contributor',FileName))
    end
    fprintf('%u GSE processed\n',Processed)
end %if Success==0


    case 'create new network'
%% CREATE NEW NETWORK
        WARNFLAG=1;
        GplRank=varargin{1};
        DataFid=varargin{2};
        cd(K.dir.geoMetadata)
        [GplName,GplDir]=uigetfile(sprintf('GPL%u.mat',GplRank),'select a GPL');
        Continue=1;
        if isnumeric(GplName)
            h=warndlg('use ''display a GPL'' to create the GPL you want to use');
            waitfor(h)
            Continue=0;
        end
        %first Continue
        if Continue
            cd(GplDir)
            load(GplName)
            %ask for geo data directory
            cd(K.dir.geoExperiments)
            GeoDir=uigetdir('*',sprintf('Select the directory containing GPL%u data',GplRank));
            if ~ischar(GeoDir)
                h=warndlg('Must have a Geo Directory. Process canceled');
                waitfor(h)
                Continue=0;
            end
            %snd Continue
            if Continue
                P.dir.geoData=GeoDir;
                %ask for algorithm type
                AlgoList={'dchip','gcrma','hook','ilm','mas5','pdnn','plier','qfarm','rdn','rma','other'};
                Selected=listdlg('promptstring','select the algorithm used','selectionmode','single','liststring',AlgoList);
                if isempty(Selected)
                    h=warndlg('Must select an algorithm. Process canceled');
                    waitfor(h)
                    Continue=0;
                end
                %third Continue
                if Continue
                    Algorithm=AlgoList{Selected};
                    if isequal(Algorithm,'other')
                        Algorithm=inputdlg('indicate the algorithm name used as an extension in result files','')
                        if isempty(Selected)
                            h=warnddlg('Must indicate an algorithm. Process canceled');
                            waitfor(h)
                            Continue=0;
                        else
                            Algorithm=Algorithm{1};
                        end
                    end
                    %fourth Continue
                    if Continue
                        %record points
                        cd(P.dir.geoData)
                        LogFid=fopen(sprintf('log_%s_%s.txt',GplName,date),'w');
                        PointOffset=0;
                        BiolOffset=0;
                        BiolNames={};
                        DataSignals=[];
                        ExpRank=0;

                        P.exp.name={};
                        P.exp.source={};
                        P.exp.reference={};
                        P.exp.pointIndex={};
                        P.exp.used=[];
                        P.exp.nb=0;

                        P.biol.name={};
                        P.biol.pointIndex={};
                        P.biol.used=[];
                        P.biol.pairs=[];
                        P.biol.nb=0;
                        
                        
                        P.point.name={};
                        P.point.expRank=[];
                        P.point.biolRank=[];
                        P.point.replicateRank=[];
                        P.point.used=[];
                        P.point.algo={};
                        P.point.minSignal=[];
                        P.point.maxSignal=[];
                        P.point.nanNb=[];
                        P.point.negNb=[];
                        P.point.nullNb=[];
                        P.point.threshNb=[];
                        P.point.threshVal=[];
                        P.point.negThreshNb=[];
                        P.point.nullDiffNb=[];
                        P.point.runNb=[];
                        P.point.negRunNb=[];
                        P.point.runVal=[];
                        P.point.diversity=[];
                        P.point.nb=0;
                        P.point.currNb=0;
                        
                        


                        ProbeSetNb=K.chip.probeSetNbs(P.chip.chipPos);
                        for GseL=1:length(Gse.gse)
                            CurrGse=Gse.gse{GseL};
                            cd(P.dir.geoData)
                            if exist(CurrGse,'dir')
                                cd(CurrGse)
                                CurrFile=sprintf('%s_%s.txt',CurrGse,Algorithm);
                                if exist(CurrFile,'file')
                                    %read the first line to recover GSM
                                    fid=fopen(CurrFile,'r');
                                    FirstLine=fgetl(fid);
                                    fclose(fid)
                                    Gsms=regexp(FirstLine,'GSM\d+','match');

                                    %control that Gsm exist in Gsm.gsm
                                    ExistGsm=zeros(1,length(Gsms));
                                    GsmPos=[];
                                    for GsmL=1:length(Gsms)
                                        CurrGsmPos=strmatch(Gsms{GsmL},Gsm.gsm,'exact');
                                        if ~isempty(CurrGsmPos)
                                            ExistGsm(GsmL)=1;
                                            GsmPos(end+1)=CurrGsmPos;
                                            fprintf(LogFid,'1\t%s\t%s\n',CurrGse,Gsms{GsmL});
                                        else
                                            if WARNFLAG
                                                h=warndlg(sprintf('%s does not exist in Gsm database',Gsms{GsmL}));
                                                waitfor(h)
                                            end
                                            fprintf(LogFid,'0\t%s\t%s\n',CurrGse,Gsms{GsmL});
                                        end
                                    end

                                    if ~isempty(GsmPos)

                                        %LOAD AND VERIFY DATA

                                        %Load algorithm output (processed signals)
                                        ExistPos=find(ExistGsm);
                                        Gsms=Gsms(ExistPos);
                                        PointNb=length(GsmPos);
                                        Output='[ProbeSet';
                                        DataSignals=[];
                                        Output=strcat(Output,',DataSignals(:,1)');
                                        DataType='%s';
                                        DataType=strcat(DataType,'%n');
                                        if PointNb>=2
                                            for PointL=2:PointNb
                                                Output=strcat(Output,sprintf(',DataSignals(:,%u)',PointL));
                                                DataType=strcat(DataType,'%n');
                                            end
                                        end
                                        Output=strcat(Output,']');
                                        eval(sprintf('%s=textread(CurrFile,''%s'',''headerlines'',1,''delimiter'',''\t'',''emptyvalue'',NaN);',Output,DataType));
                                        %keep only GSM existing in Gsm
                                        DataSignals=DataSignals(:,ExistPos);

                                        if length(ProbeSet)~=ProbeSetNb
                                            if WARNFLAG
                                                h=warndlg(sprintf('%s will not be processed (has %u probe sets instead of %u)',CurrGse,length(ProbeSet),ProbeSetNb));
                                                waitfor(h)
                                            end
                                            % gse data do not have the right size
                                            fprintf(LogFid,'0\t%s\t-2\n',CurrGse);
                                        else

                                            % detect negative values
                                            PosNb=zeros(size(DataSignals,2),1);
                                            for DataL=1:size(DataSignals,2)
                                                PosNb(DataL)=length(find(DataSignals(:,DataL)));
                                            end
                                            if length(find(PosNb~=ProbeSetNb))
                                                if WARNFLAG
                                                    h=warndlg(sprintf('%u points have negative signals',length(find(PosNb~=ProbeSetNb))));
                                                    waitfor(h)
                                                end
                                                fprintf(LogFid,'1\t%s\t%u\n',CurrGse,length(find(PosNb~=ProbeSetNb)));
                                            end

                                            % recover probe set list and probe set nb
                                            if ExpRank==0
                                                P.chip.currProbeSetNb=ProbeSetNb;
                                                P.chip.probeSetNb=ProbeSetNb;
                                                P.chip.probeSetIds=ProbeSet;
                                                ProbeSetList=ProbeSet;
                                            else

                                                %reorder eventually results
                                                if sum(cellfun(@isequal,ProbeSet,ProbeSetList))~=length(ProbeSetList)
                                                    %find SortIndex
                                                    sprintf('reordering probe sets for %s',CurrGse)
                                                    FindAllFlag=1;
                                                    [SortIndex,Temp]=find_sortorder(ProbeSetList,ProbeSet,'string',FindAllFlag);
                                                    for PointL=1:PointNb
                                                        DataSignals(:,PointL)=DataSignals(SortIndex,PointL);
                                                    end
                                                end
                                            end

                                            % CALCULATE RANKS AND SAVE DATA
                                            %detect log values
                                            for PointL=1:PointNb
                                                if max(DataSignals(:,PointL))<20
                                                    DataSignals(:,PointL)=2.^DataSignals(:,PointL);
                                                end
                                            end
                                            
                                            DataRanks=[];
                                            for PointL=1:PointNb
                                                trs_pointproperties(DataSignals(:,PointL),PointOffset+PointL);
                                                [DataRanks(:,PointL),Success]=signal2rank(DataSignals(:,PointL),0,1,P.point.threshVal(PointOffset+PointL));
                                                if ~ Success                                                
                                                    if WARNFLAG
                                                        h=warndlg(sprintf('DataRanks of %s not calculated',Gsms(PointL)));
                                                        waitfor(h)
                                                    end
                                                    fprintf(LogFid,'-2\t%s\t%s\n',CurrGse,Gsms{PointL});
                                                end
                                            end
                                            DataSignals=[];

                                            %save data
                                            
                                            cd(P.dir.data)
                                            Success=save_data(DataRanks,DataFid,[],'a','single','ieee-le');
                                            if Success==0
                                                if WARNFLAG
                                                    h=warndlg('Data not saved. Process canceled');
                                                    waitfor(h)
                                                end
                                                %gse data have not been saved
                                                fprintf(LogFid,'0\t%s\t-1\n',CurrGse);
                                            else
                                                


                                                %FILL INFORMATIONS

                                                %Experiment informations
                                                ExpRank=ExpRank+1;
                                                P.exp.name{ExpRank,1}=Gse.title{GseL};
                                                P.exp.source{ExpRank,1}='GEO';
                                                P.exp.reference{ExpRank,1}=CurrGse;
                                                P.exp.pointIndex{ExpRank,1}=[PointOffset+1:PointOffset+PointNb]';
                                                P.exp.nb=length(P.exp.name);
                                                P.exp.used=ones(P.exp.nb,1);

                                                %Point informations
                                                for PointL=1:PointNb
                                                    CurrPoint=PointOffset+PointL;
                                                    CurrGsmPos=GsmPos(PointL);
                                                    P.point.name{CurrPoint,1}=sprintf('%s_%s_R%u',Gsm.biolName{CurrGsmPos},lower(CurrGse),Gsm.replicate(CurrGsmPos));
                                                    P.point.replicateRank(CurrPoint,1)=Gsm.replicate(CurrGsmPos);
                                                    
                                                    %Biological condition information
                                                    BiolName=sprintf('%s_%s',Gsm.biolName{CurrGsmPos},lower(CurrGse));

                                                    BiolRank=strmatch(BiolName,P.biol.name,'exact');
                                                    if isempty(BiolRank)
                                                        BiolRank=length(P.biol.name)+1;
                                                        P.biol.pointIndex{BiolRank,1}=[];
                                                    end

                                                    P.biol.name{BiolRank,1}=BiolName;
                                                    P.biol.pointIndex{BiolRank,1}=[P.biol.pointIndex{BiolRank,1};PointOffset+PointL];

                                                    P.point.biolRank(CurrPoint,1)=BiolRank;
                                                    P.point.expRank(CurrPoint,1)=ExpRank;
                                                    P.point.algo{CurrPoint,1}=Algorithm;
                                                end                                                
                                                P.point.nb=length(P.point.name);
                                                P.point.currNb= P.point.nb;
                                                P.point.used=ones(P.point.nb,1);
                                                P.biol.used=ones(length(P.biol.name),1);                                                
                                                P.biol.nb=length(P.biol.name);
                                                P.biol.pairs=cell(P.biol.nb,1);                                               
                                                PointOffset= P.point.nb;
                                            end %data saved
                                        end % if length(ProbeSet)~=ProbeSetNb
                                    else %not any gsm in Gsm
                                        fprintf(LogFid,'0\t%s\t-3\n',CurrGse);
                                    end
                                else %no output data
                                    fprintf(LogFid,'0\t%s\t-4\n',CurrGse);
                                end
                            else % no gse directory
         -                       fprintf(LogFid,'0\t%s\t-5\n',CurrGse);                            
                            end 
                        end %of GseL
                        fclose(LogFid)
                    end % fourth Continue
                end %third Continue
            end %snd Continue
        end %first Continue
        varargout{1}=Continue;
end


%% FUNCTION VERIF
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

%% FUNCTION SELECT GPL
function [Species,GPL]=SELECT_GPL(SelectionFlag)
global K
cd(K.dir.geoMetadata)
Species={};
Sel=mksqlite('SELECT organism FROM gpl');
Rank=0;
for SelL=1:length(Sel)
    CurrSpecies=Sel(SelL).organism;
    if ~isempty(CurrSpecies)
        Rank=Rank+1;
        Species{Rank,1}=CurrSpecies;
    end
end
Species=unique(Species);

%SELECT ONE SPECIES
SpeciesSel=listdlg('liststring',Species,'PromptString','select one','SelectionMode','single','ListSize',[500,600]);
if ~isempty(SpeciesSel)
    %display GPL
    eval(sprintf('Sel=mksqlite(''SELECT gpl,title,manufacturer,technology FROM gpl WHERE organism=''''%s'''''');',Species{SpeciesSel}));
    Rank=0;
    Gpls={};
    Titles={};
    Technologies={};
    Manufacturers={};
    for SelL=1:length(Sel)
        CurrGpl=Sel(SelL).gpl;
        if ~isempty(CurrGpl)
            Rank=Rank+1;
            Gpls{Rank,1}=CurrGpl;
            Technologies{Rank,1}=Sel(SelL).technology;
            if length(Sel(SelL).manufacturer)>30
                Manufacturers{Rank,1}=Sel(SelL).manufacturer(1:30);
            else
                Manufacturers{Rank,1}=Sel(SelL).manufacturer;
            end
            if length(Sel(SelL).title)>50
                Titles{Rank,1}=Sel(SelL).title(1:50);
            else
                Titles{Rank,1}=Sel(SelL).title;
            end

        end
    end
    [Gpls,Index]=unique(Gpls);
    Technologies=Technologies(Index);
    Manufacturers=Manufacturers(Index);
    Titles=Titles(Index);

    %select a subset of manufacturerrs
    ManKeywords=inputdlg('enter list of the Keywords to select Manufacturers /nyou wrote down with this format/n  {'''','''',''''}','',1,{'{''Affymetrix''}'});
    if isempty(ManKeywords)
        ManPos=1:length(Gpls);
    elseif isempty(ManKeywords{1})
        ManPos=1:length(Gpls);
    else
        ManKeywords=ManKeywords{1};        
        eval(sprintf('ManKeywords=%s;',ManKeywords))
        %find the manufacturers corresponding to the selected keywords
        ManPos=[];
        for ManL=1:length(Manufacturers)
            for KeyL=1:length(ManKeywords)
                if ~isempty(findstr(upper(ManKeywords{KeyL}),upper(Manufacturers{ManL})))
                    ManPos=[ManPos;ManL];
                end
            end
        end
    end
    ManPos=unique(ManPos);
    Gpls=Gpls(ManPos);
    Technologies=Technologies(ManPos);
    Manufacturers=Manufacturers(ManPos);    
    Titles=Titles(ManPos); 
    %fill info
    Information=cell(length(Gpls),1);
    GseNb=zeros(length(Gpls),1);
    for GplL=1:length(Gpls)    
        Information{GplL}=Gpls{GplL};
        Information{GplL}=sprintf('%s - %s - %s - %s',Information{GplL},Titles{GplL},Technologies{GplL},Manufacturers{GplL});
        eval(sprintf('Sel=mksqlite(''SELECT gse FROM gse_gpl WHERE gpl=''''%s'''''');',Gpls{GplL}));
        Information{GplL}=sprintf('%s - %u GSE',Information{GplL},length(Sel));
        GseNb(GplL)=length(Sel);
    end
    [GseNb,SortIndex]=sort(GseNb);
    SortIndex=flipud(SortIndex);
    Information=Information(SortIndex);
    Gpls=Gpls(SortIndex);
    if isequal(SelectionFlag,'single')
        GplSel=listdlg('liststring',Information,'PromptString','select one','SelectionMode','single','ListSize',[1000,800]);
        GPL=Gpls{GplSel};
    else
        GplSel=listdlg('liststring',Information,'PromptString','select multiple','SelectionMode','multiple','ListSize',[1000,800]);
        GPL=Gpls(GplSel);
    end    
end
pwd
Species=Species{SpeciesSel};
