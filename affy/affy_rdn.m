% RDN scripts developped by  Marc Hulsman, 2010
%
%===========
% AFFY_RDN %
%===========
%
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The
% Netherlands
% http://bioinformatics.tudelft.nl/users/marc-hulsman
%
%
%  RDN(ACTION,EXPNAME,CELDIR,CHIPNAME,CHIPDIR,LIBDIR,SAMPLENAMES,PRINTORDER);
%
% INPUT
% Action: Action to apply:
%           'cel': load cel files and save them in rdn format for further
%                  processing.
%           'info': recover all information on chip, probe sets and
%                   targeted genes necessary to rdn algorithms
%           'rdn': rdn algorithm
%           'print': print results (probe set name, array factor)
% ExpName: Experiment name used to name saved results
% CelDir: Directory where cel files are stored. All CEL files present in
%         this directory will be loaded. Pay attention that only wanted CEL files
%         are present in this directory, and that no other file containing the
%         substring 'cel' in its name is present.
% ChipName: The name of the chip. cdf, gin, probe and target sequence files
%           used in affy_rdn('info') must have their name compatible with
%           ChipName (e.g.<ChipName>.cdf,<ChipName>.gin,<ChipName>.probe_tab,<ChipName>.target)
% MyChipName: The name of information files created (<MyChipName_cdf.mat,
%             <MyChipName_gin.mat,<MyChipName_pa.mat,<MyChipName_seq.mat).
%          if MyChipName='', ChipName is used instead.     
% ChipDir: The directory where to store files generated by
%          affy_rdn('info')
% LibDir: The directory where to read cdf,gin,probe_tab and target files.
%         If LibDir=<ChipDir>/LibFiles. If LibDir='', LibDir is assumed to
%         be <ChipDir>/LibFiles.
% SampleNames: The list of sample names used to print the results. If
%               sampleNames={}, cel file names found in CelDir are used
%               instead.
%
% PrintOrder: PrintOrder is a vector indicating the order print function must use.
%             If SampleNames is not empty, it is not reordered, and it is assumed that
%             PrintOrder is intended to print results in the order corresponding to
%             SampleNames but if SampleName is empty, the cel file names will be reordered
%             to fit PrintOrder. If PrintOrder is empty, the alphabetical order of cel
%             file names is used.
%
%
%
% OUTPUT
%   No output. Results are written in CelDir.
%
% VERSION
%   V01 - 12/07/2010 - works with RDN v0.20

% (c) - Michel Bellis, 2010
% CRBM -CNRS
% 1919 rte de Mende, 34293 Montpellier cedex5
% France

%affy_rdn('import cel files','aorteHP','/home/mbellis/sosma/raydataraw/2011_SOULLIER_ALGOS/cel')
%affy_rdn('recover chip info','aorteHP','/home/mbellis/sosma/raydataraw/2011_SOULLIER_ALGOS/cel',...
%'Rat230_2','m33','/usr/data/chips/affymetrix/m033')
%affy_rdn('do rdn analysis','aorteHP','/home/mbellis/sosma/raydataraw/2011_SOULLIER_ALGOS/cel',...
%'Rat230_2','m33','/usr/data/chips/affymetrix/m033')
%affy_rdn('print rdn signals','aorteHP','/home/mbellis/sosma/raydataraw/2011_SOULLIER_ALGOS/cel',...
%'','','','',{'NP-R1','NP-R2','NP-R3','NP-R4','HP-R1','HP-R2','HP-R3','HP-R4'})



function affy_rdn(Action,ExpName,CelDir,ChipName,MyChipName,ChipDir,LibDir,SampleNames,PrintOrder)

switch Action    
    case 'import cel files'
%% IMPORT CEL FILES        
        if nargin~=3
            h=errordlg('affy_rdn(''cel'') needs at least three parameters : Action, ExpName,CelDir');
            waitfor(h)
        end
        %IMPORT RAW CEL DATA
        %force Dir to be terminated by a separator
        if CelDir(end)~=filesep
            CelDir=[CelDir,filesep];
        end
        cd(CelDir)
        celpaths = dir_filter(CelDir,'.cel');
        cel = mt_readcel(celpaths);
        eval(sprintf('save %s_cel cel;',ExpName))

    case 'recover chip info'
%% RECOVER CHIP INFO        
        if nargin<6
            h=errordlg('affy_rdn(''info'') needs at least five parameters : Action, ExpName,CelDir, ChipName, MyChipName, ChipDir (LibDir)');
            waitfor(h)
        elseif nargin==6
            %default LibDir
            LibDir=fullfile(ChipDir,'libfiles');
        else
            %force Dir to be terminated by a separator
            if ChipDir(end)~=filesep
                ChipDir=[ChipDir,filesep];
            end
            if LibDir(end)~=filesep
                LibDir=[LibDir,filesep];
            end
        end

        %IMPORT AFFY INFORMATION FILES
        cd(LibDir)
        try
            eval(sprintf('cdf = mt_readcdf(''%s.CDF'');',ChipName))
        catch
            eval(sprintf('cdf = mt_readcdf(''%s.cdf'');',ChipName))
        end
        cd(ChipDir)
        eval(sprintf('save %s_cdf cdf;',MyChipName))

        cd(LibDir)
        try
            eval(sprintf('gin = mt_readgin(''%s.GIN'');',ChipName))
        catch
            eval(sprintf('gin = mt_readgin(''%s.gin'');',ChipName))
        end
        cd(ChipDir)
        eval(sprintf('save %s_gin gin;',MyChipName))


        cd(LibDir)
        eval(sprintf('pa = mt_readprobe_annot(''%s.probe_tab'');',ChipName))
        cd(ChipDir)
        eval(sprintf('save %s_pa pa;',MyChipName))


        cd(LibDir)
        eval(sprintf('seq = mt_readseq(''%s.target'');',ChipName))
        cd(ChipDir)
        eval(sprintf('save %s_seq seq;',MyChipName))

    case 'do rdn analysis'
%% DO RDN ANALYSIS
        if nargin<6
            h=errordlg('affy_rdn(''rdn'') needs six parameters : Action, ExpName,CelDir, ChipName, MyChipName,ChipDir');
            waitfor(h)
        else
            %force Dir to be terminated by a separator
            if ChipDir(end)~=filesep
                ChipDir=[ChipDir,filesep];
            end
        end

        cd(ChipDir)
        eval(sprintf('load %s_cdf;',MyChipName))
        eval(sprintf('load %s_gin;',MyChipName))
        eval(sprintf('load %s_pa;',MyChipName))
        eval(sprintf('load %s_seq;',MyChipName))
        cd(CelDir)
        %make probes structure
        %the initial state of probes is saved in order to start again rdn
        %if something fails in mt_normalize
        if ~exist(sprintf('%s_probes_start probes;',ExpName),'file')
            eval(sprintf('load %s_cel',ExpName))
            probes = mt_cel2probes(cel,cdf,gin,pa,seq);
            cd(CelDir)
            eval(sprintf('save %s_probes_start probes;',ExpName))
        else
            eval(sprintf('load %s_probes_start;',ExpName))
        end
        %rdn algorithm
        probes= mt_normalize(probes);
        cd(CelDir)
        eval(sprintf('save %s_probes_end probes;',ExpName))

    case 'print rdn signals'
%% PRINT RDN SIGNALS
        if nargin<3
            h=errordlg('affy_rdn(''print'') needs three parameters : Action, ExpName,CelDir');
            waitfor(h)
        else
            %force Dir to be terminated by a separator
            if CelDir(end)~=filesep
                CelDir=[CelDir,filesep];
            end
            if nargin<8
                SampleNames={};
                PrintOrder=[];
            elseif nargin <9
                PrintOrder=[];
            end
        end
        cd(CelDir)
        eval(sprintf('load %s_probes_end;',ExpName))
        OutFid=fopen(sprintf('%s_rdn.txt',ExpName),'w');

        %print header
        %adopt BioConductor format (only Sample Names separated by tabulation, no tabulation in first position)
        %take in account PrintOrder
        if isempty(PrintOrder)
            PrintOrder=1:length(probes.array_filenames);
        else
            %verify that PrintOrder is correct
            if length(PrintOrder)~=length(probes.array_filenames)
                h=errordlg('PrintOrder does not have the same size that results');
                waitfor(h)
                error('process canceled')
            end
            if length(unique(PrintOrder))~=length(probes.array_filenames)
                h=errordlg('PrintOrder has some index repeated');
                waitfor(h)
                error('process canceled')
            end
            if max(unique(PrintOrder))~=length(probes.array_filenames)
                h=errordlg('PrintOrder has some index greater than the results size');
                waitfor(h)
                error('process canceled')
            end
        end
        if isempty(SampleNames)            
            PrintNames=probes.array_filenames(PrintOrder);
        else
            %verify that SampleNames is correct
            if length(SampleNames)~=length(probes.array_filenames)
                h=errordlg('SampleNames does not have the same size that results');
                waitfor(h)
                error('process canceled')
            end
            PrintNames=SampleNames;
        end
        
        %calculate log2(RDN signal)
        Log2Signal = repmat(probes.overall_factors, size(probes.array_factors,1) ,1) + probes.array_factors;

        %write column names
        fprintf(OutFid,'%s',PrintNames{1});
        if length(PrintNames)>1
            for ChipL=2:length(PrintNames)
                fprintf(OutFid,'\t%s',PrintNames{ChipL});
            end
        end
        fprintf(OutFid,'\n');
                
        %print log2(RDN signal)                    
        for PsL=1:length(probes.name)
            fprintf(OutFid,'%s',probes.name{PsL});
            for ChipL=1:length(PrintNames)
                fprintf(OutFid,'\t%f',Log2Signal(PrintOrder(ChipL),PsL));
            end
            fprintf(OutFid,'\n');
        end
        
        fclose(OutFid);
end


