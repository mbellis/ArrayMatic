%========================%
% FUNCTION TF_GENEMATRIX %
%========================%

% TF_GENEMATRIX recover localisation of TF in a list of genes (need ModuleMaster to generate
% intermediate results)

%INPUT PARAMETERS
% 1    ChipRank: chip rank
% 2     NtsFile: file containing NetsTensor results
% 3  TargetFile: file describing the merging process used to construct the chip (with rank
%                equal to ChipRank) from the first chip indicated in TargetFile
% 4     FromTsn: the first NtsFile used if the process is repeated across several chip with
%                successive mergin schemes ( human + mouse and then + rat for exemple)
% 5    InfoFlag: recover gene names and gene ids if equal to one

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

%tf_genematrix({[2,3],[5,8,27],[6]},{'human','mouse','rat'},1)
%tf_genematrix({[2,3],[5,8,27],[6]},{'human','mouse','rat'},0)
%tf_genematrix({[6]},{'rat'},0)
%tf_genematrix({[5,8,27]},{'mouse'},0)


function tf_genematrix(ChipRanks,Species,ModuleFlag)
global K

if ModuleFlag==0
    cd(K.dir.common)
    load tflist
end

for SpeciesL=1:length(Species)
    ChipRank=ChipRanks{SpeciesL};
    GeneId={};
    for ChipL=1:length(ChipRank)       
        cd(K.dir.chip)        
        [PsId,CurrGeneId,GeneName,Class,PNb]=textread(sprintf('m%u_gene.txt',ChipRank(ChipL)),'%s%s%s%u%u','delimiter','\t');
        GeneId=unique([GeneId;CurrGeneId]);
        EnsPos=strmatch('ENS',GeneId);
        GeneId=GeneId(EnsPos);
    end
    BlocNb=ceil(length(GeneId)/100);
    if ModuleFlag                
        cd(fullfile(K.dir.wpm,Species{SpeciesL}))
        for BlocL=1:BlocNb
            fid=fopen(sprintf('Cluster%u.txt',BlocL),'w');
            for GeneL=(BlocL-1)*100+1:min(BlocL*100,length(GeneId))
                fprintf(fid,'%s\n',GeneId{GeneL});
            end
            fclose(fid)
        end
    else
        cd(fullfile(K.dir.wpm,Species{SpeciesL},'mm'))
        EnsId=[];        
        TfRank=[];
        TfStrand=[];
        TfStart=[];
        TfEnd=[];
        TfMatchScore=[];
        TfCoreScore=[];
        TfPercMatch=[];
        TfLogWx=[];
        NotFound=[];
        BlocCount=1;
        for BlocL=101:144                            
            fid=fopen(sprintf('Cluster%u.txt.MatrixScan.txt',BlocL),'r');
            if fid>0
                BlocL
                %count header line nb
                HeaderNb=-1;
                while 1
                    CurrLine=fgetl(fid);
                    HeaderNb=HeaderNb+1;
                    if ~ischar(CurrLine)
                        break
                    end
                    if ~isequal(CurrLine(1),';')&~isequal(CurrLine(1),'#')&isempty(findstr(CurrLine,'limit'))
                        break
                    end
                end
                fclose(fid);                
                [CurrEnsId,temp,temp,Name,CurrStrand,Start,End,Sequence,MatchScore,CoreScore,PercMatch,LogWx,temp]=textread(sprintf('Cluster%u.txt.MatrixScan.txt',BlocL),...
                    '%s%s%s%s%c%d%d%s%.3f%.3f%.3f%.3f%s','delimiter','\t','headerlines',HeaderNb);
                %recover ensembl rank from ensembl id
                Prefix=regexp(CurrEnsId{1},'[A-Z]*','match');
                CurrEnsId=strrep(CurrEnsId,Prefix,'');
                try
                    CurrEnsId=str2num(cell2mat(CurrEnsId));                                
                catch                    
                    fprintf('individual correction of EnsId')
                    for EnsL=1:length(CurrEnsId)
                        CurrEnsId{EnsL}=CurrEnsId{EnsL}(1:11);
                    end                    
                    CurrEnsId=str2num(cell2mat(CurrEnsId));                                
                end
                EnsId=[EnsId;CurrEnsId];                
                CurrTfRank=zeros(size(Name));
                UName=unique(Name);
                for TfL=1:length(UName)
                    NamePos=strmatch(UName{TfL},Name,'exact');
                    TfPos=strmatch(UName{TfL},upper(Tf.name),'exact');           
                    CurrTfRank(NamePos)=TfPos(1);
                end
                TfRank=[TfRank;uint16(CurrTfRank)];
                Strand=zeros(size(CurrStrand));
                Strand(find(CurrStrand=='D'))=1;
                TfStrand=[TfStrand;uint8(Strand)];                
                TfStart=[TfStart;uint16(abs(Start))];
                TfEnd=[TfEnd;uint16(abs(End))];
                TfMatchScore=[TfMatchScore;uint8(round(MatchScore*100))];
                TfCoreScore=[TfCoreScore;uint8(round(CoreScore*100))];
                TfPercMatch=[TfPercMatch;uint8(round(PercMatch*100))];
                TfLogWx=[TfLogWx;uint8(LogWx)];
            else
                NotFound(end+1,1)=BlocL;
            end
            fid=fopen('bloc.txt','w');
            fprintf(fid,'%u\n',BlocL);
            fclose(fid)
            if mod(BlocL,100)==0|BlocL==144
                BlocCount=BlocCount+1;
                cd(K.dir.wpm)
                eval(sprintf('save %s_tfpos_%u EnsId TfRank TfStrand TfStart TfEnd TfMatchScore TfCoreScore TfPercMatch TfLogWx NotFound',Species{SpeciesL},BlocCount))
                cd(fullfile(K.dir.wpm,Species{SpeciesL},'mm'))
                if BlocL<BlocNb
                    EnsId=[];
                    TfRank=[];
                    TfStrand=[];
                    TfStart=[];
                    TfEnd=[];
                    TfMatchScore=[];
                    TfCoreScore=[];
                    TfPercMatch=[];
                    TfLogWx=[];
                    NotFound=[];
                end
            end
        end
    end
end


