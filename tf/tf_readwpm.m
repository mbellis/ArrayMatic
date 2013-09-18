%=====================%
% FUNCTION TF_READWPM %
%=====================%

% TF_READWPM reads the Weight Position Matrix of cell paper :
% 'DNA-Binding Specificities of Human Transcription Factors'
% A. Jolma et al. (2013) 152: 327


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
function tf_readwpm()
global K

SpeciesName={'human Homo sapiens','mouse Mus musculus','rat Rattus norvegicus'};
cd(K.dir.wpm)
%read wpm text file
fid=fopen('Cell_2013_152_327_Jolma.txt','r');
fod=fopen('wpm_Jolma.txt','w');
Tf={};
Tf.name={};
Tf.gene={};
Tf.family={};
Tf.type={};
Tf.species={};
Tf.wpm={};
Tf.psRank={};


cd(K.dir.gene)
Format=repmat('%s',1,3);
Output='[';
for i=1:3
    Output=[Output,sprintf('EnsemblIds{%u} ',i)];
end
Output=[Output,']'];
eval(sprintf('%s=textread(''hmr_ensembl.txt'',Format,''delimiter'',''\\t'');',Output));

[EnsId{1},Desc{1},GeneName{1}]=textread('human_ensembl.txt','%s%s%s','delimiter','\t');
[EnsId{2},Desc{2},GeneName{2}]=textread('mouse_ensembl.txt','%s%s%s','delimiter','\t');
[EnsId{3},Desc{3},GeneName{3}]=textread('rat_ensembl.txt','%s%s%s','delimiter','\t');


while 1
    CurrLine=fgetl(fid);
    if ~ischar(CurrLine)
        break
    end
    CurrTabPos=regexp(CurrLine,'\t');
    if ~isequal(CurrLine(1:CurrTabPos(1)-1),'symbol')&~isequal(CurrLine(1:CurrTabPos(1)-1),'1')
        Line{1}=CurrLine;
        TabPos=CurrTabPos;
        Line{2}=fgetl(fid);
        CurrTabPos=regexp(Line{2},'\t');
        if isequal(Line{2}(1:CurrTabPos(1)-1),'A')
            %load the next three lines
            for LineL=1:3
                Line{LineL+2}=fgetl(fid);
            end
            %export matrix
            fprintf(fod,'//\n');
            Tf.gene{end+1,1}=Line{1}(1:TabPos(1)-1);
            GenePos=strmatch(Tf.gene{end},Tf.gene,'exact');
            Tf.name{end+1,1}=sprintf('%s-%02u',Tf.gene{end},length(GenePos));
            fprintf(fod,'NA  %s\n',Tf.name{end});
            fprintf(fod,'XX\n');
            Tf.family{end+1,1}=Line{1}(TabPos(1)+1:TabPos(2)-1);
            fprintf(fod,'HC  %s\n',Tf.family{end});
            fprintf(fod,'XX\n');
            if double(Line{1}(2))<91
                Tf.species{end+1,1}='human';
                Species=1;
            else
                Tf.species{end+1,1}='mouse';
                Species=2;                
            end
            GenePos=strmatch(upper(Tf.gene{end}),upper(GeneName{Species}),'exact');                
            if ~isempty(GenePos)
                for SpeciesL=1:3
                    BF{SpeciesL}={};
                end
                for GeneL=1:length(GenePos)
                    EnsPos=strmatch(EnsId{Species}{GenePos(GeneL)},EnsemblIds{Species},'exact');
                    for EnsL=1:length(EnsPos)
                        for SpeciesL=1:3                            
                            BF{SpeciesL}{end+1,1}=EnsemblIds{SpeciesL}{EnsPos(EnsL)};                           
                            EnsPos1=strmatch(BF{SpeciesL}{end},EnsId{SpeciesL},'exact');                            
                            for EnsL1=1:length(EnsPos1)
                                BF{SpeciesL}{end+1,1}=GeneName{SpeciesL}{EnsPos1(EnsL1)};
                            end
                        end
                    end
                end
                for SpeciesL=1:3
                    BF{SpeciesL}=unique(BF{SpeciesL});
                    if ~isempty(BF{SpeciesL})
                        fprintf(fod,'BF EnsemblSynonyme;');
                        if length(BF{Species})>1
                            for NameL=1:length(BF{SpeciesL})-1
                                fprintf(fod, ' %s,',BF{SpeciesL}{NameL});
                            end
                        end
                        fprintf(fod, ' %s; species: %s.\n',BF{SpeciesL}{end},SpeciesName{SpeciesL});
                    end
                end
                fprintf(fod,'XX\n');
            end
            fprintf(fod,'OS  %s\n',Tf.species{end});
            fprintf(fod,'XX\n');
            Tf.type{end+1,1}=Line{1}(TabPos(8)+1:TabPos(9)-1);
            fprintf(fod,'TY  %s\n',Tf.type{end});
            fprintf(fod,'XX\n');
            fprintf(fod,'P0 A C G T\n');
            WPM=[];
            for LineL=1:4
                Val=regexp(Line{LineL+1},'\d*','match');
                for ValL=1:length(Val)
                    WPM(ValL,LineL)=str2num(Val{ValL});
                end
            end
            Tf.wpm{end+1,1}=WPM;
            for PosL=1:size(WPM,1)
                fprintf(fod,'%02u %u %u %u %u\n',PosL,WPM(PosL,1),WPM(PosL,2),WPM(PosL,3),WPM(PosL,4));
            end
            fprintf(fod,'XX\n');
        end
    end
end
fprintf(fod,'//\n');
fclose(fid)
fclose(fod)

%find TF gene in human, mouse and rat chips
cd(K.dir.gene)
GeneFile='hmr_ensembl.txt';
SpeciesNb=3;
Format=repmat('%s',1,SpeciesNb);
Output='[';
for i=1:SpeciesNb
    Output=[Output,sprintf('EnsemblIds{%u} ',i)];
end
Output=[Output,']'];
eval(sprintf('%s=textread(GeneFile,Format,''delimiter'',''\\t'');',Output));

%recover  Human ensemblID
cd(K.dir.gene)
[EnsId,Desc,GeneName]=textread('human_ensembl.txt','%s%s%s','delimiter','\t');
GeneName=upper(GeneName);
UGeneName=unique(Tf.gene);
Tf.unique=UGeneName;
TfEnsId=cell(size(UGeneName));
FoundTf=zeros(size(UGeneName));
for GeneL=1:length(UGeneName)
    EnsPos=strmatch(upper(UGeneName{GeneL}),GeneName,'exact');
    if ~isempty(EnsPos)
        FoundTf(GeneL)=length(EnsPos);    
        for EnsL=1:length(EnsPos)
            TfEnsId{GeneL}{EnsL}=EnsId{EnsPos(EnsL)};
        end   
    else
        TfEnsId{GeneL}='';
    end
end

%find ps rank corresponding to TF
%write txt files usable with ModuleMaster
%use signal first to test
ChipRank=[2,3,5,6,8,27];
for ChipL=1:6
    %recover information on current chip
    cd(K.dir.chip)
    [PsId,EnsId,GeneName,Class,PNb]=textread(sprintf('m%u_gene.txt',ChipRank(ChipL)),'%s%s%s%u%u','delimiter','\t');
    %recover all ps ranks corresponding to Tf
    Tf.psRank{ChipRank(ChipL),1}=cell(length(TfEnsId),1);
    Tf.isOnChip{ChipRank(ChipL),1}=uint8(zeros(length(TfEnsId),1));
    TfPsRank=[];
    TfGeneName={};
    for GeneL=1:length(TfEnsId)
        CurrRank=[];
        for EnsL=1:length(TfEnsId{GeneL})
            if  ChipRank(ChipL)<=3
                EnsPos=strmatch(TfEnsId{GeneL}{EnsL},EnsId,'exact');
                if ~isempty(EnsPos)
                    CurrRank=[CurrRank;EnsPos];
                end
            else
                %find corresponding EnsId in another species
                EnsPos=strmatch(TfEnsId{GeneL}{EnsL},EnsemblIds{1},'exact');
                for EnsL1=1:length(EnsPos)
                    if ChipRank(ChipL)==6
                        EnsPos1=strmatch(EnsemblIds{3}{EnsPos(EnsL1)},EnsId,'exact');
                    else
                        EnsPos1=strmatch(EnsemblIds{2}{EnsPos(EnsL1)},EnsId,'exact');
                    end
                    if ~isempty(EnsPos1)
                        CurrRank=[CurrRank;EnsPos1];
                    end
                end
            end
        end
        Tf.psRank{ChipRank(ChipL)}{GeneL,1}=CurrRank;
        if ~isempty(CurrRank)
            Tf.isOnChip{ChipRank(ChipL)}(GeneL,1)=1;
        end
        for PsL=1:length(CurrRank)
            TfGeneName{end+1,1}=Tf.unique{GeneL};
        end
        TfPsRank=[TfPsRank;CurrRank];
    end
    
    %load DataRanks   
    cd(fullfile(K.dir.tables,sprintf('m%u',ChipRank(ChipL)),sprintf('m%u_data',ChipRank(ChipL))))
    load DataRanks

    % print DataRanks
    cd(K.dir.wpm)
    fid=fopen(sprintf('m%u_rank4wpm.txt',ChipRank(ChipL)),'w');
    fprintf(fid,'probe set id\tgene symbol\tgene name\tensembl\t');
    Format=repmat('%.4f\t',1,size(DataRanks,2)-1);
    Format=[Format,'%.4f\n'];
    for PointL=1:size(DataRanks,2)-1
        fprintf(fid,sprintf('point_%u\t',PointL));
    end
    fprintf(fid,sprintf('point_%u\n',size(DataRanks,2)));
    for PsL=1:min(length(PsId),size(DataRanks,1))
        fprintf(fid,'%s\t',PsId{PsL});
        fprintf(fid,'%s\t',GeneName{PsL});        
        fprintf(fid,'%s\t',EnsId{PsL});
        fprintf(fid,Format,DataRanks(PsL,:));
    end
    fclose(fid)   
end

%make correspondance between Tf.name and Tf.unique
Tf.toUnique=zeros(size(Tf.name));
for GeneL=1:length(Tf.unique)
    GenePos=strmatch(Tf.unique{GeneL},Tf.gene,'exact');
    Tf.toUnique(GenePos)=GeneL;
end

cd(K.dir.common)
save tflist Tf












