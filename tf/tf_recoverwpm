%========================%
% FUNCTION TF_DISPLAYCVM %
%========================%

% TF_DISPLAYCVM displays for a given cluster of genes the CVM of genes and all transcription
% factors found in the promotor region of these genes

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

%tf_displaycvm(3,86,'m3_tsn_from_m93n1ton3_corr60',40,'m93n1ton3',[1])cd(K.dir.chip)
[PsId,GeneId,GeneName,Class,PNb]=textread('m2_gene.txt','%s%s%s%u%u','delimiter','\t');
GeneId1=unique(GeneId);
[PsId,GeneId,GeneName,Class,PNb]=textread('m3_gene.txt','%s%s%s%u%u','delimiter','\t');
GeneId=unique([GeneId1;GeneId]);
EnsPos=strmatch('ENS',GeneId);
GeneId=GeneId(EnsPos);

cd(fullfile(K.dir.wpm,'human'))
BlocNb=ceil(length(GeneId)/100);
for BlocL=1:BlocNb
    fid=fopen(sprintf('Cluster%u.txt',BlocL),'w');
    for GeneL=(BlocL-1)*100+1:min(BlocL*100,length(GeneId))
        fprintf(fid,'%s\n',GeneId{GeneL});
    end
    fclose(fid)
end

cd(K.dir.chip)
[PsId,GeneId,GeneName,Class,PNb]=textread('m5_gene.txt','%s%s%s%u%u','delimiter','\t');%GeneName1=unique(GeneName);
GeneId1=unique(GeneId);
[PsId,GeneId,GeneName,Class,PNb]=textread('m8_gene.txt','%s%s%s%u%u','delimiter','\t');
GeneId1=unique([GeneId1;GeneId]);
[PsId,GeneId,GeneName,Class,PNb]=textread('m27_gene.txt','%s%s%s%u%u','delimiter','\t');
GeneId=unique([GeneId1;GeneId]);
EnsPos=strmatch('ENS',GeneId);
GeneId=GeneId(EnsPos);

cd(fullfile(K.dir.wpm,'mouse'))
BlocNb=ceil(length(GeneId)/100);
for BlocL=1:BlocNb
    fid=fopen(sprintf('Cluster%u.txt',BlocL),'w');
    for GeneL=(BlocL-1)*100+1:min(BlocL*100,length(GeneId))
        fprintf(fid,'%s\n',GeneId{GeneL});
    end
    fclose(fid)
end

cd(K.dir.chip)
[PsId,GeneId,GeneName,Class,PNb]=textread('m6_gene.txt','%s%s%s%u%u','delimiter','\t');
GeneId=unique(GeneId);
EnsPos=strmatch('ENS',GeneId);
GeneId=GeneId(EnsPos);

cd(fullfile(K.dir.wpm,'rat'))
BlocNb=ceil(length(GeneId)/100);
for BlocL=1:BlocNb
    fid=fopen(sprintf('Cluster%u.txt',BlocL),'w');
    for GeneL=(BlocL-1)*100+1:min(BlocL*100,length(GeneId))
        fprintf(fid,'%s\n',GeneId{GeneL});
    end
    fclose(fid)
end

