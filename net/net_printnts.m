%=======================%
% FUNCTION NET_PRINTNTS %
%=======================%

% NET_PRINTNTS print the gene names or gene ids for each cluster found in nts

%INPUT PARAMETERS
% 1 ChipRank: chip rank
% 2  NtsFile: file containing NetsTensor results
% 3    CluNb: number of clusters to be printed


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

%net_printnts(2,'m2_tsn_from_m94n1ton3_corr60','m94n1ton3',5,30,0)
%net_printnts(3,'m3_tsn_from_m94n1ton3_corr60','m94n1ton3',5,30,0)
%net_printnts(5,'m5_tsn_from_m94n1ton3_corr60','m94n1ton3',5,30,0)
%net_printnts(6,'m6_tsn_from_m94n1ton3_corr60','m94n1ton3',5,30,0)
%net_printnts(8,'m8_tsn_from_m94n1ton3_corr60','m94n1ton3',5,30,0)
%net_printnts(27,'m27_tsn_from_m94n1ton3_corr60','m94n1ton3',5,30,0)

% net_printnts(2,'m2_tsn_from_m93n1ton3_corr60','m93n1ton3',5,40,0)
% net_printnts(3,'m3_tsn_from_m93n1ton3_corr60','m93n1ton3',5,40,0)
% net_printnts(5,'m5_tsn_from_m93n1ton3_corr60','m93n1ton3',5,40,0)
% net_printnts(8,'m8_tsn_from_m93n1ton3_corr60','m93n1ton3',5,40,0)
% net_printnts(27,'m27_tsn_from_m93n1ton3_corr60','m93n1ton3',5,40,0)




function net_printnts(ChipRank,NtsFile,FromTsn,MinGeneNb,DensVal,CluNb)
global K


TensorDir=fullfile(K.dir.net,sprintf('m%u',ChipRank),'tsn','result');
cd(TensorDir)
if exist(sprintf('%s.mat',NtsFile))
    eval(sprintf('load %s',NtsFile))
else
  h=errordlg(sprintf('%s does \nnot exist',NtsFile));
  waitfor(h)
  error('process canceled');  
end

DensPos=find(Densities(UsedDensities)==DensVal);
if ~isempty(DensPos)
    OutputDir=fullfile(K.dir.net,sprintf('m%u',ChipRank(1)),'tsn','result',FromTsn);
    mkdir(OutputDir)
    cd(OutputDir)
    %ModuleMaster directory
    mkdir('mm')
    if CluNb==0
        CluNb=length(NtsGeneNames{DensPos});
    else
        CluNb=min(CluNb,length(NtsGeneNames{DensPos}));
    end
    for CluL=1:CluNb
        if length(NtsGeneNames{DensPos}{CluL})>=MinGeneNb
            fid=fopen(sprintf('Cluster%u.txt',CluL),'w');
            for GeneL=1:length(NtsGeneNames{DensPos}{CluL})
                fprintf(fid,'%s\n',NtsGeneNames{DensPos}{CluL}{GeneL});
            end
            fclose(fid);
        end
    end
end
