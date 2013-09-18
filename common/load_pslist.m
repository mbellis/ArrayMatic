%======================%
% FUNCTION LOAD_PSLIST %
%======================%


% LOAD_PSLIST loads a probe set list

%INPUT PARAMETERS

% 1 ChipRank : chip rank
% 2 ListRank : list rank

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



function PsList=load_pslist(ChipRank,ListRank)
global K

cd(fullfile(K.dir.net,sprintf('m%u',ChipRank),'list'))
fid=fopen(sprintf('m%u_pslist%u.u32',ChipRank,ListRank),'r');
fseek(fid,0,'eof');
ByteNb=ftell(fid);
fclose(fid);
PsNb=ByteNb/4;
PsList=load_data(sprintf('m%u_pslist%u.u32',ChipRank,ListRank),'./',PsNb,1,'uint32','ieee-le');