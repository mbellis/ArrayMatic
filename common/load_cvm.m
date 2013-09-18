%==================%
% FUNCTION LOAD_CVM %
%==================%

% LOAD_CVM load a part of a CVM

%INPUT PARAMETERS
% 1  ChipRank: chip rank
% 2   NetRank: network (CVM) rank
% 3 RowPsRank: lines of probe set ranks to be loaded
% 3 ColPsRank: columns of probe set ranks to be loaded
% 4  CorrFlag:load CORR
% 5  AntiFlag: load ANTI
% 6  DiffFlag: make difference C-A

%C=load_cvm(94,4,1:200,1:200,1,1,1)


function [C,vargout]=load_cvm(ChipRank,NetRank,RowPsRank,ColPsRank,CorrFlag,AntiFlag,DiffFlag)
global K

if DiffFlag==1&(CorrFlag==0|AntiFlag==0)
    h=errordlg('DiffFlag=1 needs t\nthat CorrFlag and AntiFlag\nequal 1');
    waitfor(h)
    error('process canceled')
end
C=[];
A=[];
ChipPos=find(K.chip.rank==ChipRank);
PsNb=K.chip.probesetNb(ChipPos);
cd(fullfile(K.dir.net,sprintf('m%u',ChipRank),sprintf('n%u',NetRank)))
if CorrFlag
    C=load_data(sprintf('c_m%u_n%u.4mat',ChipRank,NetRank),'./',PsNb,PsNb,'uint8','ieee-le',RowPsRank,ColPsRank);
end
if AntiFlag
    A=load_data(sprintf('a_m%u_n%u.4mat',ChipRank,NetRank),'./',PsNb,PsNb,'uint8','ieee-le',RowPsRank,ColPsRank);
    if nargout==2
        vargout{1}=A;
    end
end
if DiffFlag
    c=triu(C)-triu(A);
    c(find(c<0))=0;
    A=tril(A,-1)-tril(C,-1);
    A(find(A<0))=0;
    C=c+A;
end

