%TITLE NET_DESCRIPTION
%generate the information in K.net{ModelRank} of networks corresponding to
% the pairs of biological conditions groups as indicated in GRPList

%INPUT PARAMETERS
%GrpList: Indicates the pairs of biologicla conditions used to construct networks
%         [1,2;7,2;3,8]=> will make three networks, one by comparing P.bioL.grp{1} vs P.biol.grp{2},
%         the second by comparing P.bioL.grp{7} vs P.biol.grp{2}
%         and the third by comparing P.bioL.grp{3} vs P.biol.grp{8}.
%Fdr: fdr value(s)
%Sensitivity: sensitivity value(s)

%EXTERNAL FILES

%OUTPUT PARAMETERS
%net_description([[1:20]',[[2:20],1]'],[0.01,0.10],[1,1])


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

function net_description(GrpList,Fdr,Sensitivity)
global K P 
ModelRank=P.chip.chipRank;
Tempo=K.net;
if length(Tempo)<ModelRank
    Tempo{ModelRank}=[];
end
NetNb=size(GrpList,1);
for NetL=1:NetNb
    FirstGrp=GrpList(NetL,1);
    SndGrp=GrpList(NetL,2);
    if isempty(Tempo{ModelRank})
        Pos=1;
    else
        Pos=length(Tempo{ModelRank}.name)+1;
    end
    if Pos==1
        CurrRank=1;
    else
        CurrRank=setdiff([1:Pos],Tempo{ModelRank}.rank);
        CurrRank=CurrRank(1);
    end    
    Tempo{ModelRank}.name{Pos,1}=sprintf('biolgrp %u vs biolgrp %u',FirstGrp,SndGrp);
    Tempo{ModelRank}.rank(Pos,1)=CurrRank;
    Tempo{ModelRank}.biolRank{Pos,1}{1}=P.biol.grp.biolRanks{FirstGrp};
    Tempo{ModelRank}.biolRank{Pos,1}{2}=P.biol.grp.biolRanks{SndGrp};
    Tempo{ModelRank}.compNb(Pos,1)=length(P.biol.grp.biolRanks{FirstGrp})*length(P.biol.grp.biolRanks{SndGrp});
    Tempo{ModelRank}.fdr(Pos,1)=Fdr(1);
    Tempo{ModelRank}.fdr(Pos,2)=Fdr(2);
    Tempo{ModelRank}.s(Pos,1)=Sensitivity(1);
    Tempo{ModelRank}.s(Pos,2)=Sensitivity(2);
    Tempo{ModelRank}.blocNb(Pos,1)=ceil(P.chip.currProbeSetNb/100);
    Tempo{ModelRank}.blocSize(Pos,1)=100;
    Tempo{ModelRank}.comment{Pos,1}='';
    Tempo{ModelRank}.netMade(Pos,1)=1;
end
Tempo{ModelRank}.nb=length(Tempo{ModelRank}.name);
cd(K.dir.common)
save netlist Tempo
K.net=Tempo;
'stop'


function notused()
Tempo=K.net;
ModelRank=11;
NetPos=44;
NetRank=NetPos;
NetPos1=1;
Tempo{ModelRank}.name{NetPos}='Union of networs 1 to 20';
Tempo{ModelRank}.rank(NetPos)=NetRank;
Tempo{ModelRank}.biolRank{NetPos,1}{1}=[];
Tempo{ModelRank}.biolRank{NetPos,1}{2}=[];
Tempo{ModelRank}.compNb(NetPos,1)=Tempo{ModelRank}.compNb(NetPos1,1);
Tempo{ModelRank}.fdr(NetPos,1)=Tempo{ModelRank}.fdr(NetPos1,1);
Tempo{ModelRank}.fdr(NetPos,2)=Tempo{ModelRank}.fdr(NetPos1,2);
Tempo{ModelRank}.s(NetPos,1)=Tempo{ModelRank}.s(NetPos1,1);
Tempo{ModelRank}.s(NetPos,2)=Tempo{ModelRank}.s(NetPos1,2);
Tempo{ModelRank}.blocNb(NetPos,1)=Tempo{ModelRank}.blocNb(NetPos1,1);
Tempo{ModelRank}.blocSize(NetPos,1)=Tempo{ModelRank}.blocSize(NetPos1,1);
Tempo{ModelRank}.comment{NetPos,1}='';
Tempo{ModelRank}.netMade(NetPos,1)=1;
cd(K.dir.common)
save netlist Tempo
K.net=Tempo;



