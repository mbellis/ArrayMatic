% SELECT_CHIPSET - Allows to select a chip set and eventually a chip if the
% chip set is composed of several chips

% c) Michel Bellis
% arraymathics@gmail.com

% INPUT


% OUTPUT

% ChipRank: the rank of the chip
% Type: type of chip set {T: transcriptomic, G: genomic)
% ProbeSetNb: the number of probe sets (or the number of probes if probes
%             are not grouped in probe sets)
% Gpl: GEO platform name (GPLxxxx)
% CompName: the name of the compagny {affy, nimb, abu, other}
% Chromosomes: a list of chromosomes if probes are ordered on the chip
% (e.g. genomic chips from NimbleGene)
% Success: (1/0) success of the selection

% VERSIONS
%
% V01 - 2010 05 05 - First version

function [ChipRank,ChipPos,Type,ProbeSetNb,Gpl,CompName,Chromosomes,Success] =select_chipset()
global K


ProbeSetNb=0;
CompName='';
Chromosomes={};

%order chip set by species
[Temp,SortIndex]=sort(K.chipSet.name);
[Temp,SortIndex1]=sort(K.chipSet.species(SortIndex));
SortIndex=SortIndex(SortINdex1);
ChipSet.rank=[0;K.chipSet.rank(SortIndex)];
ChipSet.myName=['not listed';K.chipSet.myName(SortIndex)];
ChipSet.name=['not listed';K.chipSet.name(SortIndex)];
ChipSet.shortName=['not listed';K.chipSet.shortName(SortIndex)];
ChipSet.species=['not listed';K.chipSet.species(SortIndex)];
ChipSet.probesetNb=['UserType';K.chipSet.probesetNb(SortIndex)];
ChipSet.probeNb=['UserType';K.chipSet.probeNb(SortIndex)];
ChipSet.compName=['UserType';K.chipSet.compName(SortIndex)];
ChipSet.geoName=['UserType';K.chipSet.geoName(SortIndex)];


Success=0;
Continue=1;
while Continue
    [ChipSel,ChipOK] = listdlg('ListString',ChipSet.name,'SelectionMode','single','ListSize',[600,600],'Name','Analysis of Raw Data','PromptString','Select Chip Model');
    if ChipOK
        ChipPos=ChipSel-1;
        Continue=0;
        Type=ChipSet.type(ChipSel);
        ProbeSetNbs=ChipSet.probesetNb{ChipSel};
        CompName=ChipSet.compName{ChipSel};
        ChipRank=ChipSet.rank(ChipSel);
        ProbeSetNb=ProbeSetNbs(1);
        Gpl=ChipSet.geoName{ChipSel};
        Success=1;
    else
        Answer=questdlg('Do you want to cancel CHIP SET SELECTION ?','CHIP SET SELECTION','Yes','No','No');
        if isequal(Answer,'Yes')
            Continue=0;
        end
    end
end

