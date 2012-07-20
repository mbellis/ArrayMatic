% SELECT_CHISET - Allows to select a chip set and eventually a chip if the
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

ChipRank=0;
ChipPos=0;
Type='';
ProbeSetNb=0;
Gpl='';
CompName='';
Chromosomes={};
Success=0;

%order chip set by species
[Temp,SortIndex]=sort(K.chipSet.name);
[Temp,SortIndex1]=sort(K.chipSet.species(SortIndex));
SortIndex=SortIndex(SortIndex1);
ChipSet.rank=[0;K.chipSet.rank(SortIndex)];
ChipSet.myName=['not listed';K.chipSet.myName(SortIndex)];
ChipSet.name=['not listed';K.chipSet.name(SortIndex)];
ChipSet.shortName=['not listed';K.chipSet.shortName(SortIndex)];
ChipSet.species=['not listed';K.chipSet.species(SortIndex)];
ChipSet.probesetNb=[0;K.chipSet.probesetNb(SortIndex)];
ChipSet.probeNb=[0;K.chipSet.probeNb(SortIndex)];
ChipSet.compName=[0;K.chipSet.compName(SortIndex)];
ChipSet.geoName=[0;K.chipSet.geoName(SortIndex)];


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

