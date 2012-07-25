% SELECT_CHIP - Allows to select a chip set and eventually a chip if the
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

function [ChipRank,ChipPos,Type,ProbeSetNb,Gpl,CompName,Chromosomes,Success] =select_chip()
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
[Temp,SortIndex]=sort(K.chip.name);
[Temp,SortIndex1]=sort(K.chip.species(SortIndex));
SortIndex=SortIndex(SortIndex1);
Chip.rank=[0;K.chip.rank(SortIndex)];
Chip.myName=['not listed';K.chip.myName(SortIndex)];
Chip.name=['not listed';K.chip.name(SortIndex)];
Chip.type=['not listed';K.chip.type(SortIndex)];
Chip.shortName=['not listed';K.chip.shortName(SortIndex)];
Chip.species=['not listed';K.chip.species(SortIndex)];
Chip.probesetNb=[0;K.chip.probesetNb(SortIndex)];
Chip.probeNb=[0;K.chip.probeNb(SortIndex)];
Chip.compName=[0;K.chip.compName(SortIndex)];
Chip.geoName=[0;K.chip.geoName(SortIndex)];
ChipName=cell(length(Chip.name),1);
for ChipL=1:length(Chip.name)
    ChipName{ChipL}=sprintf('%s - %s - %s',[Chip.geoName{ChipL},repmat('_',1,10-length(Chip.geoName{ChipL}))],[Chip.myName{ChipL},repmat('_',1,4-length(Chip.myName{ChipL}))],Chip.name{ChipL});
end


Continue=1;
while Continue
    [ChipSel,ChipOK] = listdlg('ListString',ChipName,'SelectionMode','single','ListSize',[600,600],'Name','Analysis of Raw Data','PromptString','Select Chip Model');
    if ChipOK        
        Continue=0;
        ChipRank=Chip.rank(ChipSel);        
        ChipPos=find(K.chip.rank==ChipRank);
        Type=K.chip.type(ChipPos);
        ProbeSetNb=K.chip.probesetNb(ChipPos);
        CompName=K.chip.compName{ChipPos};        
        Gpl=K.chip.geoName{ChipPos};
        Success=1;
    else
        Answer=questdlg('Do you want to cancel CHIP SET SELECTION ?','CHIP SET SELECTION','Yes','No','No');
        if isequal(Answer,'Yes')
            Continue=0;
        end
    end
end

