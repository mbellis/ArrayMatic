% SELECT_CHIPSET - Allows to select a chip set and eventually a chip if the
% chip set is composed of several chips

% c) Michel Bellis
% arraymathics@gmail.com

% INPUT
% ChipFlag: (1/0) indicates if chip selection must be done.

% OUTPUT
% ChipSetRank: the rank of the chip set in the chip list
% ChipRank: the rank of the chip in the chip set
% Type: type of chip set {T: transcriptomic, G: genomic)
% ProbeSetNb: the number of probe sets (or the number of probes if probes
%             are not grouped in probe sets)
% CompName: the name of the compagny {affy, nimb, abu, other}
% Chromosomes: a list of chromosomes if probes are ordered on the chip
% (e.g. genomic chips from NimbleGene)
% Success: (1/0) success of the selection

% VERSIONS
%
% V01 - 2010 05 05 - First version

function [ChipSetRank,ChipRank,Type,ProbeSetNb,CompName,Chromosomes,Success] =select_chipset(ChipFlag)
global K



ProbeSetNb=0;
CompName='';
Chromosomes={};

%order chip set by species
[Temp,SortIndex]=sort(K.chipSet.species);
ChipSet.rank=[0;K.chipSet.rank(SortIndex)];
ChipSet.myName=['not listed';K.chipSet.myName(SortIndex)];
ChipSet.name=['not listed';K.chipSet.name(SortIndex)];
ChipSet.shortName=['not listed';K.chipSet.shortName(SortIndex)];
ChipSet.species=['not listed';K.chipSet.species(SortIndex)];
ChipSet.probesetNb=['UserType';K.chipSet.probesetNb(SortIndex)];
ChipSet.probeNb=['UserType';K.chipSet.probeNb(SortIndex)];
ChipSet.compName=['UserType';K.chipSet.compName(SortIndex)];
ChipSet.geoName=['UserType';K.chipSet.geoName(SortIndex)];
%ChipSet.gpls=['UserType';K.chipSet.gpls(SortIndex)];

Success=0;
Continue=1;
while Continue
    [ChipSel,ChipOK] = listdlg('ListString',ChipSet.name,'SelectionMode','single','ListSize',[600,600],'Name','Analysis of Raw Data','PromptString','Select Chip Model');
    if ChipOK
        Continue=0;
        Type=ChipSet.type(ChipSel);
        ProbeSetNbs=ChipSet.probesetNb{ChipSel};
        CompName=ChipSet.compName{ChipSel};
        Chips=ChipSet.chips{ChipSel};
        ChipSetRank=ChipSet.rank(ChipSel);
        %Gpls=ChipSet.gpls{ChipSel};
        if ChipFlag==0
            ChipRank=1;
            ProbeSetNb=ProbeSetNbs(1);
            %Gpl=Gpls(1);
            Success=1;
        else
            if length(Chips)==1              
                ChipRank=1;
                ProbeSetNb=ProbeSetNbs(1);
                %Gpl=Gpls(1);
                Success=1;
            else
                ChipList=cell(length(Chips),1);
                for i=1:length(Chips)
                    ChipList{i}=sprintf('Chip %c',Chips(i));
                end
                Continue=1;                
                while Continue
                    [ChipRank,OK] = listdlg('ListString',ChipList,'SelectionMode','single','ListSize',[400,100],'Name','CHIP SET SELECTION','PromptString','Select the chip to be analyzed (if you have several chip, create a project for each)');
                    if OK
                        Continue=0;
                        ProbeSetNb=ProbeSetNbs(ChipRank);                        
                        %Gpl=Gpls(ChipRank);
                        Success=1;
                    else
                        Answer=questdlg('Do you want to cancel chip selection ?','CHIP SET SELECTION','Yes','No','No');
                        if isequal(Answer,'Yes')
                            Continue=0;                            
                        end
                    end
                end
            end
        end
    else
        Answer=questdlg('Do you want to cancel CHIP SET SELECTION ?','CHIP SET SELECTION','Yes','No','No');
        if isequal(Answer,'Yes')
            Continue=0;
        end
    end
end

