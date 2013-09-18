% CHIP2CHIP translates a liste of probe sets in one chip towards another chips
%====================%
% FUNCTION CHIP2CHIP %
%====================%
%
% CHIP2CHIP translates a liste of probe sets in one chip towards another chips
%
%INPUT PARAMETERS
% 1 ChipRank: source and target chip rank
% 2 DirName: source and target directory
% 3 FileName: source and target file
% 4 CommonPsFile: source and target correspondance between the two chips
% 5 ListType: type of probe set list (e.g. name of the variable to be used), either:
%             - 'Clu' organized in Clu{Type}{List}{Network}(:,Limit)
%             - 'PsRank' a simple list of probe set ranks

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


function chip2chip(ChipRank,DirName,FileName,CommonPsFile,ListType)
%chip2chip([5,8], {'fullfile(K.dir.net,''m5'',''mcl'')','fullfile(K.dir.net,''m8'',''mcl'')'}, {'m5_mcl_2004n5','m8_mcl_m5_2004n5'}, {'m5n123_m8n228_combinedps_corr60','m8n228_m5n123_combinedps_corr60'}, 'Clu')
%chip2chip([5,27],{'fullfile(K.dir.net,''m5'',''mcl'')','fullfile(K.dir.net,''m27'',''mcl'')'},{'m5_mcl_2004n5','m27_mcl_m5_2004n5'},{'m5n123_m27n164_combinedps_corr60','m27n164_m5n123_combinedps_corr60'},'Clu')

global K

% PsNb
ChipPos=[find(K.chip.rank==ChipRank(1)),find(K.chip.rank==ChipRank(2))];
PsNb=[K.chip.probesetNb(ChipPos(1)),K.chip.probesetNb(ChipPos(2))];

% load source file
eval(sprintf('cd(%s)',DirName{1}))
eval(sprintf('load %s',FileName{1}))

% construct correspondance files necessary to translate probe set ranks
%first chip ps rank => common ps rank
cd(K.dir.chip)
eval(sprintf('load %s;',CommonPsFile{1}))
% construct correspondance between the first chip probe set rank and the common order of ps
ComPsRank{1}=zeros(PsNb(1),1);
for PsL=1:length(PsRank)
    ComPsRank{1}(PsRank{PsL})=PsL;
end

% common ps rank => second chip ps rank
eval(sprintf('load %s;',CommonPsFile{2}))
% find the largest group of merged probe sets
MaxPsNb=0;
for PsL=1:length(PsRank)
    MaxPsNb=max(MaxPsNb,length(PsRank{PsL}));
end
ComPsRank{2}=zeros(length(PsRank),MaxPsNb);
% construct correspondance between the common order of ps and the second chip
for PsL1=1:length(PsRank)
    for PsL2=1:length(PsRank{PsL1})
        ComPsRank{2}(PsL1,PsL2)=PsRank{PsL1}(PsL2);
    end
end
clear PsRank

%translate
switch ListType
    case 'Clu'
        NewClu={};
        for TypeL=1:length(Clu)
            for ListL=1:length(Clu{TypeL})
                for NetL=1:length(Clu{TypeL}{ListL})
                    NewClu{TypeL}{ListL}{NetL}=zeros(PsNb(2),size(Clu{TypeL}{ListL}{NetL},2));
                    for LimitL=1:size(Clu{TypeL}{ListL}{NetL},2)
                        %first chip ps rank => common ps rank
                        PsPos=find(Clu{TypeL}{ListL}{NetL}(:,LimitL));
                        CurrClu=zeros(size(ComPsRank{2},1),1);
                        PosPos=find(ComPsRank{1}(PsPos));
                        PsPos=PsPos(PosPos);
                        % !!! IF SEVERAL PROBE SETS ARE MERGED, IT IS THE CLUSTER INFORMATION OF THE LAST ONE WHICH IS USED !!!
                        CurrClu(ComPsRank{1}(PsPos))=Clu{TypeL}{ListL}{NetL}(PsPos,LimitL);
                        %common ps rank => second chip ps rank
                        for MaxPsL=1:MaxPsNb
                            PsPos=find(ComPsRank{2}(:,MaxPsL));
                            NewClu{TypeL}{ListL}{NetL}(ComPsRank{2}(PsPos,MaxPsL),LimitL)=CurrClu(PsPos);
                            length(find(NewClu{TypeL}{ListL}{NetL}(:,LimitL)))
                        end
                    end
                end
            end
        end
        Clu=NewClu;
        eval(sprintf('cd(%s)',DirName{2}))
        eval(sprintf('save %s Clu',FileName{2}))
        
    case 'PsRank'
        NewPsRank
        %first chip ps rank => common ps rank
        CurrPsRank=unique(ComPsRank{1}(PsRank));
        %common ps rank => second chip ps rank
        for MaxPsL=1:MaxPsNb
            PsPos=find(ComPsRank{2}(:,MaxPsL));
            NewPsRank(ComPsRank{2}(PsPos,MaxPsL))=CurrPsRank(PsPos);
        end
        PsRank=NewPsRank;        
        eval(sprintf('cd(%s)',DirName{2}))
        eval(sprintf('save %s PsRank',FileName{2}))
end