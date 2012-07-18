%==========
% AMC_INIT
%==========

% AMC_INIT - Initiate or reinitiate an arraymatic session

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

%% first initiation
clear all

global C D F K M N O P R S T X Y DataRanks
%load menu window
F.h.amcmenu=findobj('tag','amc_menu');
if ~isempty(F.h.amcmenu)
    delete(F.h.amcmenu(1:end))
end
F.h.amcmenu=openfig('amc_menu');
set(F.h.amcmenu,'visible','on');

set(gcf,'pointer','watch')

F.gh.amcmenu=guihandles(F.h.amcmenu);
guidata(F.h.amcmenu,F.gh.amcmenu)
set(0,'tag','')
%manage TRANSCRIPTOME menu
set(get(F.gh.amcmenu.TRANSCRIPTOME,'children'),'enable','off')
set(get(F.gh.amcmenu.TRANSCRIPTOME,'children'),'checked','off')
TrsChildren=get(F.gh.amcmenu.TRANSCRIPTOME,'children');
for ChildL=1:length(TrsChildren)
    set(get(TrsChildren(ChildL),'children'),'enable','off')
    set(get(TrsChildren(ChildL),'children'),'checked','off')
end
set(F.gh.amcmenu.LOAD_A_PROJECT,'enable','on')
set(F.gh.amcmenu.NEW_PROJECT,'enable','on')
set(F.gh.amcmenu.PROJECT_MANAGEMENT,'enable','on')


%recover screen pixel size
set(0,'unit','pixels');
K.par.pixelscrsz = get(0,'ScreenSize');

%indicates the current state of the application
if ~isfield(K,'flag')
    K.flag.application='initamc';
end

%recover station type
if isunix==1
    K.flag.station='unix';
elseif ispc==1
    K.flag.station='windows';
end

%recover root containing arraymatic files
if isequal(K.flag.station,'windows')
    try
        cd('e:\')
    catch
    end
    [Temp,SosDir]=uigetfile('*','select any file in the disk that contains the application files you want to use');
    K.dir.root=SosDir(1:2);
else
    K.dir.root=fullfile('/home','mbellis');
end

%initiate directories in K.dir
amc_dir

%recover current vertion of Matlab
CurrMatVer=ver;
for VerL=1:length(CurrMatVer)
    if isequal(CurrMatVer(VerL).Name,'MATLAB')
        MatVer=CurrMatVer(VerL).Version;
        break
    end
end
K.par.matver=round(str2num(MatVer));


%Load or construct Chip Set and Network Lists
try
    cd(K.dir.common)
catch
    %modify K.dir
    K.dir.common='E:\sosma\arraymatic\amcdata\common';
    K.dir.test= 'E:\sosma\arraymatic\amcdata\test';
    K.dir.metadata= 'E:\sosma\arraymatic\amcdata\metadata';
    K.dir.geoMetadata= 'E:\sosma\arraymatic\amcdata\metadata\geo';
    K.dir.affyMetadata= 'E:\sosma\arraymatic\amcdata\metadata\affy';
    K.dir.chipMetadata= 'E:\sosma\arraymatic\amcdata\metadata\chip';
    K.dir.cliques= 'E:\sosma\arraymatic\amcdata\cliques';
    K.dir.chip= 'E:\sosma\arraymatic\amcdata\chip';
end
cd(K.dir.common)

if exist('netlist.mat')==2
    load netlist
    K.net=Tempo;
    clear Tempo
else
    amc_netlist
end

cd(K.dir.common)
if exist('chiplist.mat')==2
    load chiplist
    K.chipSet=Tempo;
    clear Tempo
else
    %amc_chiplist
    cd(K.dir.chipMetadata)
    [K.chipSet.rank,K.chipSet.myName,K.chipSet.name,K.chipSet.shortName,K.chipSet.species,K.chipSet.probesetNb,K.chipSet.probeNb,K.chipSet.compName,temp,temp,K.chipSet.geoName]=textread('chip.txt','%u%s%s%s%s%u%u%s%s%s%s','delimiter','\t');
    cd(K.dir.chip)
    [K.chipSet.ref.rank,K.chipSet.ref.signal]=textread('signal_vs_rank.txt','%f%f','delimiter','\t');
    Tempo=K.chipSet;
    cd(K.dir.common)
    save chiplist Tempo
end



%% common
K.flag.application='initamc';
set(0,'unit','normalized');
figure(F.h.amcmenu)
set(gcf,'pointer','arrow')

clear Var RootUnit VerL CurrMatVer MatVer Root SosDir Temp