% CREATE_DIR - Create a diretory or use an existing one

% c) Michel Bellis
% arraymatic@gmail.com

% INPUT
% Message : the message displayed in the dialog box

% OUTPUT
% ParentDir : parent directory
% NewDir : new directory
% Success : [0/1] indicates the success of operation
% ErrMessage : if Success=0 indicates the reason of failure

% VERSIONS
%
% V01 - 2010 05 04- first version

function [NewDir,DirName,Success,ErrMessage]=create_dir(Message)
Input=1;
NewDir='';
while Input
    Success=1;
    ErrMessage='';
    ParentDir=uigetdir('*.*',Message);
    cd(ParentDir)
    DirName=inputdlg({'Give a new directory name (space will be replaced by underscore)'},'',1,{''});
    DirName=strrep(DirName{1},' ','_');
    if exist(sprintf('.%c%s',filesep,DirName),'dir')
        NameChoice=questdlg(sprintf('%s already exists. Do you want to ',DirName),'','use it','use another name','cancel','use another name');
        if isequal(NameChoice,'cancel')
            Success=0;
            ErrMessage='Canceled by user';
            Input=0;
        elseif isequal(NameChoice,'use it')
            Input=0;
            NewDir=fullfile(ParentDir,DirName);
        end
    else        
        [Success,ErrMessage]=mkdir(ParentDir,DirName);
        if Success==0
            Another=questdlg(sprintf('creation of %s failed : %s. Another try ?',DirName,ErrMessage),'','yes','yes','no');
            if isequal(Another,'no')
                Input=0;                
            end
        else
            Input=0;
            NewDir=fullfile(ParentDir,DirName);
        end
    end
end
