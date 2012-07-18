function [varargout]=loadsave_comp(Action,ResRank,Columns,Dir,ZVarFlag,FdrFlag,SensitivityFlag,PvFlag,FcFlag,varargin)

switch Action
    case 'load'
        if nargout~=ZVarFlag+FdrFlag+SensitivityFlag+PvFlag+FcFlag+1
            h=errordlg('nargout not correct');
            waitfor(h)
            error('process canceled')
        end
        if nargin~=12
            h=errordlg('nargin not correct');
            waitfor(h)
            error('process canceled')
        end

        ProbeSetNb=varargin{1};
        PointNb=varargin{2};
        EmptyFlag=varargin{3};
        cd(Dir)
        VarRank=0;


        if ZVarFlag
            Values=DATALOAD('ZVar',Dir,ResRank,Columns,ProbeSetNb,PointNb,EmptyFlag);
            VarRank=VarRank+1;
            varargout{VarRank}=Values;
        end
        if FdrFlag
            Values=DATALOAD('Fdr',Dir,ResRank,Columns,ProbeSetNb,PointNb,EmptyFlag);
            VarRank=VarRank+1;
            varargout{VarRank}=Values;
        end
        if SensitivityFlag
            Values=DATALOAD('Sensitivity',Dir,ResRank,Columns,ProbeSetNb,PointNb,EmptyFlag);
            VarRank=VarRank+1;
            varargout{VarRank}=Values;
        end
        if PvFlag
            Values=DATALOAD('Pv',Dir,ResRank,Columns,ProbeSetNb,PointNb,EmptyFlag);
            VarRank=VarRank+1;
            varargout{VarRank}=Values;
        end
        if FcFlag
            Values=DATALOAD('Fc',Dir,ResRank,Columns,ProbeSetNb,PointNb,EmptyFlag);
            VarRank=VarRank+1;
            varargout{VarRank}=Values;
        end
        VarRank=VarRank+1;
        load(sprintf('TotalVar_%02u.mat',ResRank))
        varargout{VarRank}=TotalVar;


    case 'save'

        if nargin~=ZVarFlag+FdrFlag+SensitivityFlag+PvFlag+FcFlag+11
            h=errordlg('nargin not correct');
            waitfor(h)
            error('process canceled')
        end
        ProbeSetNb=varargin{1};
        VarRank=1;        
        if ZVarFlag
            VarRank=VarRank+1;
            Values=varargin{VarRank};
            DATASAVE('ZVar',Dir,ResRank,Values,Columns,ProbeSetNb)            
        end
        if FdrFlag
            VarRank=VarRank+1;
            Values=varargin{VarRank};
            DATASAVE('Fdr',Dir,ResRank,Values,Columns,ProbeSetNb)
        end
        if SensitivityFlag
            VarRank=VarRank+1;
            Values=varargin{VarRank};
            DATASAVE('Sensitivity',Dir,ResRank,Values,Columns,ProbeSetNb)
        end
        if PvFlag
            VarRank=VarRank+1;
            Values=varargin{VarRank};
            DATASAVE('Pv',Dir,ResRank,Values,Columns,ProbeSetNb)        
        end
        if FcFlag
            VarRank=VarRank+1;
            Values=varargin{VarRank};
            DATASAVE('Fc',Dir,ResRank,Values,Columns,ProbeSetNb)        
        end
        VarRank=VarRank+1;
        Values=varargin{VarRank};
        cd(Dir)
        if exist(sprintf('TotalVar_%02u.mat',ResRank),'file')
            load(sprintf('TotalVar_%02u.mat',ResRank))
        else
            TotalVar.inc=[];
            TotalVar.dec=[];
        end
        if length(Columns)==1&Columns==0
            TotalVar.inc=Values(1);
            TotalVar.dec=Values(2);
        else
            TotalVar.inc(Columns,1)=Values(1);
            TotalVar.dec(Columns,1)=Values(2);
        end
        eval(sprintf('save TotalVar_%02u TotalVar',ResRank))        
end

function Values=DATALOAD(Type,Dir,ResRank,Columns,ProbeSetNb,PointNb,EmptyFlag)
File=sprintf('%s_%02u.float32le',Type,ResRank);
cd(Dir)
if exist(File,'file')
    if Columns==0
        Values=load_data (File,Dir,ProbeSetNb,PointNb,'single','ieee-le');
    else
        Values=load_data (File,Dir,ProbeSetNb,PointNb,'single','ieee-le',1:ProbeSetNb,Columns);
    end
else
    if EmptyFlag
        if Columns==0
            Values=single(zeros(ProbeSetNb,PointNb));
        else
            Values=single(zeros(ProbeSetNb,length(Columns)));
        end
    else
        h=errordlg('ZVar not loaded. Process canceled.');
        waitfor(h)
        error('process canceled')
    end
end

function DATASAVE(Type,Dir,ResRank,Values,Columns,ProbeSetNb)

File=sprintf('%s_%02u.float32le',Type,ResRank);
if Columns==0
    save_data(Values,File,Dir,'r+','single','ieee-le');
else
    save_data(Values,File,Dir,'r+','single','ieee-le',ProbeSetNb,[],Columns);    
end
