%===================
% FUNCTION LOAD_DATA
%===================
% 
%LOAD_DATA loads data from a binary file

%INPUT PARAMETERS
% 1      DataFile: data file name
% 2       DataDir: data file directory name
% 3        LineNb: number of line of the matrix stored in DataFile
% 4         ColNb: number of columns of the matrix stored in DataFile
%                  data are stored in column order
% 5     Precision: type of data 
%                  Precisions={'int8','int16','int32','int64','uint8','uint16',
%                 'uint32','uint64','float32','float64','double','single'};
% 6 MachineFormat: endianess
%                  MachineFormat={'ieee-be','ieee-le','b','l'}                 
% varargin:
% 7     LineIndex: the index of lines to be loaded
%                  Can be alone (=> ColIndex is then set to 1:ColNb)
% 8      ColIndex: the index of columns to be loaded
%                  In this case, LineIndex must be indicated even if not
%                  particular selection is made (=> LineIndex=1:LineNb)

%OUTPUT
%Success: indicates success/fail of the process {0,1}
%   Data: a {LineNb or length(LineIndex)} x {ColNb or length(ColIndex}) matrix 


function [Data,varargout]=load_data(DataFile,DataDir,LineNb,ColNb,Precision,MachineFormat,varargin)

Precisions={'int8','int16','int32','int64','uint8','uint16','uint32','uint64','float32','float64','double','single'};
MachineFormats={'ieee-be','ieee-le','b','l'};
Success=1;
Data=[];
if isempty(DataDir)
    KeepOpen=1;
else
    KeepOpen=0;
end
NargOut=nargout;

%verify Precision
if isempty(find(strmatch(Precision,Precisions,'exact')))
    Success=MESSAGE(0,1,NargOut,sprintf('Precision %s is not allowed (use (u)int{8,16,32,64},float{32,64},single,double)',Precision));
elseif  isequal(Precision,'single')
    Precision='float32';
end
%verify machine format
if isempty(find(strmatch(MachineFormat,MachineFormats,'exact')))
    Success=MESSAGE(0,1,NargOut,sprintf('machine format %s is not allowed (use b,l,ieee-le,ieeee-be)',MachineFormat));
end


%open file if necessary
if Success
    if ~isempty(DataDir)
        cd(DataDir)
        [Fid,Msg]=fopen(DataFile,'r',MachineFormat);
        if Fid==-1
            Success=MESSAGE(Fid,1,NargOut,sprintf('%s cannot be opened: %s',DataFile,Msg));    
        end
    else
        Fid=DataFile;
    end   
end
if Success
    if nargin==6
        %read all data
        Status=fseek(Fid,0,'bof');
        if Status==-1
                Success=MESSAGE(Fid,KeepOpen,NargOut,sprintf('cannot position at 0 in %s',DataFile));
        end
        if LineNb>0
        [Data,Count]=fread(Fid,[LineNb,ColNb],Precision,MachineFormat);
        if Count~=LineNb*ColNb
            Success=MESSAGE(Fid,KeepOpen,NargOut,sprintf('%u items read instead of %u',Count,LineNb*ColNb));
        end
        else
            Data=fread(Fid,Inf,Precision,MachineFormat);
        end
    else
        switch Precision
            case {'int8','uint8'}
                Size=1;
            case {'int16','uint16'}
                Size=2;
            case {'int32','uint32','float32','single'}
                Size=4;
            case {'int64','uint64','float64','double'};
                Size=8;
        end
        LineIndex=varargin{1};
        LineFlag=1;
        if isempty(LineIndex)
            LineIndex=1:LineNb;
            LineFlag=0;
        end
        if nargin==8
            ColIndex=varargin{2};
            if isempty(ColIndex)
                ColIndex=1:ColNb;
            end
        else
            ColIndex=1:ColNb;
        end
        Data=zeros(length(LineIndex),length(ColIndex));
        for ColL=1:length(ColIndex)
            Status=fseek(Fid,(ColIndex(ColL)-1)*LineNb*Size,'bof');
            if Status==-1
                Success=MESSAGE(Fid,KeepOpen,NargOut,sprintf('cannot position at %u in %s',(ColIndex(ColL)-1)*Size,DataFile));
                break
            end
            [CurrCol,Count]=fread(Fid,[LineNb,1],Precision);
            
            if Count==LineNb
                if LineFlag
                    Data(:,ColL)=CurrCol(LineIndex);
                else
                    Data(:,ColL)=CurrCol;
                end
            else
                Success=MESSAGE(Fid,KeepOpen,NargOut,sprintf('%u items read instead of %u',Count,LineNb));
                waitfor(h)
                break
            end
        end
    end
end
if Success
    if ~isequal(Precision,'float64')&~isequal(Precision,'double')
        if isequal(Precision,'float32')
            Precision='single';
        end
        eval(sprintf('Data=%s(Data);',Precision));
    end
end
if nargout==2
    varargout{1}=Success;        
end
if KeepOpen==0 & Fid>0 & Success
    fclose(Fid);
end


function Success=MESSAGE(Fid,KeepOpen,NargOut,Message)
if NargOut==2
    h=warndlg(Message);
    waitfor(h)
else
    h=errordlg(Message);
    waitfor(h)
    error('process canceled')
end
if KeepOpen==0 & Fid>0
    fclose(Fid);
end
Success=0;