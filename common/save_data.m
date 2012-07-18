%SAVE_DATA
%save data in a binary form

%INPUT
% Data: The data to be saved. The data are saved
% at the position given by ColumnIndex
% DataFile
% DataDir
% Permission
% Precision
% varargin
% LineNb (varargin{1}): the number of lines in the matrix
% LineIndex (varargin{2}): the position of the inserted lines(s)
% ColumnIndex (varargin{3}): the position of the inserted column(s) in the
% written matrix. If Length(ColumnIndex)==1 and there exist several columns
% in Data, then all the bloc is written at this column position. Otherwise,
% Data columns are written successively to the indicated column positions.



%VERSIONS
% V04 - 2010 07 20 - Add LineIndex. ColumnIndex may have several columns.
%                    Add MESSAGE function.
% V03 - 2010 07 05 - Pass Fid in parmeter (to prevent to open/close the same file
%                       several times)
% V02 - 2010 06 22 - Can write columns - Uses machine format
% V01 - 2010 05 14- First version

function varargout=save_data(Data,DataFile,DataDir,Permission,Precision,MachineFormat,varargin)

Permissions={'w','a','r+','w+','a+'};
Precisions={'int8','int16','int32','int64','uint8','uint16','uint32','uint64','float32','float64','double','single'};
MachineFormats={'ieee-be','ieee-le','b','l'};
Success=1;
if isempty(DataDir)
    KeepOpen=1;
else
    KeepOpen=0;
end
NargOut=nargout;
    
%verify nargin
if nargin~=6&nargin~=9    
    Success=MESSAGE(0,1,NargOut,'save_data must have 6 or 9 parameters');
end
if Success
    %verify permission
    if isempty(find(strmatch(Permission,Permissions,'exact')))
        Success=MESSAGE(0,1,NargOut,sprintf('permission %s is not allowed (use w,a,r+,w+,a+,W pr A)',Permission));        
    end
    %verify Precision
    if Success
        if isempty(find(strmatch(Precision,Precisions,'exact')))
            Success=MESSAGE(0,1,NargOut,sprintf('Precision %s is not allowed (use (u)int{8,16,32,64},float{32,64},single,double)',Precision));
        elseif  isequal(Precision,'single')
            Precision='float32';
        end
    end
    %verify machine format
    if isempty(find(strmatch(MachineFormat,MachineFormats,'exact')))
        Success=MESSAGE(0,1,NargOut,sprintf('machine format %s is not allowed (use b,l,ieee-le,ieeee-be)',MachineFormat));
    end
    %open file if necessary
    if Success
        if ~isempty(DataDir)
            cd(DataDir)
            if isequal(Permission,'r+')
                %create if not exist otherwise can't be used
                if ~exist(DataFile,'file')
                    Fid=fopen(DataFile,'w',MachineFormat);
                    fclose(Fid)
                end
            end
            [Fid,Msg]=fopen(DataFile,Permission,MachineFormat);
            if Fid==-1
                Success=MESSAGE(Fid,1,NargOut,sprintf('%s cannot be created: %s',DataFile,Msg));                
            end
        else
            Fid=DataFile;
        end
    end
    if Success
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
        if nargin==6
            %position at the end of file and write all the Data
            Status=fseek(Fid,0,'eof');
            if Status==-1
                Success=MESSAGE(Fid,KeepOpen,NargOut,sprintf('cannot position at the right position of %s ',DataFile));
            else
                Count=fwrite(Fid,Data,Precision);
                if Count~=size(Data,1)*size(Data,2)
                    Success=MESSAGE(Fid,KeepOpen,NargOut,sprintf('%u items written instead of %u',Count,size(Data,1)*size(Data,2)));
                end
            end
        elseif nargin==9
            %write Data at given positions
            LineNb=varargin{1};
            LineIndex=varargin{2};
            ColumnIndex=varargin{3};
            if isempty(LineIndex)
                LineFlag=0;
            else
                LineFlag=1;
            end
            if isempty(ColumnIndex)
                ColumnIndex=1;
            end
            if length(ColumnIndex)>1
                %verify the size of Data
                if length(ColumnIndex)~=size(Data,2)             
                    Success=MESSAGE(Fid,KeepOpen,NargOut,sprintf('Data has %u columns and ColumnIndex %u.',size(Data,2),length(ColumnIndex)));
                end
            end
            if Success
                if length(ColumnIndex)==1
                    %position at the right column position
                    Status=fseek(Fid,LineNb*(ColumnIndex(1)-1)*Size,'bof');
                    if Status==-1
                        Success=MESSAGE(Fid,KeepOpen,NargOut,sprintf('cannot position at the right position of %s ',DataFile));     
                    else
                        if LineFlag
                            %replace existing lines by Data
                            [CurrData,Count]=fread(Fid,[LineNb,size(Data,2)],Precision);
                            if Count~=LineNb*size(Data,2)
                                Success=MESSAGE(Fid,KeepOpen,NargOut,sprintf('%u items read instead of %u',Count,LineNb*size(Data,2)));
                            end
                            if ~isequal(Precision,'float64')&~isequal(Precision,'double')
                                if isequal(Precision,'float32')
                                    Precision='single';
                                end
                                eval(sprintf('CurrData=%s(CurrData);',Precision));
                            end
                            CurrData(LineIndex,:)=Data;
                            Data=CurrData;
                            %set to cursor to the right place
                            fseek(Fid,LineNb*(ColumnIndex(1)-1)*Size,'bof');
                        end
                        Count=fwrite(Fid,Data,Precision);
                        if Count~=size(Data,1)*size(Data,2)
                            Success=MESSAGE(Fid,KeepOpen,NargOut,sprintf('%u items written instead of %u',Count,size(Data,1)*size(Data,2)));
                        end
                    end
                else
                    for ColL=1:size(Data,2)
                        %position at the right column position
                        Status=fseek(Fid,LineNb*(ColumnIndex(ColL)-1)*Size,'bof');
                        if Status==-1          
                            Success=MESSAGE(Fid,KeepOpen,NargOut,sprintf('cannot position at the right position of %s ',DataFile));
                            break
                        else
                            if LineFlag
                                %replace existing lines by Data
                                [CurrData,Count]=fread(Fid,[LineNb,1],Precision);
                                if Count~=LineNb
                                    Success=MESSAGE(Fid,KeepOpen,NargOut,sprintf('%u items read instead of %u',Count,LineNb));
                                    break
                                end
                                if ~isequal(Precision,'float64')&~isequal(Precision,'double')
                                    if isequal(Precision,'float32')
                                        Precision='single';
                                    end
                                    eval(sprintf('CurrData=%s(CurrData);',Precision));
                                end
                                CurrData(LineIndex)=Data(:,ColL);                                
                                %set to cursor to the right place
                                fseek(Fid,LineNb*(ColumnIndex(ColL)-1)*Size,'bof');
                            else
                                CurrData=Data(:,ColL);
                            end
                            Count=fwrite(Fid,CurrData,Precision);
                            if Count~=size(CurrData,1)*size(CurrData,2)
                                Success=MESSAGE(Fid,KeepOpen,NargOut,sprintf('%u items written instead of %u',Count,size(Data,1)));
                                break
                            end
                        end
                    end
                end
            end
    
        end
        if nargout==1
            varargout{1}=Success;        
        end
        if KeepOpen==0 & Fid>0 & Success
            fclose(Fid);
        end
    end
end

function Success=MESSAGE(Fid,KeepOpen,NargOut,Message)
if NargOut==1
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