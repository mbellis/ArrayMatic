% ExportType: 'clique'
% CLimit
% ALimit
% CRatio : C/A inferior limit
% ARatio : A/C inferior limit
% export net in text format used as input in Jerome Program
% of graphical network representation
%net_keiko(2,80,'',50,50,1,1,1);
%net_keiko(2,80,'',0,0,1,1,1);
function net_keiko(ChipRank,NetRank,ExportType,CLimit,ALimit,CRatio,ARatio,varargin)

tic
global K P
KeikoDir='E:\Sosma\RayProg\App\Cpp\keiko';

if nargin>4
    DiffVal=varargin{1};
else
    DiffVal=0;
end
ChipPos=find(K.chip.rank==ChipRank);
PsNb=K.chip.probesetNb(ChipPos);
NetDir=fullfile(K.dir.net,sprintf('m%u',ChipRank),sprintf('n%u',NetRank));

Bindex=zeros(PsNb,1);

BlocSize=500;
BlocStart=1;
BlocNb=ceil(PsNb/BlocSize);
BlocEnd=BlocNb;


cd(NetDir)
if ~exist(fullfile(NetDir,'keiko'),'dir')
    mkdir('keiko')
end
ResDir=fullfile(NetDir,'keiko');

cd(ResDir)
   
IndexFile=sprintf('m%u_index',ChipRank);
FileName=sprintf('%s.txt',IndexFile);

Probeset={};
Continue=1;

if exist(FileName,'file')==2
    Continue=0;
end


if ~isequal(ExportType,'clique')
    if DiffVal==0
        RezoFile=sprintf('m%u_n%u_cl%u_cr%u_al%u_ar%u',ChipRank,NetRank,CLimit,CRatio,Alimit,ARatio);
    else
        RezoFile=sprintf('m%u_n%u_cl%u_cr%u_al%u_ar%u_diffval',ChipRank,NetRank,CLimit,CRatio,ALimit,ARatio);
    end
else
    RezoFile=sprintf('m%u_n%u_cl%u_cr%u_al%u_ar%u_clique',ChipRank,NetRank,CLimit,CRatio,ALimit,ARatio);
end

% if ~isequal(ExportType,'clique')
%     BatFile=sprintf('%s.bat',RezoFile);
%     BatFid = fopen(BatFile,'w');
%     fprintf(BatFid,'copy %s.txt %s\\val.txt\n',RezoFile,KeikoDir);
%     fprintf(BatFid,'cd %s\n',KeikoDir);
%     fprintf(BatFid,'keiko val def index');
%     fclose(BatFid)
% end

if isequal(ExportType,'clique')
    % selection on climit is made at the level of cpp prog
    CLimit=1;
end

FileName=sprintf('%s.txt',RezoFile);
StatFile=sprintf('%s_stat.txt',RezoFile);

Continue=1;
if exist(FileName,'file')==2
    Answer=questdlg('do you want to overwrite existing NETWORK file ?','','Yes','No','Yes');
    if isequal(Answer,'No')
        Continue=0;
    end
end

if Continue==1
    rand('twister',5489);
    % EXPORT CORR/ANTI VALUES
    RezoFid = fopen(FileName,'w');
    cd(NetDir)
    ATotalNb=0;
    CTotalNb=0;

    for BlocL=BlocStart:BlocEnd
        OffSet=BlocSize*(BlocL-1);
        FirstLine=(BlocL-1)*BlocSize+1;
        if BlocL<BlocEnd
            LastLine=FirstLine+BlocSize-1;
        else
            LastLine=PsNb;
        end
        LineNb=LastLine-FirstLine+1;
        
        %load CORR
        DataFile=sprintf('c_m%u_n%u.4mat',ChipRank,NetRank);
        C=load_data(DataFile,'./',PsNb,PsNb,'uint8','ieee-le',FirstLine:LastLine);
        MemC=C;              
        % forces the diagonal to equal 0 ( 1 not supported by jerome's prog )        
        StartCol=BlocSize*(BlocL-1);
        for PosL=1:LineNb
            C(PosL,StartCol+PosL)=0;
        end
        
        %load ANTI
        DataFile=sprintf('a_m%u_n%u.4mat',ChipRank,NetRank);
        A=load_data(DataFile,'./',PsNb,PsNb,'uint8','ieee-le',FirstLine:LastLine);
        MemA=A;
        
        %process CORR        
        %eliminate point below limits
        C=single(C);
        A=single(A);
        CIndex=find(A~=0&C./A<=CRatio);
        C(CIndex)=0;
        clear CIndex        
        %C must be >= CLimit to be eventually processed to C-A
        if DiffVal==1         
                %CIndex=find(A~=0&C./A>CRatio);
                %C(CIndex)=C(CIndex)-A(CIndex);
                C=C-A;
                %remove negative values
                %CIndex=C<0;
                %C(CIndex)=0;
        end
        CIndex=find(C<CLimit);
        C(CIndex)=0;        
        clear A
        if ~isequal(ExportType,'clique')
            save MemC C
            clear C;
        end
        clear CIndex
        
        %process ANTI
        if ~isequal(ExportType,'clique')
            %DataFile=sprintf('c_m%u_n%u.4mat',ChipRank,NetRank);
            %C=load_data(DataFile,'./',PsNb,PsNb,'uint8','ieee-le',FirstLine:LastLine);
            %DataFile=sprintf('a_m%u_n%u.4mat',ChipRank,NetRank);
            %A=load_data(DataFile,'./',PsNb,PsNb,'uint8','ieee-le',FirstLine:LastLine);
            C=single(MemC);
            clear MemC
            A=single(MemA);
            clear MemA
            AIndex=find(C~=0&A./C<=ARatio);
            A(AIndex)=0;
            clear AIndex            
            %A must be >= ALimit to be eventually processed to A-C
            if DiffVal==1
                %AIndex=find(C~=0&A./C>ARatio);
                %A(AIndex)=A(AIndex)-C(AIndex);
                A=A-C;
                %remove negative values
                %AIndex=A<0;
                %A(AIndex)=0;
            end
            AIndex=A<ALimit;
            A(AIndex)=0;
            clear C
            clear AIndex
            load MemC;
        end

        %write each line
        for LineL=1:size(C,1)

            CLine=C(LineL,:);
            %rank of current probeset minus 1 (zero based index)
            
            Bindex(LineL+OffSet)=1;
            FirstLineFlag=0;
            %random position
            Alpha=rand(1)*2*pi;
            Beta=rand(1)*2*pi;
            z=100*sin(Beta);
            x=100*cos(Alpha);
            y=100*sin(Alpha);            
            if ~isempty(find(CLine))                
                        fprintf(RezoFid,'>%d\t0\t%.3f\t%.3f\t%.3f\n',LineL+OffSet-1,x,y,z);
                FirstLineFlag=1;
                Index=find(CLine);
                CTotalNb=CTotalNb+length(Index);
                [Val SortIndex]=sort(CLine(Index));
                Index=Index(SortIndex);
                Val=fliplr(Val);
                Index=fliplr(Index);
                Val=Val*100;
                for WBloc=1:length(Index)
                    %rank of the probeset correlated to the current probest (minus 1 (zero based index))
                    fprintf(RezoFid,'%d\t%d\t%d\n',Index(WBloc)-1,Val(WBloc),0);
                    Bindex(Index(WBloc))=1;
                end
            end
            if ~isequal(ExportType,'clique')
                ALine=A(LineL,:);
                if ~isempty(find(ALine))
                    if FirstLineFlag==0
                        fprintf(RezoFid,'>%d\t0\t%.3f\t%.3f\t%.3f\n',LineL+OffSet-1,x,y,z);
                    end
                    Index=find(ALine);
                    ATotalNb=ATotalNb+length(Index);
                    [Val SortIndex]=sort(ALine(Index));
                    Index=Index(SortIndex);
                    Val=Val*100;
                    for WBloc=1:length(Index)
                        fprintf(RezoFid,'%d\t%d\t%d\n',Index(WBloc)-1,-Val(WBloc),0);
                        Bindex(Index(WBloc))=1;
                    end
                end
            end
        end
        %clear A C Index
        fprintf('%u upon %u\n',BlocL,BlocNb)
        if ~isequal(ExportType,'clique')
            ATotalNb
        end
        CTotalNb

    end
    fclose(RezoFid);
    cd(K.dir.amcResults)
    StatFid = fopen(StatFile,'w');
    fprintf(StatFid,'C Total Edge Nb : %d\n ',CTotalNb);
    if ~isequal(ExportType,'clique')
        fprintf(StatFid,'A Total Edge Nb : %d\n ',ATotalNb);
    end
    fprintf(StatFid,'A Total Ps Nb : %d\n ',length(find(Bindex)));
    fclose(StatFid);
end
toc