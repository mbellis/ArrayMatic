% export net in text format used as input in Jerome Program
% of graphical network representation
function net_keiko(ModelRank,NetRank,ExportType,CLimit,ALimit,CRatio,ARatio,varargin)

tic

global K P
KeikoDir='E:\Sosma\RayProg\App\Cpp\keiko';

if nargin>4
    DiffVal=varargin{1};
else
    DiffVal=0;
end

NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank),sprintf('n%05u',NetRank));
PsNb=P.chip.probeSetNb;
Bindex=zeros(PsNb,1);

BlocSize=1000;
BlocStart=1;
BlocNb=ceil(PsNb/BlocSize);
BlocEnd=BlocNb;
% CRatio=K.net{ModelRank}.export.cratio(ExportRank);

cd(NetDir)
if ~exist(fullfile(NetDir,'keiko'),'dir')
    mkdir('keiko')
end
ResDir=fullfile(NetDir,'keiko');

cd(ResDir)
   
IndexFile=sprintf('m%u_index',ModelRank);
FileName=sprintf('%s.txt',IndexFile);

Probeset={};
Continue=1;

if exist(FileName,'file')==2
    Continue=0;
end
if Continue==1    
    IndexFid = fopen(FileName,'w');
    PsRank=1:PsNb;
    Probeset=P.chip.probeSetIds;
    for PsL=1:length(PsRank)-1
        fprintf(IndexFid,'%s\n',Probeset{PsL});
    end
    PsL=PsL+1;
    fprintf(IndexFid,'%s',Probeset{PsL});
    fclose(IndexFid)
end

% EXPORT DEF IF NECESSARY
DefFile=sprintf('m%u_def',ModelRank);
FileName=sprintf('%s.csv',DefFile);
Continue=1;
if exist(FileName,'file')==2
    Continue=0;
end
if Continue==1
    DefFid = fopen(FileName,'w');
    for PsL=1:PsNb
        fprintf(DefFid,'"%s","%s","%s","","","%s"',Probeset{PsL},Probeset{PsL},Probeset{PsL},Probeset{PsL});
    end    
    fclose(DefFid)
end

if ~isequal(ExportType,'clique')
    if DiffVal==0
        RezoFile=sprintf('m%u_n%u_cl%u_cr%u_al%u_ar%u',ModelRank,NetRank,CLimit,CRatio,Alimit,ARatio);
    else
        RezoFile=sprintf('m%u_n%u_cl%u_cr%u_al%u_ar%u_diffval',ModelRank,NetRank,CLimit,CRatio,ALimit,ARatio);
    end
else
    RezoFile=sprintf('m%u_n%u_cl%u_cr%u_al%u_ar%u_clique',ModelRank,NetRank,CLimit,CRatio,ALimit,ARatio);
end
if ~isequal(ExportType,'clique')
    BatFile=sprintf('%s.bat',RezoFile);
    BatFid = fopen(BatFile,'w');
    fprintf(BatFid,'copy %s.txt %s\\val.txt\n',RezoFile,KeikoDir);
    fprintf(BatFid,'copy %s\\%s.csv %s\\def.csv\n',K.dir.amcResults,DefFile,KeikoDir);
    fprintf(BatFid,'copy %s\\%s.txt %s\\index.txt\n',K.dir.amcResults,IndexFile,KeikoDir);
    fprintf(BatFid,'cd %s\n',KeikoDir);
    fprintf(BatFid,'keiko val def index');
    fclose(BatFid)
end

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
        DataFile=sprintf('c_m%u_n%u.4mat',ModelRank,NetRank);
        C=load_data(DataFile,'./',PsNb,PsNb,'uint8','ieee-le',FirstLine:LastLine);
        MemC=C;
        
      
        % forces the diagonal to equal 0 ( 1 not supported by jerome's prog )        
        StartCol=BlocSize*(BlocL-1);
        for PosL=1:LineNb
            C(PosL,StartCol+PosL)=0;
        end
        
        %load ANTI
        DataFile=sprintf('a_m%u_n%u.4mat',ModelRank,NetRank);
        A=load_data(DataFile,'./',PsNb,PsNb,'uint8','ieee-le',FirstLine:LastLine);
        MemA=A;
        %process CORR
        
        %eliminate point below limits
        C=single(C);
        A=single(A);
        CIndex=find(A~=0&C./A<=CRatio);
        C(CIndex)=0;
        clear CIndex
        CIndex=find(C<CLimit);
        C(CIndex)=0;
        %C must be >= CLimit to be eventually processed to C-A
        if DiffVal==1         
                CIndex=find(A~=0&C./A>CRatio);
                C(CIndex)=C(CIndex)-A(CIndex);
        end
        clear A
        if ~isequal(ExportType,'clique')
            save MemC C
            clear C;
        end
        clear CIndex
        
        %process ANTI
        if ~isequal(ExportType,'clique')
            %DataFile=sprintf('c_m%u_n%u.4mat',ModelRank,NetRank);
            %C=load_data(DataFile,'./',PsNb,PsNb,'uint8','ieee-le',FirstLine:LastLine);
            %DataFile=sprintf('a_m%u_n%u.4mat',ModelRank,NetRank);
            %A=load_data(DataFile,'./',PsNb,PsNb,'uint8','ieee-le',FirstLine:LastLine);
            C=single(MemC);
            clear MemC
            A=single(MemA);
            clear MemA
            AIndex=find(C~=0&A./C<=ARatio);
            A(AIndex)=0;
            clear AIndex
            AIndex=A<ALimit;
            A(AIndex)=0;
            %A must be >= ALimit to be eventually processed to A-C
            if DiffVal==1
                AIndex=find(C~=0&A./C>ARatio);
                A(AIndex)=A(AIndex)-C(AIndex);
            end
            clear C
            clear AIndex
            load MemC;
        end


        %write each line

        for LineL=1:size(C,1)

            CLine=C(LineL,:);
            %rank of current probeset minus 1 (zero based index)
            fprintf(RezoFid,'>%d\t0\t0\t0\t0\n',LineL+OffSet-1);
            Bindex(LineL+OffSet)=1;
            if ~isempty(CLine)
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
                if ~isempty(ALine)
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
        BlocL
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