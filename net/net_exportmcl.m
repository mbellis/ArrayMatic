function net_exportmcl(ModelRank,NetRank,MclRank,CluRank)
%net_exportmcl([8,27],[55,24],[5,5],[1,5])
global K
cd(K.dir.affyMetadata)
if length(ModelRank)>2 |length(NetRank)>2
    h=errordlg('at most 2 networks');
    waitfor(h)
    error('process canceled')
end
if length(ModelRank)==1
    ModelRank=[ModelRank,ModelRank];
end
if length(NetRank)==1
    NetRank=[NetRank,NetRank];
end

FileName=sprintf('m%u_m%u_commonps.mat',min(ModelRank),max(ModelRank));

if exist(FileName,'file')
    load(FileName)
    if ModelRank(1)>ModelRank(2)
        Temp=ComPsRank;
        ComPsRank(:,1)=ComPsRank(:,2);
        ComPsRank(:,2)=Temp(:,1);
        clear Temp
    end
else
    h=errordlg(sprintf('no correspondance file between m%y and m%u',ModelRank(1),ModelRank(2)));
    waitfor(h)
    error('process canceled')
end

ComClu=zeros(length(ComPsRank),2);
for ModelL=1:2
    NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank(ModelL)),sprintf('n%05u',NetRank(ModelL)),'mcl');
    cd(NetDir)
    load(sprintf('m%un%u_mcl.mat',ModelRank(ModelL),NetRank(ModelL)))
    if ModelRank(1)~=ModelRank(2)
        ComClu(:,ModelL)=Clu(ComPsRank(:,ModelL),MclRank(ModelL));
    else
        ComClu(:,ModelL)=Clu(:,MclRank(ModelRank));
    end
    clear Clu
end
PsPos=find(ComClu(:,1)==CluRank(1)&ComClu(:,2)==CluRank(2));

%calculateps index
PsRank{1}=ComPsRank(PsPos,1);
PsRank{1}=ComPsRank(PsPos,2);
%export MCL region
for ModelL=1:2
%load data
NetDir=fullfile(K.dir.net,sprintf('m%03u',ModelRank(ModelL)),sprintf('n%05u',NetRank(ModelL)));
cd(NetDir)
Corr=load_data(sprintf('c_m%u_n%u.4mat',ModelRank(ModelL),NetRank(ModelL)),'./',K.chip.probeSetNbs{ModelRank(ModelL)}(1),K.chip.probeSetNbs{ModelRank(ModelL)}(1),'uint8','ieee-le',PsRank{ModelL},PsRank{ModelL});
Anti=load_data(sprintf('a_m%u_n%u.4mat',ModelRank(ModelL),NetRank(ModelL)),'./',K.chip.probeSetNbs{ModelRank(ModelL)}(1),K.chip.probeSetNbs{ModelRank(ModelL)}(1),'uint8','ieee-le',PsRank{ModelL},PsRank{ModelL});
'stop'

%output fot GAP


%output for W
end