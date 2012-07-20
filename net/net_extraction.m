function net_extraction(ModelRank,NetRank)
%net_extraction([8,27,42],[54,1])
global K

if length(ModelRank)~=3
    h=errordlg('needs 3 chips model');
    waitfor(h)
    error('process canceled')
end
if length(NetRank)~=2
    h=errordlg('needs 2 networks');
    waitfor(h)
    error('process canceled')
end

cd(K.dir.affyMetadata)
FileName=sprintf('m%u_m%u_commonps.mat',min(ModelRank(1:2)),max(ModelRank(1:2)));
if exist(FileName,'file')
    load(FileName)
    if ModelRank(1)>ModelRank(2)
        ComPsRank(:,1)=ComPsRank(:,2);
    end
    ColRanks=ComPsRank(:,1);
    clear ComPsRank
else
    h=errordlg(sprintf('no correspondance file between m%y and m%u',ModelRank(1),ModelRank(2)));
    waitfor(h)
    error('process canceled')
end

PsNb(1)=K.chip.probeSetNbs{ModelRank(1)}(1);
PsNb(2)=length(ColRanks);
BlocNb=ceil(PsNb(2)/1000);
NetDir{1}=fullfile(K.dir.net,sprintf('m%03u',ModelRank(1)),sprintf('n%05u',NetRank(1)));
NetDir{2}=fullfile(K.dir.net,sprintf('m%03u',ModelRank(3)),sprintf('n%05u',NetRank(2)));
mkdir(K.dir.net,sprintf('m%03u',ModelRank(3)))
mkdir(fullfile(K.dir.net,sprintf('m%03u',ModelRank(3))),sprintf('n%05u',NetRank(2)))
for BlocL=1:BlocNb
    if BlocL<BlocNb
        RawRanks=ColRanks((BlocL-1)*1000+1:BlocL*1000);
    else
        RawRanks=ColRanks((BlocL-1)*1000+1:PsNb(2));
    end
    Data=load_data(sprintf('c_m%u_n%u.4mat',ModelRank(1),NetRank(1)),NetDir{1},PsNb(1),PsNb(1),'uint8','ieee-le',RawRanks,ColRanks);
    save_data(Data',sprintf('c_m%u_n%u.4mat',ModelRank(3),NetRank(2)),NetDir{2},'a','uint8','ieee-le')
    Data=load_data(sprintf('a_m%u_n%u.4mat',ModelRank(1),NetRank(1)),NetDir{1},PsNb(1),PsNb(1),'uint8','ieee-le',RawRanks,ColRanks);
    save_data(Data',sprintf('a_m%u_n%u.4mat',ModelRank(3),NetRank(2)),NetDir{2},'a','uint8','ieee-le')
end
