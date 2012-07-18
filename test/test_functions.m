%test existing functions

%function affy_rdn(Action,ExpName,CelDir,ChipName,ChipDir,LibDir)

function test_functions(Action)
global K

switch Action
    case 'rdn'
        RootDir='/home/mbellis/sosma';
        ExpName='microcebus';
        CelDir=fullfile(RootDir,'raydataraw','2009_DEVAU_MICROCEBUS','cel');
        ChipName='HG-U133_Plus_2';
        ChipDir=fullfile(RootDir,'data','affy','HG-U133_Plus_2');
        LibDir=fullfile(RootDir,'data','affy','Rice','LibFiles');
        SampleNames={'AF-973a-R01','AF-896f-R02','OF-008b-R01','OF-042b-R02','OF-043b-R03','OF-906f-R04','OF-046h-R05','OF-118h-R06','OF-143h-R07','OM-916d-R08','OM-921c-R09','OM-901g-R10','YF-052b-R01','YF-986e-R02','YF-984f-R03','YF-139h-R04','YM-954a-R05','YM-985e-R06'};
        PrintOrder=[15,9,1,2,3,11,4,6,8,12,13,10,5,18,16,7,14,17];

        affy_rdn('cel',ExpName,CelDir,ChipName,ChipDir,LibDir)
        affy_rdn('info',ExpName,CelDir,ChipName,ChipDir,LibDir)
        affy_rdn('rdn',ExpName,CelDir,ChipName,ChipDir,LibDir)
        affy_rdn('print',ExpName,CelDir,ChipName,ChipDir,LibDir,SampleNames,PrintOrder)


        ExpName='RIZ2010';
        CelDir=fullfile(RootDir,'raydataraw','2010_BELLAFIORE','raw_data','cel_files');
        ChipName='Rice';
        ChipDir=fullfile(RootDir,'data','affy','Rice');
        LibDir=fullfile(RootDir,'data','affy','Rice','LibFiles');


        affy_rdn('cel',ExpName,CelDir,ChipName,ChipDir,LibDir)
        affy_rdn('info',ExpName,CelDir,ChipName,ChipDir,LibDir)
        affy_rdn('rdn',ExpName,CelDir,ChipName,ChipDir,LibDir)
        affy_rdn('print',ExpName,CelDir,ChipName,ChipDir,LibDir)


    case 'data_geo'
        %DATA_GEO
        %function data_geo(Action,varargin)
        data_geo('example','GPL198','GSE21076')

    case 'save_load'

        %SAVE_DATA - LOAD_DATA
        Precisions={'int8','int16','int32','int64','uint8','uint16','uint32','uint64','float32','float64'};

        %matrix used to generate the different saved/loaded matrix
        Rp=rand(100);
        Neg=ceil(rand(5000,1)*10000);
        Rn=Rp;
        Rn(Neg)=-Rn(Neg);

        R={};
        %int8
        R{1}=int8(Rn*100);
        %int16
        R{2}=int16(Rn*(2^15-1));
        %int32
        R{3}=int32(Rn*(2^31-1));
        %int64
        R{4}=int64(Rn*(2^63-1));
        %int8
        R{5}=uint8(abs(R{2}));
        %int16
        R{6}=uint16(abs(R{2}));
        %int32
        R{7}=uint32(abs(R{3}));
        %int64
        R{8}=uint64(Rp*(2^63-1));
        %float32
        R{9}=single(R{4});
        %float64
        R{10}=double(R{4});
        
        %nullData
        %int8
        Z{1}=int8(zeros(100));
        %int16
        Z{2}=int16(zeros(100));
        %int32
        Z{3}=int32(zeros(100));
        %int64
        Z{4}=int64(zeros(100));
        %int8
        Z{5}=Z{1};
        %int16
        Z{6}=Z{2};
        %int32
        Z{7}=Z{3};
        %int64
        Z{8}=Z{4};
        %float32
        Z{9}=single(zeros(100));
        %float64
        Z{10}=double(zeros(100));

        DataDir=K.dir.test;
        Success=1;
        for PrecL=1:length(Precisions)
            DataFile=sprintf('test_%s',Precisions{PrecL});

            %save data
            Success=save_data(R{PrecL},DataFile,DataDir,'w',Precisions{PrecL},'ieee-le');
            if Success==0
                h=warndlg(sprintf('failed to write test_%s',Precisions{PrecL}));
                waitfor(h)
                Success=0;
            end

            %load all data
            [Data,Success]=load_data(DataFile,DataDir,100,100,Precisions{PrecL},'ieee-le');
            if Success
                if length(find(Data~=R{PrecL}))~=0
                    h=warndlg(sprintf('loaded test_%s is not identical to saved test_%s',Precisions{PrecL},Precisions{PrecL}));
                    waitfor(h)
                    Success=0;
                end
            else
                h=warndlg(sprintf('failed to load test_%s',Precisions{PrecL}));
                waitfor(h)
                Success=0;
            end
            
            %replace by zeros
            %save data
            Success=save_data(Z{PrecL},DataFile,DataDir,'w',Precisions{PrecL},'ieee-le');
            if Success==0
                h=warndlg(sprintf('failed to write test_%s',Precisions{PrecL}));
                waitfor(h)
                Success=0;
            end


            for InsertL=3
                %save partial data
                switch(InsertL)
                    case 1
                        LineIndex=1:2:100;
                        ColumnIndex=[3,10,25,78];
                        Insert=R{PrecL}(LineIndex,ColumnIndex);
                    case 2
                        LineIndex=1:2:100;
                        ColumnIndex=[];
                        Insert=R{PrecL}(LineIndex,:);
                    case 3
                        LineIndex=[];
                        ColumnIndex=[3,10,25,78];
                        Insert=R{PrecL}(:,ColumnIndex);
                end
                Success=save_data(Insert,DataFile,DataDir,'r+',Precisions{PrecL},'ieee-le',100,LineIndex,ColumnIndex);
                if Success==0
                    h=warndlg(sprintf('failed to write test_%s in InsertL %u',Precisions{PrecL},InsertL));
                    waitfor(h)
                    Success=0;
                end

                %load partial data
                [Data,Success]=load_data(DataFile,DataDir,100,100,Precisions{PrecL},'ieee-le',LineIndex,ColumnIndex);
                if Success
                    if length(find(Data~=Insert))~=0
                        h=warndlg(sprintf('loaded Insert of test_%s is not identical to those saved test_%s in InsertL %u',Precisions{PrecL},Precisions{PrecL},InsertL));
                        waitfor(h)
                        Success=0;
                    end
                else
                    h=warndlg(sprintf('failed to load Insert of test_%s in InsertL %u',Precisions{PrecL},InsertL));
                    waitfor(h)
                    Success=0;
                end
            end
        end
        if Success
            h=warndlg('test for save/load are all OK');
            waitfor(h)
        end
end






