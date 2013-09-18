%TRS_NPY2PS is used to write binary files from npy format (itself writen by python program mat:
% py from matlab table
%the format is ps x comp

%INPUT PARAMETERS
%ModelRank: chip model rank
%OldTableSize: the current size
%NewTableSize: the wanted size

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
%flip table
%trs_npy2ps([2,3,5],{[10,100],[10],[10,100]},[13,25,25])
%trs_npy2ps([2,3,5],{[10],[10],[10]},[13,25,25])
%trs_npy2ps([6,8,27],{[10,100],[10,100],[10,100]},[10,46,26])
%trs_npy2ps([6,8,27],{[10],[10],[10]},[10,46,26])
%trs_npy2ps([8,27],{[10,100],[10,100]},[46,26])
%trs_npy2ps([8,27],{[1000],[1000]},[46,26])
function trs_npy2ps(ChipRank,Fdr,BlocNb)
global K

ChipNb=length(ChipRank);

for ChipL=1:ChipNb
    for FdrL=1:length(Fdr{ChipL})

        Continue=1;
        tic
        CurrChipRank=ChipRank(ChipL)
        %load point info
        K.dir.point='/home/mbellis/sosma/raydata/mlf/point';
        cd(K.dir.point)
        eval(sprintf('load m%u_point',CurrChipRank))
        BiolNb=length(P.biol.scoupleindex{1});
        RowNb=(BiolNb-1)*BiolNb/2;
        cd(fullfile(K.dir.net,sprintf('m%u',CurrChipRank),'table'))
        %delete eventualy existing files
        TableL=1;
        while exist(sprintf('t_m%u_fdr_%u.uint8',CurrChipRank,TableL),'file')
            eval(sprintf('delete t_m%u_fdr_%u.uint8',CurrChipRank,TableL))
            TableL=TableL+1;
        end
        TableL=1;
        while exist(sprintf('table_m%u_fdr%u_%u.uint8',CurrChipRank,Fdr{ChipL}(FdrL),TableL),'file')
            eval(sprintf('delete table_m%u_fdr%u_%u.uint8',CurrChipRank,Fdr{ChipL}(FdrL),TableL))
            TableL=TableL+1;
        end           
        RowNb1=0;
        RowNb2=0;
        for TableL=1:BlocNb(ChipL)        
            %npy format have a header            
            fid=fopen(sprintf('table_m%03u_A_fdr_%03u_%u.npy',CurrChipRank,TableL,Fdr{ChipL}(FdrL)),'r','ieee-le');
            if fid==-1
                fid=fopen(sprintf('table_m%u_fdr_%u_%u.npy',CurrChipRank,Fdr{ChipL}(FdrL),TableL),'r','ieee-le');
                if fid==-1
                    if Fdr{ChipL}(FdrL)==100
                        fid=fopen(sprintf('table_m%03u_A_fdr_%03u_01.npy',CurrChipRank,TableL),'r','ieee-le');
                    end
                end
            end   
            if fid>0
                Comp=fread(fid,'uint8');
            else
                'stop'
            end
            fclose(fid);
            %clear header (size coded in little endian sur 2 bytes (9 and 10))
            HeadLength=bin2dec([dec2bin(Comp(10)),dec2bin(Comp(9))]);
            Comp(1:10+HeadLength)=[];
            if mod(length(Comp),RowNb)==0
                ColNb=round(length(Comp)/RowNb);
                if TableL==1
                    %for the next step when matrix is transposed
                    RowNb1=ColNb;
                elseif TableL==BlocNb(ChipL)
                    RowNb2=ColNb;
                end
            else
                sprintf('row nb=%u, comp nb=%u mod=%',RowNb,length(Comp),mod(length(Comp),RowNb))
                Continue=0;
                break
            end
            Comp=reshape(Comp,RowNb,ColNb);
            Comp=Comp';
            save_data(Comp,sprintf('t_m%u_fdr_%u.uint8',CurrChipRank,TableL),'./','r+','uint8','ieee-le');
        end
        toc

        if Continue
            tic
            %reconstruct Table to have colomn or results
            ColNb=RowNb;
            LoopNb=ceil(RowNb/15000);
            %for BlocL=2:ceil(ColNb/15000)
            for BlocL=1:LoopNb
                ColRank=(BlocL-1)*15000+1:min(BlocL*15000,ColNb);
                Comp=[];
                %load data and concatenate columns
                for TableL=1:BlocNb(ChipL)
                    if TableL<BlocNb(ChipL)
                        RowNb=RowNb1;
                    else
                        RowNb=RowNb2;
                    end
                    Comp=[Comp;load_data(sprintf('t_m%u_fdr_%u.uint8',CurrChipRank,TableL),'./',RowNb,ColNb,'uint8','ieee-le',1:RowNb,ColRank)];
                end
                save_data(Comp,sprintf('table_m%u_fdr%u_%u.uint8',CurrChipRank,Fdr{ChipL}(FdrL),BlocL),'./','r+','uint8','ieee-le');
            end
            toc
            %delete files
            for TableL=1:BlocNb(ChipL)
                eval(sprintf('delete t_m%u_fdr_%u.uint8',CurrChipRank,TableL))
            end
        end
    end
end