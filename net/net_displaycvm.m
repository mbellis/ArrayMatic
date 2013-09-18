%=========================%
% FUNCTION NET_DISPLAYCVM %
%=========================%

%NET_DISPLAYCVM: displays a CVM eventually structured by a MCL result or/and a NTS result
% by using ImageMagick: several intermediate bmp files are first written, and finally merged
% into a single image.
%INPUT PARAMETERS
% 1 ChipRank: chip rank
% 2 NetRank: network(s) to be displayed 
% 3 PsOrder: list od ordered probe set ranks
% 4 FileName: name added to chiprank and netrank information in output file names;
% 5 ColorNb: range of colors
% 6 ImageSize: size of the final image recomposed by ImageMagik
% 7 ClearBmpFlag: erase intermediate bmp files

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


function net_displaycvm(ChipRank,NetRank,PsOrder,FileName,ColorNb,ImageSize,ClearBmpFlag)
global K

PsNb=length(PsOrder);
BlocNb=ceil(PsNb/500);
if ~isempty(FileName)
    Output=sprintf('m%u_n%u_%s',ChipRank,NetRank,FileName);
else
    Output=sprintf('m%u_n%u',ChipRank,NetRank);
end
Colors=colors(colormap,ColorNb);
NetDir=fullfile(K.dir.net,sprintf('m%u',ChipRank),sprintf('n%u',NetRank));
mkdir(NetDir)
for BlocL1=1:BlocNb
    fprintf('processing bloc %u\n',BlocL1)
    CurrPs1=PsOrder((BlocL1-1)*500+1:min(PsNb,BlocL1*500));
    %process all clu and save bmp files
    for BlocL2=BlocL1:BlocNb
        CurrPs2=PsOrder((BlocL2-1)*500+1:min(PsNb,BlocL2*500));
        [C,A]=load_cvm(ChipRank,NetRank,CurrPs1,CurrPs2,1,1,0);
        A=A{1};
        if BlocL2==BlocL1
            %make composit bloc (C and A values)
            C=triu(C,1)+tril(A,-1);
        end
        FileName=sprintf('%s_bloc%u_vs_bloc%u.bmp',Output,BlocL1,BlocL2);
        cd(NetDir)
        imwrite(C,Colors,FileName,'bmp')
        if BlocL2~=BlocL1
            FileName=sprintf('%s_bloc%u_vs_bloc%u.bmp',Output,BlocL2,BlocL1);
            imwrite(A',Colors,FileName,'bmp')
        end
    end
end
%write commands for montage (image magick)
cd(NetDir)
FileNameList={};
for BlocL1=1:BlocNb
    for BlocL2=BlocL1:BlocNb
        FileName=sprintf('%s_bloc%u_vs_bloc%u.bmp',Output,BlocL1,BlocL2);
        FileNameList{end+1,1}=FileName;
        if BlocL2~=BlocL1
            FileName=sprintf('%s_bloc%u_vs_bloc%u.bmp',Output,BlocL2,BlocL1);
            FileNameList{end+1,1}=FileName;
        end
    end
end
%reorder lines
FileOrder=zeros(BlocNb,BlocNb);
Rank=0;
for BlocL1=1:BlocNb
    for BlocL2=BlocL1:BlocNb
        Rank=Rank+1;
        FileOrder(BlocL1,BlocL2)=Rank;
        if BlocL2~=BlocL1
            Rank=Rank+1;
            FileOrder(BlocL2,BlocL1)=Rank;
        end
    end
end
FileOrder=FileOrder';

fid=fopen(sprintf('%s_imagelist.txt',Output),'w');
CurrFileOrder=FileOrder(1:BlocNb,1:BlocNb);
CurrFileNameList=FileNameList(CurrFileOrder(:));
for NameL=1:length(CurrFileNameList)
    fprintf(fid,'%s\n',CurrFileNameList{NameL});
end
fclose(fid)

fid=fopen(sprintf('%s.bat',Output),'w');
fprintf(fid,'montage @%s_imagelist.txt -tile %ux%u -geometry "10x10<" -mode Concatenate  %s.bmp',Output,BlocNb,BlocNb,Output)
fclose(fid)
!chmod a+x *.bat
eval(sprintf('!./%s.bat',Output))
eval(sprintf('!convert %s.bmp -resize %ux%u resized_%s.bmp',Output,ImageSize,ImageSize,Output))
FileName=sprintf('%s.bmp',Output);
delete(FileName)
if ClearBmpFlag
    for BlocL1=1:BlocNb
        %process all clu and save bmp files
        for BlocL2=BlocL1:BlocNb
            FileName=sprintf('%s_bloc%u_vs_bloc%u.bmp',Output,BlocL1,BlocL2);
            delete(FileName)
            if BlocL2~=BlocL1
                FileName=sprintf('%s_bloc%u_vs_bloc%u.bmp',Output,BlocL2,BlocL1);
                delete(FileName)
            end
        end
    end
end
