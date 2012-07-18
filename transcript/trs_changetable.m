%TRS_CHANGETABLE is used in large analysis to change the size of matlab table containing the fdr values of comparisons.
%see trs_maketableTwo types of tables are constructed :

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

function trs_changetable(ModelRank,OldTableSize,NewTableSize)
global K

NetDir=fullfile(K.dir.tables,sprintf('m%u',ModelRank),sprintf('m%u_data',ModelRank));
cd(NetDir)
OldTableNb=ceil(K.chipSet.probeSetNbs{ModelRank}/OldTableSize);
NewTableNb=ceil(K.chipSet.probeSetNbs{ModelRank}/NewTableSize);
if NewTableNb<OldTableNb
    if mod(NewTableSize,OldTableSize)>0
        h=errordlg('NewTableSize must be a multiple of OldTableSize');
        waitfor(h)
        error('process canceled')
    else
        for TypeL=1:2            
            MergedNb=NewTableSize/OldTableSize;
            for NewTableL=1:NewTableNb
                Table=[];
                for OldTableL=((NewTableL-1)*MergedNb)+1:min(NewTableL*MergedNb,OldTableNb)
                    load(sprintf('m%u_%u_%u',ModelRank,OldTableL,TypeL));
                    if TypeL==1
                        Table=[Table;Table1];
                    else
                        Table=[Table;Table2];
                    end
                end
                if TypeL==1
                    Table1=Table;
                else
                    Table2=Table;
                end
                eval(sprintf('save m%u_%u_%u Table%u',ModelRank,NewTableL,TypeL,TypeL));
                %erase not used tables
            end
            for TableL=NewTableNb+1:OldTableNb
                delete(sprintf('m%u_%u_%u.mat',ModelRank,TableL,TypeL))
            end            
        end
    end
end
