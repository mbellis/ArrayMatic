%TRS_NEWNETWORK

%INPUT PARAMETERS

%EXTERNAL FILES

%OUTPUT PARAMETERS


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
function trs_newnetwork()
KNetMem=K.net;
ModelRank=P.chipRank;
Continue=0;
Cancel=0;
Default='All';
NewNet=0;
while Continue==0    
    Name=inputdlg({'Give a short name for the network file'},'',1,{sprintf('NR(%u) FDR(1pm) 2010',length(P.biol.grp.grpBiolRanks{1}))});
    Name=Name{1};
    if ~isempty(findstr(upper(Name),'_PERM'))
        h=warndlg('do not use "_PERM" in the net name (added automatically for randomized data)');
        waitfor(h)
    else
        if length(K.net)>=ModelRank
            if isfield(K.net{ModelRank},'name')
                if ~isempty(strmatch(Name,K.net{ModelRank}.name,'exact'))
                    Answer=questdlg([Name,' already exist. Do you want to overwrite it ?'],'','Yes','No','No');
                    if isequal(Answer,'Yes')
                        NetPos=strmatch(Name,K.net{ModelRank}.name,'exact');
                        Continue=1;
                    else
                        Answer=questdlg('Do you want to cancel ?','','Yes','No','No');
                        if isequal(Answer,'Yes')
                            Cancel=1;
                            Continue=1;
                        else
                            Default=Name;
                        end
                    end
                else
                    NewNet=1;
                    Continue=1;
                end
            else
                NewNet=1;
                Continue=1;
            end
        else
            NewNet=1;
            Continue=1;
        end
    end
end
if Cancel~=1
    if NewNet==1
        if length(K.net)>=ModelRank
            if isfield(K.net{ModelRank},'name')
                K.net{ModelRank}.nb=length(K.net{ModelRank}.name)+1;
            else
                K.net{ModelRank}.nb=1;
            end
        else
            K.net{ModelRank}.nb=1;
        end
        NetPos=K.net{ModelRank}.nb;
    end
    K.net{ModelRank}.name{NetPos,1}=Name;
    K.net{ModelRank}.netMade(NetPos,1)=0;

    FLimit=ONELIMIT('Select the first FDR value (0< FDR <=1)','0.01');
    SLimit=ONELIMIT('Select the first sensitivity value (0< S <=1)','1');
    K.net{ModelRank}.fdr(NetPos,1)=FLimit;
    K.net{ModelRank}.s(NetPos,1)=SLimit;
    SndFdr=questdlg('Do you want to use a second fdr ?','','yes','no','yes');
    if isequal(SndFdr,'yes')
        FLimit=ONELIMIT('Select the second FDR value (0< FDR <=1)','0.10');
        SLimit=ONELIMIT('Select the second sensitivity value (0< S <=1)','1');
        K.net{ModelRank}.fdr(NetPos,2)=FLimit;
        K.net{ModelRank}.s(NetPos,2)=SLimit;
    else
        K.net{ModelRank}.fdr(NetPos,2)=0;
        K.net{ModelRank}.s(NetPos,2)=0;
    end


    % biol conditions must be indicated as absolute ranks. First table has an index of biol cond really used, which allows to
    % transform absolute ranks in relative positions
    SelNb=questdlg('Select the number of biological groups','','1','2','1')
    SelNb=str2num(SelNb);
    SelType=questdlg('Type of selection','COM_RAY','manual selection','list selection','load a file','list selection')
    for SelL=1:SelNb
        Continue=1;
        switch  SelType
            case 'manual selection'
                %K.net{ModelRank}.refchiprank{NetPos,1}(1)=RefChipRank;
                [BiolRank,BiolPos]=SELECTBIOL('select biological conditions',1)
                % Rank are used because user can change biol cond used for searching scouples
                K.net{ModelRank}.biolRank{NetPos,1}{SelL}=sort(unique(BiolRank));
            case 'list selection'
                [Sel Ok]=listdlg('promptstring','select a list of biological conditions','selectionmode','multiple','liststring',P.biol.grp{ChipRank}.name);
                if Ok==1
                    K.net{ModelRank}.biolRank{NetPos,1}{SelL}=P.biol.grp{ChipRank}.grpBiolRanks{Sel(1)};
                    if length(Sel)>1
                        for i=2:length(Sel)
                            K.net{ModelRank}.biolRank{NetPos,1}{SelL}=[K.net{ModelRank}.biolRank{NetPos,1}{SelL};P.biol.grp{ChipRank}.grpBiolRanks{Sel(i)}];
                        end
                    end
                    K.net{ModelRank}.biolRank{NetPos,1}{SelL}=sort(unique(K.net{ModelRank}.biolRank{NetPos,1}{SelL}));
                else
                    Continue=0;
                    h=warndlg('process canceled');
                    waitfor(h)
                end
            case 'load a file'
                [MatFile,MatFileDir]=uigetfile('*.txt*','Select the biological condition index file ');
                eval(['cd ',MatFileDir])
                BiolIndex=textread(MatFile,'%n');
                K.net{ModelRank}.biolRank{NetPos,1}{SelL}=sort(unique(BiolIndex));
        end
    end
    if SelNb==2
        %eliminate common biol cond
        Biol1=K.net{ModelRank}.biolRank{NetPos,1}{1};
        Biol2=K.net{ModelRank}.biolRank{NetPos,1}{2};
        Clear1=[];
        Clear2=[];
        for BiolL=1:length(Biol1)
            Pos=find(Biol2==Biol1(BiolL));
            if ~isempty(Pos)
                Clear1=[Clear1;BiolL];
                Clear2=[Clear2;Pos];
            end
        end
        if length(Clear1)>0
            h=warndlg(sprintf('%u common biol cond eliminated',length(unique(Biol1(Clear1)))));
            waitfor(h)
            Biol1(Clear1)=[];
            K.net{ModelRank}.biolRank{NetPos,1}{1}=sort(unique(Biol1));
            Biol2(Clear2)=[];
            K.net{ModelRank}.biolRank{NetPos,1}{2}=sort(unique(Biol2));
        end
        BiolNb1=length(K.net{ModelRank}.biolRank{NetPos,1}{1});
        BiolNb2=length(K.net{ModelRank}.biolRank{NetPos,1}{2});
        K.net{ModelRank}.compNb(NetPos,1)=BiolNb1*BiolNb2;
    else
        BiolNb=length(K.net{ModelRank}.biolRank{NetPos,1}{1});
        K.net{ModelRank}.compNb(NetPos,1)=BiolNb*(BiolNb-1)/2;
    end
    if Continue==1
        % while DONETFASTA is not run genenb and Blocnb are not known
        K.net{ModelRank}.blocNb(NetPos,1)=0;
        %find the first free rank
        if NetPos==1
            K.net{ModelRank}.rank(NetPos,1)=1;
        else
            RankList=K.net{ModelRank}.rank;
            MaxRank=max(RankList);
            CurrRank=MaxRank+1;
            for RankL=1:MaxRank
                if isempty(find(RankList==RankL))
                    CurrRank=RankL;
                    break
                end
            end
            K.net{ModelRank}.rank(NetPos,1)=CurrRank;
        end
        Comment=inputdlg({'write a comment'},'',1,{''});
        K.net{ModelRank}.comment{NetPos,1}=Comment{1};
    else
        K.net=KNetMem;
    end
end