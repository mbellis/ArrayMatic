%PLOT_DISTANCES

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
function trs_plotdistances(DataType)
global P

cd(P.dir.data)
if exist('Tree.mat','file')
    load Tree
else
    h=errordlg('run prepare data for tree first');
    waitfor(h)
    error('process canceled')
end
if isequal(DataType,'points')
    TRank=1;
elseif isequal(DataType,'biological conditions')
    TRank=2;
end

Distances=squareform(T{TRank}.distances);
XTick=find(diff(P.point.expRank)==1);
%determine the range of signals
Factor=round(median(T{TRank}.distances)/50);
%all point
h=figure;
set(gcf,'color',[1,1,1])
image(Distances/Factor);
if TRank==1
    set(h,'name','POINT DISTANCES')
    title('POINT DISTANCES');
else
    set(h,'name','CONDITIONS DISTANCES')
    title('CONDITIONS DISTANCES');
end
set(gca,'tickdir','out')
set(gca,'xtick',XTick)
set(gca,'xticklabel','')
%diagonal by bloc of 500 points
set_figsize('960spx')
cd(P.dir.resTree)
if TRank==1
    eval(sprintf('saveas(h,''D%u_all_%ups_%s.png'',''png'')',TRank,length(T{TRank}.psIndex),date))
elseif TRank==2
    eval(sprintf('saveas(h,''D%u_all_%ups_%s.png'',''png'')',TRank,length(T{TRank}.biolIndex),date))
end

BlocNb=ceil(P.point.nb/500);
MemPos=0;
for BlocL=1:BlocNb
    if BlocL<BlocNb
        Pos=find(XTick>BlocL*500);
        Pos=Pos(1);
        set(gca,'xlim',[MemPos,XTick(Pos)])
        set(gca,'ylim',[MemPos,XTick(Pos)])
        MemPos=XTick(Pos);
    else
        Pos=P.point.nb;
        set(gca,'xlim',[MemPos,Pos])
        set(gca,'ylim',[MemPos,Pos])

    end

    cd(P.dir.resTree)
    if TRank==1
        eval(sprintf('saveas(h,''D%u_%01u_%ups_%s.png'',''png'')',TRank,BlocL,length(T{TRank}.psIndex),date))
    elseif TRank==2
        eval(sprintf('saveas(h,''D%u_%01u_%ups_%s.png'',''png'')',TRank,BlocL,length(T{TRank}.biolIndex),date))
    end
end


for TRank=1:2
    h=figure;
    set(gcf,'color',[1,1,1])
    set(h,'name','DISTANCE DISTRIBUTION')
    if TRank==1
        set(h,'name',sprintf('DISTANCE DISTRIBUTION (%u  PROBE SETS - PV < %.02f',length(T{1}.psIndex),T{1}.parPv))
        hist(T{TRank}.distances,100)
        title(sprintf('DISTANCE DISTRIBUTION (%u  PROBE SETS - PV < %.02f',length(T{1}.psIndex),T{1}.parPv));
    elseif TRank==2
        set(h,'name',sprintf('DISTANCE DISTRIBUTION (%u  PROBE SETS - CV > %.02f',length(T{2}.biolIndex),T{2}.parCv))
        hist(T{TRank}.distances,100)
        title(sprintf('DISTANCE DISTRIBUTION (%u  PROBE SETS - CV > %.02f',length(T{2}.biolIndex),T{2}.parCv));
    end
    xlabel('distance')
    ylabel('frequency')
end



