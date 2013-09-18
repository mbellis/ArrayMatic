%not used

cd('/home/mbellis/sosma/data/cacn/')
GoBranch={'bp','mf','cc'};

%import GoTermFinder files
for TypeL=1:3
        sprintf('gotermfinder_%s_mm_only.txt',GoBranch{TypeL})
        [Go{TypeL},Def{TypeL},Pv1{TypeL},Pv2{TypeL},Found{TypeL},ListSize{TypeL},TotalFound{TypeL},TotalSize{TypeL},Fdr{TypeL},FalsePos{TypeL},GeneList{TypeL}]=textread(sprintf('gotermfinder_%s_mm_only.txt',GoBranch{TypeL}),'%s%s%s%s%u%u%u%u%s%s%s','delimiter','\t','headerlines',12,'bufsize',30000);
end

% %Transform GO Ids into numbers
% MemGo=Go;
% clear Go
% for TypeL=1:3
%     for FileL=1:6
%         Go{TypeL}=zeros(length(MemGo{TypeL}),1);
%         for GoL=1:length(MemGo{TypeL})
%             try
%             Go{TypeL}(GoL)=str2num(MemGo{TypeL}{GoL}(4:end));
%             catch
%             end
%         end
%     end
% end

%Keep power of Pv1
MemPv1=Pv1;
clear Pv1
for TypeL=1:3
    for FileL=1:6
        Pv1{TypeL}=zeros(length(MemPv1{TypeL}),1);
        for PvL=1:length(MemPv1{TypeL})            
            Pos=findstr(MemPv1{TypeL}{PvL},'-');
            if ~isempty(Pos)
            Pv1{TypeL}(PvL)=str2num(MemPv1{TypeL}{PvL}(Pos+1:end));
            else
                Pv1{TypeL}(PvL)=round(-log10(str2num(MemPv1{TypeL}{PvL})));
            end
            
        end
    end
end

for TypeL=1:3
[temp,SortOrder]=sort(Found{TypeL});
SDef{TypeL}=Def{TypeL}(flipud(SortOrder));
end