TYPE='ALL';
COMPS=[9:16];
cd(P.dir.data)
load Comp
% Col=[];
% for i=1:length(M{1}.compName)
%     if ~isempty(findstr(M{1}.compName{i},['_',TYPE]))
%         Col=[Col;i];
%     end
% end

Fdr=load_data('Fdr_01.float32le',P.dir.data,P.chip.currProbeSetNb,P.chip.currProbeSetNb,'single','ieee-le',1:P.chip.currProbeSetNb,COMPS);
ZVar=load_data('ZVar_01.float32le',P.dir.data,P.chip.currProbeSetNb,P.chip.currProbeSetNb,'single','ieee-le',1:P.chip.currProbeSetNb,COMPS);
Sensitivity=load_data('Sensitivity_01.float32le',P.dir.data,P.chip.currProbeSetNb,P.chip.currProbeSetNb,'single','ieee-le',1:P.chip.currProbeSetNb,COMPS);
FC=load_data('Fc_01.float32le',P.dir.data,P.chip.currProbeSetNb,P.chip.currProbeSetNb,'single','ieee-le',1:P.chip.currProbeSetNb,COMPS);
load TotalVar_01
TotalVar.inc=TotalVar.inc(COMPS);
TotalVar.dec=TotalVar.dec(COMPS);
cd(P.dir.results)
%fid=fopen(sprintf('aorteHP_%s_results_%s.txt',lower(TYPE),date),'w');
fid=fopen(sprintf('DRG_crush_results_%s.txt',date),'w');

CompNames=M{1}.compName(COMPS);

Header='comparisons:';
CompNb=length(CompNames);
for i=1:CompNb
Header=sprintf('%s\t%s\t\t\t',Header,CompNames{i});
end
fprintf(fid,'%s\n',Header);

Header='total inc:';
for CompL=1:CompNb
    Header=sprintf('%s\t%u\t\t\t',Header,round(TotalVar.inc(CompL)));
end
fprintf(fid,'%s\n',Header);

Header='total dec:';
for CompL=1:CompNb
    Header=sprintf('%s\t%u\t\t\t',Header,round(TotalVar.dec(CompL)));
end
fprintf(fid,'%s\n',Header);

Header='probe set';
for i=1:CompNb
Header=sprintf('%s\tZVar\tFDR\tSens\tFC',Header);
end
fprintf(fid,'%s\n',Header);


for i=1:length(Fdr)
    fprintf(fid,'%s',P.chip.probeSetIds{i});
    for j=1:CompNb
        fprintf(fid,'\t%.3f\t%.3f',ZVar(i,j),Fdr(i,j),Sensitivity(i,j),FC(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid)
        