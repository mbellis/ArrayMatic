cd(K.dir.geoMetadata)
K.tmp.geoDbId=mksqlite(0,'open','GEOmetadb.sqlite');
AffyGpls=mksqlite('SELECT * FROM GPL WHERE manufacturer="Affymetrix"');
GplNb=length(AffyGpls);
Gpl=cell(GplNb,1);
Msg=cell(GplNb,1);
Title=cell(GplNb,1);
GsmNb=zeros(GplNb,1);
GseNb=zeros(GplNb,1);
PsNb=zeros(GplNb,1);
sprintf('processing %u Affymetrix GPL',GplNb)
CdfPos=zeros(GplNb,1);
TGsmNb=0;
CdfGsmNb=0;
TGseNb=0;
CdfGseNb=0;
Dist=ones(GplNb,1)*-1;
for GplL=1:GplNb
    if isequal('non-commercial',AffyGpls(GplL).distribution)
        Dist(GplL)=0;
    elseif isequal('custom-commercial',AffyGpls(GplL).distribution)
        Dist(GplL)=1;
    elseif isequal('commercial',AffyGpls(GplL).distribution)
        Dist(GplL)=2;
    end
end
        

for GplL=1:GplNb
    GplL
    CurrGpl=AffyGpls(GplL).gpl;
    Gpl{GplL}=CurrGpl;
    CurrTitle=mksqlite(sprintf('SELECT title FROM gpl WHERE gpl="%s"',CurrGpl));
    CurrPsNb=mksqlite(sprintf('SELECT data_row_count FROM gpl WHERE gpl="%s"',CurrGpl));
    CurrGsm=mksqlite(sprintf('SELECT gsm FROM gsm WHERE gpl="%s"',CurrGpl));    
    CurrGse=mksqlite(sprintf('SELECT gse FROM gse_gpl WHERE gpl="%s"',CurrGpl));    
    PsNb(GplL)=CurrPsNb(1).data_row_count;
    Title{GplL}=CurrTitle(1).title;
    GsmNb(GplL)=length(CurrGsm);
    GseNb(GplL)=length(CurrGse);
    if ~isempty(findstr('CDF',upper(Title{GplL})))
        CdfPos(GplL)=1;
        CdfGsmNb=CdfGsmNb+GsmNb(GplL);
        CdfGseNb=CdfGseNb+GseNb(GplL);
    else
        TGsmNb=TGsmNb+GsmNb(GplL);
        TGseNb=TGseNb+GseNb(GplL);
    end
    Msg{GplL}=sprintf('%s (%u): %s (%u ps) => %u GSM => %u Gse (%u CDF) & %u Gsm (%u CDF)',...
        Gpl{GplL},GplL,Title{GplL},PsNb(GplL),GsmNb(GplL),TGseNb,CdfGseNb,...
        TGsmNb,CdfGsmNb);
    sprintf('%s',Msg{GplL})
end
fid=fopen('Affymetrix_CDF.txt','w');
for GplL=1:GplNb
    fprintf(fid,'%u\t%s\t%s\t%u\t%u\t%u\n',CdfPos(GplL),Gpl{GplL},Title{GplL},PsNb(GplL),GseNb(GplL),GsmNb(GplL));
end
fclose(fid)