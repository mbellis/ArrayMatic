function rdn_v2_riz()

%IMPORT AFFY INFORMATION FILES
RootDir='/home/mbellis/sosma';

cd(fullfile(RooDir,'data','affy','Rice','LibFiles'))
cdf = mt_readcdf('Rice.cdf',1);
cd(fullfile(RootDir,'data','affy','Rice',''))
save Rice_cdf cdf

cd(fullfile(RooDir,'data','affy','Rice','LibFiles'))
gin = mt_readgin('Rice.gin',1);
cd(fullfile(RootDir,'data','affy','Rice',''))
save Rice_gin gin

cd(fullfile(RooDir,'data','affy','Rice','LibFiles'))
pa = mt_readprobe_annot('Rice.probe_tab');
cd(fullfile(RootDir,'data','affy','Rice',''))
save Rice_pa pa

cd(fullfile(RootDir,'data','affy','Rice','LibFiles'))
seq = mt_readseq('Rice.target',1);
cd(fullfile(RootDir,'data','affy','Rice',''))
save Rice_seq seq


%IMPORT RAW CEL DATA

cd(fullfile(RootDir,'raydataraw','2010_BELLAFIORE','raw_data','cel_files'))
cel=mt_readcel_mb({'00hR1.CEL','00hR2.CEL','00hR3.CEL','00hR4.CEL','00hR5.CEL','00hR6.CEL','48hR1.CEL','48hR2.CEL','48hR3.CEL','96hR1.CEL','96hR2.CEL','96hR3.CEL','06dR4.CEL','06dR5.CEL','06dR6.CEL'},1);
save RIZ2010_cel cel


%PROBE ANNOTATION
cd(fullfile(RootDir,'data','affy','Rice','LibFiles'))
load Rice_cdf
load Rice_gin
load Rice_seq
load Rice_pa
cd(fullfile(RootDir,'raydataraw','2010_BELLAFIORE','raw_data','cel_files'))
load RIZ2010_cel
probes = mt_cel2probes(cel,cdf,gin,pa,seq);
cd(fullfile(RootDir,'raydataraw','2010_BELLAFIORE','raw_data','cel_files'))
save RIZ2010_probes_start probes

probes= mt_normalize(probes);
cd(fullfile(RootDir,'raydataraw','2010_BELLAFIORE','raw_data','cel_files'))
save RIZ2010_probes_end probes
