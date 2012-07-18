%% COMMON DIRECTORIES
if isunix
    K.dir.common=fullfile('/home','mbellis','array2','mbellis','sosma','arraymatic','amcdata','common');
    K.dir.test=fullfile('/home','mbellis','array2','mbellis','sosma','arraymatic','amcdata','test');
    K.dir.affyChipData=fullfile('/usr','data','chips','affymetrix');
    K.dir.geoExperiments=fullfile('/usr','data','experiments','geo');
    K.dir.metadata=fullfile('/home','mbellis','array2','mbellis','sosma','arraymatic','amcdata','metadata');
    K.dir.geoMetadata=fullfile('/home','mbellis','array2','mbellis','sosma','arraymatic','amcdata','metadata','geo');
    K.dir.affyMetadata=fullfile('/home','mbellis','array2','mbellis','sosma','arraymatic','amcdata','metadata','affy');
    K.dir.chipMetadata=fullfile('/home','mbellis','array2','mbellis','sosma','arraymatic','amcdata','metadata','chip');
    K.dir.cliques=fullfile('/home','mbellis','array2','mbellis','sosma','arraymatic','amcdata','cliques');
    K.dir.chip=fullfile('/home','mbellis','array2','mbellis','sosma','arraymatic','amcdata','chip');
    K.dir.gene=fullfile('/home','mbellis','array2','mbellis','sosma','arraymatic','amcdata','gene');
    K.dir.list=fullfile('/home','mbellis','array2','mbellis','sosma','arraymatic','amcdata','list');
    K.dir.net=fullfile('/home','mbellis','array1','sosma','net');
    K.dir.tables=fullfile('/usr','data','net');
    K.dir.amcResults=fullfile('/home','mbellis','sosma','amcres');

else
    K.dir.common=fullfile('X:','mbellis','sosma','arraymatic','amcdata','common');
    K.dir.test=fullfile('X:','mbellis','sosma','arraymatic','amcdata','test');
    K.dir.affyChipData=fullfile('E:','sosma','chips','affymetrix')
    K.dir.geoExperiments=fullfile('E:','sosma','data','experiments','geo');
    K.dir.metadata=fullfile('X:','mbellis','sosma','arraymatic','amcdata','metadata');
    K.dir.geoMetadata=fullfile('X:','mbellis','sosma','arraymatic','amcdata','metadata','geo');
    K.dir.affyMetadata=fullfile('X:','mbellis','sosma','arraymatic','amcdata','metadata','affy');
    K.dir.chipMetadata=fullfile('X:','mbellis','sosma','arraymatic','amcdata','metadata','chip');
    K.dir.cliques=fullfile('X:','mbellis','sosma','arraymatic','amcdata','cliques');
    K.dir.chip=fullfile('X:','mbellis','sosma','arraymatic','amcdata','chip');
    K.dir.gene=fullfile('X:','mbellis','sosma','arraymatic','amcdata','gene');
    K.dir.gene=fullfile('X:','mbellis','sosma','arraymatic','amcdata','list');
    K.dir.net=fullfile('E:','mbellis','net');
    K.dir.tables=fullfile('E:','mbellis','tables');
    K.dir.amcResults=fullfile('E:','mbellis','sosma','amcres');
end


