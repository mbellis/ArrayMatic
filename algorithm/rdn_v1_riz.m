function rdn(action)

switch action
    case 'import information'
        %IMPORT AFFY INFORMATION FILES
        RootDir=CLEAR_ALL;

        cd(fullfile(RooDir,'data','affy','Rice','CD_Rice','Full','Rice','LibFiles'))
        cdf = mt_readcdf('Rice.cdf',1);
        cd(fullfile(RootDir,'data','affy','Rice',''))
        save m060(rice)_cdf cdf
        cd(fullfile(RootDir,'data','affy','Rice','CD_Rice','Full','Rice','LibFiles'))
        gin = mt_readgin('Rice.gin',1);
        cd(fullfile(RootDir,'data','affy','Rice',''))
        save m060(rice)_gin gin
        cd(fullfile(RootDir,'data','affy','Rice',''))
        probeannot = mt_readprobe_annot('Rice.probe_tab');
        save m060(rice)_probeannot probeannot
        cd(fullfile(RootDir,'data','affy','Rice'))
        seq = mt_readseq('Rice.target',1);
        save m060(rice)_seq seq

    case 'read cel'
        %IMPORT RAW CEL DATA
        RootDir=CLEAR_ALL;
        cd(fullfile(RootDir,'raydataraw','2010_BELLAFIORE','raw_data','cel_files'))
        cel=mt_readcel_mb({'00hR1.CEL','00hR2.CEL','00hR3.CEL','00hR4.CEL','00hR5.CEL','00hR6.CEL','48hR1.CEL','48hR2.CEL','48hR3.CEL','96hR1.CEL','96hR2.CEL','96hR3.CEL','06dR4.CEL','06dR5.CEL','06dR6.CEL'},1);
        save RIZ2010_cel cel

    case 'read probes'
        %PROBE ANNOTATION
        cd(fullfile(RootDir,'data','affy','rice',''))
        load m060(rice)_cdf
        load m060(rice)_gin
        load m060(rice)_seq
        load m060(rice)_probeannot
        cd(fullfile(RootDir,'raydataraw','2010_BELLAFIORE','raw_data','cel_files'))
        load RIZ2010_cel
        probes = mt_cel2probes_mb(cel,cdf,gin,probeannot,seq,0);
        probes.pm=single(probes.pm);
        probes.mm=single(probes.mm);
        probes.ind=single(probes.ind);
        probes.indmm=single(probes.indmm);
        probes.position=single(probes.position);
        probes.pm_pos=single(probes.pm_pos);
        probes.mm_pos=single(probes.mm_pos);
        cd(fullfile(RootDir,'raydataraw','2010_BELLAFIORE','raw_data','cel_files'))
        save RIZ2010_probes_new probes

    case 'estimate background'
        %ESTIMATE OPTICAL & BACKGROUND SIGNAL
        cd(fullfile(RootDir,'raydataraw','2010_BELLAFIORE','raw_data','cel_files'))
        load RIZ2010_probes_new
        nprobe = size(probes.pm,2);
        narray = size(probes.pm,1);
        ngene = length(probes.name);
        pm=probes.pm;
        sequence=probes.sequence;
        qq_factors=[];
        image_factors=[];
        seq_correction=[];
        amp_correction=[];
        if isfield(probes,'qq_factors')
            qq_factors=probes.qq_factors;
        end
        if isfield(probes,'image_factors')
            image_factors=probes.image_factors;
        end
        if isfield(probes,'seq_correction')
            seq_correction=probes.seq_correction;
        end
        if isfield(probes,'ampt_correction')
            amp_correction=probes.amp_correction;
        end
        %nprobes = mt_bg_est(probes)
        clear probes
        [seqbg_factors,seqbg] = mt_bg_est_mb(nprobe,narray,ngene,qq_factors,image_factors,seq_correction,amp_correction);
        cd(fullfile(RootDir,'raydataraw','2010_BELLAFIORE','raw_data','cel_files'))
        load RIZ2010_probes_new
        probes.seqbg_factors=seqbg_factors;
        probes.seqbg=seqbg;
        clear seqbg seqbg_factors
        save RIZ2010_probes_new probes

    case 'estimate model'
        %ESTIMATE SEQUENCE MODEL PARAMETERS
        cd(fullfile(RootDir,'raydataraw','2010_BELLAFIORE','raw_data','cel_files'))
        load RIZ2010_probes_new
        nprobe = size(probes.pm,2);
        narray = size(probes.pm,1);
        ngene = length(probes.name);
        pm=probes.pm;
        sequence=probes.sequence;
        qq_factors=[];
        image_factors=[];
        seq_correction=[];
        amp_correction=[];
        seqbg=probes.seqbg;
        seqbg_factors=probes.seqbg_factors;
        ind=probes.ind;
        if isfield(probes,'qq_factors')
            qq_factors=probes.qq_factors;
        end
        if isfield(probes,'image_factors')
            image_factors=probes.image_factors;
        end
        if isfield(probes,'seq_correction')
            seq_correction=probes.seq_correction;
        end
        if isfield(probes,'ampt_correction')
            amp_correction=probes.amp_correction;
        end
        %with amplification processing (does not work)
        %gene_sequence=probes.gene_sequence;
        %clear probes
        %J=[];
        %[seq_norm,seq_factors,seq_correction,amp_correction] = mt_cor_hybamp_mb(nprobe,narray,ngene,qq_factors,image_factors,seq_correction,amp_correction,seqbg,seqbg_factors,ind,'m_estimation','use_amplification');
        %without amplification processing
        clear probes
        J=[];
        [seq_norm,seq_factors,seq_correction,amp_correction] = mt_cor_hybamp_mb(nprobe,narray,ngene,qq_factors,image_factors,seq_correction,amp_correction,seqbg,seqbg_factors,ind,'m_estimation');
        cd(fullfile(RootDir,'raydataraw','2010_BELLAFIORE','raw_data','cel_files'))
        load RIZ2010_probes_new
        probes.seq_norm=seq_norm;
        probes.seq_factors=seq_factors;
        probes.seq_correction=seq_correction;
        probes.amp_correction=amp_correction;
        J=[];
        sequence={};
        pm=[];
        gene_sequence={};
        save RIZ2010_probes_new probes
        clear seqbg gene_ sequence ind seqbg seqbg_factors
        clear signal se_norm seq_factors seq_correction amp_correction
        clear probes

    case 'correct image'
        %CORRECT IMAGE
        cd(fullfile(RootDir,'raydataraw','2010_BELLAFIORE','raw_data','cel_files'))
        load RIZ2010_probes_new
        nprobe = size(probes.pm,2);
        narray = size(probes.pm,1);
        ngene = length(probes.name);
        qq_factors=[];
        image_factors=[];
        seq_correction=[];
        amp_correction=[];
        if isfield(probes,'qq_factors')
            qq_factors=probes.qq_factors;
        end
        if isfield(probes,'image_factors')
            image_factors=probes.image_factors;
        end
        if isfield(probes,'seq_correction')
            seq_correction=probes.seq_correction;
        end
        if isfield(probes,'ampt_correction')
            amp_correction=probes.amp_correction;
        end
        pm=probes.pm;
        signal = mt_real_signal_mb(narray,nprobe,qq_factors,image_factors,seq_correction,amp_correction);
        pm=[];

        pmsize=size(probes.pm);
        seqbg_factors=probes.seqbg_factors;
        chipnrows=probes.nrows;
        chipncols=probes.ncols;
        pm_pos=probes.pm_pos;
        clear probes

        image_factors = mt_cor_image_mb(pmsize,signal,seqbg_factors,chipnrows,chipncols,pm_pos,9,'use_optical');
        cd(fullfile(RootDir,'raydataraw','2010_BELLAFIORE','raw_data','cel_files'))
        load RIZ2010_probes_new
        probes.image_factors=image_factors;
        save RIZ2010_probes_new probes

    case 'quantile normalization'
        %QUANTILE NORMALIZATON
        <dataset> = mt_cor_qq(<dataset>

    case 'summarize data'
        %CORRECT FOR ARRAY LOCATION EFFECTS & CALCULATE SUMMARIZED DATA
        cd(fullfile(RootDir,'raydataraw','2010_BELLAFIORE','raw_data','cel_files'))
        load RIZ2010_probes_new
        nprobe = size(probes.pm,2);
        narray = size(probes.pm,1);
        ngene = length(probes.name);

        lab_idx={};
        if(isfield(probes,'lab'))
            labels = unique(probes.lab);
            nlab = length(labels);
            lab_idx = {};
            for i = 1:nlab
                lab_idx{i} = find(probes.lab == labels(i));
            end;
        end;
        ind=probes.ind;
        pmsize=size(probes.pm);
        seqbg=probes.seqbg;
        qq_factors=[];
        image_factors=[];
        seq_correction=[];
        amp_correction=[];
        if isfield(probes,'qq_factors')
            qq_factors=probes.qq_factors;
        end
        if isfield(probes,'image_factors')
            image_factors=probes.image_factors;
        end
        if isfield(probes,'seq_correction')
            seq_correction=probes.seq_correction;
        end
        if isfield(probes,'ampt_correction')
            amp_correction=probes.amp_correction;
        end
        pm=probes.pm;
        clear probes
        signal = mt_real_signal_mb(narray,nprobe,qq_factors,image_factors,seq_correction,amp_correction);
        pm=[];
        r=mt_sum_plm_mb(nprobe,narray,ngene,lab_idx,ind,pmsize,seqbg,signal,'m_estimation','keep_probe_info');
        cd(fullfile(RootDir,'raydataraw','2010_BELLAFIORE','raw_data','cel_files'))
        load RIZ2010_probes_new
        probes.probe_factors=single(r.probe_factors);
        probes.removed_factors=single(r.removed_factors);
        probes.array_factors=single(r.array_factors);
        probes.overall_factors=single(r.overall_factors);
        probes.resids=single(r.resids);
        save RIZ2010_probes_new probes




        function RootDir=CLEAR_ALL()
        clear all
        RootDir=fullfile('/home','mbellis','sosma');


