#!/usr/bin/env snakemake

"""
Snakemake PoPs pipeline
----------------------

Snakemake pipeline for conducting PoPs on munged summary statistics
"""

import glob
import numpy as np
import snakemake.utils as snk
snk.min_version('7.0')

module_prefix = 'results'
n_batch = config['empirical_pops']['n_batch']
if not config['empirical_pops']['batch_iterations']:
    n_batch = config['empirical_pops']['n_iterations']
n_per_batch = int(np.ceil(config['empirical_pops']['n_iterations'] / n_batch))

################################################################################
# FUNCTIONS
################################################################################
def build_pops_options(
    gene_locations,
    feature_input,
    magma_input,
    chr_curr,
    control_features,
    pval_threshold,
    keep_hla_genes,
    out_prefix
):
    cmd = ''
    # Build inputs
    cmd += '--gene_annot_path {} '.format(gene_locations)
    cmd += '--feature_mat_prefix {} '.format(feature_input)
    cmd += '--magma_prefix {} '.format(magma_input)

    # Identify # feature chunks
    n_feats = len(glob.glob('{}.cols.*.txt'.format(feature_input)))
    cmd += '--num_feature_chunks {} '.format(n_feats)

    # Chromosome specification
    # Train on leave-one-out
    if chr_curr != 'all':
        train_chr = config['chromosomes'].copy()
        train_chr.remove(chr_curr)
        train_str = ' '.join(train_chr)
        cmd += '--project_out_covariates_chromosomes {} '.format(train_str)
        cmd += '--feature_selection_chromosomes {} '.format(train_str)
        cmd += '--training_chromosomes {} '.format(train_str)

    # Check for control features
    if control_features != '':
        cmd += '--control_features_path {} '.format(control_features)

    # Now specify settings
    cmd += '--use_magma_covariates '
    cmd += '--use_magma_error_cov '
    cmd += '--feature_selection_p_cutoff {} '.format(str(pval_threshold))
    cmd += '--method ridge '
    cmd += '--save_matrix_files '

    if keep_hla_genes:
        cmd += '--project_out_covariates_keep_hla '
        cmd += '--feature_selection_keep_hla '
        cmd += '--training_keep_hla '
    else:
        cmd += '--project_out_covariates_remove_hla '
        cmd += '--feature_selection_remove_hla '
        cmd += '--training_remove_hla '

    cmd += '--out_prefix {}'.format(out_prefix)
    return(cmd)

################################################################################

def get_input(wildcards):
    input_list = []
    input_list.extend(expand(
        [
            '{pref}/{study}/empirical_pops_pvals.tsv',
            '{pref}/{study}/empirical_pops_null_dist.tsv',
            '{pref}/{study}/empirical_pops_pvals.png',
            '{pref}/{study}/plots/pops_distributions_density.png',
            '{pref}/{study}/plots/pops_distributions_qq.png',
            '{pref}/{study}/plots/pops_distributions_ecdf.png',

            # cleaned files
            '{pref}/{study}/null_pops/clean.done',
            '{pref}/features__null_swarm/clean.done'
        ],
        pref = module_prefix,
        study = list(config['study'].keys())
    ))

    input_list.append('{}/merged_study_results.tsv.gz'.format(module_prefix))
    return(input_list)

rule all:
    input:
        get_input

def get_study_sumstat(wildcards):
    return([config['study'][wildcards.study]])

rule prepare_study_sumstats:
    input:
        get_study_sumstat
    output:
        '{}/{{study}}/sumstats.tsv'.format(module_prefix)
    shell:
        'zcat < {input} > {output}'

rule prepare_features:
    input:
        config['features']
    output:
        '{}/features_raw/feats.done'.format(module_prefix)
    params:
        out_prefix = '{}/features_raw'.format(module_prefix)
    shell:
        'for file in {input}; do '
        '    cp $file {params.out_prefix}/.; '
        'done; '
        'touch {output}'

################################################################################
# Step 1: Run MAGMA
################################################################################
rule magma__get_gene_locations:
    input:
        config['gene_locations']
    output:
        '{}/magma_gene_location.tsv'.format(module_prefix)
    shell:
        'cat {input} | cut -f1,3-5 > {output}'

rule magma__annotate_genes:
    input:
        varloc = '{}.bim'.format(config['magma']['score_genes']['ld_reference']),
        gloc = '{}/magma_gene_location.tsv'.format(module_prefix)
    output:
        out = '{}/genwin_{}_{}kb.genes.annot'.format(
            module_prefix,
            str(config['magma']['annotate']['window_upstream_kb']).replace('.', 'pt'),
            str(config['magma']['annotate']['window_downstream_kb']).replace('.', 'pt')
        )
    params:
        window = '{},{}'.format(
            str(config['magma']['annotate']['window_upstream_kb']),
            str(config['magma']['annotate']['window_downstream_kb'])
        ),
        out_prefix = '{}/genwin_{}_{}kb'.format(
            module_prefix,
            str(config['magma']['annotate']['window_upstream_kb']).replace('.', 'pt'),
            str(config['magma']['annotate']['window_downstream_kb']).replace('.', 'pt')
        )
    shell:
        'magma '
        '--annotate window={params.window} '
        '--snp-loc {input.varloc} '
        '--gene-loc {input.gloc} '
        '--out {params.out_prefix}'

def get_magma_run_input(wildcards):
    in_dict = {
        'ss': '{}/{}/sumstats.tsv'.format(module_prefix, wildcards.study),
        'ld_ref_bed': '{}.bed'.format(config['magma']['score_genes']['ld_reference']),
        'ld_ref_bim': '{}.bim'.format(config['magma']['score_genes']['ld_reference']),
        'ld_ref_fam': '{}.fam'.format(config['magma']['score_genes']['ld_reference'])
    }
    if config['magma']['annotate']['run_process']:
        in_dict['annot'] = '{}/genwin_{}_{}kb.genes.annot'.format(
            module_prefix,
            str(config['magma']['annotate']['window_upstream_kb']),
            str(config['magma']['annotate']['window_downstream_kb'])
        )
    else:
        in_dict['annot'] = 'data/magma.genes.annot'
    
    return(in_dict)

rule magma__run:
    """
    Run MAGMA
    """
    input:
        unpack(get_magma_run_input)
    output:
        raw = '{}/{{study}}/magma.genes.raw'.format(module_prefix),
        out = '{}/{{study}}/magma.genes.out'.format(module_prefix)
    params:
        ld_ref = config['magma']['score_genes']['ld_reference'],
        ss_descr = 'ncol={} use={},{}'.format(
            config['magma']['score_genes']['sumstat_n_col'],
            config['magma']['score_genes']['sumstat_variant_id_col'],
            config['magma']['score_genes']['sumstat_pvalue_col']
        ),
        out_prefix = lambda wildcards: '{}/{}/magma'.format(
            module_prefix,
            wildcards.study
        )
    shell: 
        'magma '
        '--bfile {params.ld_ref} '
        '--gene-annot {input.annot} '
        '--pval {input.ss} {params.ss_descr} ' 
        '--gene-model snp-wise=mean '
        '--out {params.out_prefix}'

# ################################################################################
# # Step 2: Run nominal PoPS
# ################################################################################
rule pops__munge_features:
    """
    Select features to use for PoPS
    """
    input:
        feats = '{}/features_raw/feats.done'.format(module_prefix),
        controls = config['pops']['control_features'],
        gloc = config['gene_locations']
    output:
        feats = '{}/features_processed/feats.rows.txt'.format(module_prefix)
    params:
        feat_dir = '{}/features_raw'.format(module_prefix),
        out_pref = '{}/features_processed/feats'.format(module_prefix),
        script = srcdir('scripts/pops_v0pt2__munge_feature_directory.py')
    shell:
        '{params.script} '
        '--gene_annot_path {input.gloc} '
        '--feature_dir {params.feat_dir} '
        '--control_path {input.controls} '
        '--nan_policy raise '
        '--max_cols 5000 '
        '--save_prefix {params.out_pref}'
    
rule pops__run:
    """
    Run PoPs
    """
    input:
        magma_raw = '{}/{{study}}/magma.genes.raw'.format(module_prefix),
        magma_out = '{}/{{study}}/magma.genes.out'.format(module_prefix),
        feat_rows = '{}/features_processed/feats.rows.txt'.format(module_prefix),
        gloc = config['gene_locations']
    output:
        preds = '{}/{{study}}/swarm_chr/nominal_pops.{{chr}}.preds'.format(module_prefix),
        coefs = '{}/{{study}}/swarm_chr/nominal_pops.{{chr}}.coefs'.format(module_prefix),
        margs = '{}/{{study}}/swarm_chr/nominal_pops.{{chr}}.marginals'.format(module_prefix)
    params:
        pops_options = lambda wildcards, input: build_pops_options(
            gene_locations = input.gloc,
            feature_input = '{}/features_processed/feats'.format(module_prefix),
            magma_input = '{}/{}/magma'.format(module_prefix, wildcards.study),
            chr_curr = wildcards.chr,
            control_features = config['pops']['control_features'],
            pval_threshold = 0.5,
            keep_hla_genes = config['pops']['keep_hla_genes'],
            out_prefix = '{}/{}/swarm_chr/nominal_pops.{}'.format(
                module_prefix,
                wildcards.study,
                wildcards.chr
            )
        ),
        script = srcdir('scripts/pops_v0pt2__run.py')
    shell:
        '{params.script} {params.pops_options}'

rule pops__filter_chr:
    """
    PoPs filter chromosomes
    """
    input:
        preds = '{}/{{study}}/swarm_chr/nominal_pops.{{chr}}.preds'.format(module_prefix),
        gloc = config['gene_locations']
    output:
        preds = '{}/{{study}}/swarm_chr/nominal_pops.{{chr}}.results'.format(module_prefix)
    params:
        script = srcdir('scripts/filter_pops_chrom_results.py')
    shell:
        'if [[ "{wildcards.chr}" == "all" ]]; then '
        '   cp {input.preds} {output.preds}; '
        'else '
        '   {params.script} '
        '   --pops_results {input.preds} '
        '   --gene_locations {input.gloc} '
        '   --chromosome {wildcards.chr} '
        '   --output_file {output.preds}; '
        'fi'

def get_pops_gather_input(wildcards):
    in_list = [ '{}/{{study}}/swarm_chr/nominal_pops.all.results'.format(module_prefix) ]
    if config['loco']:
        in_list = expand(
            '{pref}/{{study}}/swarm_chr/nominal_pops.{chr}.results',
            pref = module_prefix,
            chr = [ str(x) for x in config['chromosomes'] ]
        )
    return(in_list)

rule pops__gather:
    """
    Gather nominal PoPs results
    """
    input:
        get_pops_gather_input
    output:
        '{}/{{study}}/nominal_pops.tsv'.format(module_prefix)
    run:
        if len(input) == 1:
            shell('set +o pipefail; cat {input} | cut -f1,2 > {output}')
        else:
            shell((
                'set +o pipefail; '
                'cat {input[0]} | head -n 1 | cut -f1,2 > {output}; '
                'for file in {input}; do '
                '   cat $file | tail -n +2 | cut -f1,2 >> {output}; '
                'done'
            ))

################################################################################
# Step 3: Run empirical PoPS
################################################################################
rule null_pops__shuffle_features:
    """
    Shuffle features and munge
    """
    input:
        feats = '{}/features_processed/feats.rows.txt'.format(module_prefix)
    output:
        done = '{}/features__null_swarm/done.txt'.format(module_prefix)
    params:
        feat_dir = '{}/features_processed'.format(module_prefix),
        script = srcdir('scripts/null_pops__shuffle_mtx.py')
    run:
        for bat in range(1, n_batch+1):
            for i in range(1, n_per_batch+1):
                out_dir = '{}/features__null_swarm/batch{}_{}_{}'.format(
                    module_prefix,
                    bat,
                    i,
                    n_per_batch
                )
                shell('mkdir -p {out_dir}; cp {params.feat_dir}/* {out_dir}/.;')
                seed = bat*i

                # Get file - use 0 if shuffling controls, 1 if else
                init_index = 0 if config['pops']['permute_controls'] else 1
                files = ','.join([
                    '{}/feats.mat.{}.npy'.format(out_dir, x) for x in range(
                        init_index,
                        len(glob.glob('{}/feats.mat.*.npy'.format(out_dir)))
                    )
                ])
                shell((
                    '{params.script} '
                    '--feature_rows {out_dir}/feats.rows.txt '
                    '--feature_matrices {files} '
                    '--seed {seed} '
                ))
        shell('touch {output.done}')

rule null_pops__run:
    """
    Run and filter in one step for reduction in computation
    """
    input:
        magma_raw = '{}/{{study}}/magma.genes.raw'.format(module_prefix),
        magma_out = '{}/{{study}}/magma.genes.out'.format(module_prefix),
        feats = '{}/features__null_swarm/done.txt'.format(module_prefix),
        gloc = config['gene_locations']
    output:
        out = '{}/{{study}}/null_pops/swarm/batch{{bat}}/done.{{chr}}.txt'.format(module_prefix)
    params:
        pops_script = srcdir('scripts/pops_v0pt2__run.py'),
        filter_script = srcdir('scripts/filter_pops_chrom_results.py')
    run:
        for i in range(1, n_per_batch+1):
            base_dir = '{}/{}/null_pops/swarm/batch{}/perm{}_{}'.format(
                module_prefix,
                wildcards.study,
                wildcards.bat,
                i,
                n_per_batch
            )
            out_prefix = '{}/null_pops.{}'.format(base_dir, wildcards.chr)
            feat_prefix = '{}/features__null_swarm/batch{}_{}_{}'.format(
                module_prefix,
                wildcards.bat,
                i,
                n_per_batch
            )
            pops_options = build_pops_options(
                gene_locations = input.gloc,
                feature_input = '{}/feats'.format(feat_prefix),
                magma_input = '{}/{}/magma'.format(module_prefix, wildcards.study),
                chr_curr = wildcards.chr,
                control_features = config['pops']['control_features'],
                pval_threshold = 0.5,
                keep_hla_genes = config['pops']['keep_hla_genes'],
                out_prefix = out_prefix
            )
            shell((
                'mkdir -p {base_dir}; '
                '{params.pops_script} {pops_options}; '
                'if [[ "{wildcards.chr}" == "all" ]]; then '
                '   cp {out_prefix}.preds {out_prefix}.results; '
                'else '
                '   {params.filter_script} '
                '   --pops_results {out_prefix}.preds '
                '   --gene_locations {input.gloc} '
                '   --chromosome {wildcards.chr} '
                '   --output_file {out_prefix}.results; '
                'fi' 
            ))

            # Intermediate clean to reduce storage footprint
            shell((
                'rm -rf {out_prefix}.coefs; '
                'rm -rf {out_prefix}.marginals; '
                'rm -rf {out_prefix}.matdata; '
                'rm -rf {out_prefix}.preds; '
                'rm -rf {out_prefix}.traindata'
            ))

        shell('touch {output.out}')

def get_null_pops_gather_input(wildcards):
    if config['loco']:
        in_list = expand(
            '{pref}/{{study}}/null_pops/swarm/batch{bat}/done.{chr}.txt',
            pref = module_prefix,
            bat = list(range(1, n_batch+1)),
            chr = [ str(x) for x in config['chromosomes'] ]
        )
    else:
        in_list = expand(
            '{pref}/{{study}}/null_pops/swarm/batch{bat}/done.{chr}.txt',
            pref = module_prefix,
            bat = list(range(1, n_batch+1)),
            chr = 'all'
        )
    return(in_list)

rule null_pops__gather:
    """
    Gather nominal PoPs results
    """
    input:
        get_null_pops_gather_input
    output:
        '{}/{{study}}/null_pops/null_pops.tsv'.format(module_prefix)
    run:
        in_list = expand(
            '{pref}/{study}/null_pops/swarm/batch{bat}/perm{i}_{j}/null_pops.{chr}.results',
            pref = module_prefix,
            study = wildcards.study,
            bat = list(range(1, n_batch+1)),
            i = list(range(1, n_per_batch+1)),
            j = n_per_batch,
            chr = [ str(x) for x in config['chromosomes'] ] if config['loco'] else 'all'
        )
        print('Merging {} results files...'.format(len(in_list)))

        if len(input) == 1:
            shell('set +o pipefail; cat {in_list} | cut -f1,2 > {output}')
        else:
            shell((
                'set +o pipefail; '
                'cat {in_list[0]} | head -n 1 | cut -f1,2 > {output}; '
                'for file in {in_list}; do '
                '   cat $file | tail -n +2 | cut -f1,2 >> {output}; '
                'done'
            ))

rule null_pops__clean_pops:
    """
    Clean null directoriess
    """
    input:
        '{}/{{study}}/null_pops/null_pops.tsv'.format(module_prefix)
    output:
        '{}/{{study}}/null_pops/clean.done'.format(module_prefix)
    params:
        swarm = lambda wildcards: '{}/{}/null_pops/swarm'.format(
            module_prefix,
            wildcards.study
        )
    shell:
        'rm -rf {params.swarm}; '
        'touch {output}'

rule null_pops__clean_features:
    """
    Clean features
    """
    input:
        expand(
            '{pref}/{study}/null_pops/null_pops.tsv',
            pref = module_prefix,
            study = list(config['study'].keys())
        )
    output:
        '{}/features__null_swarm/clean.done'.format(module_prefix)
    params:
        feats = '{}/features__null_swarm/batch'.format(module_prefix)
    shell:
        'rm -rf {params.feats}*; '
        'touch {output}'
################################################################################
# Step 4: Calculate pvalue
################################################################################
rule pval__calculate_pval:
    input:
        nom = '{}/{{study}}/nominal_pops.tsv'.format(module_prefix),
        null = '{}/{{study}}/null_pops/null_pops.tsv'.format(module_prefix)
    output:
        pvals = '{}/{{study}}/empirical_pops_pvals.tsv'.format(module_prefix),
        null_dist = '{}/{{study}}/empirical_pops_null_dist.tsv'.format(module_prefix)
    params:
        grouping_vars = 'ENSGID',
        score_column = 'PoPS_Score',
        basis = 'grouped' if config['empirical_pops']['grouped'] else 'pooled',
        abs_val = '--use_absolute_value' if config['empirical_pops']['use_absolute_value'] else '',
        script = srcdir('scripts/emp_pops__calculate_empirical_pval.py')
    shell:
       '{params.script} '
       '--observed {input.nom} '
       '--null {input.null} '
       '--index_col {params.grouping_vars} '
       '--score_col {params.score_column} '
       '--comparison_operator ">=" '
       '--basis {params.basis} '
       '--save_merged_null {output.null_dist} '
       '--output_file {output.pvals} '
       '{params.abs_val}'

rule pval__annotate:
    """
    Annotate final PoPs results
    """
    input:
        pvals = '{}/{{study}}/empirical_pops_pvals.tsv'.format(module_prefix),
        e2s = config['gene_locations']
    output:
        tmp = temp('{}/{{study}}/temp.empirical_pops_pvals__annotated.tsv'.format(module_prefix)),
        rez = '{}/{{study}}/empirical_pops_pvals__annotated.tsv'.format(module_prefix)
    params:
        perl_script = srcdir('scripts/merge.plx'),
        pval_corr_script = srcdir('scripts/emp_pops__pval_correction.R')
    shell: 
        '{params.perl_script} -a -s 2 {input.e2s}:1 {input.pvals}:1 | tail -n +2 | '
        'awk -v FS="\\t" -v OFS="\\t" \'{{print}} '
        'BEGIN{{print "ENSGID", "symbol", "chr", "start", "end", "tss", "score", "pvalue"}}\' '
        '> {output.tmp}; '
        '{params.pval_corr_script} '
        '--in_file {output.tmp} '
        '--pvalue_col pvalue '
        '--correction_method bh '
        '--out_file {output.tmp}; '
        'cat {output.tmp} | awk -v FS="\\t" -v OFS="\\t" '
        '\'NR==1{{$(NF+1)="trait"}} NR>1{{$(NF+1)="{wildcards.study}"}}1\' '
        '> {output.rez}'

rule pval__plot_pops:
    input:
        pvals = '{}/{{study}}/empirical_pops_pvals__annotated.tsv'.format(module_prefix),
        gloc = config['gene_locations']
    output:
        png = '{}/{{study}}/empirical_pops_pvals.png'.format(module_prefix)
    params:
        script = srcdir('scripts/emp_pops__plot_pval_manhattan.R')
    shell:
        '{params.script} '
        '--in_file {input.pvals} '
        '--gene_loc {input.gloc} '
        '--output_file {output.png}'

rule pval__plot_distribution:
    input:
        nom = '{}/{{study}}/nominal_pops.tsv'.format(module_prefix),
        null_dist = '{}/{{study}}/empirical_pops_null_dist.tsv'.format(module_prefix)
    output:
        dens = '{}/{{study}}/plots/pops_distributions_density.png'.format(module_prefix),
        qq = '{}/{{study}}/plots/pops_distributions_qq.png'.format(module_prefix),
        ecdf = '{}/{{study}}/plots/pops_distributions_ecdf.png'.format(module_prefix)
    params:
        score_col = 'PoPS_Score',
        out_pref = lambda wildcards: '{}/{}/plots/pops_distributions'.format(
            module_prefix,
            wildcards.study
        ),
        script = srcdir('scripts/emp_pops__plot_distributions.R')
    shell:
        '{params.script} '
        '--observed {input.nom} '
        '--expected {input.null_dist} '
        '--score_col {params.score_col} '
        '--output {params.out_pref}'

################################################################################
# Step 5: Merge data
################################################################################
rule merge__results:
    """
    Merge files back together
    """
    input:
        expand(
            '{pref}/{study}/empirical_pops_pvals__annotated.tsv',
            pref = module_prefix,
            study = list(config['study'].keys())
        )
    output:
        rez = '{}/merged_study_results.tsv.gz'.format(module_prefix)
    params:
        input_lst = lambda wildcards, input: ','.join(input),
        script = srcdir('scripts/merge_dataframes.py')
    shell: 
        '{params.script} '
        '--in_list {params.input_lst} '
        '--output_file {output.rez}'