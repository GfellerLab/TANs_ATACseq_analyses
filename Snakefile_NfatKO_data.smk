import pandas as pd
from itertools import chain
from itertools import combinations

# Read the metadata file matching bam files, narrowPeaks files, samples, and conditions
metadata_df = pd.read_csv(config['data']['metadata'], sep = "\t", header = 0)

conditions = metadata_df['groups'].unique()
samples = metadata_df['sample_name']

# Create a dictionary to map each sample to its condition
sample_to_condition = dict(zip(metadata_df['sample_name'], metadata_df['groups']))

# Create a dictionary to hold BAM files grouped by condition
metadata_df['bam_file'] = str(config["preprocessing_path"]) + metadata_df['sample_name'] + '/aligned_mm10/' + metadata_df['sample_name'] + '_sort_dedup.bam'
metadata_df['shifted_bam_file'] = str(config["preprocessing_path"]) + metadata_df['sample_name'] + '/aligned_mm10/' + metadata_df['sample_name'] + '_sort_dedup_shifted.bam'


metadata_df2 = metadata_df[metadata_df['sample_name'] != 'N12']
print(metadata_df2)
bam_files_by_condition = {condition: metadata_df2[metadata_df2['groups'] == condition]['bam_file'].tolist() for condition in conditions}
shifted_bam_files_by_condition = {condition: metadata_df2[metadata_df2['groups'] == condition]['shifted_bam_file'].tolist() for condition in conditions}
narrowPeaks_files = metadata_df2['summits']

directions = ["up", "down"]
logFC_thr = [0, 1, 2]
comparisons = ["Nfatc1KO_VS_WT"]
genesets = ["GOBP",  "kegg_pathway", "reactome" ]
clusters = ["up","down"]

clustering_versions = ["peaks_clustering"]
condition_pairs = list(combinations(conditions, 2))

rule all:
    input:
        str(config["analysis_folder"]) + 'PCA/PCA_PC1-3.pdf',
        str(config["analysis_folder"]) + 'chromVAR/chromvar_scores.txt',
        str(config["analysis_folder"]) + 'raw_counts_N12out.txt',
        str(config["analysis_folder"]) + 'metadata_N12out.txt',
        str(config["analysis_folder"]) + '/TOBIAS_res/NFATC1_footprint_comparison_all.pdf',
        str(config["analysis_folder"]) + '/TOBIAS_res/NFATC1_footprint_comparison_WTbound.pdf',
        str(config["analysis_folder"]) + '/chromVAR/chromVar_TOBIASchange-vs-deviation_test.pdf',
        expand(str(config["analysis_folder"]) + '/TOBIAS_res/{condition}.merged.bam', condition=conditions),
        str(config["analysis_folder"]) + '/TOBIAS_res/merged_peaks.bed',
        expand(str(config["analysis_folder"]) + '/TOBIAS_res/footprinting_scores/{condition}_footprints.bw', condition=conditions),
        expand(str(config["analysis_folder"]) + '/TOBIAS_res/BINDetect_diff_res/{comp_pair}/bindetect_results.txt', comp_pair=comparisons),
        str(config["analysis_folder"]) + '/DAR_HOMER/differential_analysis_output.txt',
        expand(str(config["analysis_folder"]) + '/DAR_HOMER/{comparison}/{comparison}_{direction}_logFC{logFC}.txt', comparison=comparisons, direction=directions, logFC=logFC_thr),
        expand(str(config["analysis_folder"]) + '/DAR_HOMER/{comparison}/annotated_{comparison}_{direction}_logFC{logFC}.txt', comparison=comparisons, direction=directions, logFC=logFC_thr),
        expand(str(config["analysis_folder"]) + '/peaks_clustering/C{cluster}/cluster_peaks.txt', cluster=clusters),
        expand(str(config["analysis_folder"]) + '/pathways_enrichment_res/{geneset}_enrichment.txt', geneset = genesets),
        expand(str(config["analysis_folder"]) + '/HOMER_res/C{cluster}/homerMotifs.all.motifs', cluster=clusters),
        str(config["analysis_folder"]) + '/HOMER_res/all_denovo_TFs.pdf',
        str(config["analysis_folder"]) + '/RNA_data/TFactivity_heatmap100_v2.pdf',
        str(config["analysis_folder"]) + '/RNA_data/peaks_NFATKO-DEGsenrichment/NfatKODARs_DEGs_overlap_window100000.pdf',
        str(config["analysis_folder"]) + '/RNA_res/noWT5_diff_genes_all.xlsx',
        str(config["analysis_folder"]) + '/RNA_res/RNA_proteomics_gsea_common.xlsx',
        str(config["analysis_folder"]) + '/decoupleR_atlases_res//human_atlas/Nfat_TFs_scores.txt'

# include footprinting rules (TOBIAS) to the pipeline
include: str(config["snakemake_path"]) + "/workflow/scripts/Motifs_analyses/TOBIAS_analysis/TOBIAS_rules.smk"

rule PCA_samplesCorrelation:
    '''
    Gather raw counts and save in a count matrix file.
    '''
    input:
        counts_file = config['data']['counts'],
        metadata_file = config['data']['metadata'],
        fasta = config['data']['fasta_file'],
        pfm_file = config['data']['pfm_file'],
        KP_clusters = config['data']['KP_clusters'],
        counts_KP_peaks = config['data']['counts_KP_peaks']
    output:
        pca_res = str(config["analysis_folder"]) + 'PCA/PCA_PC1-3.pdf',
        chromVAR_res  = str(config["analysis_folder"]) + 'chromVAR/chromvar_scores.txt',
        counts_N12out = str(config["analysis_folder"]) + 'raw_counts_N12out.txt',
        metadata_N12out = str(config["analysis_folder"]) + 'metadata_N12out.txt'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['run_chromVAR']['mem_mb'],
        time = config['run_chromVAR']['time']
    threads: config['run_chromVAR']['cpus-per-task']
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/Nfactc1_KO/atac_data_analysis.r {input.metadata_file} {input.counts_file} {input.fasta} {params.output_path}/ {input.pfm_file} \
        {input.KP_clusters} {input.counts_KP_peaks}
        """

rule generate_chromVAR_related_figures:
    input:
        metadata_file = rules.PCA_samplesCorrelation.output.metadata_N12out,
        motif_annotations = config['data']['motif_annotations'],
        dev = str(config["analysis_folder"]) + '/chromVAR/dev_scores.rds',
        zscore = str(config["analysis_folder"]) + '/chromVAR/z_scores.rds',
        footprinting_res = str(config["analysis_folder"]) + '/TOBIAS_res/BINDetect_diff_res/Nfatc1KO_VS_WT/bindetect_results.txt'
    output:
        str(config["analysis_folder"]) + '/chromVAR/chromVar_TOBIASchange-vs-deviation_test.pdf'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/Nfactc1_KO/chromVAR_figures.r {input.dev} {input.zscore} {input.metadata_file} {input.motif_annotations} {input.footprinting_res} {params.output_path}/chromVAR/
        """

rule DAR_HOMER:
    '''
    Identify differential peaks using HOMER.
    '''
    input:
        counts_file = rules.PCA_samplesCorrelation.output.counts_N12out
    output:
        diff_res = str(config["analysis_folder"]) + '/DAR_HOMER/differential_analysis_output.txt'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_latest.sif"
    params:
        output_path = config["analysis_folder"]
    shell:
        """
        mkdir -p {params.output_path}/DAR_HOMER/;
        getDiffExpression.pl {input.counts_file} Nfatc1KO Nfatc1KO WT WT WT WT Nfatc1KO \
        -vst -AvsA -DESeq2 > {params.output_path}/DAR_HOMER/differential_analysis_output.txt
        """

rule extract_diff_peaks:
    '''
    Identify differential peaks using HOMER.
    '''
    input:
        diff_res = rules.DAR_HOMER.output.diff_res
    output:
        diff_peaks = str(config["analysis_folder"]) + '/DAR_HOMER/{comparison}/{comparison}_{direction}_logFC{logFC}.txt'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v4.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/Motifs_analyses/HOMER_analyses/retrieve_diff_peaks.r {input.diff_res} {wildcards.comparison} {wildcards.direction} {wildcards.logFC} {params.output_path}/DAR_HOMER/{wildcards.comparison}/
        """
rule annotate_diff_peaks:
    '''
    Identify differential peaks using HOMER.
    '''
    input:
        diff_peaks = str(config["analysis_folder"]) + '/DAR_HOMER/{comparison}/{comparison}_{direction}_logFC{logFC}.txt'
    output:
        diff_peaks = str(config["analysis_folder"]) + '/DAR_HOMER/{comparison}/annotated_{comparison}_{direction}_logFC{logFC}.txt'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v4.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    wildcard_constraints:
        logFC = "0|1|2"
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/Motifs_analyses/HOMER_analyses/annotate_diff_peaks.r {input.diff_peaks} {wildcards.comparison} {wildcards.direction} {wildcards.logFC} {params.output_path}/DAR_HOMER/{wildcards.comparison}/
        """
# peak clustering here is simply up DARs and down DARs
rule peak_clustering:
    '''
    Cluster the peaks differentially expressed.
    '''
    input:
        counts_file = rules.PCA_samplesCorrelation.output.counts_N12out,
        metadata_file = rules.PCA_samplesCorrelation.output.metadata_N12out,
        fasta = config['data']['fasta_file'],
        diff_res = rules.DAR_HOMER.output.diff_res
    output:
        str(config["analysis_folder"]) + '/peaks_clustering/cluster_peaks_annotation.txt',
        cluster_peaks = expand(str(config["analysis_folder"]) + '/peaks_clustering/C{cluster}/cluster_peaks.txt', cluster = clusters)
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/Nfactc1_KO/peaks_clustering.r {input.counts_file} {input.metadata_file} {input.diff_res} {input.fasta} {params.output_path}/peaks_clustering/
        """

rule peaks_pathway_enrichment:
    '''
    Perform pathways enrichment in each cluster.
    '''
    input:
        annotated_peaks = str(config["analysis_folder"]) + '/peaks_clustering/cluster_peaks_annotation.txt'
    output:
        str(config["analysis_folder"]) + '/pathways_enrichment_res/{geneset}_enrichment.txt'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/peak_annotation_latest.sif"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['pathway_enrichment']['mem_mb'],
        time = config['pathway_enrichment']['time']
    threads: config['pathway_enrichment']['cpus-per-task']
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"],
        gs = "{geneset}"
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/annotation/annotation_analysis.r {input.annotated_peaks} {params.output_path}/pathways_enrichment_res/ {params.gs} {threads} mm10
        """

# on each set of peaks: run homer:
rule motif_enrichment_HOMER:
    '''
    Identify enriched motifs .
    '''
    input:
        peaks_set = str(config["analysis_folder"]) + '/peaks_clustering/C{cluster}/cluster_peaks.txt',
        motif_file = config['data']['homer_motifs_file']
    output:
        homer_enrichment = str(config["analysis_folder"]) + '/HOMER_res/C{cluster}/homerMotifs.all.motifs'
    container:
        "docker://nfcore/atacseq:latest"
    params:
        homer_bin_path = "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/homer/bin/",
        output_path = config["analysis_folder"]
    threads: config['motif_enrichment_HOMER']['cpus-per-task']
    shell:
        """
        mkdir -p {params.output_path}/HOMER_res/C{wildcards.cluster}/;
        PATH=$PATH:{params.homer_bin_path};
        {params.homer_bin_path}findMotifsGenome.pl {input.peaks_set} mm10 {params.output_path}/HOMER_res/C{wildcards.cluster}/ -mset vertebrates -mknown {input.motif_file} -mcheck {input.motif_file} -mask -p {threads};
        {params.homer_bin_path}compareMotifs.pl {params.output_path}/HOMER_res/C{wildcards.cluster}/homerMotifs.all.motifs -known {input.motif_file} {params.output_path}/HOMER_res/C{wildcards.cluster}/ -cpu {threads};
        """

rule motif_enrichment_figure:
    input:
        expand(str(config["analysis_folder"]) + '/HOMER_res/C{cluster}/homerMotifs.all.motifs', cluster = clusters),
        motif_file = config['data']['motif_annotations']
    output:
        str(config["analysis_folder"]) + '/HOMER_res/all_denovo_TFs.pdf'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/Motifs_analyses/HOMER_analyses/homer_de_novo.r {params.output_path}/HOMER_res/ {input.motif_file}
        """

rule peaks_enrichment_nfatKO:
    input:
        dar_path = str(config["analysis_folder"]) + '/peaks_clustering/cluster_peaks_annotation.txt',
        kp_cluster_peaks = config['data']['KP_clusters'],
        nfatko_deg = str(config["analysis_folder"]) + 'RNA_data/noWT5_diff_genes_all.xlsx'
    output:
        str(config["analysis_folder"]) + '/RNA_data/peaks_NFATKO-DEGsenrichment/NfatKODARs_DEGs_overlap_window100000.pdf'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v4.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/Nfactc1_KO/peak_enrichment_DEGs.r \
        {input.nfatko_deg} \
        {params.output_path}/RNA_data/peaks_NFATKO-DEGsenrichment/ \
        100000 \
        {input.dar_path} \
        {input.kp_cluster_peaks}

        Rscript {params.smk_path}/workflow/scripts/Nfactc1_KO/peak_enrichment_DEGs.r \
        {input.nfatko_deg} \
        {params.output_path}/RNA_data/peaks_NFATKO-DEGsenrichment/ \
        50000 \
        {input.dar_path} \
        {input.kp_cluster_peaks}

        Rscript {params.smk_path}/workflow/scripts/Nfactc1_KO/peak_enrichment_DEGs.r \
        {input.nfatko_deg} \
        {params.output_path}/RNA_data/peaks_NFATKO-DEGsenrichment/ \
        10000 \
        {input.dar_path} \
        {input.kp_cluster_peaks}
        """


rule NFATC1_footprinting:
    '''
    Run differential footprinting.
    '''
    input:
        expand(str(config["analysis_folder"]) + '/TOBIAS_res/corrected_signals/{condition}.merged_corrected.bw', condition=conditions),
        str(config["analysis_folder"]) + 'TOBIAS_res/footprinting_scores/Nfatc1KO_bindetect_results/_NFAC1.H12CORE.0.P.B/beds/_NFAC1.H12CORE.0.P.B_all.bed'
    output:
        str(config["analysis_folder"]) + '/TOBIAS_res/NFATC1_footprint_comparison_all.pdf',
        str(config["analysis_folder"]) + '/TOBIAS_res/NFATC1_footprint_comparison_WTbound.pdf'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
    params:
        output_path = config["analysis_folder"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['diff_footprinting']['mem_mb'],
        time = config['diff_footprinting']['time']
    threads: config['diff_footprinting']['cpus-per-task']
    shell:
        """
        export MPLCONFIGDIR={params.output_path}/TOBIAS_res/;
        TOBIAS PlotAggregate --TFBS {params.output_path}TOBIAS_res/footprinting_scores/Nfatc1KO_bindetect_results/_NFAC1.H12CORE.0.P.B/beds/_NFAC1.H12CORE.0.P.B_all.bed \
        --signals {params.output_path}/TOBIAS_res/corrected_signals/WT.merged_corrected.bw {params.output_path}/TOBIAS_res/corrected_signals/Nfatc1KO.merged_corrected.bw \
        --output {params.output_path}/TOBIAS_res/NFATC1_footprint_comparison_all.pdf \
        --share_y both --plot_boundaries --signal-on-x;

        TOBIAS PlotAggregate --TFBS {params.output_path}TOBIAS_res/footprinting_scores/WT_bindetect_results/_NFAC1.H12CORE.0.P.B/beds/_NFAC1.H12CORE.0.P.B_WT_footprints_bound.bed \
        --signals {params.output_path}/TOBIAS_res/corrected_signals/WT.merged_corrected.bw {params.output_path}/TOBIAS_res/corrected_signals/Nfatc1KO.merged_corrected.bw \
        --output {params.output_path}/TOBIAS_res/NFATC1_footprint_comparison_WTbound.pdf \
        --share_y both --plot_boundaries --signal-on-x
        """


#################################### RNA rules #########################################
rule RNA_analysis:
    input:
        rna_data_path = config['data']['rna_data_path']
    output:
        degs = str(config["analysis_folder"]) + '/RNA_res/noWT5_diff_genes_all.xlsx'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v4.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/Nfactc1_KO/rna_analysis.r \
        {input.rna_data_path} \
        {params.output_path}/RNA_res/

        """

rule Proteomics_RNA_gsea:
    input:
        rna_degs = rules.RNA_analysis.output.degs,
        prot_degs = config["data"]["proteomics_diff_res"],
        geneset_path = config['input_data_path']
    output:
        str(config["analysis_folder"]) + '/results_proteomics/gsea_res.xlsx',
        str(config["analysis_folder"]) + '/RNA_res/RNA_proteomics_gsea_common.xlsx'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v4.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"],
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/proteomics/run_gsea.r \
        {input.prot_degs} \
        {input.geneset_path}/mouse_genesets_symbols.xlsx \
        {params.output_path}/results_proteomics/

        Rscript {params.smk_path}/workflow/scripts/Nfactc1_KO/rna_gsea.r \
        {input.rna_degs} \
        {input.geneset_path}/mouse_genesets_symbols.xlsx \
        {params.output_path}/RNA_res/ \
        {params.output_path}/results_proteomics/gsea_res.xlsx
        """

rule neutrophil_atlas_decoupler:
    input:
        atlases_path = config["data"]["neutrophil_atlases"]
    output:
        str(config["analysis_folder"]) + '/decoupleR_atlases_res//human_atlas/Nfat_TFs_scores.txt'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v4.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"],
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/neutrophil_atlases/decoupleR.r \
        {input.atlases_path} \
        {params.output_path}/decoupleR_atlases_res/
        """



# workflow_path=/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/snakemake_pipeline/workflow/
# snakemake --snakefile ${workflow_path}Snakefile_NfatKO_data.smk --profile ${workflow_path}config_files/slurm/ --configfile ${workflow_path}config_files/config_NfatKO.yml --dry-run
