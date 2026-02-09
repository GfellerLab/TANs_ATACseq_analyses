import pandas as pd
from itertools import chain
from itertools import combinations

# Read the metadata file matching bam files, narrowPeaks files, samples, and conditions
metadata_df = pd.read_csv(config['data']['metadata'], sep = "\t", header = 0)
conditions = metadata_df['groups'].unique()
samples = metadata_df['sample_name']
narrowPeaks_files = metadata_df['summits']

# Create a dictionary to map each sample to its condition
sample_to_condition = dict(zip(metadata_df['sample_name'], metadata_df['groups']))

# Create a dictionary to hold BAM files grouped by condition
metadata_df['bam_file'] = str(config["preprocessing_path"]) + metadata_df['sample_name'] + '/aligned_mm10/' + metadata_df['sample_name'] + '_sort_dedup.bam'
metadata_df['shifted_bam_file'] = str(config["preprocessing_path"]) + metadata_df['sample_name'] + '/aligned_mm10/' + metadata_df['sample_name'] + '_sort_dedup_shifted.bam'

bam_files_by_condition = {condition: metadata_df[metadata_df['groups'] == condition]['bam_file'].tolist() for condition in conditions}
shifted_bam_files_by_condition = {condition: metadata_df[metadata_df['groups'] == condition]['shifted_bam_file'].tolist() for condition in conditions}

directions = ["up", "down"]
logFC_thr = [0, 1, 2]
comparisons = ["KP_lung_SiglecF_high_VS_KP_lung_SiglecF_low", "KP_lung_SiglecF_high_VS_Healthy_lung", "KP_lung_SiglecF_high_VS_KP_blood", "KP_lung_SiglecF_high_VS_Healthy_blood", "KP_lung_SiglecF_low_VS_Healthy_lung", "KP_lung_SiglecF_low_VS_KP_blood", "KP_lung_SiglecF_low_VS_Healthy_blood", "Healthy_lung_VS_KP_blood", "Healthy_lung_VS_Healthy_blood", "Healthy_blood_VS_KP_blood"]
genesets = ["GOBP",  "kegg_pathway", "reactome" ]

clustering_versions = ["peaks_clustering_V4"]
condition_pairs = list(combinations(conditions, 2))

rule all:
    input:
        str(config["analysis_folder"]) + '/PCA/pca_coordinates.rds',
        str(config["analysis_folder"]) + '/chromVAR/dev_scores.rds',
        str(config["analysis_folder"]) + '/DAR_HOMER/differential_analysis_output.txt',
        expand(str(config["analysis_folder"]) + '/DAR_HOMER/{comparison}/{comparison}_{direction}_logFC{logFC}.txt', comparison=comparisons, direction=directions, logFC=logFC_thr),
        expand(str(config["analysis_folder"]) + '/DAR_HOMER/{comparison}/annotated_{comparison}_{direction}_logFC{logFC}.txt', comparison=comparisons, direction=directions, logFC=logFC_thr),
        expand(str(config["analysis_folder"]) + '/TOBIAS_res/{condition}.merged.bam', condition=conditions),
        str(config["analysis_folder"]) + '/TOBIAS_res/merged_peaks.bed',
        expand(str(config["analysis_folder"]) + '/TOBIAS_res/footprinting_scores/{condition}_footprints.bw', condition=conditions),
        expand(str(config["analysis_folder"]) + '/TOBIAS_res/BINDetect_diff_res/{comp_pair}/bindetect_results.txt', comp_pair=comparisons),
        str(config["analysis_folder"]) + '/chromVAR/chromVar_TOBIASchange-vs-deviation_SiglecFhiVSSiglecFlo.pdf',
        # str(config['input_data_path']) + '/mouse_genesets_symbols.txt',
        # str(config['input_data_path']) + '/human_genesets_symbols.txt',
        expand(str(config["analysis_folder"]) + '/peaks_clustering_V4/C{cluster}/cluster_peaks.txt', cluster = [1,2,3,4,5]),
        expand(str(config["analysis_folder"]) + '/pathways_enrichment_res_V4/{geneset}_enrichment.txt', geneset = genesets),
        expand(str(config["analysis_folder"]) + '/atlases_res_V4/C{cluster}_associated_genes/human_atlas/TSSdist1e+05_signatures_bySubtype.pdf', cluster = [1,2,3,4,5]),
        expand(str(config["analysis_folder"]) + '/HOMER_res_V4/C{cluster}/homerMotifs.all.motifs', cluster=[1,2,3,4,5]),
        str(config["analysis_folder"]) + '/HOMER_res_V4/all_denovo_TFs.pdf',
        str(config["analysis_folder"]) + '/peaks_clustering_V4/peak_clustering_annotations_withNonDARs.pdf',
        str(config["analysis_folder"]) + '/peaks_clustering_V4/upset_plot.pdf',
        expand(str(config["analysis_folder"]) + '/{peak_clustering_version}/human_atlas/TSSdist1e+05_heatmap.pdf', peak_clustering_version=clustering_versions),
        expand(str(config["analysis_folder"]) + '/{peak_clustering_version}/peaks_KP-DEGsenrichment/clusterPeaks_DEGs_overlap_window100000.pdf', peak_clustering_version=clustering_versions),
        expand(str(config["analysis_folder"]) + '/{peak_clustering_version}/peaks_NFATKO-DEGsenrichment/clusterPeaks_DEGs_overlap_window100000.pdf', peak_clustering_version=clustering_versions),
        str(config["analysis_folder"]) + '/AddChromatinModule_results/tss_samplesScores.txt'




# include footprinting rules (TOBIAS) to the pipeline
include: str(config["snakemake_path"]) + "/workflow/scripts/Motifs_analyses/TOBIAS_analysis/TOBIAS_rules.smk"

rule PCA_samplesCorrelation:
    '''
    Gather raw counts and save in a count matrix file.
    '''
    input:
        counts_file = config['data']['counts'],
        metadata_file = config['data']['metadata'],
        fasta = config['data']['fasta_file']
    output:
        pca_res = str(config["analysis_folder"]) + '/PCA/pca_coordinates.rds'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/pca_analysis_KP.r {input.metadata_file} {input.counts_file} {input.fasta} {params.output_path}/PCA/
        """

rule run_chromVAR:
    '''
    Compute deviation scores and save it in a matrix.
    '''
    input:
        counts_file = config['data']['counts'],
        metadata_file = config['data']['metadata'],
        pfm_file = config['data']['pfm_file']
    output:
        chromvar_dev_res = str(config["analysis_folder"]) + '/chromVAR/dev_scores.rds',
        chromvar_z_res = str(config["analysis_folder"]) + '/chromVAR/z_scores.rds'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['run_chromVAR']['mem_mb'],
        time = config['run_chromVAR']['time']
    threads: config['run_chromVAR']['cpus-per-task']
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/Motifs_analyses/chromVAR_analysis.r {input.metadata_file} {input.counts_file} {input.pfm_file} {params.output_path}/chromVAR/ mm10
        """

rule DAR_HOMER:
    '''
    Identify differential peaks using HOMER.
    '''
    input:
        counts_file = config['data']['counts']
    output:
        diff_res = str(config["analysis_folder"]) + '/DAR_HOMER/differential_analysis_output.txt'
    # container:
    #     "docker://nfcore/atacseq:latest"
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_latest.sif"
    params:
        output_path = config["analysis_folder"]
    shell:
        """
        mkdir -p {params.output_path}/DAR_HOMER/;
        getDiffExpression.pl {input.counts_file} KP_lung_SiglecF_high KP_lung_SiglecF_high Healthy_blood KP_blood Healthy_lung KP_lung_SiglecF_low Healthy_blood KP_lung_SiglecF_high Healthy_blood KP_blood Healthy_lung KP_lung_SiglecF_low KP_blood Healthy_lung KP_lung_SiglecF_low KP_lung_SiglecF_high Healthy_blood KP_blood Healthy_lung KP_lung_SiglecF_low \
        -batch 2 3 4 4 4 4 1 4 5 5 5 5 1 1 1 1 2 2 2 2 -vst -AvsA -DESeq2 > {params.output_path}/DAR_HOMER/differential_analysis_output.txt
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

rule generate_chromVAR_related_figures:
    input:
        metadata_file = config['data']['metadata'],
        motif_annotations = config['data']['motif_annotations'],
        dev = str(config["analysis_folder"]) + '/chromVAR/dev_scores.rds',
        zscore = str(config["analysis_folder"]) + '/chromVAR/z_scores.rds',
        footprinting_res = str(config["analysis_folder"]) + '/TOBIAS_res/BINDetect_diff_res/KP_lung_SiglecF_high_VS_KP_lung_SiglecF_low/bindetect_results.txt',
        footprinting_res2 = str(config["analysis_folder"]) + '/TOBIAS_res/BINDetect_diff_res/KP_lung_SiglecF_high_VS_Healthy_lung/bindetect_results.txt'
    output:
        str(config["analysis_folder"]) + '/chromVAR/chromVar_TOBIASchange-vs-deviation_SiglecFhiVSSiglecFlo.pdf',
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/Motifs_analyses/chromVAR_figures.r {input.dev} {input.zscore} {input.metadata_file} {input.motif_annotations} {input.footprinting_res} {input.footprinting_res2} {params.output_path}/chromVAR/
        """

rule get_genesets:
    output:
        str(config['input_data_path']) + '/mouse_genesets_symbols.txt',
        str(config['input_data_path']) + '/human_genesets_symbols.txt'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v4.sif"
    params:
        output_path = config["input_data_path"],
        senescence_genes = str(config['input_data_path']) + '/Senescence_genesets_Saul_2022_Nat_Comm.xlsx',
        mouse_mapping = str(config['input_data_path']) + '/HOM_MouseHumanSequence.rpt',
        smk_path = config["snakemake_path"]
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/input_data_generation/genesets_save.r {params.senescence_genes} {params.mouse_mapping} {params.output_path}
        """



#########################################################
# Peaks clustering V4
#########################################################
# get group specific peaks
rule peak_clustering_V4:
    '''
    Cluster the peaks differentially expressed.
    '''
    input:
        counts_file = config['data']['counts'],
        metadata_file = config['data']['metadata'],
        fasta = config['data']['fasta_file'],
        diff_res = rules.DAR_HOMER.output.diff_res
    output:
        str(config["analysis_folder"]) + '/peaks_clustering_V4/cluster_peaks_annotation.txt',
        cluster_peaks = expand(str(config["analysis_folder"]) + '/peaks_clustering_V4/C{cluster}/cluster_peaks.txt', cluster = [1,2,3,4,5])
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/Motifs_analyses/HOMER_analyses/peak_clustering.r {input.counts_file} {input.metadata_file} {input.diff_res} {input.fasta} {params.output_path}/peaks_clustering_V4/ V4
        """

rule peaks_pathway_enrichment_V4:
    '''
    Perform pathways enrichment in each cluster.
    '''
    input:
        annotated_peaks = str(config["analysis_folder"]) + '/peaks_clustering_V4/cluster_peaks_annotation.txt'
    output:
        str(config["analysis_folder"]) + '/pathways_enrichment_res_V4/{geneset}_enrichment.txt'
    # container:
    #     "docker://agabriel/peak_annotation:latest"
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
        Rscript {params.smk_path}/workflow/scripts/annotation/annotation_analysis.r {input.annotated_peaks} {params.output_path}/pathways_enrichment_res_V4/ {params.gs} {threads} mm10
        """

rule atlases_signatures_V4:
    input:
        expand(str(config["analysis_folder"]) + '/peaks_clustering_V4/C{cluster}/cluster_peaks.txt', cluster = [1,2,3,4,5])
    output:
        str(config["analysis_folder"]) + '/atlases_res_V4/C{cluster}_associated_genes/human_atlas/TSSdist1e+05_signatures_bySubtype.pdf',
        str(config["analysis_folder"]) + '/atlases_res_V4/C{cluster}_associated_genes/mouse_atlas/TSSdist1e+05_signatures_bySubtype.pdf'
    # container:
    #     "docker://agabriel/peak_annotation:latest"
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v4.sif"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['atlases_signatures']['mem_mb'],
        time = config['atlases_signatures']['time']
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"],
        cluster = lambda wildcards: wildcards.cluster
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/neutrophil_atlases/atlas_signatures.r \
        /work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/data/raw/Neutrophil_atlas/ \
        {params.output_path}/peaks_clustering_V4/C{params.cluster}/cluster_peaks.txt \
        {params.output_path}/atlases_res_V4/C{params.cluster}_associated_genes/
        """

# on each set of peaks: run homer:
rule motif_enrichment_HOMER_V4:
    '''
    Identify enriched motifs .
    '''
    input:
        peaks_set = str(config["analysis_folder"]) + '/peaks_clustering_V4/C{cluster}/cluster_peaks.txt',
        motif_file = config['data']['homer_motifs_file']
    output:
        homer_enrichment = str(config["analysis_folder"]) + '/HOMER_res_V4/C{cluster}/homerMotifs.all.motifs'
    container:
        "docker://nfcore/atacseq:latest"
    params:
        homer_bin_path = "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/homer/bin/",
        output_path = config["analysis_folder"]
    threads: config['motif_enrichment_HOMER']['cpus-per-task']
    shell:
        """
        mkdir -p {params.output_path}/HOMER_res_V4/C{wildcards.cluster}/;
        PATH=$PATH:{params.homer_bin_path};
        {params.homer_bin_path}findMotifsGenome.pl {input.peaks_set} mm10 {params.output_path}/HOMER_res_V4/C{wildcards.cluster}/ -mset vertebrates -mknown {input.motif_file} -mcheck {input.motif_file} -mask -p {threads};
        {params.homer_bin_path}compareMotifs.pl {params.output_path}/HOMER_res_V4/C{wildcards.cluster}/homerMotifs.all.motifs -known {input.motif_file} {params.output_path}/HOMER_res_V4/C{wildcards.cluster}/ -cpu {threads};
        """
rule motif_enrichment_figure_V4:
    input:
        expand(str(config["analysis_folder"]) + '/HOMER_res_V4/C{cluster}/homerMotifs.all.motifs', cluster = [1,2,3,4,5]),
        motif_file = config['data']['motif_annotations']
    output:
        str(config["analysis_folder"]) + '/HOMER_res_V4/all_denovo_TFs.pdf',
        str(config["analysis_folder"]) + '/HOMER_res_V4/homer_enrichment_results.txt'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/Motifs_analyses/HOMER_analyses/homer_de_novo.r {params.output_path}/HOMER_res_V4/ {input.motif_file}
        """

rule peak_related_figures_V4:
    '''
    Upset plot with peaks overlap info + peaks annotation barplots
    '''
    input:
        diff_res = rules.DAR_HOMER.output.diff_res,
        peak_clustering = str(config["analysis_folder"]) + '/peaks_clustering_V4/cluster_peaks_annotation.txt',
        all_called_peaks = config['data']['indiv_samples_peaks'],
        counts_file = config['data']['counts'],
        metadata_file = config['data']['metadata'],
        fasta = config['data']['fasta_file'],
        motif_annotations = config['data']['motif_annotations'],
        homer_res = str(config["analysis_folder"]) + '/HOMER_res_V4/homer_enrichment_results.txt'
    output:
        str(config["analysis_folder"]) + '/peaks_clustering_V4/peak_clustering_annotations_withNonDARs.pdf',
        str(config["analysis_folder"]) + '/peaks_clustering_V4/upset_plot.pdf'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v4.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/preprocessing/get_peakRelated_figures_KPlung.r {input.all_called_peaks} {params.output_path}/peaks_clustering_V4/ {input.peak_clustering} {input.diff_res} \
        {input.metadata_file} {input.counts_file} {input.fasta} {input.homer_res} {input.motif_annotations}
        """

rule pancreas_validation:
    input:
        annotated_peaks = str(config["analysis_folder"]) + '/{peak_clustering_version}/cluster_peaks_annotation.txt',
        counts_file = config['Pancreas_data']['counts'],
        metadata_file = config['Pancreas_data']['metadata'],
        fasta = config['data']['fasta_file']
    output:
        str(config["analysis_folder"]) + '/{peak_clustering_version}/pancreas_heatmap.pdf',
        str(config["analysis_folder"]) + '/{peak_clustering_version}/pancreas_boxplots.pdf'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v4.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/pancreas_analyses/pancreas_validation.r {input.metadata_file} {input.counts_file} {input.fasta} {params.output_path}/{wildcards.peak_clustering_version}/ {input.annotated_peaks}
        """

rule neutrophil_atlas_heatmap:
    input:
        annotated_peaks = str(config["analysis_folder"]) + '/{peak_clustering_version}/cluster_peaks_annotation.txt'
    output:
        str(config["analysis_folder"]) + '/{peak_clustering_version}/human_atlas/TSSdist1e+05_heatmap.pdf'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['atlases_signatures']['mem_mb'],
        time = config['atlases_signatures']['time']
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/neutrophil_atlases/signature_heatmaps.r \
        /work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/data/raw/Neutrophil_atlas/ \
        {input.annotated_peaks} \
        {params.output_path}/{wildcards.peak_clustering_version}/
        """

rule peaks_enrichment:
    input:
        expand(str(config["analysis_folder"]) + '/DAR_HOMER/{comparison}/{comparison}_{direction}_logFC{logFC}.txt', comparison=comparisons, direction=directions, logFC=logFC_thr),
        annotated_cluster_peaks = str(config["analysis_folder"]) + '/{peak_clustering_version}/cluster_peaks_annotation.txt',
        kp_deg = str(config["snakemake_path"]) + '/input_data/KPlung_RNA_diff_genes.xlsx'
    output:
        str(config["analysis_folder"]) + '/{peak_clustering_version}/peaks_KP-DEGsenrichment/clusterPeaks_DEGs_overlap_window100000.pdf'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v4.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"],
        dar_path = str(config["analysis_folder"]) + '/DAR_HOMER/KP_lung_SiglecF_high_VS_KP_lung_SiglecF_low/'
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/Motifs_analyses/HOMER_analyses/peak_enrichment_tests.r \
        {input.kp_deg} \
        {params.output_path}/{wildcards.peak_clustering_version}/peaks_KP-DEGsenrichment/ \
        100000 \
        {params.dar_path} \
        {input.annotated_cluster_peaks}
        """

rule peaks_enrichment_nfatKO:
    input:
        expand(str(config["analysis_folder"]) + '/DAR_HOMER/{comparison}/{comparison}_{direction}_logFC{logFC}.txt', comparison=comparisons, direction=directions, logFC=logFC_thr),
        annotated_cluster_peaks = str(config["analysis_folder"]) + '/{peak_clustering_version}/cluster_peaks_annotation.txt',
        kp_deg = str(config["snakemake_path"]) + '/input_data/NFATKO_RNA_diff_genes.xlsx'
    output:
        str(config["analysis_folder"]) + '/{peak_clustering_version}/peaks_NFATKO-DEGsenrichment/clusterPeaks_DEGs_overlap_window100000.pdf'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v4.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"],
        dar_path = str(config["analysis_folder"]) + '/DAR_HOMER/KP_lung_SiglecF_high_VS_KP_lung_SiglecF_low/'
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/Motifs_analyses/HOMER_analyses/peak_enrichment_tests.r \
        {input.kp_deg} \
        {params.output_path}/{wildcards.peak_clustering_version}/peaks_NFATKO-DEGsenrichment/ \
        100000 \
        {params.dar_path} \
        {input.annotated_cluster_peaks}
        """

rule addChromatinModule:
    input:
        counts_file = config['data']['counts'],
        metadata_file = config['data']['metadata']
    output:
        str(config["analysis_folder"]) + '/AddChromatinModule_results/tss_samplesScores.txt'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v4.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"],
        senescence_genes = str(config['input_data_path']) + '/Senescence_genesets_Saul_2022_Nat_Comm.xlsx',
        db_genesets = str(config['input_data_path']) + '/peaks_db_genesets_mouse.rds'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['addChromatinModule']['mem_mb'],
        time = config['addChromatinModule']['time']
    threads: config['addChromatinModule']['cpus-per-task']
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/module_scores/addchromatinscore_mouse.r \
        {input.metadata_file} \
        {input.counts_file} \
        {params.output_path}/AddChromatinModule_results/ \
        {params.senescence_genes} \
        {params.db_genesets}
        """


# workflow_path=/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/snakemake_pipeline/workflow/
# snakemake --snakefile ${workflow_path}Snakefile_KP_lung.smk --profile ${workflow_path}config_files/slurm/ --configfile ${workflow_path}config_files/config_KPlung.yml
