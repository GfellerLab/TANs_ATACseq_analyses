import pandas as pd
from itertools import chain
from itertools import combinations

directions = ["up", "down"]
logFC_thr = [0, 1, 2]
conditions_pairwise = "Tumor", "Adjacent_lung", "Blood"
comparisons_pairwise = ["Tumor_VS_Adjacent_lung", "Tumor_VS_Blood", "Adjacent_lung_VS_Blood"]

conditions = "Tumor", "Normal"
comparisons = ["Tumor_VS_Normal"]
genesets = ["GOBP",  "kegg_pathway", "reactome" ]

clusters = ["up","down"]
clusters_pairwise = [1,2,3,4,5]
clusters_pairwise_grouped = ["up", "down"]


rule all:
    input:
        str(config["analysis_folder"]) + 'merged_data/raw_counts_merged.txt',
        str(config["analysis_folder"]) + 'chromVAR/chromvar_scores.csv',
        str(config["analysis_folder"]) + 'DAR_HOMER/differential_analysis_output.txt',
        expand(str(config["analysis_folder"]) + 'DAR_HOMER/{comparison}/{comparison}_{direction}_logFC{logFC}.txt', comparison=comparisons, direction=directions, logFC=logFC_thr),
        expand(str(config["analysis_folder"]) + 'DAR_HOMER/{comparison}/annotated_{comparison}_{direction}_logFC{logFC}.txt', comparison=comparisons, direction=directions, logFC=logFC_thr),
        expand(str(config["analysis_folder"]) + 'peaks_clustering/C{cluster}/cluster_peaks.txt', cluster = clusters),
        expand(str(config["analysis_folder"]) + 'atlases_res/C{cluster}_associated_genes/mouse_atlas/TSSdist1e+05_signatures_bySubtype.pdf', cluster=clusters),
        expand(str(config["analysis_folder"]) + 'pathways_enrichment_res/{geneset}_enrichment.txt', geneset = genesets),
        expand(str(config["analysis_folder"]) + 'HOMER_res/C{cluster}/homerMotifs.all.motifs', cluster=clusters),
        expand(str(config["analysis_folder"]) + 'DAR_HOMER_pairwise/{comparison_pairwise}/{comparison_pairwise}_{direction}_logFC{logFC}.txt', comparison_pairwise=comparisons_pairwise, direction=directions, logFC=logFC_thr),
        expand(str(config["analysis_folder"]) + 'DAR_HOMER_pairwise/{comparison_pairwise}/annotated_{comparison_pairwise}_{direction}_logFC{logFC}.txt', comparison_pairwise=comparisons_pairwise, direction=directions, logFC=logFC_thr),
        expand(str(config["analysis_folder"]) + 'peaks_clustering_pairwise/C{cluster_pairwise}/cluster_peaks.txt', cluster_pairwise = clusters_pairwise),
        expand(str(config["analysis_folder"]) + 'HOMER_res_pairwise/C{cluster_pairwise}/homerMotifs.all.motifs', cluster_pairwise = clusters_pairwise),
        str(config["analysis_folder"]) + 'HOMER_res_pairwise/all_denovo_TFs.pdf',
        str(config["analysis_folder"]) + 'HOMER_res/all_denovo_TFs.pdf',
        expand(str(config["analysis_folder"]) + 'pathways_enrichment_res_pairwise/{geneset}_enrichment.txt', geneset = genesets),
        str(config["analysis_folder"]) + 'peaks_clustering/peak_clustering_annotations_withNonDARs.pdf',
        str(config["analysis_folder"]) + 'peaks_clustering_pairwise/homer_res_heatmap.pdf',
        str(config["analysis_folder"]) + 'AddChromatinModule_results/tss_samplesScores.txt'




rule PCA_chromVAR:
    input:
        counts_file = config['human_data']['counts'],
        metadata_file = config['human_data']['metadata'],
        fasta = config['human_data']['fasta_file'],
        pfm_file = config['human_data']['pfm_file'],
        motif_annotations = config['human_data']['motif_annotations']
    output:
        pca_res = str(config["analysis_folder"]) + 'PCA/PCA_plots_merged_samples.pdf',
        merged_counts_res = str(config["analysis_folder"]) + 'merged_data/raw_counts_merged.txt',
        metadata_merged = str(config["analysis_folder"]) + 'merged_data/metadata_merged.txt',
        merged_corrected_counts_res = str(config["analysis_folder"]) + 'merged_data/corrected_counts_merged.txt',
        chromvar_res = str(config["analysis_folder"]) + 'chromVAR/chromvar_scores.csv',
        top_TAN_peaks_up = str(config["analysis_folder"]) + 'PCA/top_pc1_up.txt',
        top_TAN_peaks_down = str(config["analysis_folder"]) + 'PCA/top_pc1_down.txt'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    resources:
        mem_mb = config["PCA_chromVAR"]["mem_mb"]
    threads: 3
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/human_data/pca_analysis.r {input.counts_file} {input.metadata_file} {input.fasta} {params.output_path}/ {input.pfm_file} {input.motif_annotations} TRUE
        """

rule DAR_HOMER:
    '''
    Identify differential peaks using HOMER.
    '''
    input:
        counts_file = str(config["analysis_folder"]) + 'merged_data/raw_counts_merged.txt'
    output:
        diff_res = str(config["analysis_folder"]) + 'DAR_HOMER/differential_analysis_output.txt'
    params:
        output_path = config["analysis_folder"]
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_latest.sif"
    shell:
        """
        mkdir -p {params.output_path}/DAR_HOMER/;
        getDiffExpression.pl {input.counts_file} Normal Normal Tumor Normal Normal Tumor Normal Normal Tumor Normal Normal Tumor \
        -batch 1 1 1 1 1 1 2 2 2 2 1 2 -vst -AvsA -DESeq2 > {params.output_path}/DAR_HOMER/differential_analysis_output.txt
        """

rule extract_diff_peaks:
    '''
    Identify differential peaks using HOMER.
    '''
    input:
        diff_res = rules.DAR_HOMER.output.diff_res
    output:
        diff_peaks = str(config["analysis_folder"]) + 'DAR_HOMER/{comparison}/{comparison}_{direction}_logFC{logFC}.txt'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
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
        diff_peaks = str(config["analysis_folder"]) + 'DAR_HOMER/{comparison}/{comparison}_{direction}_logFC{logFC}.txt'
    output:
        diff_peaks = str(config["analysis_folder"]) + 'DAR_HOMER/{comparison}/annotated_{comparison}_{direction}_logFC{logFC}.txt'
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


rule peak_clustering:
    '''
    Cluster the peaks differentially expressed.
    '''
    input:
        counts_file = str(config["analysis_folder"]) + 'merged_data/raw_counts_merged.txt',
        metadata_file = str(config["analysis_folder"]) + 'merged_data/metadata_merged.txt',
        fasta = config['human_data']['fasta_file'],
        diff_res = rules.DAR_HOMER.output.diff_res
    output:
        str(config["analysis_folder"]) + 'peaks_clustering/cluster_peaks_annotation.txt',
        cluster_peaks = expand(str(config["analysis_folder"]) + 'peaks_clustering/C{cluster}/cluster_peaks.txt', cluster = clusters)
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/human_data/peak_clustering.r {input.counts_file} {input.metadata_file} {input.diff_res} {input.fasta} {params.output_path}/peaks_clustering/ TRUE
        """

rule peaks_pathway_enrichment:
    '''
    Perform pathways enrichment in each cluster.
    '''
    input:
        annotated_peaks = str(config["analysis_folder"]) + 'peaks_clustering/cluster_peaks_annotation.txt'
    output:
        str(config["analysis_folder"]) + 'pathways_enrichment_res/{geneset}_enrichment.txt'
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
        Rscript {params.smk_path}/workflow/scripts/annotation/annotation_analysis.r {input.annotated_peaks} {params.output_path}/pathways_enrichment_res/ {params.gs} {threads} hg38
        """

# on each set of peaks: run homer:
rule motif_enrichment_HOMER:
    '''
    Identify enriched motifs .
    '''
    input:
        peaks_set = str(config["analysis_folder"]) + 'peaks_clustering/C{cluster}/cluster_peaks.txt',
        motif_file = config['human_data']['homer_motifs_file']
    output:
        homer_enrichment = str(config["analysis_folder"]) + 'HOMER_res/C{cluster}/homerMotifs.all.motifs'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_latest.sif"
    params:
        homer_bin_path = "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/homer/bin/",
        output_path = config["analysis_folder"]
    threads: config['motif_enrichment_HOMER']['cpus-per-task']
    shell:
        """
        mkdir -p {params.output_path}/HOMER_res/C{wildcards.cluster}/;
        PATH=$PATH:{params.homer_bin_path};
        {params.homer_bin_path}findMotifsGenome.pl {input.peaks_set} hg38 {params.output_path}/HOMER_res/C{wildcards.cluster}/ -mset vertebrates -mknown {input.motif_file} -mcheck {input.motif_file} -mask -p {threads};
        {params.homer_bin_path}compareMotifs.pl {params.output_path}/HOMER_res/C{wildcards.cluster}/homerMotifs.all.motifs -known {input.motif_file} {params.output_path}/HOMER_res/C{wildcards.cluster}/ -cpu {threads};
        """

rule motif_enrichment_figure:
    input:
        expand(str(config["analysis_folder"]) + 'HOMER_res/C{cluster}/homerMotifs.all.motifs', cluster = clusters),
        motif_file = config['human_data']['motif_annotations']
    output:
        str(config["analysis_folder"]) + 'HOMER_res/all_denovo_TFs.pdf'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/Motifs_analyses/HOMER_analyses/homer_de_novo.r {params.output_path}/HOMER_res/ {input.motif_file}
        """

rule atlases_signatures:
    input:
        expand(str(config["analysis_folder"]) + 'peaks_clustering/C{cluster}/cluster_peaks.txt', cluster = clusters),
        diff_res = rules.DAR_HOMER.output.diff_res
    output:
        str(config["analysis_folder"]) + 'atlases_res/C{cluster}_associated_genes/human_atlas/TSSdist1e+05_signatures_bySubtype.pdf',
        str(config["analysis_folder"]) + 'atlases_res/C{cluster}_associated_genes/mouse_atlas/TSSdist1e+05_signatures_bySubtype.pdf'
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
        Rscript {params.smk_path}/workflow/scripts/human_data/atlas_signatures.r \
        /work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/data/raw/Neutrophil_atlas/ \
        {params.output_path}/peaks_clustering/C{params.cluster}/cluster_peaks.txt \
        {params.output_path}/atlases_res/C{params.cluster}_associated_genes/ {input.diff_res}
        """

rule gsva_run:
    input:
        counts_file = str(config["analysis_folder"]) + 'merged_data/corrected_counts_merged.txt',
        metadata_file = str(config["analysis_folder"]) + 'merged_data/metadata_merged.txt',
        fasta = config['human_data']['fasta_file'],
        genesets_info = config['input_data_path'],
        gtf = config['human_data']['gtf']
    output: str(config['analysis_folder']) + 'GSVA/gsva_scores.txt'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v4.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/annotation/gsva_analysis.r {input.counts_file} {input.metadata_file} {input.fasta} \
        {params.output_path}/GSVA/ {input.genesets_info}/human_genesets_symbols.xlsx {input.gtf}
        """


# Run pairwise differential analysis_folder
rule DAR_HOMER_pairwise:
    '''
    Identify differential peaks using HOMER.
    '''
    input:
        counts_file = str(config["analysis_folder"]) + 'merged_data/raw_counts_merged.txt'
    output:
        diff_res = str(config["analysis_folder"]) + 'DAR_HOMER_pairwise/differential_analysis_output.txt'
    params:
        output_path = config["analysis_folder"]
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_latest.sif"
    shell:
        """
        mkdir -p {params.output_path}/DAR_HOMER_pairwise/;
        getDiffExpression.pl {input.counts_file} Blood Adjacent_lung Tumor Blood Adjacent_lung Tumor Blood Adjacent_lung Tumor Blood Adjacent_lung Tumor \
        -batch 1 1 1 1 1 1 2 2 2 2 1 2 -vst -AvsA -DESeq2 > {params.output_path}/DAR_HOMER_pairwise/differential_analysis_output.txt
        """

rule extract_diff_peaks_pairwise:
    '''
    Identify differential peaks using HOMER.
    '''
    input:
        diff_res = rules.DAR_HOMER_pairwise.output.diff_res
    output:
        diff_peaks = str(config["analysis_folder"]) + 'DAR_HOMER_pairwise/{comparison_pairwise}/{comparison_pairwise}_{direction}_logFC{logFC}.txt'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/Motifs_analyses/HOMER_analyses/retrieve_diff_peaks.r {input.diff_res} {wildcards.comparison_pairwise} {wildcards.direction} {wildcards.logFC} {params.output_path}/DAR_HOMER_pairwise/{wildcards.comparison_pairwise}/
        """

rule annotate_diff_peaks_pairwise:
    '''
    Identify differential peaks using HOMER.
    '''
    input:
        diff_peaks = str(config["analysis_folder"]) + 'DAR_HOMER_pairwise/{comparison_pairwise}/{comparison_pairwise}_{direction}_logFC{logFC}.txt'
    output:
        diff_peaks = str(config["analysis_folder"]) + 'DAR_HOMER_pairwise/{comparison_pairwise}/annotated_{comparison_pairwise}_{direction}_logFC{logFC}.txt'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v4.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    wildcard_constraints:
        logFC = "0|1|2"
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/Motifs_analyses/HOMER_analyses/annotate_diff_peaks.r {input.diff_peaks} {wildcards.comparison_pairwise} {wildcards.direction} {wildcards.logFC} {params.output_path}/DAR_HOMER_pairwise/{wildcards.comparison_pairwise}/
        """

rule peak_clustering_pairwise:
    '''
    Cluster the peaks differentially expressed.
    '''
    input:
        counts_file = str(config["analysis_folder"]) + 'merged_data/raw_counts_merged.txt',
        metadata_file = str(config["analysis_folder"]) + 'merged_data/metadata_merged.txt',
        fasta = config['human_data']['fasta_file'],
        diff_res = rules.DAR_HOMER_pairwise.output.diff_res
    output:
        str(config["analysis_folder"]) + 'peaks_clustering_pairwise/cluster_peaks_annotation.txt',
        cluster_peaks = expand(str(config["analysis_folder"]) + 'peaks_clustering_pairwise/C{cluster_pairwise}/cluster_peaks.txt', cluster_pairwise = clusters_pairwise)
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/human_data/peak_clustering_pairwise.r {input.counts_file} {input.metadata_file} {input.diff_res} {input.fasta} {params.output_path}/peaks_clustering_pairwise/ TRUE
        """

rule peaks_pathway_enrichment_pairwise:
    '''
    Perform pathways enrichment in each cluster.
    '''
    input:
        annotated_peaks = str(config["analysis_folder"]) + 'peaks_clustering_pairwise/cluster_peaks_annotation.txt'
    output:
        str(config["analysis_folder"]) + 'pathways_enrichment_res_pairwise/{geneset}_enrichment.txt'
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
        Rscript {params.smk_path}/workflow/scripts/annotation/annotation_analysis.r {input.annotated_peaks} {params.output_path}/pathways_enrichment_res_pairwise/ {params.gs} {threads} hg38
        """

rule motif_enrichment_HOMER_pairwise:
    '''
    Identify enriched motifs .
    '''
    input:
        peaks_set = str(config["analysis_folder"]) + 'peaks_clustering_pairwise/C{cluster_pairwise}/cluster_peaks.txt',
        motif_file = config['human_data']['homer_motifs_file']
    output:
        homer_enrichment = str(config["analysis_folder"]) + 'HOMER_res_pairwise/C{cluster_pairwise}/homerMotifs.all.motifs'
    container:
        "docker://nfcore/atacseq:latest"
    params:
        homer_bin_path = "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/homer/bin/",
        output_path = config["analysis_folder"]
    threads: config['motif_enrichment_HOMER']['cpus-per-task']
    shell:
        """
        mkdir -p {params.output_path}/HOMER_res_pairwise/C{wildcards.cluster_pairwise}/;
        PATH=$PATH:{params.homer_bin_path};
        {params.homer_bin_path}findMotifsGenome.pl {input.peaks_set} hg38 {params.output_path}/HOMER_res_pairwise/C{wildcards.cluster_pairwise}/ -mset vertebrates -mknown {input.motif_file} -mcheck {input.motif_file} -mask -p {threads};
        {params.homer_bin_path}compareMotifs.pl {params.output_path}/HOMER_res_pairwise/C{wildcards.cluster_pairwise}/homerMotifs.all.motifs -known {input.motif_file} {params.output_path}/HOMER_res_pairwise/C{wildcards.cluster_pairwise}/ -cpu {threads};
        """

rule motif_enrichment_figure_pairwise:
    input:
        expand(str(config["analysis_folder"]) + 'HOMER_res_pairwise/C{cluster_pairwise}/homerMotifs.all.motifs', cluster_pairwise = clusters_pairwise),
        motif_file = config['human_data']['motif_annotations']
    output:
        str(config["analysis_folder"]) + 'HOMER_res_pairwise/all_denovo_TFs.pdf',
        str(config["analysis_folder"]) + 'HOMER_res_pairwise/homer_enrichment_results.txt'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/Motifs_analyses/HOMER_analyses/homer_de_novo.r {params.output_path}/HOMER_res_pairwise/ {input.motif_file}
        """


rule peak_related_figures:
    '''
    Upset plot with peaks overlap info + peaks annotation barplots
    '''
    input:
        diff_res = rules.DAR_HOMER.output.diff_res,
        peak_clustering = str(config["analysis_folder"]) + 'peaks_clustering/cluster_peaks_annotation.txt',
        peak_clustering_pairwise = str(config["analysis_folder"]) + 'peaks_clustering_pairwise/cluster_peaks_annotation.txt',
        counts_file = str(config["analysis_folder"]) + 'merged_data/raw_counts_merged.txt',
        metadata_file = str(config["analysis_folder"]) + 'merged_data/metadata_merged.txt',
        fasta = config['human_data']['fasta_file'],
        motif_annotations = config['human_data']['motif_annotations'],
        # all_called_peaks = config['human_data']['indiv_samples_peaks'],
        homer_res = str(config["analysis_folder"]) + 'HOMER_res_pairwise/homer_enrichment_results.txt'
    output:
        str(config["analysis_folder"]) + 'peaks_clustering/peak_clustering_annotations_withNonDARs.pdf',
        str(config["analysis_folder"]) + 'peaks_clustering_pairwise/homer_res_heatmap.pdf'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v4.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/human_data/get_peakRelated_figures.r {input.peak_clustering} {input.peak_clustering_pairwise} {input.diff_res} \
         {input.metadata_file} {input.counts_file} {input.fasta} {input.homer_res} {input.motif_annotations} {params.output_path}/peaks_clustering/
        """

rule addChromatinModule:
    input:
        counts_file = str(config["analysis_folder"]) + 'merged_data/raw_counts_merged.txt',
        metadata_file = str(config["analysis_folder"]) + 'merged_data/metadata_merged.txt'
    output:
        str(config["analysis_folder"]) + '/AddChromatinModule_results/tss_samplesScores.txt'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v4.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"],
        senescence_genes = str(config['input_data_path']) + '/Senescence_genesets_Saul_2022_Nat_Comm.xlsx',
        db_genesets = str(config['input_data_path']) + '/peaks_db_genesets_human.rds',
        mouse_mapping = str(config['input_data_path']) + '/HOM_MouseHumanSequence.rpt'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['addChromatinModule']['mem_mb'],
        time = config['addChromatinModule']['time']
    threads: config['addChromatinModule']['cpus-per-task']
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/module_scores/addchromatinscore_human.r \
        {input.metadata_file} \
        {input.counts_file} \
        {params.output_path}/AddChromatinModule_results/ \
        {params.senescence_genes} \
        {params.db_genesets} \
        {params.mouse_mapping}
        """



# mamba activate /users/agabrie4/mambaforge/envs/snakemake-8.6.0
# workflow_path=/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/snakemake_pipeline/workflow/
# snakemake --snakefile ${workflow_path}Snakefile_human.smk --profile ${workflow_path}config_files/slurm/ --configfile ${workflow_path}config_files/config_HumanData.yml --dry-run
