import pandas as pd
from itertools import chain
from itertools import combinations

# Read the metadata file 
metadata_df = pd.read_csv(config['data']['metadata'], sep = "\t", header = 0)

conditions = metadata_df['groups'].unique()
samples = metadata_df['sample_name']
narrowPeaks_files = metadata_df['summits']

# Create a dictionary to map each sample to its condition
sample_to_condition = dict(zip(metadata_df['sample_name'], metadata_df['groups']))

rule all:
    input:
        str(config["analysis_folder"]) + '/chromVAR/dev_scores.rds'

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

rule generate_chromVAR_related_figures:
    input:
        metadata_file = config['data']['metadata'],
        motif_annotations = config['data']['motif_annotations'],
        dev = str(config["analysis_folder"]) + '/chromVAR/dev_scores.rds',
        zscore = str(config["analysis_folder"]) + '/chromVAR/z_scores.rds'
    output:
        str(config["analysis_folder"]) + '/chromVAR/chromVar_TOBIASchange-vs-deviation_SiglecFhiVSSiglecFlo.pdf',
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
    params:
        output_path = config["analysis_folder"],
        smk_path = config["snakemake_path"]
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/pancreas_analyses/get_chromvar_figures.r {input.dev} {input.zscore} {input.metadata_file} {input.motif_annotations} {params.output_path}/chromVAR/
        """
# workflow_path=/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/snakemake_pipeline/workflow/
# snakemake --snakefile ${workflow_path}Snakefile_Pancreas_data.smk --profile ${workflow_path}config_files/slurm/ --configfile ${workflow_path}config_files/config_Pancreas.yml --dry-run
