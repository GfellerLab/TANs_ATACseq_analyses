from pathlib import Path
import os
import re
import pandas as pd
import random
SEED = 448

print(config)

prostateSamples = list(config["prostateHanSamples"])
cleanProstateSamples = list(config["cleanProstateHanSamples"])
pathCellRangerProstate = config["pathCellRangerProstateData"]
celltype_values = ["Luminal", "Basal", "Neuroendocrine", "Neuron", "seminalVesicle", "Mesenchymal1", "Mesenchymal2","Endothelial", "Neutrophil", "Tcell", "Bcell"]
wdir = os.getcwd()
rule all:
  input:
    expand(str(config["analysis_folder"]) + "output/{sample}/metacells.rds",sample = prostateSamples),
    expand(str(config["analysis_folder"]) + "output/{clean_sample}/metacells.rds",clean_sample = cleanProstateSamples),
    str(config["analysis_folder"]) + "output_signac_q0/combined.metacells.rds",
    str(config["analysis_folder"]) +"output_signac_q0/neutrophils.rds",
    str(config["analysis_folder"]) +"output_signac_q0/neutrophils_withKPpeaks_chromVAR.rds"


rule metacellPerSample:
  input: fragment = str(config["raw_data_path"]) + "/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/{sample}/atac_fragments.tsv.gz",
         index = str(config["raw_data_path"]) + "/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/{sample}/atac_fragments.tsv.gz.tbi",
         h5 = str(config["raw_data_path"]) + "/sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/{sample}/filtered_feature_bc_matrix.h5",
         consensusPeaks = str(config["consensus_peaks"])
  output: str(config["analysis_folder"]) +"output/{sample}/metacells.rds"
  threads: 10
  params:
      raw_data_path = config["raw_data_path"],
      output_path = config["analysis_folder"],
      smk_path = config["snakemake_path"]
  singularity: str(config["snakemake_path"]) + "workflow/config_files/supercell_multiomics_v2.sif"
  resources:
      mem_mb = config["MC_build"]["mem_mb"],
      time = config["MC_build"]["time"]
  shell: "Rscript {params.smk_path}/workflow/scripts/prostate_data/atlas_generation/metacellPerSampleCL.R -i {input.h5} \
         -f {input.fragment} \
         -o {params.output_path}/output/{wildcards.sample}/ \
         -c {input.consensusPeaks} \
         -g 10"

rule Integrating_metacells_RNA_ATAC_prostateCancer:
  input: expand(str(config["analysis_folder"]) + "output/{clean_sample}/metacells.rds",clean_sample = cleanProstateSamples)
  output: str(config["analysis_folder"]) + "output_signac_q0/combined.metacells.rds"
  singularity: str(config["snakemake_path"]) + "workflow/config_files/supercell_multiomics_v2.sif"
  params:
      raw_data_path = config["raw_data_path"],
      output_path = config["analysis_folder"],
      smk_path = config["snakemake_path"]
  resources:
      mem_mb = config["Integrating_metacells_RNA_ATAC_prostateCancer"]["mem_mb"],
      time = config["Integrating_metacells_RNA_ATAC_prostateCancer"]["time"]
  shell: "Rscript {params.smk_path}/workflow/scripts/prostate_data/atlas_generation/Integrating_metacells_RNA_ATAC_prostateCancer.R \
  -o {params.output_path}/output_signac_q0/ -i signac -t q0"

rule annotate_metacell_cluster:
  input: rds = str(config["analysis_folder"]) +"output_signac_q0/combined.metacells.rds"
  output: str(config["analysis_folder"]) +"output_signac_q0/neutrophils.rds"
  singularity: str(config["snakemake_path"]) + "workflow/config_files/supercell_multiomics_v2.sif"
  params:
      raw_data_path = config["raw_data_path"],
      output_path = config["analysis_folder"],
      smk_path = config["snakemake_path"]
  resources:
      mem_mb = config["Integrating_metacells_RNA_ATAC_prostateCancer"]["mem_mb"],
      time = config["Integrating_metacells_RNA_ATAC_prostateCancer"]["time"]
  shell: "Rscript {params.smk_path}/workflow/scripts/prostate_data/atlas_generation/cell_annotation.r -i {input.rds} -o {params.output_path}output_signac_q0/"

rule neutrophil_analysis:
  input: str(config["analysis_folder"]) +"output_signac_q0/neutrophils.rds"
  output: str(config["analysis_folder"]) +"output_signac_q0/neutrophils_withKPpeaks_chromVAR.rds"
  singularity: str(config["snakemake_path"]) + "workflow/config_files/atacseq_pipeline_v4.sif"
  params:
      raw_data_path = config["raw_data_path"],
      output_path = config["analysis_folder"],
      smk_path = config["snakemake_path"]
  threads: 10
  resources:
      mem_mb = config["Integrating_metacells_RNA_ATAC_prostateCancer"]["mem_mb"],
      time = config["Integrating_metacells_RNA_ATAC_prostateCancer"]["time"]
  shell: "Rscript {params.smk_path}/workflow/scripts/prostate_data/atlas_generation/neutrophil_analysis.r"

# snakemake --snakefile ${workflow_path}Snakefile_prostate_atlas.smk --profile ${workflow_path}config_files/slurm/ --configfile ${workflow_path}config_files/config_ProstateData.yml
