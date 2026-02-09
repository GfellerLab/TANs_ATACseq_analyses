This repository contains the code used to perform the analyses shown in Kiss, Gabriel et al..

The raw ATAC-Seq data were processed using the snakemake pipeline pepatac_processing.smk using the following command line:

```
snakemake --snakefile ${workflow_path}pepatac_processing.smk --profile ${workflow_path}config_files/slurm/ --configfile ${workflow_path}config_files/config_pepatac_XX.yml
# config_pepatac_XX.yml was replaced by the conig file corresponding to each dataset (config_pepatac_KPData.yml, config_files/config_pepatac_PancreasData.yml, config_files/config_pepatac_NfatKO.yml,config_files/config_pepatac_humanData.yml)
```

After ATAC-Seq pre-processing, the KP lung ATAC-Seq data were analyzed using the snakemake pipeline Snakefile_KP_lung.smk using the following command line:

```
snakemake --snakefile ${workflow_path}Snakefile_KP_lung.smk --profile ${workflow_path}config_files/slurm/ --configfile ${workflow_path}config_files/config_KPlung.yml
```

After ATAC-Seq pre-processing, the Human ATAC-Seq data were analyzed using the snakemake pipeline Snakefile_human.smk using the following command line:

```
snakemake --snakefile ${workflow_path}Snakefile_human.smk --profile ${workflow_path}config_files/slurm/ --configfile ${workflow_path}config_files/config_HumanData.yml 
```


After ATAC-Seq pre-processing, the NfatKO ATAC-Seq data were analyzed using the snakemake pipeline Snakefile_NfatKO_data.smk using the following command line:

```
snakemake --snakefile ${workflow_path}Snakefile_NfatKO_data.smk --profile ${workflow_path}config_files/slurm/ --configfile ${workflow_path}config_files/config_NfatKO.yml 
```

After ATAC-Seq pre-processing, the pancreas data (public data re-analyzed) were analyzed using the snakemake pipeline Snakefile_Pancreas_data.smk using the following command line:

```
snakemake --snakefile ${workflow_path}Snakefile_Pancreas_data.smk --profile ${workflow_path}config_files/slurm/ --configfile ${workflow_path}config_files/config_Pancreas.yml
```

The prostate 10xmultiome data (public data re-analyzed) were analyzed using the snakemake pipeline Snakefile_prostate_atlas.smk using the following command line:

```
snakemake --snakefile ${workflow_path}Snakefile_prostate_atlas.smk --profile ${workflow_path}config_files/slurm/ --configfile ${workflow_path}config_files/config_ProstateData.yml
```
