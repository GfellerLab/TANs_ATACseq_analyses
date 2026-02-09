from Parse_config_pepatac import *
from pathlib import Path

cfg = ParsedConfig(config)

rule all:
    input:
        cfg.get_pepatac_output_files(config),
        str(cfg.OUTDIR) + f"/consensus_peaks/raw_counts.txt"


# Download all necessary files to run PEPATAC
# rule PEPATAC_inputs_download:
#     # singularity: str(config["snakemake_path"]) + "config/atac_bulk_processing_pepatac2_latest.sif"
#     singularity: str(config["snakemake_path"]) + "/workflow/config_files/pepatac-longer-reads_latest.sif"
#     output:
#         str(cfg.genome_folder / "alias" / cfg.genome_folder / "ensembl_rb/default" / f"{cfg.genome}.gff.gz")
#     shell:
#         """
#         genome={cfg.genome}
#         refgenie init -c {cfg.genome_folder}genome_config.yaml
#         export REFGENIE={cfg.genome_folder}genome_config.yaml
#
#         if [ $genome = "hg38" ]
#         then
#             refgenie pull hg38/fasta hg38/bowtie2_index hg38/refgene_anno hg38/ensembl_gtf hg38/ensembl_rb
#             refgenie pull rCRSd/fasta rCRSd/bowtie2_index human_repeats/fasta human_repeats/bowtie2_index
#         elif [ $genome = 'mm10' ]
#         then
#             refgenie pull mm10/fasta mm10/bowtie2_index mm10/refgene_anno mm10/ensembl_gtf mm10/ensembl_rb
#             refgenie pull mouse_chrM2x/fasta mouse_chrM2x/bowtie2_index
#         fi
#         """

rule SRA_download:
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/pepatac-longer-reads_latest.sif"
    output:
        str(cfg.fastq_path / "{sample}_1.fastq"),
        str(cfg.fastq_path / "{sample}_2.fastq")
    params:
        download_path = cfg.fastq_path,
        ncbi_settings = config["ncbi_config"]
    shell:
        """
        FILE={params.download_path}/{wildcards.sample}_1.fastq
        echo $FILE
        if [ -f "$FILE" ]; then
            echo "already downloaded"
        else
            export NCBI_SETTINGS={params.ncbi_settings}
            prefetch -v {wildcards.sample} -O {params.download_path} --max-size 50000000000
            fasterq-dump -f -O {params.download_path} {params.download_path}/{wildcards.sample}/{wildcards.sample}.sra
            vdb-validate -x {params.download_path}/{wildcards.sample}/{wildcards.sample}.sra
        fi
        """


rule run_PEPATAC:
    '''
    Run PEPATAC preprocessing
    '''
    input:
        fastq1 = str(cfg.fastq_path / f"{{sample}}{cfg.fastq_suffix1}"),
        fastq2 = str(cfg.fastq_path / f"{{sample}}{cfg.fastq_suffix2}")
    output:
        str(cfg.OUTDIR) + f"/{{sample}}/aligned_{cfg.genome}/{{sample}}_sort_dedup.bam",
        str(cfg.OUTDIR) + f"/{{sample}}/peak_calling_{cfg.genome}/{{sample}}_peaks_normalized.narrowPeak",
        str(cfg.OUTDIR) + f"/{{sample}}/peak_calling_{cfg.genome}/{{sample}}_peaks_rmBlacklist.narrowPeak"
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/pepatac-longer-reads_latest.sif"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['PEPATAC_params']['mem_mb'],
        time = config['PEPATAC_params']['time']
    threads: config['PEPATAC_params']['cpus-per-task']
    params:
        genome_version = config['PEPATAC_params']['genome_version'],
        genome_folder = config['PEPATAC_params']['genome_folder'],
        genome_size = config['PEPATAC_params']['genome_size'],
        blacklist_file = config['PEPATAC_params']['blacklist_file'],
        yaml_file = cfg.pepatac_config,
        output_path = cfg.OUTDIR

    shell:
        """
        genome={params.genome_version}

        if [ "$genome" = "hg38" ] || [ "$genome" = "hg19" ]; then
            prealignment_params="--prealignment-index rCRSd={params.genome_folder}alias/rCRSd/bowtie2_index/default/rCRSd human_repeats={params.genome_folder}alias/human_repeats/bowtie2_index/default/human_repeats"
        elif [ "$genome" = "mm10" ]; then
            prealignment_params="--prealignment-index mouse_chrM2x={params.genome_folder}alias/mouse_chrM2x/bowtie2_index/default/mouse_chrM2x"
        else
            prealignment_params=""
        fi

        pepatac.py --single-or-paired paired  -C {params.yaml_file} \
        $prealignment_params \
        --genome {params.genome_version} --genome-index {params.genome_folder}alias/{params.genome_version}/bowtie2_index/default/. \
        -gs {params.genome_size} \
        --chrom-sizes {params.genome_folder}alias/{params.genome_version}/fasta/default/{params.genome_version}.chrom.sizes \
        --sample-name {wildcards.sample} --input {input.fastq1} --input2 {input.fastq2}  -O {params.output_path} --peak-caller macs2 \
        --trimmer trimmomatic \
        --aligner bowtie2 \
        --deduplicator samtools \
        -P {threads} \
        --TSS-name {params.genome_folder}alias/{params.genome_version}/refgene_anno/default/{params.genome_version}_TSS.bed \
        --blacklist {params.blacklist_file}

        prefix={params.output_path}/{wildcards.sample}/aligned_{params.genome_version}/{wildcards.sample}_sort
        samtools flagstat $prefix.bam > $prefix.bam.flagstat
        samtools idxstats $prefix.bam > $prefix.bam.idxstats
        samtools stats $prefix.bam > $prefix.bam.stats

        prefix={params.output_path}/{wildcards.sample}/aligned_{params.genome_version}/{wildcards.sample}_sort_dedup
        samtools flagstat $prefix.bam > $prefix.bam.flagstat
        samtools idxstats $prefix.bam > $prefix.bam.idxstats
        samtools stats $prefix.bam > $prefix.bam.stats

        prefix={params.output_path}/{wildcards.sample}/prealignments/
        if [ $genome = "hg38" -o $genome = "hg19" ]
        then
            fastqc --noextract --outdir $prefix $prefix*human_repeats_unmap_R1.fq.gz
            fastqc --noextract --outdir $prefix $prefix*human_repeats_unmap_R2.fq.gz
        elif [ $genome = 'mm10' ]
        then
            fastqc --noextract --outdir $prefix $prefix*mouse_chrM2x_unmap_R1.fq.gz
            fastqc --noextract --outdir $prefix $prefix*mouse_chrM2x_unmap_R2.fq.gz
        fi
        rm -r {params.output_path}/{wildcards.sample}/prealignments/
        rm -r {params.output_path}/{wildcards.sample}/fastq/
        """

rule run_fastqc:
    '''
    Run PEPATAC preprocessing
    '''
    input:
        fastq1 = str(cfg.fastq_path / f"{{sample}}{cfg.fastq_suffix1}"),
        fastq2 = str(cfg.fastq_path / f"{{sample}}{cfg.fastq_suffix2}")
    output:
        str(cfg.OUTDIR) + f"/fastqc_results/{{sample}}_fastqc.html"
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/pepatac-longer-reads_latest.sif"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['PEPATAC_params']['mem_mb'],
        time = '1:00:00'
    params:
        output_path = cfg.OUTDIR

    shell:
        """
        mkdir -p {params.output_path}fastqc_results/
        fastqc {input.fastq1} -o {params.output_path}fastqc_results/
        fastqc {input.fastq2} -o {params.output_path}fastqc_results/
        """

rule run_multiQC:
    input:
        path ('fastqc/*') from fastqc_res.collect()
        path ('alignment/filtered/*') from filtered_bam_flagstat_mqc.collect()
        path ('alignment/dedup/*') from dedup_bam_flagstat_mqc.collect()
        path ('preseq/*') from preseq_mqc.collect()
        path ('macs2/*') from macs2_mqc.collect()
        path ('fastqc_filtered/*') from fq_filtered_files.collect()
        file mqc_config from Channel.fromPath(params.multiqc_config).collect()
    output:
      file("multiqc_report.html") into multiqc_post
      file("multiqc_report_data") into multiqc_post_data

    shell:
    """
    multiqc -v . -n multiqc_report.html --comment "ATAC-Seq QC report" --config multiqc_config.yaml
    """

rule run_bam_shift:
    input:
        bam = str(cfg.OUTDIR) + f"/{{sample}}/aligned_{cfg.genome}/{{sample}}_sort_dedup.bam"
    resources:
        time = config['run_bam_shift']['time'],
        mem_mb = lambda wildcards, attempt: attempt * config['run_bam_shift']['mem_mb']
    threads: config['run_bam_shift']['cpus-per-task']
    output:
        str(cfg.OUTDIR) + f"/{{sample}}/aligned_{cfg.genome}/{{sample}}_sort_dedup_shifted.bam"

    params:
        genome_version = config['PEPATAC_params']['genome_version'],
        output_path = cfg.OUTDIR,
        smk_path = config["snakemake_path"]
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_latest.sif"
    shell:
        """
        alignmentSieve -b {input.bam} \
        -o {params.output_path}/{wildcards.sample}/aligned_{params.genome_version}/{wildcards.sample}_shifted.bam \
         --minMappingQuality 10 -p 10 --ATACshift --filterMetrics {params.output_path}/{wildcards.sample}/aligned_{params.genome_version}/shifting_log.txt
        samtools sort {params.output_path}/{wildcards.sample}/aligned_{params.genome_version}/{wildcards.sample}_shifted.bam -o {params.output_path}/{wildcards.sample}/aligned_{params.genome_version}/{wildcards.sample}_sort_dedup_shifted.bam
        rm {params.output_path}/{wildcards.sample}/aligned_{params.genome_version}/{wildcards.sample}_shifted.bam
        samtools index {params.output_path}/{wildcards.sample}/aligned_{params.genome_version}/{wildcards.sample}_sort_dedup_shifted.bam
        """


rule get_consensus_peaks:
    input:
        expand("{dataPath}/{sample}/peak_calling_{genomeV}/{sample}_peaks_normalized.narrowPeak", sample=cfg.Sample_IDs, dataPath=cfg.OUTDIR, genomeV=cfg.genome),
    resources:
        time = config['get_consensus_peaks']['time']
    threads: config['get_consensus_peaks']['cpus-per-task']
    output:
        str(cfg.OUTDIR) + f"/consensus_peaks/peaks.saf"
    params:
        genome_version = config['PEPATAC_params']['genome_version'],
        output_path = cfg.OUTDIR,
        genome_folder = config['PEPATAC_params']['genome_folder'],
        blacklist_file = config['PEPATAC_params']['blacklist_file'],
        chrom_sizes = config['PEPATAC_params']['chrom_sizes'],
        config = config['PEPATAC_params']['config'],
        sampleList = config[cfg.DATASET]["Sample_IDs_list"],
        data_path = cfg.OUTDIR,
        smk_path = config["snakemake_path"],
        dataset = cfg.DATASET,
        KP_lung_path = config["KP_path"]
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atac_bulk_processing_pepatac2_latest.sif"
    shell:
        """
        dataset={params.dataset}
        if [ "$dataset" = "Human_bulk_atac" ]; then
            Rscript {params.smk_path}/workflow/scripts/human_data/human_peak_consensus.r 3 human {params.config} {params.chrom_sizes} {params.blacklist_file} {params.data_path} {params.smk_path}/{params.sampleList} {params.output_path}/consensus_peaks/ bed
        elif [ "$dataset" = "KP_lung" ]; then
            Rscript {params.smk_path}/workflow/scripts/preprocessing/get_consensus_peaks_KPlung.r 2 mouse {params.config} {params.chrom_sizes} {params.blacklist_file} {params.data_path} {params.smk_path}/{params.sampleList} {params.output_path}/consensus_peaks/ bed
        elif [ "$dataset" = "NfatKO" ]; then
            Rscript {params.smk_path}/workflow/scripts/Nfactc1_KO/get_consensus_peaks.r 2 mouse {params.config} {params.chrom_sizes} {params.blacklist_file} {params.data_path} {params.smk_path}/{params.sampleList} {params.output_path}/consensus_peaks/ bed
            mkdir -p {params.output_path}/consensus_peaks_KP/
            cp {params.KP_lung_path}/consensus_peaks/peaks.saf {params.output_path}/consensus_peaks_KP/
        else
            Rscript {params.smk_path}/workflow/scripts/pancreas_analyses/get_consensus_peaks.r 2 mouse {params.config} {params.chrom_sizes} {params.blacklist_file} {params.data_path} {params.smk_path}/{params.sampleList} {params.output_path}/consensus_peaks/ bed
            mkdir -p {params.output_path}/consensus_peaks_KP/
            cp {params.KP_lung_path}/consensus_peaks/peaks.saf {params.output_path}/consensus_peaks_KP/
        fi
        """

rule get_counts:
    input:
        expand("{dataPath}/{sample}/aligned_{genomeV}/{sample}_sort_dedup.bam", sample=cfg.Sample_IDs, dataPath=cfg.OUTDIR, genomeV=cfg.genome),
        saf = str(cfg.OUTDIR) + f"/consensus_peaks/peaks.saf"
    resources:
        time = config['get_counts_params']['time']
    threads: config['get_counts_params']['cpus-per-task']

    output:
        str(cfg.OUTDIR) + f"/consensus_peaks/featureCounts.txt",
        str(cfg.OUTDIR) + f"/consensus_peaks/featureCounts.txt.summary"
    params:
        genome_version = config['PEPATAC_params']['genome_version'],
        output_path = cfg.OUTDIR
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atac_bulk_processing_pepatac2_latest.sif"
    shell:
        """
        bam_files={params.output_path}/*/aligned_{params.genome_version}/*_sort_dedup.bam
        echo $bam_files
        featureCounts -F SAF -O --fracOverlap 0.2 -T 5 -p -a {input.saf} -o {params.output_path}/consensus_peaks/featureCounts.txt $bam_files
        """

rule save_counts_matrix:
    input:
        counts = str(cfg.OUTDIR) + f"/consensus_peaks/featureCounts.txt"
    output:
        str(cfg.OUTDIR) + f"/consensus_peaks/raw_counts.txt"
    params:
        smk_path = config["snakemake_path"],
        output_path = cfg.OUTDIR
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atac_bulk_processing_pepatac2_latest.sif"
    shell:
        """
        Rscript {params.smk_path}/workflow/scripts/preprocessing/get_counts_matrix.r {input.counts} {params.output_path}/consensus_peaks/
        """


# Snakemake command:
# workflow_path=/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/collaborations/Mate_Pittet_lab/snakemake_pipeline/workflow/
# snakemake --snakefile ${workflow_path}pepatac_processing.smk --profile ${workflow_path}config_files/slurm/ --configfile ${workflow_path}config_files/config_pepatac_XX.yml
