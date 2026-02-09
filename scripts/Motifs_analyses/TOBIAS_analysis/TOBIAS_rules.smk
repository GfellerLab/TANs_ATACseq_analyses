
rule merge_bams:
    '''
    Merge bam files in order to run footprinting in each condition.
    '''
    input:
        lambda wildcards: shifted_bam_files_by_condition[wildcards.condition]
    output:
        str(config["analysis_folder"]) + '/TOBIAS_res/{condition}.merged.bam'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
    params:
        output_path = config["analysis_folder"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['footprinting']['mem_mb'],
        time = config['footprinting']['time']
    threads: config['footprinting']['cpus-per-task']
    shell:
        """
        samtools merge --threads {threads} -f {output} {input};
        samtools index {output}
        """

rule merge_peaks:
    '''
    Merge peaks called in each sample using MACS2 narrowPeaks files.
    '''
    input:
        narrowPeaks_files
    output:
        merged_peaks = str(config["analysis_folder"]) +'/TOBIAS_res/merged_peaks.bed'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
    params:
        output_path = config["analysis_folder"]
    shell:
        """
        cat {input} | sort -k1,1 -k2,2n | bedtools merge > {params.output_path}/TOBIAS_res/merged_peaks.bed
        """

rule footprinting:
    '''
    Run footprinting.
    '''
    input:
        merged_bam_files = str(config["analysis_folder"]) +'/TOBIAS_res/{condition}.merged.bam',
        merged_peaks = str(config["analysis_folder"]) +'/TOBIAS_res/merged_peaks.bed',
        genome = config['data']['tobias_fasta'],
        peak_file = config['data']['peak_file'],
        motifs_file = config['data']['tobias_motifs_file']
    output:
        footprinting_scores = str(config["analysis_folder"]) +'/TOBIAS_res/footprinting_scores/{condition}_footprints.bw'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
    params:
        output_path = config["analysis_folder"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['footprinting']['mem_mb'],
        time = config['footprinting']['time']
    threads: config['footprinting']['cpus-per-task']
    shell:
        """
        export MPLCONFIGDIR={params.output_path}/TOBIAS_res/;
        TOBIAS ATACorrect --bam {input.merged_bam_files} --genome {input.genome} --peaks {input.merged_peaks} --outdir {params.output_path}/TOBIAS_res/corrected_signals/ --cores {threads}
        TOBIAS FootprintScores --signal {params.output_path}/TOBIAS_res/corrected_signals/{wildcards.condition}.merged_corrected.bw --regions {input.merged_peaks} --output {params.output_path}/TOBIAS_res/footprinting_scores/{wildcards.condition}_footprints.bw --cores {threads};
        mkdir -p {params.output_path}/TOBIAS_res/footprinting_scores/{wildcards.condition}_bindetect_results/;
        TOBIAS BINDetect --signals {params.output_path}/TOBIAS_res/footprinting_scores/{wildcards.condition}_footprints.bw --motifs {input.motifs_file} --genome {input.genome} --peaks {input.peak_file} --outdir {params.output_path}/TOBIAS_res/footprinting_scores/{wildcards.condition}_bindetect_results/ --cores {threads}
        """

rule diff_footprinting:
    '''
    Run differential footprinting.
    '''
    input:
        expand(str(config["analysis_folder"]) + '/TOBIAS_res/footprinting_scores/{condition}_footprints.bw', condition=conditions),
        genome = config['data']['tobias_fasta'],
        peak_file = config['data']['peak_file'],
        motifs_file = config['data']['tobias_motifs_file']
    output:
        str(config["analysis_folder"]) + '/TOBIAS_res/BINDetect_diff_res/{comp_pair}/bindetect_results.txt'
    singularity: str(config["snakemake_path"]) + "/workflow/config_files/atacseq_pipeline_v3.sif"
    params:
        output_path = config["analysis_folder"],
        cond1 = lambda wildcards: wildcards.comp_pair.split('_VS_')[0],
        cond2 = lambda wildcards: wildcards.comp_pair.split('_VS_')[1],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['diff_footprinting']['mem_mb'],
        time = config['diff_footprinting']['time']
    threads: config['diff_footprinting']['cpus-per-task']
    shell:
        """
        export MPLCONFIGDIR={params.output_path}/TOBIAS_res/;
        TOBIAS BINDetect --signals {params.output_path}/TOBIAS_res/footprinting_scores/{params.cond1}_footprints.bw \
        {params.output_path}/TOBIAS_res/footprinting_scores/{params.cond2}_footprints.bw \
        --cond_names {params.cond1} {params.cond2} --motifs {input.motifs_file} \
        --peaks {input.peak_file} --genome {input.genome} \
        --outdir {params.output_path}/TOBIAS_res/BINDetect_diff_res/{wildcards.comp_pair} --cores {threads}
        """

# uropa --bed {input.peak_file} --gtf {params.gtf_file} --show_attributes gene_id gene_name --feature_anchor start --distance 20000 10000 --feature gene
# cut -f 1-6,16-17 peaks_finalhits.txt | tail -n +2 > peaks_annotated.bed
