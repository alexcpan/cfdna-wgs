# Read trimming per NCI
rule trimmomatic:
    input:
        read1 = config["data_dir"] + "/fastq/raw/{read_id}_R1.fastq.gz",
        read2 = config["data_dir"] + "/fastq/raw/{read_id}_R2.fastq.gz",
    params:
        adapter_fasta = config["adapter_fastq"],
    output:
        read1 = config["data_dir"] + "/fastq/processed/{read_id}_proc_R1.fastq.gz",
        read1_unpr = config["data_dir"] + "/fastq/unpaired/{read_id}_unpr_R1.fastq.gz",
        read2 = config["data_dir"] + "/fastq/processed/{read_id}_proc_R2.fastq.gz",
        read2_unpr = config["data_dir"] + "/fastq/unpaired/{read_id}_unpr_R2.fastq.gz",
    log:
        int = config["data_dir"] + "/logs/trimmomatic_trimlog_{read_id}.log",
        main = config["data_dir"] + "/logs/trimmomatic_{read_id}.log",
    shell:
        """
        trimmomatic PE \
                    -threads {config[threads]} \
                    -trimlog {log.int} \
                    {input.read1} {input.read2} \
                    {output.read1} {output.read1_unpr} \
                    {output.read2} {output.read2_unpr} \
                    ILLUMINACLIP:{params.adapter_fasta}:2:30:10 \
                    LEADING:10 TRAILING:10 MAXINFO:50:0.97 MINLEN:20 &> {log.main}
        """

# BWA alignment
rule align:
    input:
        read1 = config["data_dir"] + "/fastq/processed/{read_id}_proc_R1.fastq.gz",
        read2 = config["data_dir"] + "/fastq/processed/{read_id}_proc_R2.fastq.gz",
    output:
        config["data_dir"] + "/bam/{read_id}.sam",
    log:
        config["data_dir"] + "/logs/align_{read_id}.log"
    shell:
        """
        bwa mem -M -t 4 {config[bwa_index]} {input.read1} {input.read2} > {output}
	"""

# FastQC
rule fastqc:
    input:
        raw =  config["data_dir"] + "/fastq/raw/{read_id}_{read}.fastq.gz",
        proc = config["data_dir"] + "/fastq/processed/{read_id}_proc_{read}.fastq.gz",
    params:
        out_dir = config["data_dir"] + "/qc",
    output:
        raw_html = config["data_dir"] + "/qc/{read_id}_{read}_fastqc.html",
        proc_html = config["data_dir"] + "/qc/{read_id}_proc_{read}_fastqc.html",
    log:
        raw = config["data_dir"] + "/logs/fastqc_raw_{read_id}_{read}.log",
        proc = config["data_dir"] + "/logs/fastqc_proc_{read_id}_{read}.log",
    shell:
        """
        fastqc --outdir {params.out_dir} \
        --quiet \
        --threads {config[threads]} {input.raw} &> {log}
        fastqc --outdir {params.out_dir} \
        --quiet \
        --threads {config[threads]} {input.proc} &> {log}
        """

# Alignment deduplication and sorting
rule alignment_processing:
    input:
        config["data_dir"] + "/bam/{read_id}.sam",
    output:
        bam = config["data_dir"] + "/bam/{read_id}_raw.bam",
        dedup = temp(config["data_dir"] + "/bam/{read_id}_dedup_unsort.bam"),
        sort = config["data_dir"] + "/bam/{read_id}_dedup.bam",
        index = config["data_dir"] + "/bam/{read_id}_dedup.bam.bai",
    log:
        config["data_dir"] + "/logs/alignment_processing_{read_id}.log"
    shell:
        """
        sambamba view -t {config[threads]} -S -f bam {input} > {output.bam}
        sambamba markdup -r -t {config[threads]} {output.bam} {output.dedup}
        sambamba sort -t {config[threads]} {output.dedup} -o {output.sort}
        sambamba index -t {config[threads]} {output.sort}
        """

# Alignment samtools QC
rule alignment_qc:
    input:
        config["data_dir"] + "/bam/{read_id}_{bam_step}.bam",
    output:
        samstat = config["data_dir"] + "/qc/{read_id}_{bam_step}_samstats.txt",
        flagstat = config["data_dir"] + "/qc/{read_id}_{bam_step}_flagstat.txt",
    log:
        config["data_dir"] + "/logs/alignment_qc_{read_id}_{bam_step}.err",
    shell:
        """
        samtools stats {input} > {output.samstat} 2>{log}
        samtools flagstat {input} > {output.flagstat} 2>>{log}
        """

# Alignment downsampling
rule downsample_bams:
    input:
        config["data_dir"] + "/bam/{read_id}_dedup.bam",
    output:
        config["data_dir"] + "/bam/{read_id}_ds{milreads}.bam",
    log:
        config["data_dir"] + "/logs/downsample_bam_{read_id}_{milreads}.err"
    shell:
        """
        {config[cfdna_wgs_script_dir]}/downsample_bam.sh {input} {wildcards.milreads}000000 {output} 2>{log}
        """
