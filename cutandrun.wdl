version 1.0

workflow cutandrun {
	input{
		File samplesheet
		String outdir = "./results"
		String? email
		String? multiqc_title
		Boolean? save_reference
		Boolean? save_merged_fastq
		Boolean? save_trimmed
		Boolean? save_spikein_aligned
		Boolean? save_unaligned
		Boolean? save_align_intermed
		String aligner = "bowtie2"
		Int? clip_r1
		Int? clip_r2
		Int? three_prime_clip_r1
		Int? three_prime_clip_r2
		Int? trim_nextseq
		String? genome
		String? bowtie2
		String? gtf
		String? blacklist
		String spikein_genome = "K12-MG1655"
		String? spikein_bowtie2
		String? spikein_fasta
		File? fasta
		String igenomes_base = "s3://ngi-igenomes/igenomes"
		Boolean? igenomes_ignore
		Int? minimum_alignment_q_score
		Boolean? dedup_target_reads
		Int normalisation_c = 10000
		Boolean igg_control = true
		Float peak_threshold = 0.05
		String? gene_bed
		Int replicate_threshold = 1
		Boolean? only_genome
		Boolean? only_input
		Boolean? only_preqc
		Boolean? only_alignment
		Boolean? only_filtering
		Boolean? only_peak_calling
		Boolean? skip_fastqc
		Boolean? skip_trimming
		Boolean? skip_removeduplicates
		Boolean? skip_scale
		Boolean? skip_reporting
		Boolean? skip_igv
		Boolean? skip_heatmaps
		Boolean? skip_multiqc
		Boolean? skip_upset_plots
		String custom_config_version = "master"
		String custom_config_base = "https://raw.githubusercontent.com/nf-core/configs/master"
		String? config_profile_name
		String? config_profile_description
		String? config_profile_contact
		String? config_profile_url
		Int max_cpus = 16
		String max_memory = "128.GB"
		String max_time = "240.h"
		Boolean? help
		String publish_dir_mode = "copy"
		String? email_on_fail
		Boolean? plaintext_email
		String max_multiqc_email_size = "25.MB"
		Boolean? monochrome_logs
		String? multiqc_config
		String tracedir = "./results/pipeline_info"
		Boolean validate_params = true
		Boolean? show_hidden_params
		Boolean? enable_conda
		Boolean? singularity_pull_docker_container

	}

	call make_uuid as mkuuid {}
	call touch_uuid as thuuid {
		input:
			outputbucket = mkuuid.uuid
	}
	call run_nfcoretask as nfcoretask {
		input:
			samplesheet = samplesheet,
			outdir = outdir,
			email = email,
			multiqc_title = multiqc_title,
			save_reference = save_reference,
			save_merged_fastq = save_merged_fastq,
			save_trimmed = save_trimmed,
			save_spikein_aligned = save_spikein_aligned,
			save_unaligned = save_unaligned,
			save_align_intermed = save_align_intermed,
			aligner = aligner,
			clip_r1 = clip_r1,
			clip_r2 = clip_r2,
			three_prime_clip_r1 = three_prime_clip_r1,
			three_prime_clip_r2 = three_prime_clip_r2,
			trim_nextseq = trim_nextseq,
			genome = genome,
			bowtie2 = bowtie2,
			gtf = gtf,
			blacklist = blacklist,
			spikein_genome = spikein_genome,
			spikein_bowtie2 = spikein_bowtie2,
			spikein_fasta = spikein_fasta,
			fasta = fasta,
			igenomes_base = igenomes_base,
			igenomes_ignore = igenomes_ignore,
			minimum_alignment_q_score = minimum_alignment_q_score,
			dedup_target_reads = dedup_target_reads,
			normalisation_c = normalisation_c,
			igg_control = igg_control,
			peak_threshold = peak_threshold,
			gene_bed = gene_bed,
			replicate_threshold = replicate_threshold,
			only_genome = only_genome,
			only_input = only_input,
			only_preqc = only_preqc,
			only_alignment = only_alignment,
			only_filtering = only_filtering,
			only_peak_calling = only_peak_calling,
			skip_fastqc = skip_fastqc,
			skip_trimming = skip_trimming,
			skip_removeduplicates = skip_removeduplicates,
			skip_scale = skip_scale,
			skip_reporting = skip_reporting,
			skip_igv = skip_igv,
			skip_heatmaps = skip_heatmaps,
			skip_multiqc = skip_multiqc,
			skip_upset_plots = skip_upset_plots,
			custom_config_version = custom_config_version,
			custom_config_base = custom_config_base,
			config_profile_name = config_profile_name,
			config_profile_description = config_profile_description,
			config_profile_contact = config_profile_contact,
			config_profile_url = config_profile_url,
			max_cpus = max_cpus,
			max_memory = max_memory,
			max_time = max_time,
			help = help,
			publish_dir_mode = publish_dir_mode,
			email_on_fail = email_on_fail,
			plaintext_email = plaintext_email,
			max_multiqc_email_size = max_multiqc_email_size,
			monochrome_logs = monochrome_logs,
			multiqc_config = multiqc_config,
			tracedir = tracedir,
			validate_params = validate_params,
			show_hidden_params = show_hidden_params,
			enable_conda = enable_conda,
			singularity_pull_docker_container = singularity_pull_docker_container,
			outputbucket = thuuid.touchedbucket
            }
		output {
			Array[File] results = nfcoretask.results
		}
	}
task make_uuid {
	meta {
		volatile: true
}

command <<<
        python <<CODE
        import uuid
        print("gs://truwl-internal-inputs/nf-cutandrun/{}".format(str(uuid.uuid4())))
        CODE
>>>

  output {
    String uuid = read_string(stdout())
  }
  
  runtime {
    docker: "python:3.8.12-buster"
  }
}

task touch_uuid {
    input {
        String outputbucket
    }

    command <<<
        echo "sentinel" > sentinelfile
        gsutil cp sentinelfile ~{outputbucket}/sentinelfile
    >>>

    output {
        String touchedbucket = outputbucket
    }

    runtime {
        docker: "google/cloud-sdk:latest"
    }
}

task fetch_results {
    input {
        String outputbucket
        File execution_trace
    }
    command <<<
        cat ~{execution_trace}
        echo ~{outputbucket}
        mkdir -p ./resultsdir
        gsutil cp -R ~{outputbucket} ./resultsdir
    >>>
    output {
        Array[File] results = glob("resultsdir/*")
    }
    runtime {
        docker: "google/cloud-sdk:latest"
    }
}

task run_nfcoretask {
    input {
        String outputbucket
		File samplesheet
		String outdir = "./results"
		String? email
		String? multiqc_title
		Boolean? save_reference
		Boolean? save_merged_fastq
		Boolean? save_trimmed
		Boolean? save_spikein_aligned
		Boolean? save_unaligned
		Boolean? save_align_intermed
		String aligner = "bowtie2"
		Int? clip_r1
		Int? clip_r2
		Int? three_prime_clip_r1
		Int? three_prime_clip_r2
		Int? trim_nextseq
		String? genome
		String? bowtie2
		String? gtf
		String? blacklist
		String spikein_genome = "K12-MG1655"
		String? spikein_bowtie2
		String? spikein_fasta
		File? fasta
		String igenomes_base = "s3://ngi-igenomes/igenomes"
		Boolean? igenomes_ignore
		Int? minimum_alignment_q_score
		Boolean? dedup_target_reads
		Int normalisation_c = 10000
		Boolean igg_control = true
		Float peak_threshold = 0.05
		String? gene_bed
		Int replicate_threshold = 1
		Boolean? only_genome
		Boolean? only_input
		Boolean? only_preqc
		Boolean? only_alignment
		Boolean? only_filtering
		Boolean? only_peak_calling
		Boolean? skip_fastqc
		Boolean? skip_trimming
		Boolean? skip_removeduplicates
		Boolean? skip_scale
		Boolean? skip_reporting
		Boolean? skip_igv
		Boolean? skip_heatmaps
		Boolean? skip_multiqc
		Boolean? skip_upset_plots
		String custom_config_version = "master"
		String custom_config_base = "https://raw.githubusercontent.com/nf-core/configs/master"
		String? config_profile_name
		String? config_profile_description
		String? config_profile_contact
		String? config_profile_url
		Int max_cpus = 16
		String max_memory = "128.GB"
		String max_time = "240.h"
		Boolean? help
		String publish_dir_mode = "copy"
		String? email_on_fail
		Boolean? plaintext_email
		String max_multiqc_email_size = "25.MB"
		Boolean? monochrome_logs
		String? multiqc_config
		String tracedir = "./results/pipeline_info"
		Boolean validate_params = true
		Boolean? show_hidden_params
		Boolean? enable_conda
		Boolean? singularity_pull_docker_container

	}
	command <<<
		export NXF_VER=21.10.5
		export NXF_MODE=google
		echo ~{outputbucket}
		/nextflow -c /truwl.nf.config run /cutandrun-1.1  -profile truwl,nfcore-cutandrun  --input ~{samplesheet} 	~{"--samplesheet '" + samplesheet + "'"}	~{"--outdir '" + outdir + "'"}	~{"--email '" + email + "'"}	~{"--multiqc_title '" + multiqc_title + "'"}	~{true="--save_reference  " false="" save_reference}	~{true="--save_merged_fastq  " false="" save_merged_fastq}	~{true="--save_trimmed  " false="" save_trimmed}	~{true="--save_spikein_aligned  " false="" save_spikein_aligned}	~{true="--save_unaligned  " false="" save_unaligned}	~{true="--save_align_intermed  " false="" save_align_intermed}	~{"--aligner '" + aligner + "'"}	~{"--clip_r1 " + clip_r1}	~{"--clip_r2 " + clip_r2}	~{"--three_prime_clip_r1 " + three_prime_clip_r1}	~{"--three_prime_clip_r2 " + three_prime_clip_r2}	~{"--trim_nextseq " + trim_nextseq}	~{"--genome '" + genome + "'"}	~{"--bowtie2 '" + bowtie2 + "'"}	~{"--gtf '" + gtf + "'"}	~{"--blacklist '" + blacklist + "'"}	~{"--spikein_genome '" + spikein_genome + "'"}	~{"--spikein_bowtie2 '" + spikein_bowtie2 + "'"}	~{"--spikein_fasta '" + spikein_fasta + "'"}	~{"--fasta '" + fasta + "'"}	~{"--igenomes_base '" + igenomes_base + "'"}	~{true="--igenomes_ignore  " false="" igenomes_ignore}	~{"--minimum_alignment_q_score " + minimum_alignment_q_score}	~{true="--dedup_target_reads  " false="" dedup_target_reads}	~{"--normalisation_c " + normalisation_c}	~{true="--igg_control  " false="" igg_control}	~{"--peak_threshold " + peak_threshold}	~{"--gene_bed '" + gene_bed + "'"}	~{"--replicate_threshold " + replicate_threshold}	~{true="--only_genome  " false="" only_genome}	~{true="--only_input  " false="" only_input}	~{true="--only_preqc  " false="" only_preqc}	~{true="--only_alignment  " false="" only_alignment}	~{true="--only_filtering  " false="" only_filtering}	~{true="--only_peak_calling  " false="" only_peak_calling}	~{true="--skip_fastqc  " false="" skip_fastqc}	~{true="--skip_trimming  " false="" skip_trimming}	~{true="--skip_removeduplicates  " false="" skip_removeduplicates}	~{true="--skip_scale  " false="" skip_scale}	~{true="--skip_reporting  " false="" skip_reporting}	~{true="--skip_igv  " false="" skip_igv}	~{true="--skip_heatmaps  " false="" skip_heatmaps}	~{true="--skip_multiqc  " false="" skip_multiqc}	~{true="--skip_upset_plots  " false="" skip_upset_plots}	~{"--custom_config_version '" + custom_config_version + "'"}	~{"--custom_config_base '" + custom_config_base + "'"}	~{"--config_profile_name '" + config_profile_name + "'"}	~{"--config_profile_description '" + config_profile_description + "'"}	~{"--config_profile_contact '" + config_profile_contact + "'"}	~{"--config_profile_url '" + config_profile_url + "'"}	~{"--max_cpus " + max_cpus}	~{"--max_memory '" + max_memory + "'"}	~{"--max_time '" + max_time + "'"}	~{true="--help  " false="" help}	~{"--publish_dir_mode '" + publish_dir_mode + "'"}	~{"--email_on_fail '" + email_on_fail + "'"}	~{true="--plaintext_email  " false="" plaintext_email}	~{"--max_multiqc_email_size '" + max_multiqc_email_size + "'"}	~{true="--monochrome_logs  " false="" monochrome_logs}	~{"--multiqc_config '" + multiqc_config + "'"}	~{"--tracedir '" + tracedir + "'"}	~{true="--validate_params  " false="" validate_params}	~{true="--show_hidden_params  " false="" show_hidden_params}	~{true="--enable_conda  " false="" enable_conda}	~{true="--singularity_pull_docker_container  " false="" singularity_pull_docker_container}	-w ~{outputbucket}
	>>>
        
    output {
        File execution_trace = "pipeline_execution_trace.txt"
        Array[File] results = glob("results/*/*html")
    }
    runtime {
        docker: "truwl/nfcore-cutandrun:1.1_0.1.0"
        memory: "2 GB"
        cpu: 1
    }
}
    