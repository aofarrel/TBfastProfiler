version 1.0

workflow TBfastProfiler {
    input {
        File fastq1
        File fastq2
        Int average_qual = 30
        Boolean disable_adapter_trimming = true
        Boolean output_fastps_cleaned_fastqs = false
        Float q30_cutoff
    }
    
    parameter_meta {
        average_qual: "If one read's average quality score < avg_qual, then this read/pair is discarded. 0 means no requirement. Independent of q30_cutoff."
        disable_adapter_trimming: "Disable trimming adapters; use this if your reads already went through trimmomatic."
        output_fastps_cleaned_fastqs: "[WDL only] If true, output fastps' cleaned fastqs, otherwise ignore them to save on storage and delocalization costs."
        q30_cutoff: "If a sample's average quality score < q30_cutoff, then this sample is considered failing (did_this_sample_pass output will be false)."
    }
    
    call main {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            average_qual = average_qual,
            disable_adapter_trimming = disable_adapter_trimming,
            output_fastps_cleaned_fastqs = output_fastps_cleaned_fastqs,
    }
    
    # silly workaround to avoid errors with WDL getting made about declaring the same variable multiple times
    # (when what we actually want to do is change it if something is true)
    Boolean always_false = false
    if(main.percent_above_q30 > q30_cutoff) {
        Boolean passed_q30 = true
    }
    Boolean outcome = select_first([passed_q30, always_false])
    
    output {
        File? cleaned_fastq1 = main.very_clean_fastq1
        File? cleaned_fastq2 = main.very_clean_fastq2
        Boolean did_this_sample_pass = outcome
        File fastp_txt = main.fastp_txt
        String samp_resistance = main.samp_resistance
        String samp_strain = main.samp_strain
        File tbprofiler_json = main.tbprofiler_json
        File tbprofiler_txt = main.tbprofiler_txt
    }
}

task main {
    input {
        File fastq1
        File fastq2
        
        # fastp options
        Int average_qual
        Boolean disable_adapter_trimming
        Boolean output_fastps_cleaned_fastqs
        
        # compute setup
        Int addldisk = 15
        Int cpu = 2
        Int memory = 4
        Int preempt = 1
        Boolean ssd = false
    }
    
    parameter_meta {
        average_qual: "If one read's average quality score < avg_qual, then this read/pair is discarded. 0 means no requirement"
        disable_adapter_trimming: "Disable trimming adapters; use this if your reads already went through trimmomatic"
        output_fastps_cleaned_fastqs: "[WDL only] If true, output fastps' cleaned fastqs, otherwise ignore them. fastp will generate cleaned fastqs no matter what, so setting this to false will only save you on storage and delocalization costs."
    }
    
    Int diskSize = addldisk + ceil(2*size(fastq1, "GB"))
    String diskType = if((ssd)) then " SSD" else " HDD"
    
    # fastp arguments
    String arg_adapter_trimming = if(disable_adapter_trimming) then "--disable_adapter_trimming" else ""
    
    # This needs to be to handle inputs like sample+run+num (ERS457530_ERR551697_1.fastq)
    # or inputs like sample+num (ERS457530_1.fastq). In both cases, we want to convert to just
	# sample name (ERS457530).
	String read_file_basename = basename(fastq1)
	String sample_name = sub(read_file_basename, "_.*", "")
    
    command <<<    
    set -eux pipefail
    
    # fastp
    start=$SECONDS
    fastp --in1 "~{fastq1}" --in2 "~{fastq2}" --out1 "~{sample_name}_fastp_1.fq" --out2 "~{sample_name}_fastp_2.fq" \
        --average_qual ~{average_qual} "~{arg_adapter_trimming}" \
        --html "~{sample_name}_fastp.html" --json "~{sample_name}_fastp.json"
    
    # parse fastp outputs from JSON
    python3 << CODE
    import os
    import json
    with open("~{sample_name}_fastp.json", "r") as fastpJSON:
        fastp = json.load(fastpJSON)
    with open("fastp_summary.txt", "w") as outfile:
        for keys, values in fastp["summary"]["before_filtering"].items():
            outfile.write(f"{keys}\t{values}\n")
    with open("q30.txt", "w") as q30_rate: q30_rate.write(str(fastp["summary"]["before_filtering"]["q30_rate"]))
    with open("total_reads.txt", "w") as read_count: read_count.write(str(fastp["summary"]["before_filtering"]["total_reads"]))              
    
    # delete fastp cleaned fastqs if we dont want them to save on delocalization time
    if "~{output_fastps_cleaned_fastqs}" == "false":
        os.remove("~{sample_name}_fastp_1.fq")
        os.remove("~{sample_name}_fastp_2.fq")
    
    CODE
    echo "Finished fastp and parsing its outputs in $(( SECONDS - start ))"

    # tb profiler
    start=$SECONDS
    tb-profiler profile --read1 ~{fastq1} --read2 ~{fastq2} --prefix ~{sample_name} --txt

    # parse tbprofiler output (from the textfile outut, since JSON parsing pains my soul)
    sed -n '11p' results/~{sample_name}.results.txt | sed -r 's/^Strain: //' >> strain.txt
    sed -n '12p' results/~{sample_name}.results.txt | sed -r 's/^Drug-resistance: //' >> resistance.txt
    sed -n '13p' results/~{sample_name}.results.txt | sed -r 's/^Median Depth: //' >> depth.txt
    echo "Finished running TBProfiler and parsing its outputs $(( SECONDS - start ))"
    >>>
    
    runtime {
        cpu: cpu
		disks: "local-disk " + diskSize + diskType
		docker: "ashedpotatoes/tbfastprofiler:0.0.1"
		memory: "${memory} GB"
		preemptible: "${preempt}"
    }
    output {
        # fastqs
        File? very_clean_fastq1 = "~{sample_name}_fastp_1.fq"
        File? very_clean_fastq2 = "~{sample_name}_fastp_2.fq"
        
        # reports as files
        File fastp_html = glob("*_fastp.html")[0]
        File fastp_json = glob("*_fastp.json")[0]
        File fastp_txt  = "fastp_summary.txt"
        File tbprofiler_json = "results/~{sample_name}.results.json"
        File tbprofiler_txt = "results/~{sample_name}.results.txt"
        
        # important metrics
        Int    median_depth = read_int("depth.txt")
        Float  percent_above_q30 = read_float("q30.txt")
        String resistance = read_string("resistance.txt")
        String strain = read_string("strain.txt")
        Int    total_reads = read_int("total_reads.txt")
        
        # important metrics as string with sample name (useful if concatenating results from many samples)
        # terra doesn't support outputs passed on other outputs, so we have to read the files twice
        String samp_median_depth = "~{sample_name}\t" + read_string("depth.txt")
        String samp_percent_above_q30 = "~{sample_name}\t" + read_string("q30.txt")
        String samp_resistance = "~{sample_name}\t" + read_string("resistance.txt")
        String samp_strain = "~{sample_name}\t" + read_string("strain.txt")
        String samp_total_reads = "~{sample_name}\t" + read_string("total_reads.txt")
    }
}