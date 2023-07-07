version 1.0

workflow TBfastProfiler {
    input {
        File fastq1
        File fastq2
    }
    
    call main {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2
    }
}

task main {
    input {
        File fastq1
        File fastq2
        
        # fastp options
        Int average_qual = 30
        #Boolean disable_adaptor_trimming = true
        #Boolean output_fastps_cleaned_fastqs = false
        
        # compute setup
        Int addldisk = 15
        Int cpu = 2
        Int memory = 4
        Int preempt = 1
        Boolean ssd = false
    }
    
    parameter_meta {
        average_qual: "If one read's average quality score < avg_qual, then this read/pair is discarded. 0 means no requirement"
        #disable_adaptor_trimming: "Disable trimming adaptors; use this if your reads already went through trimmomatic"
        #output_fastps_cleaned_fastqs: "[WDL only] If true, output fastps' cleaned fastqs, otherwise ignore them. fastp will generate cleaned fastqs no matter what, so setting this to false will only save you on storage and delocalization costs."
    }
    
    Int diskSize = addldisk + ceil(2*size(fastq1, "GB"))
    String diskType = if((ssd)) then " SSD" else " HDD"
    
    # fastp arguments
    #String arg_adaptor_trimming = "--disable_adaptor_trimming" if disable_adaptor_trimming is true else ""
    
    # This needs to be to handle inputs like sample+run+num (ERS457530_ERR551697_1.fastq)
    # or inputs like sample+num (ERS457530_1.fastq). In both cases, we want to convert to just
	# sample name (ERS457530).
	String read_file_basename = basename(fastq1) # used to calculate sample name + outfile_sam
	String sample_name = sub(read_file_basename, "_.*", "")
    
    command <<<    
    # fastp
    start=$SECONDS
    fastp --in1 "~{fastq1}" --in2 "~{fastq2}" \
        --average_qual ~{average_qual} \
        --html "~{sample_name}_fastp.html" --json "~{sample_name}_fastp.json"
    
    # parse fastp outputs from JSON
    # It is very tempting to redirect fastp's stderr to a file and then parse that, since it
    # dumps much of the information we need to sterr. But if we somehow get an error in fastp,
    # that approach could be problematic.
    python3 << CODE
    import json
    with open("~{sample_name}_fastp.json", "r") as fastpJSON:
        fastp = json.load(fastpJSON)
    with open("fastp_summary.txt", "a") as outfile:
        for keys, values in fastp["summary"]["before_filtering"].items():
            outfile.write(f"{keys}\t{values}\n")
    with open("q30.txt", "w") as q30_rate: q30_rate.write(str(fastp["summary"]["before_filtering"]["q30_rate"]))
    with open("total_reads.txt", "w") as read_count: read_count.write(str(fastp["summary"]["before_filtering"]["total_reads"]))              
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
        # as file
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
        # terra doesn't support outputs pased on other outputs, so we have to read the files twice
        String samp_median_depth = "~{sample_name}\t" + read_string("depth.txt")
        String samp_percent_above_q30 = "~{sample_name}\t" + read_string("q30.txt")
        String samp_resistance = "~{sample_name}\t" + read_string("resistance.txt")
        String samp_strain = "~{sample_name}\t" + read_string("strain.txt")
        String samp_total_reads = "~{sample_name}\t" + read_string("total_reads.txt")
    }
}