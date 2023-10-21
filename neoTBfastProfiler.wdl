version 1.0
import "https://raw.githubusercontent.com/aofarrel/fastp-wdl/0.0.1/fastp_tasks.wdl" as fashtp
import "https://raw.githubusercontent.com/aofarrel/public_health_bioinformatics/smw-tbprofiler-dev/tasks/species_typing/task_tbprofiler.wdl" as tbprof
import "https://raw.githubusercontent.com/aofarrel/public_health_bioinformatics/smw-tbprofiler-dev/tasks/species_typing/task_tbp_parser.wdl" as tbprof_parser

workflow TBfastProfiler {
    input {
        File fastq1
        File fastq2
        String? operator
        
        # trimming/cleaning reads
        Int average_qual = 30
        Boolean disable_adapter_trimming = true
        
        # qc cutoffs
        Float q30_cutoff = 30
        Float pct_mapped_cutoff = 0.98
        
        # what should we do the bad stuff?
        Boolean soft_all_qc = false
        Boolean soft_q30 = true
        Boolean soft_pct_mapped = false
        Boolean use_fastps_cleaned_fastqs = true
        Int warn_if_below_this_depth = 10
        
    }
    
    parameter_meta {
        fastq1: "This sample's forward read"
        fastq2: "This sample's reverse read"
        average_qual: "If one read's average quality score < avg_qual, then this read/pair (NOT the whole sample) is discarded. 0 means no requirement. Independent of q30_cutoff."
        disable_adapter_trimming: "Disable trimming adapters; use this if your fastqs already went through trimmomatic."
        use_fastps_cleaned_fastqs: "If true, use fastps' cleaned fastqs for TBProfiler and output those cleaned fastqs as task-level outputs. If false, cleaned fastqs will be thrown out and TBProfiler will run on the fastqs you input."
        q30_cutoff: "If a sample's average quality score < q30_cutoff, then this sample is considered a failure. Independent of average_qual."
        soft_all_qc: "If true, pass_or_errorcode will always return PASS. Effectively sets soft_q30 and soft_pct_mapped to true."
        soft_q30: "If true, pass_or_errorcode will return PASS even if this sample failed the Q30 check. pass_or_warnings will still be set."
        soft_pct_mapped: "If true, pass_or_errorcode will return PASS even if this sample failed the percent mapped check. pass_or_warnings will still be set."
        warn_if_below_this_depth: "Mutations below this depth will be flagged as low-depth in the Laboratorian report. Does not affect TBProfiler JSON nor any cleaning of FASTQs."
    }
    
    call fashtp.fastp_and_parse as fastp {
        input:
            fastq_1 = fastq1,
            fastq_2 = fastq2,
            average_qual = average_qual,
            disable_adaptor_trimming = disable_adapter_trimming,
            output_cleaned_fastqs = use_fastps_cleaned_fastqs
    }
    
    call tbprof.tbprofiler as profiler {
        input:
            read1 = select_first([fastp.very_clean_fastq1, fastq1]),
            read2 = select_first([fastp.very_clean_fastq2, fastq2]),
            samplename = fastp.sample_name
    }
    
    call tbprof_parser.tbp_parser as csv_maker {
        input:
            tbprofiler_bam = profiler.tbprofiler_output_bam,
            tbprofiler_bai = profiler.tbprofiler_output_bai,
            tbprofiler_json = profiler.tbprofiler_output_json,
            samplename = fastp.sample_name,
            sequencing_method = "WGS",
            operator = select_first([operator, "operator_not_filled_in"]),
            min_depth = warn_if_below_this_depth
    }
    
    if(soft_all_qc) {
        String override = "PASS"
    }
    
    if(!(fastp.percent_above_q30 > q30_cutoff)) {
        String warning_q30 = "EARLYQC_" + (fastp.percent_above_q30*100) + "_PCT_ABOVE_Q30_(MIN_" + q30_cutoff + ")" #!StringCoercion
        if(!(soft_q30)) {
            String failed_q30 = "EARLYQC_" + (fastp.percent_above_q30*100) + "_PCT_ABOVE_Q30_(MIN_" + q30_cutoff + ")" #!StringCoercion
        }
    }
    if(!(profiler.tbprofiler_pct_reads_mapped > (100 - pct_mapped_cutoff))) {
        String warning_mapping = "EARLYQC_" + profiler.tbprofiler_pct_reads_mapped + "_PCT_MAPPED_(MIN_" + (100-pct_mapped_cutoff) + ")" #!StringCoercion
        if(!(soft_pct_mapped)) {
            String failed_mapping = "EARLYQC_" + profiler.tbprofiler_pct_reads_mapped + "_PCT_MAPPED_(MIN_" + (100-pct_mapped_cutoff) + ")" #!StringCoercion
        }
    }
    String error_or_pass = select_first([override, failed_q30, failed_mapping, "PASS"])
    Array[String] warnings = select_all([warning_q30, warning_mapping])
    if(length(warnings) < 0) {
        Array[String] no_warnings = ["PASS"]
    }
    
    output {
        File? cleaned_fastq1 = fastp.very_clean_fastq1
        File? cleaned_fastq2 = fastp.very_clean_fastq2
        
        # stats
        String status_code = error_or_pass
        Array[String] warning_codes = select_first([no_warnings, warnings])
        String resistance = profiler.tbprofiler_dr_type
        String strain = profiler.tbprofiler_sub_lineage
        Float reads_mapped = profiler.tbprofiler_pct_reads_mapped
        Int median_coverage = profiler.tbprofiler_median_coverage
        Float percent_coverage = csv_maker.tbp_parser_genome_percent_coverage
        
        # reports
        File fastp_txt = fastp.short_report
        File tbprofiler_json = profiler.tbprofiler_output_json 
        #File tbprofiler_txt = profiler.tbprofiler_txt # not in thiagen TBProfiler task
        
        # CSVs for other tools
        File tbprofiler_looker_csv = csv_maker.tbp_parser_looker_report_csv
        File tbprofiler_laboratorian_report_csv = csv_maker.tbp_parser_laboratorian_report_csv
        File tbprofiler_lims_report_csv = csv_maker.tbp_parser_lims_report_csv
        File tbprofiler_coverage_report_csv = csv_maker.tbp_parser_coverage_report
        
    }
}

