# TBfastProfiler
Run [TBProfiler](https://github.com/jodyphelan/TBProfiler) and [fastp](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6129281/) on some tuberculosis fastqs in a single WDL task. Make sure to use `--copy-input-files` if running with miniwdl.

If you only want to run fastp, use [my standalone WDL](https://github.com/aofarrel/fastp-wdl) or [Thiagen's implementation](https://github.com/theiagen/public_health_viral_genomics/blob/d75e99bd471413ed9315fb31183dcff934d79204/tasks/task_read_clean.wdl#L250).  
If you only want to run TBProfiler, use [my standalone WDL](https://github.com/aofarrel/tb_profiler).

## Justification
Every WDL task spins up a Docker image, and spinning up a Docker image has some monetary and temporal overhead if you're running on the Cloud. For very large runs, it is economical to combine very short-lived tasks. In particular, some versions of [myco](github.com/aofarrel/myco) use both TBProfiler and fastp as filtering steps.

## Caveat
*This only applies to people who want to reuse the Docker image for something else. If you're using this WDL, this caveat does not apply!*  
This WDL is designed to run TBProfiler in fastq mode, not BAM mode. If you're more interested in the Docker image than the WDL, be warned that TBfastProfiler's Docker image is based on [my dockerization of TBProfiler](https://github.com/aofarrel/tb_profiler/blob/main/Dockerfile), which explictly declares the TB reference genome to be named [NC_000962.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000962.3). There is a chance you'll run into issues if you try to use this Docker to run TBProfiler in BAM mode against BAMs created against something not called `NC_000962.3`, even if it's the same actual sequence. [See TBProfiler's documentation for more details](https://github.com/jodyphelan/TBProfiler/tree/v4.4.2#running-with-an-existing-bam-file).
