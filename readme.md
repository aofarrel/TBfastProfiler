# TBfastProfiler (aka earlyQC)
Use TBProfiler and fastp for your WDL analysis, and output Looker and other CSVs for further analysis.

## TBfastProfiler.wdl
Run [TBProfiler](https://github.com/jodyphelan/TBProfiler) and [fastp](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6129281/) on some tuberculosis fastqs in a single WDL task. Tested on both Cromwell and miniwdl. Make sure to use `--copy-input-files` if running with miniwdl.

**Justification** 
Every WDL task spins up a Docker image, and spinning up a Docker image has some monetary and temporal overhead if you're running on the Cloud. For very large runs, it is economical to combine very short-lived tasks. In particular, some versions of [myco](github.com/aofarrel/myco) use both TBProfiler and fastp as filtering steps. Combining them into one task just makes sense.

## neoTBfastProfiler.wdl
Use Thiagen's fork of TBProfiler and fastp together in the same WDL workflow. 

**Justification** 
Although this lacks some of the Docker spinup savings of the other workflow in this repo, it leverages new changes to TBProfiler to improve usages for local health jurisidictions.

### Alternatives
* If you only want to run fastp, use [my standalone WDL](https://github.com/aofarrel/fastp-wdl) or [Thiagen's implementation](https://github.com/theiagen/public_health_viral_genomics/blob/d75e99bd471413ed9315fb31183dcff934d79204/tasks/task_read_clean.wdl#L250)
* If you only want to run TBProfiler, use [my standalone WDL](https://github.com/aofarrel/tb_profiler), which can work on fastqs or bams
* If you want to run this in a genomics pipeline, consider using [my tuberculosis pipeline](https://github.com/aofarrel/myco)
* If you only care about having a Docker image with fastp and TBProfiler in it, be warned that TBfastProfiler's Docker image is based on [my dockerization of TBProfiler](https://github.com/aofarrel/tb_profiler/blob/main/Dockerfile), which explictly declares the TB reference genome to be named [NC_000962.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000962.3). There is a chance you'll run into issues if you try to use this Docker to run TBProfiler in BAM mode against BAMs created against something not called `NC_000962.3`, even if it's the same actual sequence. [See TBProfiler's documentation for more details](https://github.com/jodyphelan/TBProfiler/tree/v4.4.2#running-with-an-existing-bam-file).
 
