msir
====

Tandem repeat analyzer for microsatellite instability detection by DNA-seq

[![wercker status](https://app.wercker.com/status/f7d34d6279d11f821ed9a829b75f13a5/s/master "wercker status")](https://app.wercker.com/project/byKey/f7d34d6279d11f821ed9a829b75f13a5)

Installation
------------

```sh
$ pip install -U https://github.com/dceoy/msir/archive/master.tar.gz
```

Usage
-----

```sh
$ msir --help
Tandem repeat analyzer for microsatellite instability detection by DNA-seq

Usage:
    msir id [--debug] [--max-unit-len=<int>] [--min-rep-times=<int>]
            [--ex-region-len=<int>] [--processes=<int>] [--unit-tsv=<path>]
            <bed> <fasta>
    msir count [--debug] [--unit-tsv=<path>] [--out-dir=<path>] [--index-bam]
               [--samtools=<path>] [--cut-end-len=<int>] [--csv]
               [--processes=<int>] <bam>...
    msir -h|--help
    msir -v|--version

Options:
    -h, --help              Print help and exit
    -v, --version           Print version and exit
    --debug                 Execute a command with debug messages
    --max-unit-len=<int>    Set a maximum length for repeat units [default: 10]
    --min-rep-times=<int>   Set a minimum repeat times [default: 3]
    --ex-region-len=<int>   Search around extra regions [default: 0]
    --processes=<int>       Limit max cores for multiprocessing
    --unit-tsv=<path>       Set a TSV file of repeat units [default: ru.tsv]
    --out-dir=<path>        Pass an output directory [default: .]
    --index-bam             Index BAM or CRAM files if required
    --samtools=<path>       Set a path to samtools command
    --cut-end-len=<int>     Ignore repeats on ends of reads [default: 10]
    --csv                   Write results with CSV instead of TSV

Arguments:
    <bed>                   Path to a BED file of repetitive regions
    <fasta>                 Path to a reference genome fasta file
    <bam>                   Path to an input BAM/CRAM file

Commands:
    id                      Indentify repeat units from reference sequences
    count                   Extract and count tandem repeats in read sequences
```
