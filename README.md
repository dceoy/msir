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

1.  Identify tandem repeat units on a BED files from a reference FASTA file and write them into a TSV file.

    ```sh
    $ msir id --unit-tsv=./repeat_units.tsv ./microsatellite.bed ./hg38.fa
    ```

2.  Extract and count tandem repeats within read sequences in BAM files and write them into a TSV file.

    ```sh
    $ msir hist --hist-tsv=./repeat_counts.tsv --unit-tsv=./repeat_units.tsv sample1.bam sample2.bam
    ```

Run `msir --help` for more information about options.
