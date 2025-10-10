# Caller catalog format & choosing a catalog

## Caller catalog format

For the `--loci` argument, `strkit call` takes a list of loci in a modified BED / TSV format,
similar to methods like Straglr/Tandem-genotypes/GangSTR.

The file must be structured with a row per locus, where each row looks like one of the following:

```
chr4    3074876    3074966    [...]    CAG
chrX    147912050  147912110  [...]    ID=FMR1;MOTIF=CGG
```

The important requirements here are:

  * The fields are tab-separated
  * The rows are sorted by contig, and then by starting position
  * Locus coordinates are 0-based and half-open (start is inclusive, end is exclusive)
  * The locus motif (either plain, or as part of a [more complex locus description](#specifying-a-custom-locus-id)) 
    must come **last** in the row, but *any number of fields* can separate the end position and the motif.

As a result, STRkit can take myrid different TSV-type catalog formats as input, including
those produced from the TRF UCSC browser track, or for GangSTR, or for Straglr.

Here are a few notes on catalogs:

  * Coordinates are used to locate the STR locus in the reference genome, but may be slightly 
    expanded to better encompass the entire locus.
  * Be wary of using Tandem Repeats Finder output directly as a catalog, as it can output multiple
    rows for the same locus, or define motifs in a "compound" fashion, e.g., `ATATAT` instead of `AT`.
  * Some disease expansions can contain multiple different motifs, 
    which may be not present in the reference genome at all (for example: 
    [CANVAS](https://pubmed.ncbi.nlm.nih.gov/31230722/), [BAFME2](https://www.nature.com/articles/s41467-019-12671-y)).
    As such, we provide a mechanism to specify motifs using any 
    [IUPAC code](https://www.bioinformatics.org/sms/iupac.html). 
    Thus, the CANVAS and BAFME2 motifs can be represented as `AARRG` and `AAAWK`, respectively.
    We also add in a non-IUPAC code, `X`, which behaves like `N` in that it represents any base, 
    but instead of giving a reward of `+2` it neither penalizes nor rewards alignment, 
    and penalizes a gap. We use this internally to represent low-confidence base calls.
  * Related to the above, this can be important for diseases such as SCA37, where the motif composition 
    (rather than the actual copy number) is associated with disease 
    ([Seixas *et al.* 2017](https://doi.org/10.1016%2Fj.ajhg.2017.06.007)). Here, STRkit's motif-sized k-mer counting
    function can be used during calling with the `--count-kmers` flag. See the 
    [advanced usage](https://github.com/davidlougheed/strkit/blob/master/docs/caller_usage.md#all-optional-flags) page 
    for more.


### Specifying a custom locus ID

The last column of a STRkit-compatible locus catalog can be either just a motif, or a set of key-value pairs. The 
following are all valid values for the last column of the catalog BED file:

* `CAG` results in `ID=locus#;MOTIF=CAG` where `#` is the line index of the locus (1-indexed).
* `MOTIF=CAG` results in the same as above.
* `ID=HTT;MOTIF=CAG` specifies a custom ID alongside the `CAG` motif.

The locus ID is included in the JSON and VCF outputs, as well as logging output.


## Choosing an existing catalog

Other researchers have done extensive work in identifying and cataloguing loci for genotyping:

  * The Tandem Repeats Finder track for the UCSC browser, available as a 
    [downloadable BED file](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.trf.bed.gz),
    with the caveat that this file includes **overlapping entries**, and TRs may not always be represented in 
    their most 'essential' form (e.g., using the motif `TATATATA` instead of just `TA`). Thus, some work may be
    required to create a desirable locus catalog.
  * The researchers behind the [GangSTR](https://github.com/gymreklab/GangSTR) short-read STR genotyping method
    have prepared [several extensive STR catalogs](https://github.com/gymreklab/GangSTR#gangstr-reference-files) 
    for different human reference genomes, containing motifs up to 20bp in length. However, **these files use
    1-based closed-interval coordinates**, and should be adjusted (subtracting 1 from all start coordinates) to 
    transform them into the 0-based half-open interval coordinates when using them with STRkit.
  * We have prepared a [catalog of disease-causing or disease-associated loci](../catalogs/pathogenic_assoc.hg38.tsv) 
    for the `hg38` reference genome, partially based on the review research done by Gall-Duncan *et al.* (2022), as well
    as entries from the [STRipy database](https://stripy.org/database) 
    (DOI: [10.1002/humu.24382](https://doi.org/10.1002/humu.24382)) and our own reading of other articles.
