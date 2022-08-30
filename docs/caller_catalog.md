# Caller catalog format & choosing a catalog

## Caller catalog format

For the `--loci` argument, `strkit call` takes a list of loci in a modified BED / TSV format,
similar to methods like Straglr/Tandem-genotypes/GangSTR.

The file must be structured with a row per locus, where each row looks like:

```
chr#    10000    10101    [...]    AC
```

The important requirements here are:

  * The fields are tab-separated
  * Locus coordinates are 0-based and half-open (start is inclusive, end is exclusive)
  * The locus motif must come **last** in the row, but *any number of fields* can separate
    the end position and the motif.

As a result, STRkit can take myrid different TSV-type catalog formats as input, including
those produced from the TRF UCSC browser track, or for GangSTR, or for Straglr.

Here are a few notes on catalogs:

  * TODO: talk about expanding coordinates, TRF multiple for same TR, ...
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
    for the `hg38` reference genome, based on the review research done by. It will be available soon (TO COME).
