# Caller catalog format & choosing a catalog

## Caller catalog format

For the `--loci` argument, `strkit call` takes a list of loci in a modified BED / TSV format,
similar to methods like Straglr/Tandem-genotypes/GangSTR.

The file must be structured with a row per locus, where each row looks like:

```
chr#    10000    10100    [...]    AC
```

The important requirements here are:

  * The fields are tab-separated
  * Locus coordinates are 0-based and half-open (start is inclusive, end is exclusive)
  * The locus motif must come **last** in the row, but *any number of fields* can separate
    the end position and the motif.

As a result, STRkit can take myrid different TSV-type catalog formats as input, including
those produced from the TRF UCSC browser track, or for GangSTR, or for Straglr.

Here are a few important notes on catalogs:

  * TODO: talk about expanding coordinates, TRF multiple for same TR, ...
