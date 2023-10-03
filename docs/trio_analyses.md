# Trio analyses with STRkit

Trio datasets include genomic sequence data for a child, their mother, and their father (the "trio"). These data
can be used to discover de novo mutation (and incidental genotyping errors).

STRkit includes a Mendelian inheritance (MI) analysis tool, under the sub-command `strkit mi`.
After genotyping the trio with `strkit call`, this command can be used to discover loci which:

1. Do not respect exact MI
2. Do not respect MI allowing for a ±1 repeat unit difference 
   (Note: most true mutation occurs in 1-repeat-unit changes too! 
   See [Ellegren, 2004](https://www.nature.com/articles/nrg1348).)
3. Do not respect MI under the 95% locus confidence intervals
4. Look like de novo mutation at a read count distribution level, via a Mann-Whitney *U* test (with tie correction).
   The alternative hypothesis can be specified as either two-sided or looking for expansion in the offspring. 
   *The requirements for this test are invalidated in cases of mosaicism.*
5. Look like de novo mutation at a read count distribution level, via a chi-squared independence test,
   where the contingency table looks like the following:

| Read distribution \ Copy number | 11   | 12   | 13   |
|---------------------------------|------|------|------|
| Parent reads (best peak fit)    | 20   | 10   | 0    |
| Child reads                     | 2    | 20   | 10   |


## Trio-level

At a trio level, the chi-squared test gives (optionally multiple testing-corrected) loci with a significant
chance of containing a de novo mutation.

## Cohort-level

At a cohort level, multiple downstream analyses are possible from a collection of trio mutation analyses,
such as:

  1. Case-control analysis looking for frequency of de novo mutations in specific loci
  2. Case-control analysis looking at the incidence rate of de novo mutation

Currently, tools to automatically perform these analyses are not available in STRkit.
