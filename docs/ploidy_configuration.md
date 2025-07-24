# Ploidy configuration

> **NOTE:** Currently, STRkit only supports calling/analyzing haploid/diploid genomic regions in most parts of the code.
> Polyploid genomes may work with some parts of the code, but are not officially supported. 

STRkit uses a "ploidy configuration file" to determine the number of alleles to call, passed via `--ploidy` or 
`--sex-chr`. Several are included in the package (see [`strkit/data/ploidy_configs`](../strkit/data/ploidy_configs)),
and can be referenced by alias (`diploid_autosomes`, `XX`/`diploid_xx`, `XY`/`diploid_xy`, `haploid`) but custom 
configurations can also be specified as a path to a JSON file formatted like the files linked.
