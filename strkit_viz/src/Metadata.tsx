import {ReactNode} from "react";

const EXAMPLE_DATA = {
  "caller": {
    "name": "strkit",
    "version": "0.23.0-dev"
  },
  "contigs": [
    "chr4"
  ],
  "parameters": {
    "consensus": true,
    "count_kmers": "none",
    "flank_size": 70,
    "gm_filter_factor": 3,
    "hq": true,
    "log_level": 20,
    "log_progress_interval": 5,
    "max_rcn_iters": 80,
    "max_reads": 250,
    "min_allele_reads": 2,
    "min_avg_phred": 13,
    "min_read_align_score": 0.9,
    "min_reads": 4,
    "num_bootstrap": 100,
    "ploidy": "diploid_autosomes",
    "processes": 1,
    "read_file": "HG002.SequelII.ccs.phased.40x.chr4.bam",
    "realign": false,
    "reference_file": "/Users/davidlougheed/git/gt-poc/hg38.analysisSet.fa.gz",
    "respect_ref": false,
    // @ts-ignore
    "sample_id": null,
    // @ts-ignore
    "seed": null,
    "skip_secondary": false,
    "skip_supplementary": false,
    "snv_min_base_qual": 20,
    "snv_vcf": "00-common_all.vcf.gz",
    "targeted": false,
    "use_hp": false,
    "vcf_anchor_size": 5,
    "verbose": false
  },
  "runtime": 90.56821920804214,
  "sample_id": "HG002"
};

type MetadataRec = Record<string, null | boolean | string | number>;

const MetadataObject = ({ heading, rec }: { heading: string; rec: MetadataRec }) => {
  return <>
    <tr><th colSpan={2}>{heading}</th></tr>
    {Object.entries(rec).map(([k, v]) => (
      <tr key={k}>
        <th style={{ textAlign: "left" }}>{k}</th>
        <td>{JSON.stringify(v)}</td>
      </tr>
    ))}
  </>;
};

const Metadata = () => {
  let objs: ReactNode[] = [];
  let other: Record<string, null | boolean | string | number> = {};

  Object.entries(EXAMPLE_DATA).forEach(([k, v]) => {
    if (v instanceof Object) {
      objs.push(<MetadataObject heading={k} rec={v as MetadataRec} />)
    } else {
      other[k] = v;
    }
  })

  return <table>
    {objs}
    <MetadataObject heading="Other" rec={other} />
  </table>;
};

export default Metadata;
