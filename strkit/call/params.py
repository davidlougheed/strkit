import logging
import pathlib

from pysam import AlignmentFile

from .gmm import GMMParams
from .repeats import RepeatCountParams
from ..logger import log_levels

__all__ = ["CallParams"]


DEFAULT_RC_LOCAL_SEARCH_RANGE: int = 3
DEFAULT_RC_STEP_SIZE: int = 1


class CallParams:
    def __init__(
        self,

        logger: logging.Logger,

        read_file: str,
        reference_file: str,
        loci_file: str,
        sample_id: str | None,
        min_reads: int = 4,
        min_allele_reads: int = 2,
        max_reads: int = 250,
        min_avg_phred: int = 13,
        min_read_align_score: float = 0.9,
        max_rcn_iters: int = 50,
        num_bootstrap: int = 100,
        gm_filter_factor: int = 3,
        flank_size: int = 70,
        skip_supplementary: bool = False,
        skip_secondary: bool = False,
        sex_chroms: str | None = None,
        realign: bool = False,
        hq: bool = False,
        use_hp: bool = False,
        snv_vcf: pathlib.Path | None = None,
        snv_min_base_qual: int = 20,
        targeted: bool = False,
        respect_ref: bool = False,
        count_kmers: str = "none",  # "none" | "peak" | "read"
        consensus: bool = False,
        vcf_anchor_size: int = 5,
        # ---
        log_level: int = logging.WARNING,
        verbose: bool = False,
        seed: int | None = None,
        processes: int = 1,
    ):
        self.read_file: str = read_file
        self.reference_file: str = reference_file
        self.loci_file: str = loci_file
        self.min_reads: int = min_reads
        self.min_allele_reads: int = min_allele_reads
        self.max_reads: int = max_reads
        self.min_avg_phred: int = min_avg_phred
        self.min_read_align_score: float = min_read_align_score
        self.max_rcn_iters: int = max_rcn_iters
        self.num_bootstrap: int = num_bootstrap
        self.gm_filter_factor: int = gm_filter_factor
        self.flank_size: int = flank_size
        self.skip_supplementary: bool = skip_supplementary
        self.skip_secondary: bool = skip_secondary
        self.sex_chroms: str | None = sex_chroms
        self.realign: bool = realign
        self.hq: bool = hq
        self.use_hp: bool = use_hp
        self.snv_vcf: pathlib.Path | None = snv_vcf
        self.snv_min_base_qual: int = snv_min_base_qual
        self.targeted: bool = targeted
        self.respect_ref: bool = respect_ref
        self.count_kmers: str = count_kmers
        self.consensus: bool = consensus
        self.vcf_anchor_size: int = vcf_anchor_size
        # ---
        self.log_level: int = log_level
        self.verbose: bool = verbose
        self.seed: int | None = seed
        self.processes: int = processes

        bf = AlignmentFile(read_file, reference_filename=reference_file)

        # noinspection PyTypeChecker
        bfh = bf.header.to_dict()

        sns: set[str] = {e.get("SM") for e in bfh.get("RG", ()) if e.get("SM")}
        bam_sample_id: str | None = None

        if len(sns) > 1:
            # Error or warning or what?
            sns_str = "', '".join(sns)
            logger.warning(f"Found more than one sample ID in BAM file(s): '{sns_str}'")
        elif not sns:
            if not sample_id:
                logger.warning("Could not find sample ID in BAM file(s); sample ID can be set manually via --sample-id")
        else:
            bam_sample_id = sns.pop()

        self._sample_id_orig: str | None = sample_id
        self.sample_id = sample_id or bam_sample_id

        self._rc_params: RepeatCountParams = RepeatCountParams(
            max_iters=max_rcn_iters,
            # TODO: user params for more of this
            initial_local_search_range=DEFAULT_RC_LOCAL_SEARCH_RANGE,
            initial_step_size=DEFAULT_RC_STEP_SIZE,
        )

        # TODO: user params for more of this
        self._gmm_params: GMMParams = GMMParams(
            init_params_method="k-means++", n_init=3, pre_filter_factor=5, filter_factor=self.gm_filter_factor
        )

    @property
    def rc_params(self) -> RepeatCountParams:
        return self._rc_params

    @property
    def gmm_params(self) -> GMMParams:
        return self._gmm_params

    @classmethod
    def from_args(cls, logger: logging.Logger, p_args):
        return cls(
            logger,
            p_args.read_file,
            p_args.ref,
            p_args.loci,
            sample_id=p_args.sample_id,
            min_reads=p_args.min_reads,
            min_allele_reads=p_args.min_allele_reads,
            max_reads=p_args.max_reads,
            min_avg_phred=p_args.min_avg_phred,
            min_read_align_score=p_args.min_read_align_score,
            max_rcn_iters=p_args.max_rcn_iters,
            num_bootstrap=p_args.num_bootstrap,
            gm_filter_factor=p_args.gm_filter_factor,
            flank_size=p_args.flank_size,
            skip_supplementary=p_args.skip_supplementary,
            skip_secondary=p_args.skip_secondary,
            sex_chroms=p_args.sex_chr,
            realign=p_args.realign,
            hq=p_args.hq,
            use_hp=p_args.use_hp,
            snv_vcf=p_args.incorporate_snvs,
            snv_min_base_qual=p_args.snv_min_base_qual,
            targeted=p_args.targeted,
            respect_ref=p_args.respect_ref,
            count_kmers=p_args.count_kmers,
            consensus=p_args.consensus or not (not p_args.vcf),  # Consensus calculation is required for VCF output.
            vcf_anchor_size=min(max(p_args.vcf_anchor_size, 1), p_args.flank_size),
            # ---
            log_level=log_levels[p_args.log_level],
            verbose=p_args.verbose,
            seed=p_args.seed,
            processes=p_args.processes,
        )

    def to_dict(self, as_inputted: bool = False):
        return {
            "read_file": self.read_file,
            "reference_file": self.reference_file,
            "min_reads": self.min_reads,
            "min_allele_reads": self.min_allele_reads,
            "max_reads": self.max_reads,
            "min_avg_phred": self.min_avg_phred,
            "min_read_align_score": self.min_read_align_score,
            "max_rcn_iters": self.max_rcn_iters,
            "num_bootstrap": self.num_bootstrap,
            "gm_filter_factor": self.gm_filter_factor,
            "flank_size": self.flank_size,
            "skip_supplementary": self.skip_supplementary,
            "skip_secondary": self.skip_secondary,
            "sample_id": self._sample_id_orig if as_inputted else self.sample_id,
            "realign": self.realign,
            "hq": self.hq,
            "use_hp": self.use_hp,
            "snv_vcf": str(self.snv_vcf) if self.snv_vcf else None,
            "snv_min_base_qual": self.snv_min_base_qual,
            "targeted": self.targeted,
            "respect_ref": self.respect_ref,
            "count_kmers": self.count_kmers,
            "consensus": self.consensus,
            "vcf_anchor_size": self.vcf_anchor_size,
            "log_level": self.log_level,
            "verbose": self.verbose,
            "seed": self.seed,
            "processes": self.processes,
        }
