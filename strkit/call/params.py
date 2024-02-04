import logging
import pathlib

from pysam import AlignmentFile
from typing import Optional

from ..logger import log_levels

__all__ = ["CallParams"]


class CallParams:
    def __init__(
        self,

        logger: logging.Logger,

        read_files: tuple[str, ...],
        reference_file: str,
        loci_file: str,
        sample_id: Optional[str],
        min_reads: int = 4,
        min_allele_reads: int = 2,
        min_avg_phred: int = 13,
        num_bootstrap: int = 100,
        flank_size: int = 70,
        sex_chroms: Optional[str] = None,
        realign: bool = False,
        hq: bool = False,
        use_hp: bool = False,
        snv_vcf: Optional[pathlib.Path] = None,
        targeted: bool = False,
        fractional: bool = False,
        respect_ref: bool = False,
        count_kmers: str = "none",  # "none" | "peak" | "read"
        consensus: bool = False,
        log_level: int = logging.WARNING,
        seed: Optional[int] = None,
        processes: int = 1,
    ):
        self.read_files: tuple[str, ...] = read_files
        self.reference_file: str = reference_file
        self.loci_file: str = loci_file
        self.min_reads: int = min_reads
        self.min_allele_reads: int = min_allele_reads
        self.min_avg_phred: int = min_avg_phred
        self.num_bootstrap: int = num_bootstrap
        self.flank_size: int = flank_size
        self.sex_chroms: Optional[str] = sex_chroms
        self.realign: bool = realign
        self.hq: bool = hq
        self.use_hp: bool = use_hp
        self.snv_vcf: Optional[pathlib.Path] = snv_vcf
        self.targeted: bool = targeted
        self.fractional: bool = fractional
        self.respect_ref: bool = respect_ref
        self.count_kmers: str = count_kmers
        self.consensus: bool = consensus
        self.log_level: int = log_level
        self.seed: Optional[int] = seed
        self.processes: int = processes

        bfs = tuple(AlignmentFile(rf, reference_filename=reference_file) for rf in read_files)

        # noinspection PyTypeChecker
        bfhs = [bf.header.to_dict() for bf in bfs]

        sns: set[str] = {e.get("SM") for bfh in bfhs for e in bfh.get("RG", ()) if e.get("SM")}
        bam_sample_id: Optional[str] = None

        if len(sns) > 1:
            # Error or warning or what?
            sns_str = "', '".join(sns)
            logger.warning(f"Found more than one sample ID in BAM file(s): '{sns_str}'")
        elif not sns:
            if not sample_id:
                logger.warning("Could not find sample ID in BAM file(s); sample ID can be set manually via --sample-id")
        else:
            bam_sample_id = sns.pop()

        self._sample_id_orig: Optional[str] = sample_id
        self.sample_id = sample_id or bam_sample_id

    @classmethod
    def from_args(cls, logger: logging.Logger, p_args):
        return cls(
            logger,
            tuple(p_args.read_files),
            p_args.ref,
            p_args.loci,
            sample_id=p_args.sample_id,
            min_reads=p_args.min_reads,
            min_allele_reads=p_args.min_allele_reads,
            min_avg_phred=p_args.min_avg_phred,
            num_bootstrap=p_args.num_bootstrap,
            flank_size=p_args.flank_size,
            sex_chroms=p_args.sex_chr,
            realign=p_args.realign,
            hq=p_args.hq,
            use_hp=p_args.use_hp,
            # incorporate_snvs=p_args.incorporate_snvs,
            snv_vcf=p_args.incorporate_snvs,
            targeted=p_args.targeted,
            fractional=p_args.fractional,
            respect_ref=p_args.respect_ref,
            count_kmers=p_args.count_kmers,
            consensus=p_args.consensus or not (not p_args.vcf),  # Consensus calculation is required for VCF output.
            log_level=log_levels[p_args.log_level],
            processes=p_args.processes,
        )

    def to_dict(self, as_inputted: bool = False):
        return {
            "read_files": self.read_files,
            "reference_file": self.reference_file,
            "min_reads": self.min_reads,
            "min_allele_reads": self.min_allele_reads,
            "min_avg_phred": self.min_avg_phred,
            "num_bootstrap": self.num_bootstrap,
            "flank_size": self.flank_size,
            "sample_id": self._sample_id_orig if as_inputted else self.sample_id,
            "realign": self.realign,
            "hq": self.hq,
            "use_hp": self.use_hp,
            "snv_vcf": str(self.snv_vcf) if self.snv_vcf else None,
            "targeted": self.targeted,
            "fractional": self.fractional,
            "respect_ref": self.respect_ref,
            "count_kmers": self.count_kmers,
            "consensus": self.consensus,
            "log_level": self.log_level,
            "processes": self.processes,
        }
