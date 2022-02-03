__all__ = [
    "MIResult",
]

from typing import List, Optional, Union


class MIResult:
    def __init__(self, mi_value: float, mi_value_ci: Optional[float], non_matching: List[tuple]):
        self.mi_value = mi_value
        self.mi_value_ci = mi_value_ci
        self.non_matching = non_matching

    def as_csv_row(self, sep=","):
        return f"{self.mi_value}{sep}{self.mi_value_ci}\n"

    @staticmethod
    def _res_str(res: Union[float, int]) -> str:
        return f"{res:.1f}" if isinstance(res, float) else str(res)

    def non_matching_tsv(self, sep="\t") -> str:
        res = ""

        for nm in self.non_matching:
            res += sep.join((
                *nm[:3],  # Contig/Start/End
                nm[3],  # Motif

                # Child genotype + CIs if available
                "|".join(map(self._res_str, nm[4])),
                "|".join(",".join(map(self._res_str, ci)) for ci in nm[5]) if nm[5] else "",

                # Mother genotype + CIs if available
                "|".join(map(self._res_str, nm[6])),  # Mother genotype
                "|".join(",".join(map(self._res_str, ci)) for ci in nm[7]) if nm[7] else "",

                # Father genotype + CIs if available
                "|".join(map(self._res_str, nm[8])),
                "|".join(",".join(map(self._res_str, ci)) for ci in nm[9]) if nm[9] else "",

                nm[10],  # Reference number of copies (or blank depending on caller)
            )) + "\n"
        return res

    def __str__(self):
        if self.mi_value_ci:
            return f"MI %\tMI % (CI)\n{self.mi_value*100:.2f}\t{self.mi_value_ci*100:.2f}"
        else:
            return f"MI %\n{self.mi_value*100:.2f}"
