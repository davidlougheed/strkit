from dataclasses import dataclass
from typing import Literal, TypeAlias

__all__ = ["RepeatCountMethod", "RepeatCountParams", "get_reference_rc_params"]

RepeatCountMethod: TypeAlias = Literal["repalign", "comp"]


@dataclass(frozen=True)
class RepeatCountParams:
    method: RepeatCountMethod
    max_iters: int
    initial_local_search_range: int  # Initial value; can be narrowed within the get_repeat_count fn
    initial_step_size: int


def get_reference_rc_params(
    method: RepeatCountMethod, ref_est_cn: int, default_ref_max_iters: int
) -> RepeatCountParams:
    ref_max_iters = default_ref_max_iters
    ref_step_size = 1
    ref_local_search_range = 3

    # search less with large repeat counts, but in bigger steps, because each alignment takes a long time.
    if ref_est_cn >= 200:
        if ref_est_cn < 1000:
            ref_step_size = 3
            ref_max_iters = 200
        elif 1000 <= ref_est_cn < 2000:
            ref_step_size = 5
            ref_max_iters = 150
        elif ref_est_cn >= 2000:  # ref_cn >= 2000
            ref_step_size = 15
            ref_max_iters = 50
            ref_local_search_range = 1

    return RepeatCountParams(
        method=method,
        max_iters=ref_max_iters,
        initial_step_size=ref_step_size,
        initial_local_search_range=ref_local_search_range,
    )
