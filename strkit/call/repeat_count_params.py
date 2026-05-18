from dataclasses import dataclass
from typing import Literal, TypeAlias

__all__ = ["RepeatCountMethod", "RepeatCountParams"]

RepeatCountMethod: TypeAlias = Literal["repalign", "comp"]


@dataclass(frozen=True)
class RepeatCountParams:
    method: RepeatCountMethod
    max_iters: int
    initial_local_search_range: int  # Initial value; can be narrowed within the get_repeat_count fn
    initial_step_size: int
