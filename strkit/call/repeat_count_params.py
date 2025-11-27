from dataclasses import dataclass

__all__ = ["RepeatCountParams"]


@dataclass(frozen=True)
class RepeatCountParams:
    max_iters: int
    initial_local_search_range: int  # Initial value; can be narrowed within the get_repeat_count fn
    initial_step_size: int
