from dataclasses import dataclass

__all__ = ["RepeatCountParams"]


@dataclass
class RepeatCountParams:
    max_iters: int
    initial_local_search_range: int  # Initial value; can be narrowed within the get_repeat_count fn
    initial_step_size: int

    def __hash__(self):
        return hash((self.max_iters, self.initial_local_search_range, self.initial_step_size))
