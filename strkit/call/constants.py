import numpy as np

__all__ = [
    "NP_EMPTY_ARRAY_INT32",
    "NP_EMPTY_ARRAY_FLOAT64",
    "NP_FLOAT_32_EPSILON",
]

NP_EMPTY_ARRAY_INT32 = np.array([], dtype=np.int32)
NP_EMPTY_ARRAY_FLOAT64 = np.array([], dtype=np.float64)
NP_FLOAT_32_EPSILON = np.finfo(np.float32).eps
