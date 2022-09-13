import logging
import sys

__all__ = [
    "logger",
    "attach_stream_handler",
]

fmt = logging.Formatter(fmt="%(name)s: [%(levelname)s] %(message)s")

logger = logging.getLogger("strkit")
logger.setLevel(logging.DEBUG)


def attach_stream_handler(level: int):
    ch = logging.StreamHandler(sys.stderr)
    ch.setLevel(level)
    ch.setFormatter(fmt)
    logger.addHandler(ch)
