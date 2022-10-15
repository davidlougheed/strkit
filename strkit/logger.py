import logging
import sys

__all__ = [
    "logger",
    "attach_stream_handler",
    "create_process_logger",
]

fmt = logging.Formatter(fmt="%(name)s:\t[%(levelname)s]\t%(message)s")

logger = logging.getLogger("strkit-main")
logger.setLevel(logging.DEBUG)


def attach_stream_handler(level: int, logger_=logger):
    ch = logging.StreamHandler(sys.stderr)
    ch.setLevel(level)
    ch.setFormatter(fmt)
    logger_.addHandler(ch)


def create_process_logger(pid: int, level: int):
    lg = logging.getLogger(f"strkit-{pid}")
    lg.setLevel(logging.DEBUG)
    attach_stream_handler(level, logger_=lg)
    return lg
