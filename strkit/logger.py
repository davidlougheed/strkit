import logging
import sys

__all__ = [
    "get_main_logger",
    "attach_stream_handler",
    "create_process_logger",
    "log_levels",
]

fmt = logging.Formatter(fmt="%(name)s:\t[%(levelname)s]\t%(message)s")


def get_main_logger(level: int = logging.DEBUG):
    logger = logging.getLogger("strkit-main")
    logger.setLevel(level)
    return logger


def attach_stream_handler(level: int, logger_=None):
    ch = logging.StreamHandler(sys.stderr)
    ch.setLevel(level)
    ch.setFormatter(fmt)
    logger_.addHandler(ch)


def create_process_logger(pid: int, level: int):
    lg = logging.getLogger(f"strkit-{pid}")
    lg.setLevel(level)
    if not lg.handlers:
        attach_stream_handler(level, logger_=lg)
    return lg


log_levels = {
    "debug": logging.DEBUG,
    "info": logging.INFO,
    "warning": logging.WARNING,
    "error": logging.ERROR,
}
