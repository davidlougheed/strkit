import multiprocessing as mp

__all__ = [
    "NonDaemonicPool",
]


# Need a pool which itself can spawn realignment processes - see https://stackoverflow.com/a/53180921


class NonDaemonicProcess(mp.Process):
    @property
    def daemon(self):
        return False

    @daemon.setter
    def daemon(self, value):
        pass


class NonDaemonicContext(type(mp.get_context())):
    Process = NonDaemonicProcess


class NonDaemonicPool(mp.pool.Pool):
    # noinspection PyArgumentList
    def __init__(self, *args, **kwargs):
        kwargs["context"] = NonDaemonicContext()
        super().__init__(*args, **kwargs)
