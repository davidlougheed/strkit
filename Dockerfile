FROM python:3.13-trixie

WORKDIR /strkit

COPY LICENSE .
COPY MANIFEST.in .
COPY pyproject.toml .
COPY README.md .
COPY uv.lock .
COPY strkit strkit

RUN curl https://sh.rustup.rs -sSf > rustup-init.sh
RUN sh ./rustup-init.sh -y
ENV PATH="/root/.cargo/bin:${PATH}"

COPY --from=ghcr.io/astral-sh/uv:0.12.1 /uv /bin/
ENV UV_COMPILE_BYTECODE=1
ENV UV_NO_DEV=1

RUN uv sync --locked

ENV UV_NO_SYNC=1

CMD [ "uv", "run", "strkit" ]
