FROM python:3.12-bullseye

WORKDIR /strkit

COPY LICENSE .
COPY MANIFEST.in .
COPY pyproject.toml .
COPY README.md .
COPY setup.py .
COPY strkit strkit

RUN curl https://sh.rustup.rs -sSf > rustup-init.sh
RUN sh ./rustup-init.sh -y
ENV PATH="/root/.cargo/bin:${PATH}"

RUN pip install -U pip
RUN pip install --no-cache-dir .

CMD [ "strkit" ]
