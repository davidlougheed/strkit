FROM python:3.12-bullseye

WORKDIR /strkit

COPY LICENSE .
COPY MANIFEST.in .
COPY pyproject.toml .
COPY README.md .
COPY setup.py .
COPY strkit strkit

RUN pip install --no-cache-dir .

CMD [ "strkit" ]
