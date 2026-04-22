# Requires login to dhi.io container registry with free user creds
# docker login dhi.io --username <your-username>
FROM dhi.io/debian-base:trixie-debian13-dev

RUN apt-get update && \
    apt-get install -y cdo=2.5.1-1 \
        unzip=6.0-29 \
        nco=5.3.3-1 && \
    rm -rf /var/lib/apt/lists/*
ENV USER=nonroot \
        GROUP=nonroot
# install uv
COPY --chown=$USER:$GROUP --from=ghcr.io/astral-sh/uv:0.11.7 /uv /uvx /bin/

WORKDIR /home/${USER}
USER $USER

COPY --chown=$USER:$GROUP uv.lock uv.lock
COPY --chown=$USER:$GROUP pyproject.toml pyproject.toml

RUN touch README.md && uv sync --no-dev 

COPY --chown=$USER:$GROUP mkforcing/*.py mkforcing/
COPY --chown=$USER:$GROUP mkforcing/*.sh mkforcing/

COPY --chown=$USER:$GROUP run_get_forcing.sh run_get_forcing.sh

ENTRYPOINT [ "./run_get_forcing.sh" ]
