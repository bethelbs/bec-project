# Container image for the BMEC regulatory-element analysis.
#
# Two stages: the build stage installs the Pixi environment from pixi.lock;
# the final image keeps only that environment on plain Ubuntu (no Pixi).
#
# Build (from the project root):
#   docker build --platform linux/amd64 -f containers/Dockerfile \
#     --build-arg ENV_NAME=pipeline -t bec-project:1.0.0 .
# Run the whole analysis on the project folder:
#   docker run --rm --platform linux/amd64 -v "$PWD":/work -w /work bec-project:1.0.0 \
#     bash -c "python scripts/download_data.py && snakemake --cores all && quarto render reports/bec-report.qmd"

ARG ENV_NAME=pipeline

# ---- Build stage: install the locked environment -------------------------
FROM ghcr.io/prefix-dev/pixi:0.83.0 AS build
ARG ENV_NAME
WORKDIR /app
COPY pixi.toml pixi.lock ./
# --locked: fail if pixi.lock does not match pixi.toml
RUN pixi install --locked -e ${ENV_NAME}
# Write the environment's activation commands into an entrypoint script
RUN printf '#!/bin/bash\n' > /entrypoint.sh \
 && pixi shell-hook --locked -e ${ENV_NAME} -s bash >> /entrypoint.sh \
 && printf 'exec "$@"\n' >> /entrypoint.sh

# ---- Final image: only the environment ------------------------------------
FROM ubuntu:24.04
ARG ENV_NAME
# The environment keeps its original path, so its activation script still works
COPY --from=build /app/.pixi/envs/${ENV_NAME} /app/.pixi/envs/${ENV_NAME}
COPY --from=build --chmod=0755 /entrypoint.sh /entrypoint.sh
WORKDIR /work
ENTRYPOINT ["/entrypoint.sh"]
CMD ["bash"]
