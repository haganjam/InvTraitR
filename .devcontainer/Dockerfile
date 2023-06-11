FROM blacksnowflake/r-base-dev:4.2.3@sha256:331aff7bcfca79285a9a05197c9ed1cc1aee72d96a8a533e94ab4cc3a67ef9e4

# gdal required for terra (raster)
# libudunits2-dev for bdc (units)
RUN apt update \
    && apt install -y --no-install-recommends \
        libgdal-dev \
        libudunits2-dev \
        cmake \
    && :

# install dev dependencies
RUN R -q -e 'install.packages("renv")'

WORKDIR /invertr

COPY . .
