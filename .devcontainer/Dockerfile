FROM blacksnowflake/r-base-dev:4.2.2@sha256:658759e58f57c21739184439a9c2bfc034429c95a24090fd27a4c2860ed1bbd7

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
