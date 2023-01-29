FROM blacksnowflake/r-base-dev:4.2.2@sha256:658759e58f57c21739184439a9c2bfc034429c95a24090fd27a4c2860ed1bbd7

# gdal required for terra (raster)
RUN apt update && \
    apt install -y --no-install-recommends libgdal-dev && \
    :

WORKDIR /invertr

COPY . .
