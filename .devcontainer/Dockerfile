FROM blacksnowflake/r-base-dev:4.3.0@sha256:f15172f64fe5b3a61c6cd2641d1a2ef70f5adf8f85060290b142839a31550258

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
