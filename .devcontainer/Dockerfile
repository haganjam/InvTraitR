FROM blacksnowflake/r-base-dev:4.2.3@sha256:6bcd5e08cfe7a021dc3c80035f4f177683e0dd3279d04e567b5f747842bd30b6

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
