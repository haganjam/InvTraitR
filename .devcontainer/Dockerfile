FROM blacksnowflake/r-base-dev:4.3.1@sha256:78701eb5d2f35059d385785d118642beab784811bedc2e3bd09f4e9c9ef40cc8

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
