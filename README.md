# FW_invert_biomass_allometry

[![R-CMD-check](https://github.com/haganjam/FW_invert_biomass_allometry/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/haganjam/FW_invert_biomass_allometry/actions/workflows/R-CMD-check.yaml)

Pipeline to assign body length-dry biomass allometry equations and other functional traits to a taxonomic name based on taxonomic hierarchy in freshwater invertebrates. However, the name-matching pipeline can technically be used for any taxonomic group.

## Installation

The package is WIP. It will, hopefully, be available on CRAN soon.  
You can try to install it from GitHub but as it's a WIP, it may or may not
work depending on the tides, temperature, the color of your socks and what
you had for dinner the day before yesterday.

So, if you feel lucky, try:

```r
install_github("haganjam/FW_invert_biomass_allometry")
```

## Usage

The package exports a single function `get_trait_from_taxon`:

```r
get_trait_from_taxon(
    data,                   # data.frame with at least five columns: target taxon, life stage, latitude (dd), longitude (dd) and body size (mm) if trait == "equation"
    target_taxon,           # character string with the column name containing the taxon names
    life_stage,             # character string with the column name containing the life stages
    latitude_dd,            # character string with the column name containing the latitude in decimal degrees
    longitude_dd,           # character string with the column name containing the longitude in decimal degrees
    body_size,              # character string with the column name containing the body size data if trait = "equation"
    workflow = "workflow2", # options are "workflow1" or "workflow2" (default = "workflow2)
    max_tax_dist = 3,       # maximum taxonomic distance acceptable between the target and the taxa in the database (default = 3)
    trait = "equation",     # trait to be searched for (default = "equation")
    gen_sp_dist = 0.5       # taxonomic distance between a genus and a species(default = 0.5)
)
```

See the docs for more details.

## Companion Scripts

`companion_scripts` contains all the scripts used to create access and analyse the database. The different folders hold scripts for different tasks. The numbers of the folders and the numbers of the scripts within the folders indicate in which order the scripts should be run.

There are two scripts that hold general functions:

* `special_names.R`: Function to generate the so-called *special names* that the database uses.
* `02_function_plotting_theme.R`: Customised ggplot2 plotting function

### 01_data_cleaning

The data cleaning folder holds scripts that are used to clean the raw data that was compiled in excel files. The raw data files are then saved as .rds files and stored in the *database* folder.

### 02_create_database

There are four scripts in this folder. The first two are scripts that detect whether the correct packages are installed and describe the taxonomic distance between taxonomic ranks like families and genera:

* `01_version_package_warnings.R`
* `02_taxon_matrix_function.R`

Then, there is a script that used to create the higher-level taxonomic graphs. This works by first harmonising all taxon names in the equation database to three different taxonomic backbones: COL, GBIF and ITIS. Once the names are harmonised, we extract either the family or order of each taxon name. Descendent taxa from all unique families and orders are then extracted and compiled into igraph objects that describe how the different taxon names i.e. species, genera, families, orders etc. relate to each other. These igraph objects are exported as .rds files and stored in the *database* folder.

`03_create_taxon_database.R`

The final script is used to add biogeographical realm, major habitat type and ecoregion information to each equation in the database using the latitude and longitude data associated with each equation and Abell et al.'s (2008) global ecoregion map.

`04_set_freshwater_ecoregion_data.R`

### 05_accuracy_analysis

This folder contains the script where we test the accuracy of our method for matching names to appropriate equations. Specifically, we compare the biomass generated by selecting equations in the database to actual measured biomass that we compiled from the literature along with biomass generated from equations selected by experts.

* `01_accuracy_test_script.R`

### 06_database_characteristics

These scripts are used to examine the taxonomic and geographical coverage of the equations in our database.

## Development

There's a [devcontainer](https://containers.dev/) setup included. If you use
VSC you should be prompted to open the project in a container automatically.

`devtools` are bundled with the devcontainer. Load `library(devtools)` and you
have `load_all()`, `test()` and `check()` ready at hand.

We use `renv` to provide reproducibility as far as it gets with R.
Use `renv::snapshot()` after changing dependencies, `renv::install()` to install declared versions
of the dependencies and `renv::update()` to update to latest CRAN versions (before pushin to CRAN).

The database files will be put into an appdata dir (given by `rappdirs`) when executing
tests or when people load the actual package. If you made changes to the DB files and need
to update the files in the appdata dir there's the utility function `update_user_db()`.
