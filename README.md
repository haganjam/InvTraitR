
# FW_invert_biomass_allometry

Pipeline to assign biomass-length allometry equations to a taxonomic name based on taxonomic hierarchy and geographic (or environmental) proximity

## scripts: functions 

> taxonomic_query_functions.R

This script contains a function to pull the upstream and downstream taxa from an entry in the database from three different taxonomic backbones available via the taxize package:

+ GBIF

+ ITIS

+ BOLD

> taxonomic_database_search_functions.R

This script contains a function to search the taxonomic databases and find the most suitable allometric equation in our database given a taxon name.

## scripts: implementation

> taxonomic query

This script uses the function in 'taxonomic_query_function.R' to link up equations from our database with the taxonomic backbones and output a database with linked taxonomic information.

