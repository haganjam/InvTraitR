
# FW_invert_biomass_allometry

Pipeline to assign biomass-length allometry equations and other traits to a taxonomic name based on taxonomic hierarchy.

## next steps:

The next step is to input all the equation and length data and create the database. Then, make a test script to robustly test these different functions.

Finally, we will need to integrate the gbif database as well. This should be simply about creating gbif versions of the different databases but it will be more complicated... Probably.

## scripts: functions

This folder contains functions that are used in multiple parts of the pipeline.

> 01_get_taxon_id_function.R

To access the taxonomic backbones from taxize, a given taxon name needs to be converted into a numerical ID. This function takes a taxon name and the database and outputs the numerical ID required to use the other taxize functions.

The functions uses two other functions: (1) get_gbifid2(). It implements automatic checks on the output rather than a dynamic, user interaction and (2) extract_genus() which extracts the first name from a binomial.

## scripts: create_database

This folder contains functions and scripts that are used to create a database containing the relevant allometric equation data and the taxonomic database which contains relatedness information for taxa in the equation database.

> 01_version_package_warnings.R

This script is sourced from all the other functions. It is designed to notify the user which packages were used to implement the functions. Moreover, it tests whether the user has the same version of the packages that were originally used to write the code and outputs the R-version used to write the code.

> 02_get_downstream_function.R

The downstream function from taxize gets taxon names downstream from a given numerical ID and database. However, it can sometimes get stuck due to a bad internet connection etc. Thus, this is a slight modification of the downstream that implements multiple tries if required.

> 03_gbif_downstream_function.R

Function to get all downstream taxa from a given order to the genus level using the gbif database. The downstream taxa are then processed into a taxon list and a distance matrix with the taxonomic relatedness between all taxon names.

> 04_itis_downstream_function.R

Function to get all downstream from a given order to the genus level using the itis database. The downstream taxa are then processed into a taxon list and a distance matrix with the taxonomic relatedness between all taxon names.

> 05_taxon_matrix_itis_function.R

The Itis database uses a complex taxononomic hierarchy structure. This structure can be found on pg. 12 at the following link: 

- https://www.itis.gov/pdf/ITIS_ConceptualModelEntityDefinition.pdf

This script implements a numerical conversion of this taxonomic hierarchy structure.

> 06_get_taxonomic_information_function.R

This scripts contains two functions that combine the previous functions to package the relevant information for a given taxon name.

> 07_create_equation_database.R

This script pulls together the equation databases and links them into a list which is outputted as .rds file.

> 08_create_taxon_database.R

This function combines all these functions to get the taxon information from order to genus for all taxon names in the equation database. The unique implementation in this script is to only get the taxonomic information from unique orders. These orders are then matched to the orders of the focal taxon names in the database which saves considerable computing time and memory.

Finally, the database is exported as a .rds file which will then be searchable.

## scripts: use_database

This folder contains the scripts that allow users to extract suitable length-mass allometric equations and traits for a given taxon name given taxonomic relatedness.

> 01_search_taxon_database.R

Coming soon...


