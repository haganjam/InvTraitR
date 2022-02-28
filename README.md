
# FW_invert_biomass_allometry

Pipeline to assign biomass-length allometry equations to a taxonomic name based on taxonomic hierarchy.

## scripts: functions 

> 01_version_package_warnings.R

This script is sourced from all the other functions. It is designed to notify the user which packages were used to implement the functions. Moreover, it tests whether the user has the same version of the packages that were originally used to write the code and outputs the R-version used to write the code.

> 02_get_gbifid2_function.R

This script contains a modified version of the taxize function: get_gbifid(). It implements automatic checks on the output rather than a dynamic, user interaction.

> 03_get_taxon_id_function.R

To access the taxonomic backbones from taxize, a given taxon name needs to be converted into a numerical ID. This function takes a taxon name and the database and outputs the numerical ID required to use the other taxize functions.

> 04_get_downstream_function.R

The downstream function from taxize gets taxon names downstream from a given numerical ID and database. However, it can sometimes get stuck due to a bad internet connection etc. Thus, this is a slight modification of the downstream that implements multiple tries if required.

> 05_gbif_downstream_function.R

Function to get all downstream taxa from a given order to the genus level using the gbif database. The downstream taxa are then processed into a taxon list and a distance matrix with the taxonomic relatedness between all taxon names.

> 06_itis_downstream_function.R

Function to get all downstream from a given order to the genus level using the itis database. The downstream taxa are then processed into a taxon list and a distance matrix with the taxonomic relatedness between all taxon names.

> 07_taxon_matrix_itis_function.R

The Itis database uses a complex taxononomic hierarchy structure. This structure can be found on pg. 12 at the following link: 

- https://www.itis.gov/pdf/ITIS_ConceptualModelEntityDefinition.pdf

This script implements a numerical conversion of this taxonomic hierarchy structure.

> 08_get_taxonomic_information_function.R

This scripts contains two functions that combine the previous functions to package the relevant information for a given taxon name.


## scripts: implementation

> create_taxon_database.R

This function combines all these functions to get the taxon information from order to genus for all taxon names in the database. The unique implementation in this script is to only get the taxonomic information from unique orders. These orders are then matched to the orders of the focal taxon names in the database which saves considerable computing time.

Finally, the database is exported as a .rds file which will then be searchable.

> search_taxon_database.R

Coming soon...


