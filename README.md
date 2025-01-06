# ANSA: an automated workflow to standardise taxon names for South African alien species lists
The automated workflow Alien Names South Africa (ANSA) was developed to standardise the names of taxa for South African alien species lists. It was developed specifically using the taxonomic backbones implemented in South Africa. Although initially developed within the context of the species list (http://dx.doi.org/10.5281/zenodo.8217197) for the report "The Status of Biological Invasions and their Management in South Africa" (https://dx.doi.org/10.5281/zenodo.8217182), the workflow will be useful for South African alien species lists, more generally. 

The workflow makes use of three taxonomic backbones: that of the Global Biodiversity Information Facility (GBIF); that of the Botanical Database of Southern Africa (BODATSA); and that of The World Checklist of Vascular Plants (WCVP), accessed via Plants of the World Online (POWO) (https://powo.science.kew.org/).

The taxonomic backbone used depends on the taxonomic group:
1. The GBIF backbone is used for non-plant taxa (e.g., animals, fungi, and chromista)
2. The BODATSA backbone is used for alien plants recorded outside of captivity (i.e., those included in BODATSA)
3. The WCVP backbone is used for alien plants not found outside of captivity (i.e., those not included in BODATSA)

In addition to standardising a list of taxon names based on the three taxonomic backbones, the workflow obtains canonical names and higher taxonomic information for the taxa from GBIF, and flags issues and uncertainties that arise during the taxonomic standardisation. 

The workflow has been created in R software, version 4.4.0 (R Core Team 2024). 
