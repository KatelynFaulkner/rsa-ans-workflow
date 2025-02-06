# An automated workflow to standardise taxon names for South African alien species lists
This automated workflow was developed to standardise the names of alien taxa for South African alien species lists. It can be used to standardise digitised lists of taxon names that are in one of two formats: (1) full scientific name with authorship and date information (hereafter ‘scientific name’); and (2) scientific name without authorship and date information (hereafter ‘canonical name’). While the workflow can be used to taxonomically standardise canonical names, it is recommended that information on authorship and date is retained as this decreases errors and makes it easier to resolve issues. The workflow implements taxonomic backbones that are relevant to the South African context.

Three taxonomic backbones are implemented: that of the Global Biodiversity Information Facility (GBIF); that of the Botanical Database of Southern Africa (BODATSA); and that of The World Checklist of Vascular Plants (WCVP), accessed via Plants of the World Online (POWO) (https://powo.science.kew.org/).

The taxonomic backbone used depends on the taxonomic group:

1. The GBIF backbone is used for non-plant taxa (e.g., animals, fungi, and chromista)

2. The BODATSA backbone is used for alien plants recorded outside of captivity or cultivation in South Africa (i.e., those taxa included in BODATSA)

3. The WCVP backbone is used for alien plants not found outside of captivity or cultivation in South Africa (i.e., those taxa not included in BODATSA)

In addition to standardising a list of taxon names based on the three taxonomic backbones, the workflow obtains canonical names and higher taxonomic information for the taxa from GBIF, and detects, corrects, and flags potential issues and errors. 

The workflow has been created in R software, version 4.4.0 (R Core Team 2024). 
