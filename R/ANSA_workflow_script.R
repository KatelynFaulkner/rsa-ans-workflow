
# Preparation of the R environment

## Step 1: Install and load R cran packages

packages = c("tidyr", "dplyr",
             "rgbif", "stringdist", "gtools", "stringr", "sjmisc", "rWCVP") # eight packages are required

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

## Step 2: Install and load remote package from R cran, and install packages from GitHub (those not on the R cran)

remotes.check<-"remotes" %in% rownames(installed.packages()) # one package required from github
if(remotes.check == FALSE){
  install.packages("remotes")
  remotes::install_github('matildabrown/rWCVPdata')
}
require('rWCVPdata')

# Load input data, provide required information, and prepare these data

## Step 3: Load dataset from BODATSA

BODATSA<-read.csv("data\\inputs\\BODATSA.csv")

## Step 4: Provide required information

### Specify the date that BODATSA was downloaded

BODATSADate<-"02 December 2024" # put the date in the format shown here

## Step 5: Prepare BODATSA data

### View header of data 

head(BODATSA) 

### Assign NA to records where data are missing data

BODATSA[BODATSA == ""] <- NA 

### Add a column with the canonical names to BODATSA dataset

BODATSA <- BODATSA %>% unite("canonicalName", c(Genus, Sp1, Rank1, Sp2, 
                                                Rank2, Sp3), na.rm = TRUE, sep = " ", remove = FALSE)

### Add a column with the scientific names to BODATSA dataset

BODATSA <- BODATSA %>% unite("scientificName", c(Genus, Sp1, Author1, Rank1, Sp2, Author2, Rank2, Sp3, Author3), na.rm = TRUE, sep = " ", remove = FALSE)

## Step 6: Load list of names for standardisation

NamesDat<-read.csv("data\\inputs\\OriginalNames.csv")

## Step 7: Provide required information

### Specify whether the names are 'canonical' (without authorities) or 'scientific' (with authorities)

NameType<-"canonical" # put 'canonical' or 'scientific' here

## Step 8: Prepare of list of names for standardisation

### View header of data and get the number of taxa

head(NamesDat) # Look at the first five rows of the data. See below for an example of the output

length(NamesDat$verbatimScientificName) # how many taxa in the dataset

### Standardise the spacing formats 

NamesDat<-data.frame(verbatimScientificName = str_replace_all(
  NamesDat$verbatimScientificName, "\\s", " ")) # standardise spacing formats

NamesDat$verbatimScientificName<-str_squish(NamesDat$verbatimScientificName) #remove double spaces

## Step 9: Get canonical name and higher taxonomic information for the taxa from GBIF

GBIFCnames<-c("canonicalName", "kingdom",       
              "phylum", "class", "order","family") # create a list of ideal column names for the GBIF data

HighTax<-{} # create empty dataframe for the results
m<-length(NamesDat$verbatimScientificName)
for(j in 1:m){
  x<-data.frame(name_backbone(NamesDat$verbatimScientificName[j]))
  x[setdiff(GBIFCnames, names(x))] <- NA
  x<-x[GBIFCnames]
  z<-x[, c('canonicalName', 'kingdom', 'phylum', 'class', 'order','family')]
  HighTax<-rbind(HighTax, z) 
} # Loop that gets the canonical name and higher taxonomic information for each taxon

NamesDat<-cbind(NamesDat, HighTax) # bind the taxon names with the data from GBIF

## Step 10: Separate dataset into plants and other organisms

unique(NamesDat$kingdom) # list of kingdoms in the dataset. See below for an example of output

PlantDat<-NamesDat %>% filter(kingdom == "Plantae") # create a plants subset 

OtherOrgDat<-NamesDat %>% filter(kingdom != "Plantae"|is.na(kingdom)) # create a non-plants subset 

# Taxonomic standardisation for non-plant taxa using GBIF backbone

## Step 11: Check names non-plants against the GBIF backbone 

GBIFCnames<-c("usageKey","scientificName", "canonicalName", "rank",          
              "status", "confidence", "matchType", "kingdom",       
              "phylum", "order", "family", "genus",         
              "species", "kingdomKey", "phylumKey", "classKey",      
              "orderKey", "familyKey", "genusKey", "speciesKey",    
              "synonym", "class") # create a list of ideal column names for the GBIF data

GBIFMatch<-{} # create empty dataframe for the results
m<-length(OtherOrgDat$verbatimScientificName)
for(j in 1:m){
  x<-data.frame(name_backbone(OtherOrgDat$verbatimScientificName[j], kingdom = OtherOrgDat$kingdom[j]))
  x[setdiff(GBIFCnames, names(x))] <- NA
  x<-x[GBIFCnames]
  if(is.na(x$family)){
    u<-cbind(OtherOrgDat[j,], x$scientificName, x$synonym, x$rank, x$status, x$matchType, x$usageKey)
  }
  else if(x$status == "ACCEPTED"|x$status == "DOUBTFUL"){
    u<-cbind(OtherOrgDat[j,], x$scientificName, x$synonym, x$rank, x$status, x$matchType, x$usageKey)
  }
  else if(x$status == "SYNONYM"){
    a<-name_backbone(x$species, kingdom = x$kingdom)
    u<-cbind(OtherOrgDat[j,], a$scientificName, x$synonym, x$rank, x$status, x$matchType, x$usageKey)
    }
  colnames(u)[c(ncol(u)-5, ncol(u)-4, ncol(u)-3,ncol(u)-2, ncol(u)-1,
                ncol(u))]<-c("scientificName", "synonym", "rank", "status",
                             "matchtype", "usageKey")
  GBIFMatch<-rbind(GBIFMatch,u)
} # Loop that checks each taxon name against the GBIF backbone
 
 ### Write interim output to file 

GBIFInterim<-GBIFMatch

write.csv(GBIFInterim, "data//interim//GBIFInterim.csv")

## Step 12: Look at the results and flag issues

### Look at what the results look like. See below for an example of the output.

head(GBIFMatch) # look at the first 5 rows of the dataset

### First look at taxonomic issues

GBIFMatch %>%
  group_by(status) %>%
  count()

GBIFMatch %>%
  group_by(matchtype) %>%
  count()

GBIFMatch %>%
  group_by(rank) %>%
  count()

### Flag taxonomic issues

### Flag synonyms

GBIFMatch <- GBIFMatch %>% 
  rename("synonym.flag"= "synonym")

### Flag taxon names that were not exact matches (i.e., fuzzy matches)

GBIFMatch = GBIFMatch %>% 
  mutate(uncertain.flag=case_when((matchtype =="FUZZY")~TRUE,
                                  TRUE~FALSE))

### Flag doubtful taxon names

GBIFMatch = GBIFMatch %>% 
  group_by(status) %>% 
  mutate(doubtful.flag = ifelse(status =="DOUBTFUL" & !is.na(status), TRUE, FALSE))

### Flag unrecognised taxa-those only recognised at higher levels than the taxon (i.e., genus or above level for a species, or a species level for a sub-specific entity)

GBIFMatch = GBIFMatch %>% 
  mutate(unrecognised.flag=case_when((matchtype =="HIGHERRANK")~TRUE,
                                     (rank == "GENUS")~TRUE,
                                     (matchtype == "NONE") ~TRUE,
                                     TRUE~FALSE))

### Flag if the result should be checked for a possible error or issue

GBIFMatch<-GBIFMatch %>% 
  mutate(possibleError.flag = ifelse(uncertain.flag =="TRUE" | doubtful.flag == "TRUE" | unrecognised.flag == "TRUE", TRUE, FALSE))

### Flag sub-specific entities

GBIFMatch = GBIFMatch %>% 
  group_by(rank) %>% 
  mutate(subspecific.flag = ifelse(rank =="SPECIES"|rank == "GENUS"|is.na(rank)|matchtype == "HIGHERRANK", FALSE, TRUE))

### Flag duplicates in returned scientific names (exclude duplicates at higher levels)

GBIFMatch = GBIFMatch %>% 
  group_by(scientificName) %>% 
  mutate(scientificNameduplicate.flag = ifelse(n() > 1 & rank !="GENUS" & rank !="KINGDOM" & rank !="PHYLUM" & rank !="CLASS" & rank != "ORDER" & rank != "FAMILY" & !is.na(rank),
                                               TRUE, FALSE))

### Flag duplicates in verbatim names 

GBIFMatch = GBIFMatch %>% 
  group_by(verbatimScientificName) %>% 
  mutate(verbatimNameduplicate.flag = n() > 1)

### Flag duplicates in canonical names (exclude duplicates at higher levels)

GBIFMatch = GBIFMatch %>% 
  group_by(canonicalName) %>% 
  mutate(canonicalNameduplicate.flag = ifelse(n() > 1 & rank !="GENUS" & rank !="KINGDOM" & rank !="PHYLUM" & rank !="CLASS" & rank != "ORDER" & rank != "FAMILY" & !is.na(rank),
                                                TRUE, FALSE))
## Step 13: Post-processing

### Add column with the source of the scientific name

GBIFMatch$scientificNameSource<-ifelse(GBIFMatch$unrecognised.flag == TRUE, 
                                       "NA", paste0("GBIF", " ", "(", format(Sys.Date(), "%d %B %Y"),")", " ","https://www.gbif.org/species/", GBIFMatch$usageKey))

### Remove unwanted columns

GBIFMatch<-GBIFMatch %>% select(!(rank:usageKey)) # only keep the relevant columns
names(GBIFMatch) # names of the columns. See below for an example of the new output

### Write interim output to file 

write.csv(GBIFMatch, "data//interim//NonPlantOrganisms_GBIF_standardisation.csv", row.names = FALSE)

# Taxonomic standardisation for plant taxa using BODATSA backbone

## Step 14: Match plant names with those in BODATSA

### Merge the two datasets to identify matching names. Based on canonical or scientific name, depending on inputs

if(NameType == "canonical"){
BODATSAMatch<-
  left_join(PlantDat, BODATSA, by = join_by("verbatimScientificName" == "canonicalName"), relationship = "many-to-many")
}

if(NameType == "scientific"){
  BODATSA<- BODATSA %>% select(!(canonicalName)) # drop canonicalName from BODATSA
  BODATSAMatch<-
    left_join(PlantDat, BODATSA, by = join_by("verbatimScientificName" == "scientificName"), keep = TRUE, suffix = c('.1', '.2'))
    }

### Replace synonyms in scientificName with accepted names

BODATSAMatch<-BODATSAMatch %>%
  mutate(scientificName = replace(scientificName, is.na(Accepted)==FALSE, Accepted[is.na(Accepted)==FALSE]))

### Identify taxa with multiple matches

BODATSAMatch = BODATSAMatch %>% 
  group_by(verbatimScientificName) %>% 
  mutate(multipleMatches = n() > 1)

### Write interim output to file 

BODATSAInterim<-BODATSAMatch

write.csv(BODATSAInterim, "data//interim//BODATSAInterim.csv")

### Resolve, if possible, taxa with multiple matches

#### Take the match with an accepted name if there is one

matchesMult<-BODATSAMatch%>% filter(multipleMatches == TRUE & TaxStat == "acc")
matchesFilt<-BODATSAMatch%>% filter(!verbatimScientificName %in% matchesMult$verbatimScientificName)
BODATSAMatch<-rbind(matchesMult, matchesFilt)

#### Deal with taxa with multiple matches that could not be resolved

#### Flag these taxa 

BODATSAMatch = BODATSAMatch %>% 
  group_by(verbatimScientificName) %>% 
  mutate(multipleMatchesUnresolved = n() > 1)

#### Concatenate the results for the taxa and separate with pipe delimiter '|'

BODATSAMatch<-BODATSAMatch %>%
  group_by(verbatimScientificName) %>%
  summarise(across(everything(), ~paste0(unique(.), collapse = "|")))

## Step 15: Look at the results and flag issues

### Look at what results look like

head(BODATSAMatch) # Look at the first 5 rows of the data. See below for example output

### First look at taxonomic issues

BODATSAMatch %>%
  group_by(TaxStat) %>%
  count()

BODATSAMatch %>%
  group_by(multipleMatches) %>%
  count()

BODATSAMatch %>%
  group_by(multipleMatchesUnresolved) %>%
  count()

### Flag taxonomic issues

### Flag synonyms

BODATSAMatch = BODATSAMatch %>% 
  mutate(synonym.flag = case_when((TaxStat =="syn") ~ TRUE, TRUE ~ FALSE))

### Flag taxon names that were unresolved multiple matches

BODATSAMatch = BODATSAMatch %>% 
  mutate(uncertain.flag=case_when((multipleMatchesUnresolved =="TRUE")~TRUE,
                                  TRUE~FALSE))

### Flag doubtful taxa (uncertain, missapplied)

BODATSAMatch <- BODATSAMatch %>% mutate(doubtful.flag = str_detect(TaxStat, "unc|misap"))

#### Flag unrecognised taxa

BODATSAMatch = BODATSAMatch %>% 
  group_by(Genus) %>% 
  mutate(unrecognised.flag = ifelse(Genus =="NA", TRUE, FALSE))

### Flag if the result should be checked for a possible error or issue

BODATSAMatch<-BODATSAMatch %>% 
  mutate(possibleError.flag = ifelse(uncertain.flag =="TRUE" | doubtful.flag == "TRUE" | unrecognised.flag == "TRUE", TRUE, FALSE))

### Flag sub-specific entities

BODATSAMatch<-BODATSAMatch %>% 
  mutate(sppType = str_split(canonicalName, " ") %>% 
           lengths)

BODATSAMatch<-BODATSAMatch %>% 
  mutate(subspecific.flag = ifelse(sppType > 2, TRUE, FALSE)) 

### Flag duplicates in returned scientific names (exclude NAs)

BODATSAMatch = BODATSAMatch %>% 
  group_by(scientificName) %>% 
  mutate(scientificNameduplicate.flag = ifelse(n() > 1 & scientificName !="NA", 
                                        TRUE , FALSE))

### Flag duplicates in verbatim names

BODATSAMatch = BODATSAMatch %>% 
  group_by(verbatimScientificName) %>% 
  mutate(verbatimNameduplicate.flag = n() > 1)

### Flag duplicates in canonical names

BODATSAMatch = BODATSAMatch %>% 
  group_by(canonicalName) %>% 
  mutate(canonicalNameduplicate.flag = ifelse(n() > 1 & sppType !=1, 
                                               TRUE , FALSE))

## Step 16: Post-processing

### Add column with the source of the scientific name

BODATSAMatch$scientificNameSource<-ifelse(BODATSAMatch$unrecognised.flag == "TRUE", 
                   "NA", paste0("BODATSA", " ", "(", BODATSADate,")"," ", "http://posa.sanbi.org/"))

### Remove unwanted columns

BODATSAMatch <- BODATSAMatch %>% ungroup()
BODATSAMatch<-BODATSAMatch %>% select(!(c(Family,Genus:multipleMatchesUnresolved, sppType))) # only keep the relevant columns
names(BODATSAMatch) # names of the columns. See below for an example of the new output

### Write interim output to file 

write.csv(BODATSAMatch, "data//interim//Plants_BODATSA_standardisation.csv", row.names = FALSE)

# Taxonomic standardisation for plant taxa using WCVP backbone

## Step 17: Match plant names unrecognised by BODATSA with those in WCVP

### Prepare list of names for matching

PlantUnRec<-BODATSAMatch %>% filter(unrecognised.flag == "TRUE") # subset of the unrecognised plants
PlantUnRec<-PlantUnRec %>% select(verbatimScientificName:family) # remove unwanted columns

### Match list of names with WCVP

if(NameType == "canonical"){
matches <- wcvp_match_names(PlantUnRec,
                            name_col="verbatimScientificName",
                            author_col=NULL,
                            fuzzy=TRUE,
                            progress_bar=FALSE)
}

if(NameType == "scientific"){
  matches <- wcvp_match_names(PlantUnRec,
                              name_col="canonicalName",
                              author_col=NULL,
                              fuzzy=TRUE,
                              progress_bar=FALSE)
   }

### Get the accepted scientific name for each taxon
u<-which(wcvp_names$plant_name_id %in% matches$wcvp_accepted_id)
v <- wcvp_names[u,] %>% unite("scientificName", c(taxon_name, taxon_authors), na.rm = TRUE, sep = " ", remove = FALSE)
matches<- left_join(matches, v, by = join_by("wcvp_accepted_id" == "plant_name_id"))
matches<-matches%>% select((c(verbatimScientificName:wcvp_accepted_id, scientificName, powo_id))) # select useful columns
colnames(matches)[colnames(matches) == "family.x"] <- "family"

### Write interim output to file 

WCVPInterim<- matches

write.csv(WCVPInterim, "data//interim//WCVPInterim.csv")

### Resolve, if possible, for taxa with multiple matches

#### If the scientific name with authorities is available take the match with that authority

if(NameType == "scientific"){
  matchesAuth<-matches%>% filter(multiple_matches == TRUE & match_type == "Exact (without author)")
  matchesAuth<-matchesAuth %>%
    filter(str_detect(verbatimScientificName, wcvp_authors))
  matchesAuth<-matchesAuth %>%
    mutate(match_type=replace(match_type, match_type=="Exact (without author)", "Exact"))
  matchesAuth<-matchesAuth %>%
    mutate(multiple_matches=replace(multiple_matches, multiple_matches=="TRUE", "FALSE"))
    matches<-matches %>% 
    filter(!verbatimScientificName %in% matchesAuth$verbatimScientificName)
  matches<-rbind(matches, matchesAuth)
    }

#### For remaining taxa with multiple matches, take the match with an accepted name if there is one

matchesMult<-matches%>% filter(multiple_matches == TRUE & wcvp_status == "Accepted")
matchesFilt<-matches%>% filter(!verbatimScientificName %in% matchesMult$verbatimScientificName)
matches<-rbind(matchesMult, matchesFilt)

### Deal with taxa with multiple matches that could not be resolved

#### Flag these taxa 

matches = matches %>% 
  group_by(verbatimScientificName) %>% 
  mutate(multiple_matchesUnresolved = n() > 1)

#### Concatenate the results for the taxa and separate with pipe delimiter '|'

WCVPMatch<-matches %>%
  group_by(verbatimScientificName) %>%
  summarise(across(everything(), ~paste0(unique(.), collapse = "|")))

## Step 18: Look at the results and flag issues

### Look at what results look like

head(WCVPMatch) # Look at the first 5 rows of the data. See below for example output

### First look at taxonomic issues

WCVPMatch %>%
  group_by(wcvp_status) %>%
  count()

WCVPMatch %>%
  group_by(match_type) %>%
  count()

WCVPMatch %>%
  group_by(wcvp_rank) %>%
  count()

WCVPMatch %>%
  group_by(multiple_matches) %>%
  count()

WCVPMatch %>%
  group_by(multiple_matchesUnresolved) %>%
  count()

### Flag taxonomic issues

### Flag synonyms

WCVPMatch = WCVPMatch %>% 
  group_by(wcvp_status) %>% 
  mutate(synonym.flag = ifelse(wcvp_status =="Synonym", TRUE, FALSE))

### Flag taxon names that were not exact matches, multiple matches not resolved, and not at species level

WCVPMatch<-WCVPMatch %>% 
  mutate(uncertain.flag = case_when((match_type =="NA")~FALSE,
                                    (multiple_matchesUnresolved == TRUE)~TRUE,
                                    (match_type == "Exact")~FALSE,
                                    (match_type =="Exact (without author)") ~FALSE,
                                    TRUE~TRUE))

### Flag doubtful taxon names (e.g., invalid, illegitimate, missapplied)

WCVPMatch<-WCVPMatch %>% 
  mutate(doubtful.flag = 
           ifelse(wcvp_status == "Accepted"|wcvp_status =="Synonym"|wcvp_status == "NA", "FALSE","TRUE"))

### Flag unrecognised taxa

WCVPMatch = WCVPMatch %>% 
  mutate(unrecognised.flag=ifelse(match_type =="NA", TRUE, FALSE))

### Flag if the result should be checked for a possible error or issue

WCVPMatch<-WCVPMatch %>% 
  mutate(possibleError.flag = ifelse(uncertain.flag =="TRUE" | doubtful.flag == "TRUE" | unrecognised.flag == "TRUE", TRUE, FALSE))

### Flag sub-specific entities

WCVPMatch = WCVPMatch %>% 
  group_by(wcvp_rank) %>% 
  mutate(subspecific.flag = ifelse(wcvp_rank =="Species"|wcvp_rank == "NA", FALSE, TRUE))


### Flag duplicates in returned scientific names (exclude NAs)

WCVPMatch = WCVPMatch %>% 
  group_by(scientificName) %>% 
  mutate(scientificNameduplicate.flag = ifelse(n() > 1 & scientificName !="NA", 
                                               TRUE , FALSE))
### Flag duplicates in verbatim names

WCVPMatch = WCVPMatch %>% 
  group_by(verbatimScientificName) %>% 
  mutate(verbatimNameduplicate.flag = n() > 1)

### Flag duplicates in canonical names (wcvp_name used in this case)

WCVPMatch = WCVPMatch %>% 
  group_by(wcvp_name) %>% 
  mutate(canonicalNameduplicate.flag = ifelse(n() > 1 & scientificName !="NA", 
                                               TRUE , FALSE))

## Step 19: Post-processing

### Add column with the source of the scientific name

WCVPMatch$scientificNameSource<-ifelse(WCVPMatch$unrecognised.flag == "TRUE", 
                                       "NA", paste0("POWO", " ", "(", format(Sys.Date(), "%d %B %Y"),")", " ","https://powo.science.kew.org/", WCVPMatch$wcvp_ipni_id))

### Remove unwanted columns

WCVPMatch <- WCVPMatch %>% ungroup()
WCVPMatch <- WCVPMatch %>% select(!(c(match_type:wcvp_accepted_id, powo_id:multiple_matchesUnresolved))) # only keep the relevant columns
names(WCVPMatch) # names of the columns. See below for an example of the new output


### Write interim output to file 

write.csv(WCVPMatch, "data//interim//Plants_WCVP_standardisation.csv", row.names = FALSE)

# Merge results into one output file

## Step 20: Create overall output dataset

### Format datasets from the standardisation processes so that they are the same and can be joined

GBIFMatch<-as.data.frame(GBIFMatch)
BODATSAMatch<-as.data.frame(BODATSAMatch)
WCVPMatch<-as.data.frame(WCVPMatch)

GBIFMatch = GBIFMatch %>%
  mutate(across(everything(), as.character))
BODATSAMatch = BODATSAMatch %>%
  mutate(across(everything(), as.character))
WCVPMatch = WCVPMatch %>%
  mutate(across(everything(), as.character))

### Remove unrecognised taxa from BODATSA match

BODATSAMatch<-BODATSAMatch%>% filter(unrecognised.flag == FALSE)

### Revise flagged duplicate plant names (verbatim and scientific) so that they are across the BODATSA and WCVP matches

PlantMatch<-rbind(BODATSAMatch, WCVPMatch) # Bind the datasets together

#### Flag duplicates in returned scientific names (exclude NAs)

PlantMatch = PlantMatch %>% 
  group_by(scientificName) %>% 
  mutate(scientificNameduplicate.flag = ifelse(n() > 1 & scientificName !="NA", 
                                                TRUE , FALSE))
### Flag duplicates in verbatim names

PlantMatch = PlantMatch %>% 
  group_by(verbatimScientificName) %>% 
  mutate(verbatimNameduplicate.flag = n() > 1)

### Bind the datasets 

AllMatch<-rbind(GBIFMatch, PlantMatch) # Bind all the datasets together in one dataset

### Remove flag for canonical name duplicates

AllMatch<-AllMatch %>% select(!(canonicalNameduplicate.flag)) # remove column

### Check if the number of taxa in the output dataset is that same as in the original dataset

length(AllMatch$verbatimScientificName) == length(NamesDat$verbatimScientificName)

### Reorder the taxa in the output dataset, so that they are in the same order as they were in the initial file

AllMatch<-AllMatch[match(NamesDat$verbatimScientificName, AllMatch$verbatimScientificName),]

# Step 21: Write file

write.csv(AllMatch, "data//outputs//Taxon_names_standardised.csv", row.names = FALSE)


