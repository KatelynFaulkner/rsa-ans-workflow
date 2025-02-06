
# Step 1: Preparation of the R environment

## Install and load R CRAN packages

packages = c("tidyr", "dplyr",
             "rgbif", "stringr", "stringdist", "rWCVP", "purrr") # seven packages are required

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

## Install and load "remote" package from R CRAN, and install and load "rWCVPdata" package from GitHub

remotes.check<-"remotes" %in% rownames(installed.packages()) 
if(remotes.check == FALSE){
  install.packages("remotes")
  remotes::install_github('matildabrown/rWCVPdata') # one package required from GitHub
}
require('rWCVPdata') # load package from GitHub

# Step 2: Preparation of data

## Step 2a: Preparation of BODATSA dataset

### Load dataset from BODATSA

BODATSA<-read.csv("data\\inputs\\BODATSA.csv") # file must be named "BODATSA.csv"

### Provide required information on BODATSA

BODATSADate<-"02 December 2024" # stipulate the date the BODATSA data were downloaded in the format shown here

### View header of data 

head(BODATSA) 

### Assign NA to records where data are missing

BODATSA[BODATSA == ""] <- NA 

### In "Rank1" column replace "ssp." with "subsp."

BODATSA<-BODATSA %>%
  mutate(Rank1=replace(Rank1, Rank1=="ssp.", "subsp."))

### Add a column with the canonical names of the taxa to BODATSA dataset

BODATSA <- BODATSA %>% unite("canonicalName", c(Genus, Sp1, Rank1, Sp2, 
                                                Rank2, Sp3), na.rm = TRUE, sep = " ", remove = FALSE)

### Add a column with the scientific names of the taxa to BODATSA dataset

BODATSA <- BODATSA %>% unite("scientificName", c(Genus, Sp1, Author1, Rank1, Sp2, Author2, Rank2, Sp3, Author3), na.rm = TRUE, sep = " ", remove = FALSE)

### Correct and flag incorrectly formatted names in the "Accepted" column

#### Identify the names in the "Accepted" column that are incorrectly formatted

FormatIssue<-BODATSA%>% filter(!Accepted %in% BODATSA$scientificName & !is.na(Accepted)) 

#### Identify the accepted name in the correct format using fuzzy matching

fullScientificName<-{} # create empty dataset
for(i in 1:length(FormatIssue$Accepted)){ # loop through names with format issues
  keywords<-unlist(str_split(FormatIssue$Accepted[i], " ")) # split name
  y<-paste(keywords[1:2], collapse = " ") # create string that only contains the genus and species of the name
  z <-  BODATSA %>%
    filter(str_detect(scientificName, y)) # subset BODATSA to only include records for that genus and species
  a<-stringdist(FormatIssue$Accepted[i], z$scientificName, method="jaccard") # fuzzy match the name to those of the genus and species with jaccard method
  b<-which(a==min(a)) # get the closest match
  if(length(b)==1){
  c<-z$scientificName[b] # take most similar match if there is only one
  }
  if(length(b)>1){ # if more than one closest match
   issues<-data.frame(z$scientificName[b]) 
  names(issues)<-"scientificName"
   issues<-issues %>%
      mutate(text = str_split(scientificName, '\\s+'), # split the names of the potential matches
             words_correct = map_dbl(text, ~sum(.x %in% keywords)))
    d<-which(issues$words_correct == max(issues$words_correct)) # determine number of components in the potential matches that are equal to those in the taxon name
    if(length(d) == 1){
    c<-issues$scientificName[d] # take the one with the highest number of equal components
  }
  else if(length(d)> 1){ # if more than one option with the highest number of equal components
    g<-stringdist(FormatIssue$Accepted[i], issues$scientificName[d], method="lv") # fuzzy match using lv method
    h<-which(g==min(g)) # identify the closest match
    c<-issues$scientificName[d][h] # take the one that is the closest match
  }
  }
  fullScientificName<-rbind(fullScientificName, c) # bind the name identified through fuzzy matching to result data frame
}

FormatIssue<-cbind(FormatIssue,fullScientificName) # bind the matched names to those with format issues

#### Replace incorrectly formatted names in "Accepted" column with those identified through fuzzy matching, and flag

FormatIssueAcc<-FormatIssue %>% select(Accepted, fullScientificName) # remove unnecessary columns
FormatIssueAcc<-unique(FormatIssueAcc) # remove duplicated rows

BODATSA <- left_join(BODATSA, FormatIssueAcc, by = c("Accepted" = "Accepted")) # join the matched names to BODATSA

BODATSA <- BODATSA %>% 
  rename("AcceptedIssues"= "Accepted") # rename the original "Accepted" column "AcceptedIssues"

#### Create a new "Accepted" column with the matched names 

BODATSA<-BODATSA %>%
  mutate(Accepted = replace(AcceptedIssues, is.na(fullScientificName)==FALSE, fullScientificName[is.na(fullScientificName)==FALSE]))

#### Create a column indicating where changes were made

BODATSA<-BODATSA %>%
  mutate(synAccIssue = ifelse(!is.na(fullScientificName), TRUE, FALSE)) # column is named "synAccIssues"

#### Delete unnecessary column

BODATSA$fullScientificName<-NULL

#### Write to file

write.csv(BODATSA, "data//interim//BODATSAPrepared.csv", row.names = FALSE)

## Step 2b: Preparation of list of taxon names

### Load list of names for standardisation

NamesDat<-read.csv("data\\inputs\\OriginalNames.csv") # file must be named "OriginalNames.csv"

### Provide required information on list of names for standardisation

NameType<-"canonical" # here stipulate format of taxon names: 'canonical' or 'scientific'

### View header of data and get the number of taxa

head(NamesDat) # Look at the first five rows of the data.

length(NamesDat$verbatimScientificName) # how many taxa are in the dataset

### Standardise the spacing formats 

NamesDat<-data.frame(verbatimScientificName = str_replace_all(
  NamesDat$verbatimScientificName, "\\s", " ")) # standardise spacing formats

NamesDat$verbatimScientificName<-str_squish(NamesDat$verbatimScientificName) # remove double spaces

# Step 3: Match all taxon names to GBIF backbone

## Match all names to GBIF backbone, and get canonical names and higher taxonomic information

GBIFCnames<-c("usageKey","scientificName", "canonicalName", "rank",          
              "status", "confidence", "matchType", "kingdom",       
              "phylum", "order", "family", "genus",         
              "species", "kingdomKey", "phylumKey", "classKey",      
              "orderKey", "familyKey", "genusKey", "speciesKey",    
              "synonym", "class") # create a list of ideal column names for the GBIF data

GBIFMatch<-{} # create empty dataframe for the results
m<-length(NamesDat$verbatimScientificName)
for(j in 1:m){ # loop through the taxon names and check each taxon name against the GBIF backbone
  x<-data.frame(name_backbone(NamesDat$verbatimScientificName[j])) # get information from GBIF for that name
  x[setdiff(GBIFCnames, names(x))] <- NA 
  x<-x[GBIFCnames]
  if(is.na(x$family)){ # if the name was not matched, bind the required data
    u<-cbind(NamesDat[j,], x$canonicalName, x$kingdom, x$phylum, x$class, x$order,x$family, x$scientificName, x$synonym, x$rank, x$status, x$matchType, x$usageKey)
  }
  else if(x$status == "ACCEPTED"|x$status == "DOUBTFUL"){ # if the name is accepted or doubtful, bind the required data
    u<-cbind(NamesDat[j,], x$canonicalName, x$kingdom, x$phylum, x$class, x$order,x$family, x$scientificName, x$synonym, x$rank, x$status, x$matchType, x$usageKey)
  }
  else if(x$status == "SYNONYM"){ # if name is a synonym, fetch the accepted scientific name, and bind required columns
    a<-name_backbone(x$species, kingdom = x$kingdom)
    u<-cbind(NamesDat[j,], x$canonicalName, x$kingdom, x$phylum, x$class, x$order,x$family, a$scientificName, x$synonym, x$rank, x$status, x$matchType, x$usageKey)
  }
  colnames(u)[c(ncol(u)-12, ncol(u)-11, ncol(u)-10, ncol(u)-9, ncol(u)-8, ncol(u)-7, ncol(u)-6, ncol(u)-5, ncol(u)-4, ncol(u)-3,ncol(u)-2, ncol(u)-1,
                ncol(u))]<-c("verbatimScientificName", "canonicalName", "kingdom", "phylum", "class", "order","family", "scientificName", "synonym", "rank", "status","matchtype", "usageKey") # name the columns
  GBIFMatch<-as.data.frame(rbind(GBIFMatch,u)) # bind data to result dataframe
} 
 
## Write interim output to file 

GBIFInterim<-GBIFMatch

write.csv(GBIFInterim, "data//interim//GBIFInterim.csv", row.names = FALSE)

# Step 4: Extract plants

unique(GBIFMatch$kingdom) # list of kingdoms in the output

## Create a plants subset

PlantDat<-GBIFMatch %>% filter(kingdom == "Plantae")  # create subset of plants

PlantDat<-PlantDat %>% select(verbatimScientificName:family) # only keep the relevant columns
names(PlantDat) # names of the columns

## Create a non-plants subset

GBIFMatch<-GBIFMatch %>% filter(kingdom != "Plantae"|is.na(kingdom)) # create a non-plants subset 

# Step 5 Assess and flag issues for non-plant taxa

## Assess issues

head(GBIFMatch) # look at the first 5 rows of the dataset

### Look at how many have each status (e.g., accepted, synonym etc)
GBIFMatch %>%
  group_by(status) %>%
  count()

### Look at how many had a specific match type (e.g., exact, fuzzy etc)
GBIFMatch %>%
  group_by(matchtype) %>%
  count()

### Look at how many fall into different ranks (e.g., species, subspecies etc)
GBIFMatch %>%
  group_by(rank) %>%
  count()

## Flag taxonomic issues

### Flag synonyms

GBIFMatch <- GBIFMatch %>% 
  rename("synonym.flag"= "synonym") # give TRUE if a synonym

### Flag taxon names that were not exact matches (i.e., fuzzy matches)

GBIFMatch = GBIFMatch %>% 
  mutate(uncertain.flag=case_when((matchtype =="FUZZY")~TRUE,
                                  TRUE~FALSE)) # give TRUE if a fuzzy match

### Flag doubtful taxon names (e.g., taxonomic homonyms, missapplied taxon names)

GBIFMatch = GBIFMatch %>% 
  group_by(status) %>% 
  mutate(doubtful.flag = ifelse(status =="DOUBTFUL" & !is.na(status), TRUE, FALSE)) # give TRUE if doubtful

### Flag unmatched taxa (those not matched to an accepted name, or matched at higher levels (i.e., genus or above level for a species))

GBIFMatch = GBIFMatch %>% 
  mutate(unmatched.flag=case_when((matchtype =="HIGHERRANK")~TRUE,
                                     (rank == "GENUS")~TRUE,
                                     (matchtype == "NONE") ~TRUE,
                                     TRUE~FALSE)) # give TRUE if not matched

### Flag if the result should be checked for a possible error or issue (unmatched, doubtful, and uncertain matches)

GBIFMatch<-GBIFMatch %>% 
  mutate(possibleError.flag = ifelse(uncertain.flag =="TRUE" | doubtful.flag == "TRUE" | unmatched.flag == "TRUE", TRUE, FALSE)) # TRUE if there could be an issue

### Flag sub-specific entities

GBIFMatch = GBIFMatch %>% 
  group_by(rank) %>% 
  mutate(subspecific.flag = ifelse(rank =="SPECIES"|rank == "GENUS"|is.na(rank)|matchtype == "HIGHERRANK", FALSE, TRUE)) # TRUE if sub-specific

### Flag duplicates in returned scientific names (exclude duplicates if name is at higher taxonomic levels (e.g., genus))

GBIFMatch = GBIFMatch %>% 
  group_by(scientificName) %>% 
  mutate(scientificNameDuplicate.flag = ifelse(n() > 1 & rank !="GENUS" & rank !="KINGDOM" & rank !="PHYLUM" & rank !="CLASS" & rank != "ORDER" & rank != "FAMILY" & !is.na(rank),
                                               TRUE, FALSE)) # TRUE if a duplicate

### Flag duplicates in verbatim names 

GBIFMatch = GBIFMatch %>% 
  group_by(verbatimScientificName) %>% 
  mutate(verbatimNameDuplicate.flag = n() > 1) # TRUE if a duplicate

# Step 6: Preparation of output for non-plant taxa matched to GBIF

### Add column with the source of the scientific name, the date accessed, and web-link

GBIFMatch$scientificNameSource<-ifelse(GBIFMatch$unmatched.flag == TRUE, 
                                       "NA", paste0("GBIF", " ", "(", format(Sys.Date(), "%d %B %Y"),")", " ","https://www.gbif.org/species/", GBIFMatch$usageKey))
### Remove unwanted columns

GBIFMatch<-GBIFMatch %>% select(!(rank:usageKey)) # only keep the relevant columns
names(GBIFMatch) # names of the columns

# Step 7: Match plant names to BODATSA backbone

## Merge plant dataset and BODATSA

### Merge based on canonical name

if(NameType == "canonical"){
BODATSAMatch<-
  left_join(PlantDat, BODATSA, by = join_by("verbatimScientificName" == "canonicalName"), relationship = "many-to-many")
}

### Merge based on scientific name (with author information)

if(NameType == "scientific"){
  BODATSA<- BODATSA %>% select(!(canonicalName)) # drop canonicalName from BODATSA
  BODATSAMatch<-
    left_join(PlantDat, BODATSA, by = join_by("verbatimScientificName" == "scientificName"), keep = TRUE, suffix = c('.1', '.2'))
    }

### Replace synonyms in "scientificName" column with accepted names from "Accepted" column

BODATSAMatch<-BODATSAMatch %>%
  mutate(scientificName = replace(scientificName, is.na(Accepted)==FALSE, Accepted[is.na(Accepted)==FALSE]))

### Identify taxa with multiple matches

BODATSAMatch = BODATSAMatch %>% 
  group_by(verbatimScientificName) %>% 
  mutate(multipleMatches = n() > 1) # TRUE if there were multiple matches

### Resolve, if possible, names for taxa with multiple matches

#### Take the match that is an accepted scientific name, if there is one

matchesMult<-BODATSAMatch%>% filter(multipleMatches == TRUE & TaxStat == "acc") # create subset with multiple matches with accepted scientific names
matchesFilt<-BODATSAMatch%>% filter(!verbatimScientificName %in% matchesMult$verbatimScientificName) # create a subset without the above names
BODATSAMatch<-rbind(matchesMult, matchesFilt) # bind the subsets together

#### Address taxa with multiple matches that could not be resolved

#### Flag these taxa 

BODATSAMatch = BODATSAMatch %>% 
  group_by(verbatimScientificName) %>% 
  mutate(multipleMatchesUnresolved = n() > 1) # TRUE if multiple matches are unresolved

#### Concatenate the multiple results for these taxa and separate with a pipe delimiter '|'

BODATSAMatch<-BODATSAMatch %>%
  group_by(verbatimScientificName) %>%
  summarise(across(everything(), ~paste0(unique(.), collapse = "|")))

### Write interim output to file 

BODATSAInterim<-BODATSAMatch

write.csv(BODATSAInterim, "data//interim//BODATSAInterim.csv", row.names = FALSE)

# Step 8: Assess and flag issues for plants

## Assess issues

### Look at what results look like

head(BODATSAMatch) # Look at the first 5 rows of the data

### Look at how many have each status (e.g., accepted, synonym etc)

BODATSAMatch %>%
  group_by(TaxStat) %>%
  count()

### Look at how many were assigned accepted scientific names that had been corrected through fuzzy matching (Step 2a)  

BODATSAMatch %>%
  group_by(synAccIssue) %>%
  count() # TRUE for taxa with this issue

### Look at how many had multiple matches
BODATSAMatch %>%
  group_by(multipleMatches) %>%
  count() # TRUE for taxa with multiple matches

### Look at how many had multiple matches that were not resolved
BODATSAMatch %>%
  group_by(multipleMatchesUnresolved) %>%
  count()  # TRUE for taxa with unresolved multiple matches

## Flag taxonomic issues

### Flag synonyms

BODATSAMatch = BODATSAMatch %>% 
  mutate(synonym.flag = case_when((TaxStat =="syn") ~ TRUE, TRUE ~ FALSE)) # TRUE if a synonym

### Flag taxon names that were not exact matches (i.e., accepted names corrected through fuzzy matching)

BODATSAMatch = BODATSAMatch %>% 
  mutate(uncertain.flag=case_when((synAccIssue =="TRUE")~TRUE,
                                  TRUE~FALSE)) # TRUE if fuzzy matching involved

### Flag doubtful taxa (uncertain, missapplied, and those which have multiple matches that are not resolved)

BODATSAMatch <- BODATSAMatch %>%  mutate(doubtful.flag=case_when((multipleMatchesUnresolved =="TRUE")~TRUE,
                                                                  TRUE~FALSE)) %>% 
  mutate(doubtful.flag = str_detect(TaxStat, "unc|misap")) # TRUE if doubtful

#### Flag unmatched taxa

BODATSAMatch = BODATSAMatch %>% 
  group_by(Genus) %>% 
  mutate(unmatched.flag = ifelse(Genus =="NA", TRUE, FALSE)) # TRUE if unmatched

### Flag if the result should be checked for a possible error or issue (unmatched, doubtful, and uncertain matches)

BODATSAMatch<-BODATSAMatch %>% 
  mutate(possibleError.flag = ifelse(uncertain.flag =="TRUE" | doubtful.flag == "TRUE" | unmatched.flag == "TRUE", TRUE, FALSE)) # TRUE if possible error

### Flag sub-specific entities

BODATSAMatch <- BODATSAMatch %>% 
  group_by(Rank1) %>% 
  mutate(subspecific.flag = ifelse(Rank1 =="NA", FALSE, TRUE)) # TRUE if sub-specific

### Flag duplicates in returned scientific names (exclude NAs)

BODATSAMatch = BODATSAMatch %>% 
  group_by(scientificName) %>% 
  mutate(scientificNameDuplicate.flag = ifelse(n() > 1 & scientificName !="NA", 
                                        TRUE , FALSE)) # TRUE if a duplicate

### Flag duplicates in verbatim names

BODATSAMatch = BODATSAMatch %>% 
  group_by(verbatimScientificName) %>% 
  mutate(verbatimNameDuplicate.flag = n() > 1) # TRUE if a duplicate

# Step 9: Extract unmatched plants

### Create a subset of unmatched plants
PlantUnMatched<-BODATSAMatch %>% filter(unmatched.flag == "TRUE") # subset of the unmatched plants
PlantUnMatched<-PlantUnMatched %>% select(verbatimScientificName:family) # remove unwanted columns

### Create a subset of matched plants
BODATSAMatch<-BODATSAMatch%>% filter(unmatched.flag == FALSE) # subset of matched plants

# Step 10: Preparation of output for plants matched to BODATSA 

### Add column with the source of the scientific name, the date that BODATSA was downloaded, and a web-link

BODATSAMatch$scientificNameSource<-paste0("BODATSA", " ", "(", BODATSADate,")"," ", "http://posa.sanbi.org/")

### Remove unwanted columns

BODATSAMatch <- BODATSAMatch %>% ungroup() # ungroup columns
BODATSAMatch<-BODATSAMatch %>% select(!(c(Family,Genus:multipleMatchesUnresolved))) # only keep the relevant columns
names(BODATSAMatch) # names of the columns

# Step 11: Match unmatched plant names to WCVP backbone

### Match if verbatim names are canonical (without author information)
if(NameType == "canonical"){
matches <- wcvp_match_names(PlantUnMatched,
                            name_col="verbatimScientificName",
                            author_col=NULL,
                            fuzzy=TRUE,
                            progress_bar=FALSE)
}

### Match if verbatim names are scientific (with author information)
if(NameType == "scientific"){
  matches <- wcvp_match_names(PlantUnMatched,
                              name_col="canonicalName",
                              author_col=NULL,
                              fuzzy=TRUE,
                              progress_bar=FALSE)
   }

### Get the accepted scientific name for each taxon
u<-which(wcvp_names$plant_name_id %in% matches$wcvp_accepted_id) # which rows in the WCVP dataset have the ids of the accepted taxon names
v <- wcvp_names[u,] %>% unite("scientificName", c(taxon_name, taxon_authors), na.rm = TRUE, sep = " ", remove = FALSE) # create a dataset of those names
matches<- left_join(matches, v, by = join_by("wcvp_accepted_id" == "plant_name_id")) # join that dataset to the matches dataset
matches<-matches%>% select((c(verbatimScientificName:wcvp_accepted_id, scientificName, powo_id))) # select useful columns
colnames(matches)[colnames(matches) == "family.x"] <- "family" # rename column

### Resolve, if possible, taxon names with multiple matches

#### If the scientific name with author information is available take the match with that author

if(NameType == "scientific"){
  matchesAuth<-matches%>% filter(multiple_matches == TRUE & match_type == "Exact (without author)") # multiple matches with exact matching
  matchesAuth<-matchesAuth %>%
    filter(str_detect(verbatimScientificName, fixed(wcvp_authors))) # detect authors in the verbatim name
  matchesAuth<-matchesAuth %>%
    mutate(match_type=replace(match_type, match_type=="Exact (without author)", "Exact")) # replace the match type with "Exact"
  matchesAuth<-matchesAuth %>%
    mutate(multiple_matches=replace(multiple_matches, multiple_matches==TRUE, FALSE)) # change multiple_matches to FALSE
    matches<-matches %>% 
    filter(!verbatimScientificName %in% matchesAuth$verbatimScientificName) # get a subset that doesn't contain resolved matches 
  matches<-rbind(matches, matchesAuth) # bind resolved matches to the subset
    }

#### For remaining taxa with multiple matches, select the match with an accepted scientific name, if there is one

matchesMult<-matches%>% filter(multiple_matches == TRUE & wcvp_status == "Accepted") # subset of multiple matches, where there is an accepted taxon name
matchesFilt<-matches%>% filter(!verbatimScientificName %in% matchesMult$verbatimScientificName) # remove these taxa from the matches dataset
matches<-rbind(matchesMult, matchesFilt) # bind the two datasets together

### Address taxa with multiple matches that could not be resolved

#### Flag these taxa 

matches = matches %>% 
  group_by(verbatimScientificName) %>% 
  mutate(multiple_matchesUnresolved = n() > 1) # TRUE if multiple match could not be resolved

#### Concatenate the multiple results for these taxa and separate with pipe delimiter '|'

WCVPMatch<-matches %>%
  group_by(verbatimScientificName) %>%
  summarise(across(everything(), ~paste0(unique(.), collapse = "|")))

### Write interim output to file 

WCVPInterim<- WCVPMatch

write.csv(WCVPInterim, "data//interim//WCVPInterim.csv", row.names = FALSE)

# Step 12: Assess and flag issues for plants matched to WCVP

## Assess issues

head(WCVPMatch) # Look at the first 5 rows of the data

### Look at how many have each status (e.g., accepted, synonym etc)

WCVPMatch %>%
  group_by(wcvp_status) %>%
  count()

### Look at how many had a specific match type (e.g., exact, fuzzy etc)

WCVPMatch %>%
  group_by(match_type) %>%
  count()

### Look at how many fall into different ranks (e.g., species, subspecies etc)

WCVPMatch %>%
  group_by(wcvp_rank) %>%
  count()

### Look at how many had multiple matches

WCVPMatch %>%
  group_by(multiple_matches) %>%
  count()

### Look at how many had multiple matches that could not be resolved

WCVPMatch %>%
  group_by(multiple_matchesUnresolved) %>%
  count()

## Flag taxonomic issues

### Flag synonyms

WCVPMatch = WCVPMatch %>% 
  group_by(wcvp_status) %>% 
  mutate(synonym.flag = ifelse(wcvp_status =="Synonym", TRUE, FALSE)) # TRUE if a synonym

### Flag taxon names that were not exact matches

WCVPMatch<-WCVPMatch %>% 
  mutate(uncertain.flag = case_when((match_type =="NA")~FALSE,
                                    (match_type == "Exact")~FALSE,
                                    (match_type =="Exact (without author)") ~FALSE,
                                    TRUE~TRUE)) # TRUE if uncertain

### Flag doubtful taxon names (e.g., invalid, illegitimate, missapplied, and those with unresolved multiple matches)

WCVPMatch<-WCVPMatch %>% 
  mutate(doubtful.flag = case_when((multiple_matchesUnresolved == TRUE)~TRUE,
                                   (wcvp_status == "Accepted")~FALSE,
                                    (wcvp_status =="Synonym")~FALSE,
                                    (wcvp_status == "NA") ~FALSE,
                                    TRUE~TRUE)) # TRUE if doubtful

### Flag unmatched taxa

WCVPMatch = WCVPMatch %>% 
  mutate(unmatched.flag=ifelse(match_type =="NA", TRUE, FALSE)) # TRUE if not matched

### Flag if the result should be checked for a possible error or issue (unmatched, doubtful, and uncertain matches)

WCVPMatch<-WCVPMatch %>% 
  mutate(possibleError.flag = ifelse(uncertain.flag =="TRUE" | doubtful.flag == "TRUE" | unmatched.flag == "TRUE", TRUE, FALSE)) # TRUE if possible error or issue

### Flag sub-specific entities

WCVPMatch = WCVPMatch %>% 
  group_by(wcvp_rank) %>% 
  mutate(subspecific.flag = ifelse(wcvp_rank =="Species"|wcvp_rank == "NA", FALSE, TRUE)) # TRUE if sub-specific entity

### Flag duplicates in returned scientific names (exclude NAs)

WCVPMatch = WCVPMatch %>% 
  group_by(scientificName) %>% 
  mutate(scientificNameDuplicate.flag = ifelse(n() > 1 & scientificName !="NA", 
                                               TRUE , FALSE)) # TRUE if duplicate

### Flag duplicates in verbatim names

WCVPMatch = WCVPMatch %>% 
  group_by(verbatimScientificName) %>% 
  mutate(verbatimNameDuplicate.flag = n() > 1) # TRUE if duplicate

# Step 13: Preparation of output for plants matched to WCVP

### Add column with the source of the scientific name, date accessed, and web-link

WCVPMatch$scientificNameSource<-ifelse(WCVPMatch$unmatched.flag == "TRUE", 
                                       "NA", paste0("POWO", " ", "(", format(Sys.Date(), "%d %B %Y"),")", " ","https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:", WCVPMatch$wcvp_ipni_id))

### Remove unwanted columns

WCVPMatch <- WCVPMatch %>% ungroup() # ungroup columns
WCVPMatch <- WCVPMatch %>% select(!(c(match_type:wcvp_accepted_id, powo_id:multiple_matchesUnresolved))) # only keep the relevant columns
names(WCVPMatch) # names of the columns

# Step 14: Merge

## Prepare the datasets

### Standardise formating

GBIFMatch<-as.data.frame(GBIFMatch)
BODATSAMatch<-as.data.frame(BODATSAMatch)
WCVPMatch<-as.data.frame(WCVPMatch)

GBIFMatch = GBIFMatch %>%
  mutate(across(everything(), as.character))
BODATSAMatch = BODATSAMatch %>%
  mutate(across(everything(), as.character))
WCVPMatch = WCVPMatch %>%
  mutate(across(everything(), as.character))

### Revise flagged duplicate plant names (verbatim and scientific) so that they are across the BODATSA and WCVP matches

PlantMatch<-rbind(BODATSAMatch, WCVPMatch) # Bind the datasets together

#### Flag duplicates in returned scientific names (exclude NAs)

PlantMatch = PlantMatch %>% 
  group_by(scientificName) %>% 
  mutate(scientificNameDuplicate.flag = ifelse(n() > 1 & scientificName !="NA", 
                                                TRUE , FALSE)) # TRUE for duplicates
### Flag duplicates in verbatim names

PlantMatch = PlantMatch %>% 
  group_by(verbatimScientificName) %>% 
  mutate(verbatimNameDuplicate.flag = n() > 1) # TRUE for duplicates

## Bind the datasets 

AllMatch<-rbind(GBIFMatch, PlantMatch) # Bind all the datasets together

### Check if the number of taxa in the output dataset is that same as in the original dataset

length(AllMatch$verbatimScientificName) == length(NamesDat$verbatimScientificName) # should be TRUE

### Reorder dataset, so that the taxa are in the same order as they were in the input file

AllMatch<-AllMatch[match(NamesDat$verbatimScientificName, AllMatch$verbatimScientificName),]

## Write final output to file

write.csv(AllMatch, "data//outputs//Taxon_names_standardised.csv", row.names = FALSE)


