##### Can we use floral nutrient data from Pol-NIC database?
### April 2, 2026
### J Melanson


### Load in packages
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)

bombus_path = "/Users/jenna1/Documents/UBC/bombus_project/"

### Load in floral survey data
specimenData2022 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2022specimendata.csv"), sep = ",", header = T))
specimenData2023 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2023specimendata.csv"), sep = ",", header = T))
veg2022 = read.csv(paste0(bombus_path, "raw_data/2022vegetationdata.csv"))
veg2023 = read.csv(paste0(bombus_path, "raw_data/2023vegetationdata.csv"))
plantlist = read.csv(paste0(bombus_path, "raw_data/plant_list.csv"), header = FALSE)
colnames(plantlist) = c("Plants_Latin_Name", "common_name", "plantcode", "notes")

### Load in nutrient data
nutrients = read.csv("pollendata/floralnutrients.csv")
roulston = read_excel("pollendata/Roulston_2000_Pollen_Protein_Data.xlsx")

### Pull out just proteins and lipids
PL = nutrients[c("Plants_Latin_Name", "Plants_Family", "Proteins_Total_g100g", "Lipids_Total_g100g")]
PL$Plants_Latin_Name[PL$Plants_Latin_Name == "Vaccinum spp "] = "Vaccinium spp"
PL$genus = stringr::word(PL$Plants_Latin_Name, 1)
PL$species = stringr::word(PL$Plants_Latin_Name, 2)

# Get average protein and lipid content per species
PL = PL %>%
  mutate(avg_protein = sapply(Proteins_Total_g100g, function(x) {
    x <- trimws(gsub("\\*", "", x))
    if (grepl("±", x)) {
      # Take the first value (mean)
      as.numeric(trimws(strsplit(x, "±")[[1]][1]))
    } else if (grepl("-", x)) {
      # Average the two values
      parts <- as.numeric(trimws(strsplit(x, "-")[[1]]))
      mean(parts)
    } else {
      # Single value
      as.numeric(trimws(x))
    }
  }))

PL = PL %>%
  mutate(avg_lipid = sapply(Lipids_Total_g100g, function(x) {
    x <- trimws(gsub("\\*", "", x))
    if (grepl("±", x)) {
      # Take the first value (mean)
      as.numeric(trimws(strsplit(x, "±")[[1]][1]))
    } else if (grepl("-", x)) {
      # Average the two values
      parts <- as.numeric(trimws(strsplit(x, "-")[[1]]))
      mean(parts)
    } else {
      # Single value
      as.numeric(trimws(x))
    }
  }))

# Clean plant list
plantlist$Plants_Latin_Name[plantlist$Plants_Latin_Name == "Aescules x carnea"] = "Aescules carnea"
plantlist$Plants_Latin_Name[plantlist$Plants_Latin_Name == "Convolvulus vulgaris?"] = "Convolvulus vulgaris"
plantlist$Plants_Latin_Name[plantlist$Plants_Latin_Name == "Crocosmia × crocosmiiflora"] = "Crocosmia crocosmiiflora"
plantlist$Plants_Latin_Name[plantlist$Plants_Latin_Name == "Geranium x oxonianum"] = "Geranium oxonianum"
plantlist$Plants_Latin_Name[plantlist$Plants_Latin_Name == "Glandularia x hybrida"] = "Glandularia hybrida"
plantlist$Plants_Latin_Name[plantlist$Plants_Latin_Name == "Hemerocallis x hybridis"] = "Hemerocallis hybridis"
plantlist$Plants_Latin_Name[plantlist$Plants_Latin_Name == "Iris x germanica"] = "Iris germanica"
plantlist$Plants_Latin_Name[plantlist$Plants_Latin_Name == "Lavandula x intermedia"] = "Lavandula intermedia"
plantlist$Plants_Latin_Name[plantlist$Plants_Latin_Name == "Leucanthemum x superbum"] = "Leucanthemum superbum"
plantlist$Plants_Latin_Name[plantlist$Plants_Latin_Name == "Magnolia × soulangeana"] = "Magnolia soulangeana"
plantlist$Plants_Latin_Name[plantlist$Plants_Latin_Name == "Mentha × piperita"] = "Mentha piperita"
plantlist$Plants_Latin_Name[plantlist$Plants_Latin_Name == "Petunia x atkinsiana"] = "Petunia atkinsiana"
plantlist$Plants_Latin_Name[plantlist$Plants_Latin_Name == "Pelargonium x hortorum"] = "Pelargonium hortorum"
plantlist$Plants_Latin_Name[plantlist$Plants_Latin_Name == "Pelargonium x hybridum"] = "Pelargonium hybridum"
plantlist$Plants_Latin_Name[plantlist$Plants_Latin_Name == "Photinia x fraseri"] = "Photinia fraseri"
plantlist$Plants_Latin_Name[plantlist$Plants_Latin_Name == "Rosa × damascena"] = "Rosa damascena"
plantlist$Plants_Latin_Name[plantlist$Plants_Latin_Name == "Aescules carnea"] = "Aesculus carnea"
plantlist$Plants_Latin_Name[plantlist$Plants_Latin_Name == "Calibrchoa hybrida"] = "Calibrachoa hybrida"

beeflowers = unique(c(specimenData2022$active_flower, specimenData2023$active_flower))
beeflowers = beeflowers[str_detect(beeflowers, "nest|N/A|flying|dipu|road|dead|grass|leaf|ground", negate = T)]

plantlist = plantlist[!plantlist$plantcode %in% c("dead", "flying", "Monocotyledon", "ALMA"),]

plantlist$genus = stringr::word(plantlist$Plants_Latin_Name, 1)
plantlist$species = stringr::word(plantlist$Plants_Latin_Name, 2)

# Try a random join?
compiled = left_join(plantlist, PL, by = c("genus", "species"))
specieslevel = compiled[!is.na(compiled$Proteins_Total_g100g),]

# Now average by genus
plantlist_remaining = plantlist[!plantlist$plantcode %in% specieslevel$plantcode,]
genus_average = PL %>%
  group_by(genus) %>%
  summarize(mean_protein = mean(avg_protein, na.rm = TRUE),
            mean_lipid = mean(avg_lipid, na.rm = TRUE))
genuscompiled = left_join(plantlist_remaining, genus_average, by = "genus")
genuslevel = genuscompiled[!is.na(genuscompiled$mean_protein),]


# Now what's left?
plantlist_remaining = plantlist[(!plantlist$plantcode %in% specieslevel$plantcode) & (!plantlist$plantcode %in% genuslevel$plantcode),]


# Try with Roulston data?

roulston$genus = stringr::word(roulston$Species, 1)
roulston$species = stringr::word(roulston$Species, 2)


roulston_joined = left_join(plantlist, roulston, by = c("Plants_Latin_Name" = "Species"))
roulston_spp = roulston_joined[!is.na(roulston_joined$`Protein (%)`),]
roulston_remaining = plantlist[!plantlist$Plants_Latin_Name %in% roulston_spp$Plants_Latin_Name,]

roulston_genus = roulston %>% 
  group_by(genus) %>%
  summarize(avg_protein = mean(`Protein (%)`))
roulston_joined_genus = left_join(roulston_remaining, roulston_genus, by = "genus")
