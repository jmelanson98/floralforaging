##### Can we use floral nutrient data from Pol-NIC database?
### April 2, 2026
### J Melanson


### Load in packages
library(taxize)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(rgbif)

bombus_path = "/Users/jenna1/Documents/UBC/bombus_project/"
bombus_path = "/Users/jenna1/Documents/bombus_project/"

###################################################
### Load in survey data
###################################################
specimenData2022 = as.data.frame(read.csv("cleandata/fielddata/2022specimendata.csv"), sep = ",", header = T)
specimenData2023 = as.data.frame(read.csv("cleandata/fielddata/2023specimendata.csv"), sep = ",", header = T)
veg2022 = read.csv("cleandata/fielddata/2022vegetationdata.csv")
veg2023 = read.csv("cleandata/fielddata/2023vegetationdata.csv")
observedplants = read.csv("cleandata/plantdata/observedplants.csv")

###################################################
### Load in nutrient data
###################################################
tissier1 = read.csv("pollendata/tissier1.csv")
tissier2 = read.csv("pollendata/tissier2.csv")
roulston = read.csv("pollendata/roulston.csv")
pamminger = read.csv("pollendata/pamminger2019.csv")
# just use Maurizio, Odeoux, Requier, Szczsna, Tasier
other = read.csv("pollendata/other.csv")


###################################################
### Get full taxonomic classification of plants
###################################################
observedplants = observedplants[c("scientific_name", "common_name", "plant_code")]

# Get canonical names
resolved = gna_verifier(observedplants$scientific_name)

resolved_df = tibble(
  input_plants     = observedplants$scientific_name,
  canonicalname    = resolved$matchedCanonicalFull)
resolved_df = resolved_df[!is.na(resolved_df$canonicalname), ]

# Get full classification from GBIF
# function to get taxonomy

get_classification = function(input_name, canonical_name) {
  tryCatch({
    
    match = name_backbone(name = canonical_name)
    
    # handle NONE match by falling back to name_suggest
    if (match$matchType == "NONE") {
      print(paste0(canonical_name, " not found. Trying name lookup."))
      suggestion = name_lookup(q = canonical_name, 
                               higherTaxonKey = 6, limit = 1)$data
      if (is.null(suggestion) || nrow(suggestion) == 0) stop("No match")
      key = suggestion$key[1]
      rank_out = "GENUS"
      canonical_out = canonical_name  # just use what we passed in
      rank_out     = suggestion$rank[1]
      status_out   = NA_character_
      authority_out = NA_character_
    } else {
      key           = match$usageKey
      canonical_out = match$canonicalName
      rank_out      = match$rank
      status_out    = match$status
      authority_out = match$authorship
    }
    
    if (is.null(key) || is.na(key)) stop("No usable key")
    
    cls = classification(key, db = "gbif")[[1]]
    if (is.null(cls) || nrow(cls) == 0) stop("No classification returned")
    
    ranks = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
    tax = setNames(rep(NA_character_, length(ranks)), ranks)
    for (r in ranks) {
      val = cls$name[cls$rank == r]
      if (length(val) > 0) tax[r] = val
    }
    
    tibble(
      input_plants  = input_name,
      canonicalname = canonical_out,
      taxa_rank     = rank_out,
      taxa_status   = status_out,
      synonym_name  = ifelse(!is.na(status_out) & status_out == "SYNONYM", 
                             canonical_name, NA_character_),
      authority     = authority_out,
      kingdom       = tax["kingdom"],
      phylum        = tax["phylum"],
      class         = tax["class"],
      order         = tax["order"],
      family        = tax["family"],
      genus         = tax["genus"],
      species       = tax["species"]
    )
    
  }, error = function(e) {
    message("Failed for: ", canonical_name, " — ", e$message)
    tibble(
      input_plants  = input_name,
      canonicalname = NA_character_,
      taxa_rank     = NA_character_,
      taxa_status   = NA_character_,
      synonym_name  = NA_character_,
      authority     = NA_character_,
      kingdom       = NA_character_,
      phylum        = NA_character_,
      class         = NA_character_,
      order         = NA_character_,
      family        = NA_character_,
      genus         = NA_character_,
      species       = NA_character_
    )
  })
}


# map to input plants
taxonomy_df = bind_rows(
  Map(get_classification, resolved_df$input_plants, resolved_df$canonicalname))

# add tayberry because it's a weirdo
tayberry_row = tibble(
  input_plants  = "Rubus fruticosus x Rubus idaeus",
  canonicalname = "Rubus fruticosus x Rubus idaeus",  # hybrid cross, not really recognized
  taxa_rank     = "UNRANKED",
  taxa_status   = "UNRESOLVED",
  synonym_name  = NA,
  kingdom       = "Plantae",
  phylum        = "Tracheophyta",
  class         = "Magnoliopsida",
  order         = "Rosales",
  family        = "Rosaceae",
  genus         = "Rubus",
  species       = NA)   # no species epithet --- hybrid
taxonomy_df = rbind(taxonomy_df, tayberry_row)
taxonomy_df = distinct(taxonomy_df)

# join back to original data
observedplants = observedplants %>%
  left_join(taxonomy_df, by = c("scientific_name" = "input_plants"))


###################################################
### Clean and standardize protein database
###################################################

# tissier
tissier1_prot = tissier1[c("Plants_Latin_Name", "Proteins_Total_g100g")]
colnames(tissier1_prot) = c("scientific_name", "protein_percent")
tissier1_prot$scientific_name[tissier1_prot$scientific_name == "Vaccinum spp "] = "Vaccinium spp"
tissier1_prot$scientific_name[tissier1_prot$scientific_name == "Curcubita pepo Thunb"] = "Cucurbita pepo"
tissier1_prot = filter(tissier1_prot, !is.na(scientific_name) & scientific_name != "")
tissier1_prot = tissier1_prot[str_detect(tissier1_prot$scientific_name, "%", negate = T),]
tissier1_prot = tissier1_prot[str_detect(tissier1_prot$scientific_name, "Mix", negate = T),]
tissier1_resolved = gna_verifier(tissier1_prot$scientific_name)
t1_resolved_df = tibble(
  input_plants     = tissier1_prot$scientific_name,
  canonicalname    = tissier1_resolved$matchedCanonicalFull)
t1_resolved_df = t1_resolved_df[!is.na(t1_resolved_df$canonicalname), ]
tissier1_taxonomy = bind_rows(
  Map(get_classification, t1_resolved_df$input_plants, t1_resolved_df$canonicalname))
t1keep = tissier1_taxonomy[c("input_plants", "family")]
t1keep = distinct(t1keep)
tissier1_prot = left_join(tissier1_prot, t1keep, by = c("scientific_name" = "input_plants"))

tissier2_prot = tissier2[c("Pollen_species", "Proteins_g100g")]
colnames(tissier2_prot) = c("scientific_name", "protein_percent")
tissier2_resolved = gna_verifier(tissier2_prot$scientific_name)
t2_resolved_df = tibble(
  input_plants     = tissier2_prot$scientific_name,
  canonicalname    = tissier2_resolved$matchedCanonicalFull)
t2_resolved_df = t2_resolved_df[!is.na(t2_resolved_df$canonicalname), ]
tissier2_taxonomy = bind_rows(
  Map(get_classification, t2_resolved_df$input_plants, t2_resolved_df$canonicalname))
t2keep = tissier2_taxonomy[c("input_plants", "family")]
t2keep = distinct(t2keep)
tissier2_prot = left_join(tissier2_prot, t2keep, by = c("scientific_name" = "input_plants"))

# roulston
roulston_prot = roulston[c("Species", "Protein....")]
colnames(roulston_prot) = c("scientific_name", "protein_percent")
roulston_prot$scientific_name[roulston_prot$scientific_name == "Anemopaegna sp."] = "Anemopaegma sp."
roulston_resolved = gna_verifier(roulston_prot$scientific_name)
r_resolved_df = tibble(
  input_plants     = roulston_prot$scientific_name,
  canonicalname    = roulston_resolved$matchedCanonicalFull)
r_resolved_df = r_resolved_df[!is.na(r_resolved_df$canonicalname), ]
roulston_taxonomy = bind_rows(
  Map(get_classification, r_resolved_df$input_plants, r_resolved_df$canonicalname))
rkeep = roulston_taxonomy[c("input_plants", "family")]
rkeep = distinct(rkeep)
roulston_prot = left_join(roulston_prot, rkeep, by = c("scientific_name" = "input_plants"))


# pamminger
pammingerfilt = pamminger %>% filter(str_detect(citation, "Maurizio") |
                                   str_detect(citation, "Odeoux") |
                                 str_detect(citation, "Requier") |
                                 str_detect(citation, "Szcz") |
                                 str_detect(citation, "Tasier"))
pamminger_prot = pammingerfilt[c("sci_name", "crude_protein...")]
colnames(pamminger_prot) = c("scientific_name", "protein_percent")
pamminger_prot$scientific_name[pamminger_prot$scientific_name == "Majoranum sp."] = "Origanum majorana"
pamminger_prot$scientific_name[pamminger_prot$scientific_name == "Helianthus anuus"] = "Helianthus annus"
pamminger_prot$scientific_name[pamminger_prot$scientific_name == "Vitacea spp."] = "Vitaceae spp."
pamminger_resolved = gna_verifier(pamminger_prot$scientific_name)
p_resolved_df = tibble(
  input_plants     = pamminger_prot$scientific_name,
  canonicalname    = pamminger_resolved$matchedCanonicalFull)
p_resolved_df = p_resolved_df[!is.na(p_resolved_df$canonicalname), ]
pamminger_taxonomy = bind_rows(
  Map(get_classification, p_resolved_df$input_plants, p_resolved_df$canonicalname))
pkeep = pamminger_taxonomy[c("input_plants", "family")]
pkeep = distinct(pkeep)
pamminger_prot = left_join(pamminger_prot, pkeep, by = c("scientific_name" = "input_plants"))


# other
other_prot = other[c("scientific_name", "protein_percentage")]
colnames(other_prot) = c("scientific_name", "protein_percent")
other_prot$scientific_name = paste(stringr::word(other_prot$scientific_name, 1), stringr::word(other_prot$scientific_name, 2), sep = " ")
other_resolved = gna_verifier(other_prot$scientific_name)
o_resolved_df = tibble(
  input_plants     = other_prot$scientific_name,
  canonicalname    = other_resolved$matchedCanonicalFull)
o_resolved_df = o_resolved_df[!is.na(o_resolved_df$canonicalname), ]
other_taxonomy = bind_rows(
  Map(get_classification, o_resolved_df$input_plants, o_resolved_df$canonicalname))
okeep = other_taxonomy[c("input_plants", "family")]
okeep = distinct(okeep)
other_prot = left_join(other_prot, okeep, by = c("scientific_name" = "input_plants"))


###################################################
### Get average protein content per observation
###################################################
tissier1_avgprot = tissier1_prot %>%
  mutate(avg_protein = sapply(protein_percent, function(x) {
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

tissier2_avgprot = tissier2_prot[tissier2_prot$protein_percent != "NA**",]
tissier2_avgprot$avg_protein = tissier2_avgprot$protein_percent

tissier_data = rbind(tissier1_avgprot, tissier2_avgprot)
tissier_data = tissier_data[!is.na(tissier_data$avg_protein),]
tissier_data$avg_protein = as.numeric(tissier_data$avg_protein)
tissier_data$dataset = "tissier"
tissier_data$genus = stringr::word(tissier_data$scientific_name, 1)

roulston_prot$avg_protein = roulston_prot$protein_percent
roulston_prot$dataset = "roulston"
roulston_prot$genus = stringr::word(roulston_prot$scientific_name, 1)

pamminger_prot$avg_protein = pamminger_prot$protein_percent
pamminger_prot$dataset = "pamminger"
pamminger_prot$genus = stringr::word(pamminger_prot$scientific_name, 1)

other_prot$avg_protein = stringr::word(other_prot$protein_percent, 1)
other_prot$dataset = "other"
other_prot$genus = stringr::word(other_prot$scientific_name, 1)

alldata = rbind(tissier_data, roulston_prot, other_prot, pamminger_prot)
alldata$avg_protein = as.numeric(alldata$avg_protein)

protbyfamily = ggplot(alldata, aes(y = reorder(family, avg_protein, mean), x = avg_protein, colour = dataset)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
  scale_colour_manual(values = c("coral", "darkgreen", "navy", "orange")) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
  ylab("Family") +
  xlab("Percent protein") +
  theme_bw()
ggsave("figures/protein_by_family.png", protbyfamily, height = 4000, width = 3000, units = "px")



###################################################
### Check coverage of both datasets at family
###################################################
tissier_family = tissier_data %>% 
  group_by(family) %>%
  summarize(prot = mean(avg_protein))
roulston_family = roulston_prot %>% 
  group_by(family) %>%
  summarize(prot = mean(avg_protein))
all_family = alldata %>% 
  group_by(family) %>%
  summarize(prot = mean(avg_protein))

compiled_tissier = left_join(observedplants, tissier_family, by = "family")
# like 3/4 coverage, taxonomically speaking

compiled_roulston = left_join(observedplants, roulston_family, by = "family")
# also about 3/4 coverage

compiled_all = left_join(observedplants, all_family, by = "family")
# 83% coverage

# what if I remove unvisited plants?
beeflowers = unique(c(specimenData2022$active_flower, specimenData2023$active_flower))
beeflowers = beeflowers[str_detect(beeflowers, "nest|N/A|flying|dipu|road|dead|grass|leaf|ground", negate = T)]

visitedplants = observedplants[observedplants$plant_code %in% beeflowers,]

compiled_tissier = left_join(visitedplants, tissier_family, by = "family")
# 83% coverage

compiled_roulston = left_join(visitedplants, roulston_family, by = "family")
# 85% coverage

compiled_all = left_join(visitedplants, all_family, by = "family")
# 91% coverage



###################################################
### Check visitation percentage per family
###################################################

specdata = left_join(specimenData2022, compiled_all, by = c("active_flower" = "plant_code"))
specdata = specdata[!is.na(specdata$family),]
specdata = specdata[specdata$final_id != "",]

specdata23 = left_join(specimenData2023, compiled_all, by = c("active_flower" = "plant_code"))
specdata23 = specdata23[!is.na(specdata23$family),]
specdata23 = specdata23[str_detect(specdata23$notes, "queen", negate = T),]
specdata23 = specdata23[specdata23$final_id != "",]


bubble_data = specdata %>%
  count(final_id, family, prot) %>%
  group_by(final_id) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

bubble_data23 = specdata23 %>%
  count(final_id, family, prot) %>%
  group_by(final_id) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()


ggplot(bubble_data, aes(x = final_id, y = family, size = prop, colour = prot)) +
  geom_point() +
  scale_size_continuous(range = c(1, 10)) +
  labs(x = "Bee species", y = "Plant family", size = "Frequency") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(bubble_data23, aes(x = final_id, y = family, size = prop, colour = prot)) +
  geom_point() +
  scale_size_continuous(range = c(1, 10)) +
  labs(x = "Bee species", y = "Plant family", size = "Frequency") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
