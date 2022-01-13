# Motivation: I want to have a quick fix so that load_kegg_KAAS.pl runs properly. Basically it assumes that http://rest.kegg.jp/list/genome provides genome and taxid information. But this is not the case anymore. Therefore, I will download the taxid info fom NCBI and merge to genome info provided by KEGG. My fix is appropriate for most animals, but may cause microbial strains to receive the taxid of the species level. This file will be then provided to the load_kegg_KAAS.pl instead of the version from KEGG, which is not compatible anymore.

# Libraries

library(tidyverse)
r# Set output directory
dir_out <- getwd()

# Download NCBI names file file
system("wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz")
system("tar -zxvf taxdump.tar.gz names.dmp")


# Download genomes files

#gn:T00014	pho; Pyrococcus horikoshii OT3
#gn:T00015	mtu; Mycobacterium tuberculosis H37Rv, laboratory strain
system("wget -O genomes.txt http://rest.kegg.jp/list/genome")


# Read files in R
# I need to combine and get them in the final format or similar
#genome:T00006	mpn, MYCPN, 272634; Mycoplasma pneumoniae M129

# tax ids
taxid <- read.delim("names.dmp", 
                    sep = "\t",
                    stringsAsFactors = F,
                    header = F) %>% 
  select(V1, V3) %>% 
  unique()
colnames(taxid) <- c("taxid", "species")
#taxid    species
#1        all
#1       root
#2   Bacteria
#2   bacteria


#Genomes
genomes <- read.delim("genomes.txt",
                      sep = "\t",
                      stringsAsFactors = F,
                      header = F)
colnames(genomes) <- c("genome", "species")

#head
#genome                                          species
#gn:T00001 hin; Haemophilus influenzae Rd KW20 (serotype d)

# Because the names are not necessarily compatible between the two tables, I will have to combine them in stages

# Use lower case for joining the tables

genomes <- genomes %>% 
  mutate(genome = str_replace(genome, "gn:", "genome:")) %>% 
  separate(species, c("code", "species"), sep = "; ") %>% 
  mutate(species_mod = str_to_lower(species)) %>% 
  rename(species_genome = species)
taxid <- taxid %>% 
  mutate(species_mod = str_to_lower(species)) %>% 
  rename(species_taxid = species)

# Get only that is matching according to full name
annotated <- inner_join(genomes, taxid)


# Get species that match depending on the names without parenthesis info

rest <- anti_join(genomes, annotated, by = "genome") 

annotated <- rest %>% 
  mutate(species_mod = species_mod %>% str_remove_all(" \\(.*?\\)")) %>%
  inner_join(.,taxid) %>%  
  bind_rows(annotated, .)

rest <- anti_join(genomes, annotated, by = "genome") 

# Everything else will be annotated to the Genus + species level. This may move one level up for a lot of bacterial strains, but my main interest here is animals. So, I will move along.

annotated <- rest %>% 
  separate(species_mod, c("Genus", "Species"), sep = " ") %>%
  mutate(species_mod = paste(Genus, Species)) %>% 
  inner_join(.,taxid) %>%  
  bind_rows(annotated, .)

rest <- anti_join(genomes, annotated, by = "genome")

# Combine

annotated <- bind_rows(annotated, rest)

if (!identical(annotated %>% pull(genome) %>% sort,
               genomes %>% pull(genome) %>% sort)){
  print("Annotated data do not match with KEGG genomes")
}

# Format to save
#genome:T00006	mpn, MYCPN, 272634; Mycoplasma pneumoniae M129

annotated_format <- 
  annotated %>% 
  mutate(codes = paste0(code, ", ", code, ", ", taxid)) %>% 
  mutate(col2 = paste0(codes, "; ", species_genome)) %>% 
  select(genome, col2) %>% 
  arrange(genome)

write.table(annotated_format,
            "list.genomes.txt", 
            quote = F, 
            sep = "\t",
            row.names = F,
            col.names = F)


