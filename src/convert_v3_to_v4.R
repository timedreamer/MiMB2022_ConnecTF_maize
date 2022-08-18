# This script is to convert the GRN I built from v3 to v4 ids.
# The top 1M network was downloaded from *FSU Maize GRN* website.
# Then cut the top 100k rows `head -n 100001 ll_sam.txt > ll_sam_100k.txt`

# Author: Ji Huanng
# Date: 2022-06-06


# 0. Prep -----------------------------------------------------------------

library(here)
library(tidyverse)


# 1. Function -------------------------------------------------------------

## Convert v3 network to v4.

convertv3_v4 <- function(input_v3_file, outputv4_file){
    
    stopifnot(file.exists(here("data", "maize_network", input_v3_file)))
    
    ## Read data 
    
    v3ntwk <- read_tsv(here("data", "maize_network", input_v3_file))
    
    
    v3tov4 <- read_tsv(here("data", "maize.v3TOv4.geneIDhistory.txt"), 
                       col_names = c("v3id", "v4id", "seqChangeOrNot", 
                                     "method", "type"), skip = 1)
    
    v3tov4 <- v3tov4 %>% select(v3id, v4id) %>% distinct()
    
    ## Convert to V4ids 
    
    v4ntwk <- v3ntwk %>% left_join(v3tov4, by = c("regulator" = "v3id")) %>% 
        rename(regulatorV4 = v4id) %>% 
        left_join(v3tov4, by = c("target" = "v3id")) %>% 
        rename(targetV4 = v4id) %>% 
        select(regulatorV4, targetV4, weight) %>% 
        distinct() %>% 
        filter(grepl("^(A|G|Z)", .$regulatorV4)) %>% 
        filter(grepl("^(A|G|Z)", .$targetV4)) %>% 
        arrange(desc(weight)) %>% 
        distinct(regulatorV4, targetV4, .keep_all = TRUE)
    
    ## Save the output
    
    write_tsv(v4ntwk, here("result", outputv4_file))
}


# 3. Conversion -----------------------------------------------------------

## The outputs are in "/result" folder.

convertv3_v4(input_v3_file = "ll_leaf_100k.txt",
             outputv4_file = "maize_leaf_GRN.tsv")

convertv3_v4(input_v3_file = "ll_sam_100k.txt",
             outputv4_file = "maize_SAM_GRN.tsv")

convertv3_v4(input_v3_file = "ll_root_100k.txt",
             outputv4_file = "maize_root_GRN.tsv")

convertv3_v4(input_v3_file = "ll_seed_100k.txt",
             outputv4_file = "maize_seed_GRN.tsv")


