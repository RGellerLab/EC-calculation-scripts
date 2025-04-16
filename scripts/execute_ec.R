
##### Written by Ron Geller. Version 2 , April 2, 2025 ####################
##############################################################################
# These scripts will read your template file where all the information is for
# what is in your 96 well plates and collapse them into long format (Step 1).
# It will then combine the information from all the plate (STEP 2) and obtain
# the desired model fit (STEP 3), producing a table with all the values and 
# a graph showing the fit. IT IS IMPORTANT TO CHECK THE GRAPHS TO SEE THE FIT
# MAKES SENSE. Finally, if you have fit data of >1 EC value or model
# it will collapse it in (STEP 4; optional)
##############################################################################



# Load necessary libraries
library(tidyverse)
library(readxl)

# Source collapsing and EC functions
source("./functions_ec_calculations.R")


############## STEP 1. Collapse into long format #############################


# Read the single conditions file
conditions <- read_xlsx("../conditions_file.xlsx", sheet = 1)

# Collapse data for each row of the conditions file
lapply(1:nrow(conditions), function(x) collapse_data(condition.file = conditions, x))




############## STEP 2. OPTIONAL: Combine collapsed files #####################
############## use if you have >1 plate                  #####################

# Get directory based on conditions file
collapsed_data_dir <- paste0(unique(conditions$out.dir),"/collapsed_data/")

# List all CSV files
csv_files <- list.files(path = collapsed_data_dir, pattern = "\\.csv$", full.names = TRUE)

# Read and combine all CSV files into one data frame
df<- csv_files %>%
  set_names() %>%  # Optional: keep file names as names if you want to track source file
  map_dfr(read_csv, .id = "source_file")  # Add column with source file name

# write combined collapsed data
dir.create("../results/combined_collapsed/",recursive = T)
write_csv(df, "../results/combined_collapsed/combined_collapsed_files.csv")

# clean up
rm(csv_files, collapsed_data_dir,df)



############ STEP 3. Calculate EC values #####################################


### set general parameters ############################
# additional parameters available in function if needed

# label for X-axis of graphs
x.name="Ab"

# do you want the reciprocal EC values (e.g. for neutralization assays)?
reciprocal.ec = T # choose T/F

# how many columns and rows of graphs do you want per page?
n.row=4 
n.col=4 

# Path for combined collapsed data from above
data.file="../results/combined_collapsed/combined_collapsed_files.csv"
######################################################


# Read data and filter 0 concentrations/controls
df=read_csv(data.file) %>% filter(concentration>0) %>% 
  filter(!treatment %in% unique(conditions$no.treatment.control,
                                conditions$no.virus.control))

# get EC and plot if desired
get_ec(df,
       ec.level = 50, # inhibitory concentration to calculate (e.g. 50, 90)
       plot = T, # do you want to graph?
       model = "LL.3()", # normally LL.3(), LL.2()
       xlabel = x.name, # axis label
       reciprocal.ec = reciprocal.ec, # reciprocal?
       nrow = n.row, 
       ncol = n.col 
)

get_ec(df,
       ec.level = 90,
       plot = F,
       model = "LL.3()",
       xlabel = x.name,
       reciprocal.ec = reciprocal.ec,
       nrow = n.row, 
       ncol = n.col 
)


# clean up
rm(n.row,n.col, reciprocal.ec,x.name, collapse_data,conditions,df, data.file, get_ec)


#########  OPTIONAL: STEP 4. COMBINE DATA ######################################
#########  if more than 1 EC value obtained (e.g. EC50, EC90) ##################
#########  and/or if more than 1 model tested (e.g. testing LL.2() and LL.3()) #

# Get list of files
files <- list.files("../results/ec_data/", full.names = TRUE)

# Extract EC levels from filenames
ec.levels <- unique(unlist(regmatches(basename(files), gregexpr("EC[0-9]+", basename(files)))))

# Read and bind files per EC level into a named list
ec_data_list <- set_names(ec.levels) %>%
  map(function(lvl) {
    matched <- files[grepl(lvl, files)]
    map_dfr(matched, read_csv)
  })

# Iteratively full join all data frames on common keys
df_combined <- reduce(ec_data_list, full_join, 
                      by = c("virus", "treatment", "assay", "model")) %>%
  relocate(model, .before = everything())

# write output
write_csv(df_combined,"../results/combined_ec_table.csv")

# clean up
rm(files,ec_data_list,df_combined, ec.levels)


