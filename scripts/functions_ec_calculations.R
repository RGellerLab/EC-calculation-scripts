## Ron Geller, April 2nd, 2025
collapse_data <- function(condition.file, x) {
  # Load required libraries
  library(tidyverse)
  library(readxl)
  
  # ---- Helper to process sheet name or index ----
  process_sheet_name_or_index <- function(x) {
    if (!is.na(as.numeric(x))) {
      as.numeric(x)  # Convert numeric-like strings to numeric
    } else {
      as.character(x)  # Keep as character if not numeric
    }
  }
  
  # ---- Extract parameters ----
  template <- condition.file$template.file[x] 
  page <- process_sheet_name_or_index(condition.file$template.page[x])  # Safe page handling
  
  # Data information
  datafile <- paste0(condition.file$data.path[x], "/", condition.file$datafile.name[x])
  data.page <- process_sheet_name_or_index(condition.file$data.page[x])  # Safe page handling
  data.format <- condition.file$data.format[x]
  data.range <- condition.file$data.location[x]
  
  # Output directories
  out.dir <- condition.file$out.dir[x]
  out.dir.raw.data <- paste0(out.dir, "/raw_data/")
  out.dir.collapsed.data <- paste0(out.dir, "/collapsed_data/")
  
  # Controls and misc.
  no.virus.control <- condition.file$no.virus.control[x]
  no.treatment.control <- condition.file$no.treatment.control[x]
  subtract.background <- as.logical(condition.file$subtract.background[x])  # Ensure logical
  
  # ---- Error checks ----
  if (!file.exists(template)) stop(paste("Error: Template file not found:", template))
  if (!file.exists(datafile)) stop(paste("Error: Data file not found:", datafile))
  if (!(data.format %in% c("one.line", "plate.layout"))) stop("Error: 'data.format' must be 'one.line' or 'plate.layout'.")
  if (is.na(data.page)) stop("Error: 'data.page' is missing or improperly formatted.")
  if (is.na(data.range) || data.range == "") stop("Error: 'data.range' is missing or empty.")
  if (is.na(subtract.background)) stop("Error: 'subtract.background' must be TRUE or FALSE (check formatting).")
  
  # ---- Function to read and format template ranges ----
  read_template_range <- function(range, value_col_name) {
    read_xlsx(template, sheet = page, range = range, col_types = rep("text", 13)) %>%
      pivot_longer(cols = 2:13, values_to = value_col_name, names_to = "column") %>%
      rename(row = 1) %>%
      mutate(row = toupper(row), column = as.numeric(column)) %>%
      unite("id", row, column, sep = "")
  }
  
  # ---- Read template sections ----
  treatment.name <- read_template_range("A2:M10", "treatment")
  concentrations <- read_template_range("A13:M21", "concentration")
  replicate <- read_template_range("A24:M32", "replicate")
  virus <- read_template_range("A35:M43", "virus")
  
  # ---- Read experimental data ----
  if (data.format == "one.line") {
    data <- read_xlsx(datafile, sheet = data.page, range = data.range) %>%
      mutate(across(everything(), as.character)) %>%
      pivot_longer(cols = everything(), values_to = "fluo", names_to = "id")
  } else if (data.format == "plate.layout") {
    data <- read_xlsx(datafile, sheet = data.page, range = data.range, col_types = rep("text", 13)) %>%
      pivot_longer(cols = 2:13, values_to = "fluo", names_to = "column") %>%
      rename(row = 1) %>%
      mutate(row = toupper(row), column = as.numeric(column)) %>%
      unite("id", row, column, sep = "")
  }
  
  # ---- Merge data ----
  df <- virus %>%
    left_join(treatment.name, by = "id") %>%
    left_join(concentrations, by = "id") %>%
    left_join(replicate, by = "id") %>%
    left_join(data, by = "id") %>%
    mutate(fluo = as.numeric(fluo)) %>%
    filter(treatment != "NA", fluo != "NA", virus != "NA", concentration != "NA", replicate != "NA") %>% 
    filter(!is.na(treatment), !is.na(fluo), !is.na(virus), !is.na(concentration), !is.na(replicate)) %>%  # also remove if NA is a strin in excel
    mutate(assay = condition.file$assay[x], file.name = condition.file$out.name[x])
  
  # ---- Save raw merged data ----
  dir.create(out.dir.raw.data, recursive = TRUE, showWarnings = FALSE)
  write_csv(df, paste0(out.dir.raw.data, condition.file$out.name[x], ".csv"))
  
  # ---- Background subtraction if needed ----
  if (subtract.background) {
    if (!(no.virus.control %in% df$treatment)) {
      stop(paste("Error: 'no.virus.control' (", no.virus.control, ") not found in data."))
    }
    df <- df %>%
      group_by(virus) %>%
      mutate(fluo = fluo - mean(fluo[treatment == no.virus.control], na.rm = TRUE)) %>%
      ungroup()
  }
  
  # ---- Normalize to no-treatment control ----
  if (!(no.treatment.control %in% df$treatment)) {
    stop(paste("Error: 'no.treatment.control' (", no.treatment.control, ") not found in data."))
  }
  df <- df %>%
    group_by(virus) %>%
    mutate(remaining_signal = fluo / mean(fluo[treatment == no.treatment.control], na.rm = TRUE)) %>%
    ungroup()
  
  # ---- Remove controls ----
  df <- df %>% filter(!treatment %in% c(no.virus.control))
  
  # ---- Save final collapsed data ----
  dir.create(out.dir.collapsed.data, recursive = TRUE, showWarnings = FALSE)
  write_csv(df, paste0(out.dir.collapsed.data, condition.file$out.name[x], ".csv"))
}


## Ron Geller, April 2nd, 2025
get_ec <- function(df,
                   ec.level = 50,
                   reciprocal.ec = F,
                   out.dir = "../results",
                   out.name = "res",
                   minimum.data.value = 1e-5,
                   call.ec.lower.limit = 0.4,
                   call.ec.upper.limit = 0.5,
                   xlabel,
                   model = "LL.3()",
                   plot = TRUE,
                   nrow = 4,
                   ncol = 4,
                   title_size = 10,   # Added parameter for title size
                   annotation_size = 3) {  # Added parameter for annotation size
  
  # Load required libraries
  library(tidyverse)
  library(drc)
  library(gridExtra)
  library(scales)
  
  # ---- Prepare output directories ----
  out.dir.ecdata <- paste0(out.dir, "/ec_data/")
  out.dir.graphs <- paste0(out.dir, "/ec_graphs/")
  dir.create(out.dir.ecdata, recursive = TRUE, showWarnings = FALSE)
  dir.create(out.dir.graphs, recursive = TRUE, showWarnings = FALSE)
  
  # ---- Data Preparation ----
  df$remaining_signal[df$remaining_signal < minimum.data.value] <- minimum.data.value
  df <- df %>% filter(!is.na(remaining_signal))
  
  # Remove rows where key columns are "NA" as string
  df <- df %>%
    filter(treatment != "NA", fluo != "NA", virus != "NA", concentration != "NA", replicate != "NA") %>%
    filter(!is.na(treatment), !is.na(fluo), !is.na(virus), !is.na(concentration), !is.na(replicate))
  
  # Convert concentration to numeric if needed
  df$concentration <- as.numeric(df$concentration)
  
  # Reciprocal transformation if requested
  if (reciprocal.ec) {
    if (any(df$concentration == 0, na.rm = TRUE)) stop("Error: Cannot compute reciprocal for zero concentrations.")
    df$concentration <- 1 / df$concentration
  }
  
  # Remove zero/negative concentrations for fitting
  df <- df %>% filter(concentration > 0)
  
  # Create unique ID
  df <- df %>% unite("id", virus, treatment, assay, remove = FALSE)
  
  # ---- Initialize results ----
  model_name <- gsub("\\(\\)", "", model)
  res.ec <- tibble(
    id = unique(df$id),
    reach = NA_character_,
    ec = NA_real_,
    ec_error = NA_real_,
    model = model_name
  )
  
  plot.list <- list()
  mdl <- eval(parse(text = model))
  
  # ---- Loop through each unique ID ----
  for (i in unique(res.ec$id)) {
    temp.df <- df %>% filter(id == i)
    
    # Dynamically extract assay for title and file name
    current_assay <- unique(temp.df$assay)
    
    # Check for fitting conditions
    signal_min <- min(temp.df$remaining_signal, na.rm = TRUE)
    signal_max <- max(temp.df$remaining_signal, na.rm = TRUE)
    annotation_text <- "EC not determined"
    
    # --- Fit model if appropriate ---
    if (nrow(temp.df) >= 1 && signal_min < call.ec.lower.limit && signal_max >= call.ec.upper.limit) {
      modl <- tryCatch({
        drm(data = temp.df, remaining_signal ~ concentration, fct = mdl)
      }, error = function(e) NULL)
      
      if (!is.null(modl)) {
        eds <- tryCatch({
          ED(modl, c(ec.level))
        }, error = function(e) NULL)
        
        if (!is.null(eds)) {
          ec_value <- eds[1, "Estimate"]
          ec_se <- eds[1, "Std. Error"]
          annotation_text <- paste0(formatC(ec_value, format = "e", digits = 1), " ± ", formatC(ec_se, format = "e", digits = 1))
          res.ec[res.ec$id == i, ] <- list(i, "yes", ec_value, ec_se, model_name)
          
          # Attach fitted values for plotting
          temp.df$fitted <- fitted(modl)
          
          # Prepare label dataframe for annotation
          label.df <- tibble(
            label = annotation_text,
            predictor = i  # Helps in distinguishing conditions if needed
          )
          
          # format title
          title.lab=unlist(strsplit(i, split = "_"))
          title.lab=paste(title.lab[1], title.lab[2], "\n", 
                          title.lab[3], model_name,sep=" ")
          
          
          # Updated ggplot call
          plot.list[[i]] <- ggplot(temp.df, aes(x = concentration, y = remaining_signal)) +
            geom_point() +
            geom_line(aes(y = fitted), color = "blue") +
            geom_hline(yintercept = ec.level/100, linetype = "dashed", color = "red") +
            scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
            scale_y_continuous(limits = c(0, 2)) +
            ggtitle(title.lab) +
            ylab("Relative vs. mock") +
            xlab(xlabel) +
            geom_text(data = label.df, aes(x = Inf, y = Inf, label = label),
                      hjust = "right", vjust = "top", size = annotation_size) +  # Apply annotation size
            theme(plot.title = element_text(size = title_size))  # Apply title size
          
          next
        }
      }
    }
    
    # --- EC not determined: fallback plot and reason ---
    reason <- if (signal_min >= call.ec.lower.limit) {
      paste0("Weak effect (>", call.ec.lower.limit, ")")
    } else if (signal_max < call.ec.upper.limit) {
      paste0("Too reduced (<", call.ec.upper.limit, ")")
    } else {
      "Model fit failed"
    }
    
    res.ec[res.ec$id == i, ] <- list(i, reason, NA_real_, NA_real_, model_name)
    
    # format label
    label.df <- tibble(label = reason, x = Inf, y = Inf)
    
    # format title
    title.lab=unlist(strsplit(i, split = "_"))
    title.lab=paste(title.lab[1], title.lab[2], "\n", 
                    title.lab[3], model_name,sep = " ")
    plot.list[[i]] <- ggplot(temp.df, aes(x = concentration, y = remaining_signal)) +
      geom_point() +
      scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
      scale_y_continuous(limits = c(0, 2)) +
      ggtitle(title.lab) +
      geom_hline(yintercept = ec.level/100, linetype = "dashed", color = "red") +
      ylab("Relative vs. mock") +
      xlab(xlabel) +
      geom_text(data = label.df, aes(x = Inf, y = Inf, label = label),
                hjust = "right", vjust = "top", size = annotation_size) +
      theme(plot.title = element_text(size = title_size))  # Apply title size
  }
  
  # ---- Save plots if requested ----
  if (plot) {
    plots_per_page <- nrow * ncol
    
    # Helper to split list
    split_list <- function(lst, size) {
      split(lst, ceiling(seq_along(lst) / size))
    }
    
    plot_pages <- split_list(plot.list, plots_per_page)
    
    pdf(file = paste0(out.dir.graphs, out.name, "_EC", ec.level, "_", model_name, ".pdf"), onefile = TRUE)
    for (page in plot_pages) {
      grid.arrange(grobs = page, nrow = nrow, ncol = ncol)
    }
    dev.off()
  }
  
  # ---- Save EC results ----
  res.ec <- res.ec %>%
    separate(id, into = c("virus", "treatment", "assay"), sep = "_") %>%
    mutate(assay = assay) %>%
    rename(
      !!paste0("reach_EC", ec.level) := reach, # rename reach column
      !!paste0("EC", ec.level) := ec,  # Rename EC column
      !!paste0("EC", ec.level, "_error") := ec_error  # Rename EC error column
    )
  
  write_csv(res.ec, paste0(out.dir.ecdata, out.name, "_EC", ec.level,"_",model_name, ".csv"))
}
