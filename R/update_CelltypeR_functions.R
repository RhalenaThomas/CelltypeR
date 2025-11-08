
# CelltypeR functions rewrite in baseR when possible to reduce dependencies

# density plots
 
# with correct order

library(ggplot2)
library(dplyr)

plotdensity_flowset <- function(flowset, sample_include = "all", select_channels = "all",
                                scale_factor = 50, max_name_length = 20, facet_text_size = 6) {
  
  # Convert frames environment to list
  frames_list <- as.list(flowset@frames)
  
  # Get column display names from desc (or fallback to colnames)
  get_channel_names <- function(ff) {
    desc <- ff@parameters@data$desc
    desc[is.na(desc) | desc == ""] <- colnames(ff@exprs)[is.na(desc) | desc == ""]
    desc <- substr(desc, 1, max_name_length)
    names(desc) <- colnames(ff@exprs)
    desc
  }
  
  # Use the first frame to get the channel mapping
  channel_map <- get_channel_names(frames_list[[1]])
  
  # Extract exprs and add Sample
  fs_list <- lapply(names(frames_list), function(nm) {
    ff <- frames_list[[nm]]
    df <- as.data.frame(ff@exprs)
    df$Sample <- nm
    df
  })
  
  # Combine all samples
  df <- bind_rows(fs_list)
  
  # Only keep numeric columns
  numeric_cols <- sapply(df, is.numeric)
  df <- df[, c(which(numeric_cols), ncol(df))]
  
  # Filter channels
  all_channels <- colnames(df)[colnames(df) != "Sample"]
  if (!("all" %in% select_channels)) {
    all_channels <- intersect(select_channels, all_channels)
    df <- df[, c(all_channels, "Sample")]
  }
  
  # Filter samples
  if (!("all" %in% sample_include)) {
    df <- df %>% filter(Sample %in% sample_include)
  }
  
  # Reshape to long format
  df_long <- df %>%
    tidyr::pivot_longer(
      cols = -Sample,
      names_to = "Channel",
      values_to = "Value"
    )
  
  # Compute densities manually
  df_density <- do.call(rbind, lapply(split(df_long, list(df_long$Channel, df_long$Sample)), function(subdf) {
    d <- density(subdf$Value, n = 512)
    data.frame(
      Channel = unique(subdf$Channel),
      Sample = unique(subdf$Sample),
      x = d$x,
      y = d$y
    )
  }))
  
  # Offset y
  df_density <- df_density %>%
    group_by(Channel) %>%
    mutate(
      SampleIndex = as.numeric(factor(Sample, levels = unique(Sample))),
      y_offset = y * scale_factor + SampleIndex
    ) %>%
    ungroup()
  
  # Keep facet order according to channel_map
  df_density$Channel <- factor(df_density$Channel, levels = names(channel_map))
  
  # Plot
  ggplot(df_density, aes(x = x, y = y_offset, fill = Sample)) +
    geom_ribbon(aes(ymin = SampleIndex, ymax = y_offset), alpha = 0.4) +
    facet_wrap(~Channel, scales = "free", labeller = labeller(Channel = channel_map)) +
    theme_light() +
    guides(fill = FALSE) +
    labs(y = "Sample", x = "Value") +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text = element_text(size = facet_text_size)
    )
}






 
# harmonize function - no changes

harmonize <-  function(flowset, processing = 'retro',
                       two_peaks = c(9:length(colnames(transformed_flowset))),
                       one_peak = c(1:8), threshold = 0.01) {
  # biexp transform conditions
  biexp  <- biexponentialTransform("biexp transform",a = 0.5, b = 1, c = 0.5, d = 1, f = 0, w = 0)
  # run biexp transform function
  transformed_flowset <- transform(flowset, transformList(colnames(flowset), biexp))
  # if we want to see the biexp tranform
  if (processing == 'biexp') {
    return(transformed_flowset)
  } else if (processing == 'align') {
    normtr=gaussNorm(transformed_flowset,colnames(transformed_flowset)[two_peaks],max.lms = 2,peak.density.thr = threshold) #Detects and align 2 peaks on the marker 3,5,6,9...14.
    expbe_norm2=normtr$flowset
    normtr=gaussNorm(expbe_norm2,colnames(expbe_norm2)[one_peak],max.lms = 1,peak.density.thr = threshold)#Detects and align 1 peak
    aligned_transformed_flowset=normtr$flowset
    return(aligned_transformed_flowset)
  } else{
    normtr=gaussNorm(transformed_flowset,colnames(transformed_flowset)[two_peaks],max.lms = 2,peak.density.thr = threshold) #Detects and align 2 peaks on the marker 3,5,6,9...14.
    expbe_norm2=normtr$flowset
    normtr=gaussNorm(expbe_norm2,colnames(expbe_norm2)[one_peak],max.lms = 1,peak.density.thr = threshold)#Detects and align 1 peak
    aligned_transformed_flowset=normtr$flowset
    retrotransformed_flowset <- inversebiexponentialTransform(aligned_transformed_flowset)
    return(retrotransformed_flowset)
  }
}


 
######

 
#' Convert flowSet to a single dataframe with optional selection of measurement type
#' @param flowset A flowSet object
#' @param output_path Path to save CSV (if save.csv = TRUE)
#' @param save.csv Logical, whether to save CSV
#' @param select_measure One of c("all", "Area", "Height", "Width"). Default "Area"
#' @export
#' @examples
#' \dontrun{
#' flowset_to_csv(flowset)
#' flowset_to_csv(flowset, output_path = "path/to/location/", save.csv = TRUE)
#' }
#### just AB and channel
flowset_to_csv <- function(flowset, output_path = NULL, save.csv = FALSE,
                           select_measure = "Area") {
  
  measure_suffix <- c(Area = "-A", Height = "-H", Width = "-W")
  
  if (!(select_measure %in% c("all", "Area", "Height", "Width"))) {
    stop("select_measure must be one of 'all', 'Area', 'Height', 'Width'")
  }
  
  rename_cols <- function(ff) {
    desc <- ff@parameters@data$desc
    name <- ff@parameters@data$name
    
    # Replace empty desc with name
    desc[is.na(desc) | desc == ""] <- name[is.na(desc) | desc == ""]
    
    # Extract antibody only (first part before '-')
    antibody <- sapply(strsplit(desc, "-"), `[`, 1)
    
    # Determine measurement suffix
    meas_type <- substr(name, nchar(name), nchar(name))
    
    # If select_measure is "all", append suffix; otherwise just keep antibody
    if (select_measure == "all") {
      new_names <- paste0(antibody, "-", meas_type)
    } else {
      sel_suffix <- measure_suffix[select_measure]
      keep_cols <- grep(sel_suffix, name)
      ff <- ff[, keep_cols]
      new_names <- antibody[keep_cols]  # no suffix when only one measurement
    }
    
    colnames(ff@exprs) <- new_names
    ff
  }
  
  # Apply renaming
  frames_list <- lapply(names(flowset@frames), function(nm) {
    ff <- flowset@frames[[nm]]
    ff <- rename_cols(ff)
    df <- as.data.frame(ff@exprs)
    df$Sample <- nm
    df
  })
  
  df_combined <- dplyr::bind_rows(frames_list)
  
  # Add cell IDs per sample
  df_combined$cell <- as.factor(unlist(lapply(seq_along(flowset), function(i) {
    seq_len(nrow(flowset[[i]]@exprs))
  })))
  
  # Save CSV
  if (save.csv && !is.null(output_path)) {
    write.csv(df_combined, file = file.path(output_path, paste0(deparse(substitute(flowset)), ".csv")),
              row.names = FALSE)
  }
  
  return(df_combined)
}


 
##### updata for seurat 5 ######
 
# df_to_seurat
 
# creates a seurat object from the expression matrix and adds in meta data

 
#' Creates a seurat object from the Flow Cytometry expression data frame.
 
#'
#' Takes in a dataframe created from Flow Cytometry data using the fsc_to_fs and fsc_to_csv
#' functions. This function creates a Seurat object from the Flow Cytometry expression data
#' in a data frame format where cells are rows and Markers are columns.  A column indicating
#' the starting fsc file "Sample" is required and added as meta data into the Seurat object.
#' This now accounts for negative values that become NaN
#' @export
#' @examples
#' \dontrun{
#' make_seu(df = flow_dataframe, AB_vector = markers_names)
#' }
#' @importFrom Seurat CreateSeuratObject AddMetaData
# updated function
make_seu5 <- function(df, AB_vector) {
  # Select the marker values
  df2 <- df[, AB_vector, drop = FALSE]
  
  # Convert to matrix and transpose
  m <- as.matrix(df2)
  tm <- t(m)
  rownames(tm) <- colnames(df2)
  colnames(tm) <- rownames(df2)
  
  # Create Seurat object
  seu <- CreateSeuratObject(counts = tm)
  
  # Add metadata if available
  if ("Sample" %in% colnames(df)) {
    seu <- AddMetaData(object = seu, metadata = df$Sample, col.name = 'Sample')
  }
  
  # Instead of NormalizeData, copy counts to data
  seu@assays$RNA$data <- seu@assays$RNA$counts
  # seu <- NormalizeData(seu, normalization.method = "RC") # this also works but isn't a real normalization
  seu <- ScaleData(seu)
  # seu@assays$RNA$data@x[is.na(seu@assays$RNA$data@x)] <- 0

  seu <- Seurat::RunPCA(seu, npcs = length(AB_vector), features = AB_vector, verbose = FALSE)
  return(seu)
}







