#install.packages("LDlinkR")
#install.packages("usethis")

#### the first time you use this, run the line below first, it will open a .Renviron file. 
## usethis::edit_r_environ()
#### add the line below to the .Renviron file and save it, then restart Rstudio. 
## LDLINK_TOKEN=704ee243a529

Sys.getenv("LDLINK_TOKEN")

############ use this if website is being flagged as unsafe, otherwise ignore
library(httr)
library(LDlinkR)

# 1. Force the httr package to bypass SSL peer and host verification
httr::set_config(config(ssl_verifypeer = 0L, ssl_verifyhost = 0L))

############


library(LDlinkR)
library(reshape2)
library(data.table)
library(ggplot2)
library(ggnewscale)


# ---- Helper functions --------------------------------------------------------

# Get an LDlink matrix (r2 or D') as a numeric matrix with SNP names on both axes
get_ld_matrix <- function(snps, r2d, pop, genome_build, token) {
  ld <- as.data.frame(
    LDlinkR::LDmatrix(snps = snps, pop = pop, r2d = r2d,
                      token = token, genome_build = genome_build),
    check.names = FALSE
  )
  m <- matrix(as.numeric(as.matrix(ld[, -1])), nrow = nrow(ld),
              dimnames = list(ld[[1]], colnames(ld)[-1]))
  m
}

# Look up FORGEdb scores for the given SNPs
get_forge_scores <- function(forge_file, snps) {
  scores <- fread(forge_file)
  
  # rsID column: mostly values like "rs123"
  looks_like_rsid <- sapply(scores, function(x)
    mean(grepl("^rs[0-9]+$", as.character(head(x, 1000)))) > 0.5)
  # score column: numeric, between 0 and 10
  looks_like_score <- sapply(scores, function(x)
    is.numeric(x) && min(x, na.rm = TRUE) >= 0 && max(x, na.rm = TRUE) <= 10)
  
  rs_col    <- names(scores)[which(looks_like_rsid)[1]]
  score_col <- names(scores)[which(looks_like_score)[1]]

  forge <- data.frame(SNP = scores[[rs_col]], FORGEdb = scores[[score_col]])
  unique(forge[forge$SNP %in% snps, ])
}

# Group SNPs into LD blocks: SNPs with the same set of linked SNPs (r2 > threshold)
find_ld_blocks <- function(r2_matrix, threshold) {
  snps <- rownames(r2_matrix)
  
  # Each SNP's linked set = itself + every SNP with r2 > threshold
  linked <- sapply(snps, function(s) {
    paste(sort(union(s, snps[r2_matrix[s, ] > threshold])), collapse = ", ")
  })
  groups <- split(snps, linked)
  
  blocks <- data.frame(
    SNPs_in_group = sapply(groups, paste, collapse = ", "),
    n_SNPs        = lengths(groups),
    n_linked      = lengths(strsplit(names(groups), ", ")),
    row.names     = NULL
  )
  blocks <- blocks[order(-blocks$n_linked, blocks$SNPs_in_group), ]
  
  # Singleton = the linked set is just the SNP itself
  is_block <- blocks$n_linked > 1
  blocks$Block <- "Singleton"
  blocks$Block[is_block] <- paste0("Block", seq_len(sum(is_block)))
  
  blocks[, c("Block", "n_SNPs", "SNPs_in_group")]
}


# ---- Main function -----------------------------------------------------------

#' Plot LD (r2 / D') between SNPs with FORGEdb scores on the diagonal,
#' and print the LD blocks.

plot_ld_matrix <- function(snps,
                            gene_name,
                            disorder_name,
                            forge_file     = "rsid.scores.s.txt",
                            link_threshold = 0.95,
                            pop            = "EUR",
                            genome_build   = "grch38",
                            token          = Sys.getenv("LDLINK_TOKEN")) {
  
  snps <- unique(as.character(unlist(snps)))
  
  # 1. LD matrices (same SNPs, same order in both)
  r2_matrix <- get_ld_matrix(snps, "r2", pop, genome_build, token)
  d_matrix  <- get_ld_matrix(snps, "d",  pop, genome_build, token)
  
  snp_order <- intersect(rownames(r2_matrix), rownames(d_matrix))
  r2_matrix <- r2_matrix[snp_order, snp_order]
  d_matrix  <- d_matrix[snp_order, snp_order]
  
  dropped <- setdiff(snps, snp_order)
  if (length(dropped) > 0) message("SNPs dropped by LDlink: ", toString(dropped))
  
  # 2. Long format; rows run top to bottom, columns left to right
  ld_long <- merge(
    reshape2::melt(r2_matrix, varnames = c("SNP1", "SNP2"), value.name = "r2"),
    reshape2::melt(d_matrix,  varnames = c("SNP1", "SNP2"), value.name = "Dprime")
  )
  row_i <- match(as.character(ld_long$SNP1), snp_order)
  col_j <- match(as.character(ld_long$SNP2), snp_order)
  ld_long$SNP1 <- factor(ld_long$SNP1, levels = rev(snp_order))
  ld_long$SNP2 <- factor(ld_long$SNP2, levels = snp_order)
  
  lower_df <- ld_long[row_i >  col_j, ]   # bottom-left -> r2
  upper_df <- ld_long[row_i <  col_j, ]   # top-right   -> D'
  diag_df  <- ld_long[row_i == col_j, ]   # diagonal    -> FORGEdb
  
  # 3. FORGEdb scores on the diagonal
  forge <- get_forge_scores(forge_file, snp_order)
  missing_forge <- setdiff(snp_order, forge$SNP)
  if (length(missing_forge) > 0) {
    message("SNPs with no FORGEdb score: ", toString(missing_forge))
  }
  diag_df$FORGEdb <- forge$FORGEdb[match(as.character(diag_df$SNP1), forge$SNP)]
  diag_df$label   <- ifelse(is.na(diag_df$FORGEdb), "NA", as.character(diag_df$FORGEdb))
  
  # 4. LD blocks
  blocks <- find_ld_blocks(r2_matrix, link_threshold)
  cat("\nLD blocks (r2 >", link_threshold, "):\n")
  print(blocks, row.names = FALSE)
  
  # 5. Plot
  ggplot() +
    geom_tile(data = lower_df, aes(x = SNP2, y = SNP1, fill = r2),
              color = "white", linewidth = 0.2) +
    scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, 1),
                        name = expression(r^2), na.value = "grey90") +
    new_scale_fill() +
    geom_tile(data = upper_df, aes(x = SNP2, y = SNP1, fill = Dprime),
              color = "white", linewidth = 0.2) +
    scale_fill_gradient(low = "white", high = "steelblue4", limits = c(0, 1),
                        name = "D'", na.value = "grey90") +
    geom_tile(data = diag_df, aes(x = SNP2, y = SNP1),
              fill = "#C3B1E1", color = "black", linewidth = 0.2) +
    geom_text(data = diag_df, aes(x = SNP2, y = SNP1, label = label), size = 4) +
    coord_equal() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      plot.title.position = "plot"
    ) +
    labs(
      title = paste0("Linkage Disequilibrium of ", gene_name,
                     " SNPs Linked to ", disorder_name),
      subtitle = "Upper triangle: D'  |  Lower triangle: r\u00b2  |  Diagonal: FORGEdb score",
      x = NULL, y = NULL
    )
}

# ---- these are just the lists of snps, nothing to see here--------------------

# CACNA1C SNPs from mood disorders
cacna1c_snps_mdd_bd <- c(
  "rs1006737", 
  "rs4765905",
  "rs1024582",
  "rs10848632",
  "rs10848635",
  "rs4765913",
  "rs4765914",
  "rs10774037",
  "rs7295590",
  "rs7979389",
  "rs541349080", 
  "rs10774053",
  "rs723672",
  "rs769087",
  "rs2159100",
  "rs2370413",
  "rs7297582",
  "rs758170",
  "rs1860002",
  "rs1051375"
)

cacna1c_snps_scz <- c(
  "rs1006737", 
  "rs4765905",
  "rs1024582",
  "rs2007044",
  "rs2283291",
  "rs7297582",
  "rs216016",
  "rs216017",
  "rs12305678",
  "rs10848680",
  "rs216047",
  "rs1894979",
  "rs2238072",
  "rs115773686",
  "rs2239079",
  "rs17801265",
  "rs2239104",
  "rs2238032",
  "rs2299661",
  "rs2239015"
)

cacna1i_snps_scz <- c(
  "rs5757717",
  "rs9611198",
  "rs3788567",
  "rs132575",
  "rs136805",
  "rs713860",
  "rs738168",
  "rs5757760",
  "rs5750871",
  "rs4522708",
  "rs3788568",
  "rs5750862",
  "rs9607658",
  "rs3788568",
  "rs5750862" 
)

cacna1c_snps_substance_abuse <- c(
  "rs17801265", #AF 0.1389
  "rs7297582", #AF 0.3728
  "rs12305678", #AF 0.0910
  "rs11062086", #AF 0.18
  "rs723672", #AF 0.4687
  "rs4765687", #AF0.4325
  "rs11062277", #AF0.1761
  "rs10848680", #0.6996
  "rs2239131"
)

cacna1d_snps_cocaine <- c(
  "rs4687735",
  "rs2926559",
  "rs898407"
)

cacna1c_snps_ASD_ADHD <- c(
  "rs1006737",
  "rs1024582",
  "rs216016",
  "rs4765913",
  "rs4765905",
  "rs2370413",
  "rs10848653"
)

cacna1g_snps_ASD <- c(
  "rs198547",
  "rs198550",
  "rs198555",
  "rs757415",
  "rs12603112",
  "rs198538",
  "rs198545",
  "rs12946808"
)

cacna1i_snps_ASD <-c(
  "rs738168",
  "rs3788569",
  "rs136805",
  "rs12170452",
  "rs5750860"
)

cacna1h_snps_ASD <-c(
  "rs28365117",
  "rs370675810",
  "rs372453886",
  "rs757713867"
)

# ---- finally some action -----------------------------------------------------


print(plot_ld_matrix(
   snps = cacna1c_snps_mdd_bd,
   gene_name     = "CACNA1C",
   disorder_name = "Mood Disorders",
   forge_file    = "rsid.scores.s.txt",
   link_threshold = 0.95
 ))
 
print(plot_ld_matrix(
  snps = cacna1c_snps_scz,
  gene_name     = "CACNA1C",
  disorder_name = "Schizophrenia",
  forge_file    = "rsid.scores.s.txt",
  link_threshold = 0.95
))

print(plot_ld_matrix(
  snps = cacna1i_snps_scz,
  gene_name     = "CACNA1I",
  disorder_name = "Schizophrenia",
  forge_file    = "rsid.scores.s.txt",
  link_threshold = 0.95
))

print(plot_ld_matrix(
  snps = cacna1c_snps_substance_abuse,
  gene_name     = "CACNA1C",
  disorder_name = "Substance Abuse Disorders",
  forge_file    = "rsid.scores.s.txt",
  link_threshold = 0.95
))

print(plot_ld_matrix(
  snps = cacna1d_snps_cocaine,
  gene_name     = "CACNA1D",
  disorder_name = "Substance Abuse Disorders",
  forge_file    = "rsid.scores.s.txt",
  link_threshold = 0.95
))

print(plot_ld_matrix(
  snps = cacna1c_snps_ASD_ADHD,
  gene_name     = "CACNA1C",
  disorder_name = "ASD and ADHD",
  forge_file    = "rsid.scores.s.txt",
  link_threshold = 0.95
))

print(plot_ld_matrix(
  snps = cacna1g_snps_ASD,
  gene_name     = "CACNA1G",
  disorder_name = "Autism Spectrum Disorder",
  forge_file    = "rsid.scores.s.txt",
  link_threshold = 0.95
))

print(plot_ld_matrix(
  snps = cacna1i_snps_ASD,
  gene_name     = "CACNA1I",
  disorder_name = "Autism Spectrum Disorder",
  forge_file    = "rsid.scores.s.txt",
  link_threshold = 0.95
))

print(plot_ld_matrix(
  snps = cacna1h_snps_ASD,
  gene_name     = "CACNA1H",
  disorder_name = "Autism Spectrum Disorder",
  forge_file    = "rsid.scores.s.txt",
  link_threshold = 0.95
))






# ---- this is to make the big figure ------------------------------------------

library(ggplot2)
library(ggnewscale)



# ---- 1. Disorder -> SNP lists (CACNA1C only) ----------------------------------
# Names become the column labels on the right of the heatmap
disorder_lists <- list(
  "Mood disorders"    = cacna1c_snps_mdd_bd,
  "Schizophrenia"     = cacna1c_snps_scz,
  "Substance abuse"   = cacna1c_snps_substance_abuse,
  "ASD / ADHD"        = cacna1c_snps_ASD_ADHD
)
disorder_lists <- lapply(disorder_lists, function(x) unique(as.character(unlist(x))))

disorder_cols <- c(
  "Mood disorders"  = "#E69F00",
  "Schizophrenia"   = "#009E73",
  "Substance abuse" = "#CC79A7",
  "ASD / ADHD"      = "#0072B2"
)

all_cacna1c_snps <- unique(unlist(disorder_lists))

# ---- 2. Same heatmap as before, but for the combined SNP set --------------------
p <- plot_ld_matrix(
  snps           = all_cacna1c_snps,
  gene_name      = "CACNA1C",
  disorder_name  = "Any Disorder",
  forge_file     = "rsid.scores.s.txt",
  link_threshold = 0.95
)

# ---- 3. Build the annotation grid (SNP x disorder) ------------------------------
# Use the SNP order/levels the heatmap actually used (LDlink may drop some SNPs)
snp_levels <- levels(p$layers[[1]]$data$SNP1)
n_snps     <- length(snp_levels)

annot <- expand.grid(SNP = snp_levels, disorder = names(disorder_lists),
                     stringsAsFactors = FALSE)
annot$present <- mapply(function(s, d) s %in% disorder_lists[[d]],
                        annot$SNP, annot$disorder)
annot$annot_fill <- ifelse(annot$present, disorder_cols[annot$disorder], "grey95")
annot$x   <- n_snps + 0.5 + match(annot$disorder, names(disorder_lists))  # right of the heatmap
annot$SNP <- factor(annot$SNP, levels = snp_levels)

label_df <- data.frame(
  x     = n_snps + 0.5 + seq_along(disorder_lists),
  label = names(disorder_lists)
)

# ---- 4. Add the annotation columns to the heatmap -------------------------------
p_combined <- p +
  new_scale_fill() +
  geom_tile(data = annot, aes(x = x, y = SNP, fill = annot_fill),
            color = "white", linewidth = 0.2) +
  scale_fill_identity() +
  geom_text(data = label_df, aes(x = x, y = 0.2, label = label),
            angle = 45, hjust = 1, vjust = 1, fontface = "bold", size = 4.2) +
  coord_equal(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 40, b = 90, l = 10))

print(p_combined)






