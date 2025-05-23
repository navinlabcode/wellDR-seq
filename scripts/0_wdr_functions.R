findNormalCell <- function(scCNA,assay = "segment_ratios",resolution = "auto",
                          remove_XY = TRUE,simul = TRUE,seed = 17) {
  #scCNA  <- varbin_mtx_all_log_knn 
  # bindings for NSE (non-standard evaluation)
  is_aneuploid <- NULL
  
  if (remove_XY == FALSE & simul == TRUE) {
    stop("Argument simul can't be used if remove_XY == FALSE.")
  }
  
  if (resolution != "auto" & !is.numeric(resolution)) {
    stop("Resolution must be of class numeric")
  }
  
  # retrieving data
  rg <- as.data.frame(SummarizedExperiment::rowRanges(scCNA))
  seg <- SummarizedExperiment::assay(scCNA, assay)
  ncells <- ncol(scCNA)
  
  if (remove_XY == TRUE) {
    rg <- rg %>%
      dplyr::filter(
        !grepl("X", seqnames),
        !grepl("Y", seqnames)
      )
    
    seg <- seg[1:nrow(rg), ]
  }
  
  # calculating the coefficient of variation
  cv <- vapply(
    seg, function(z) {
      sd(z) / mean(z)
    },
    numeric(1)
  )
  
  if (simul == TRUE) {
    withr::with_seed(seed,
                     cv_simul <- rnorm(ncells,
                                       mean = 0,
                                       sd = 0.01
                     )
    )
    
    names(cv_simul) <- paste0("simul", 1:length(cv_simul))
    
    cv <- c(cv_simul, cv)
  }
  
  if (resolution == "auto") {
    fit <- tryCatch(
      withr::with_seed(seed, mixtools::normalmixEM(cv)),
      error = function(e) {
        message("Could not identify aneuploid cells in the dataset.")
        message("Marking all cells as diploid.")
        message("Check colData(scCNA)$find_normal_cv.")
        return("error")
      }
    )
    
    # determining resolution
    if (length(fit) > 1) {
      resolution <- fit$mu[1] + 5 * fit$sigma[1]
      #resolution <- fit$mu[1] + 3 * fit$sigma[1]
    } else {
      resolution <- 1
    }
  }
  
  if (simul == TRUE) {
    cv <- cv[!grepl("simul", names(cv))]
  }
  
  cv_df <- data.frame(sample = names(cv),
                      CV = cv)
  
  cv_df_low_cv <- cv_df %>%
    dplyr::mutate(is_aneuploid = case_when(
      CV > resolution ~ TRUE,
      TRUE ~ FALSE
    ))
  
  message(
    "Copykit detected ",
    nrow(cv_df_low_cv %>%
           dplyr::filter(is_aneuploid == FALSE)),
    " that are possibly diploid cells using a resolution of: ",
    round(resolution, 3)
  )
  
  # reordering info to add to metadata
  info <-
    cv_df_low_cv[match(
      SummarizedExperiment::colData(scCNA)$sample,
      cv_df_low_cv$sample
    ), ]
  
  # SummarizedExperiment::colData(scCNA)$is_aneuploid <- info$is_aneuploid
  SummarizedExperiment::colData(scCNA)$is_normal <- !(info$is_aneuploid)
  SummarizedExperiment::colData(scCNA)$find_normal_cv <- round(info$CV, 2)
  
  message("Added information to colData(CopyKit).")
  
  return(scCNA)
}
readVarbin <- function(dir, remove_Y = FALSE, genome_version = c('hg19'), 
                          bin_size = c('200k', '100k'), clean_names = TRUE) {
  # Reads a copy number directory and produces
  # a scCNA object as output
  # dir <- met_path
  # remove_Y = FALSE
  # clean_names = FALSE
  # genome_version = 'hg19'
  # bin_size = c('200k')
  # checks
  # dir <- paste0("./map_seg_output/", pro_name1)
  # remove_Y = T
  
  if (fs::file_exists(fs::dir_ls(
    path = dir,
    recurse = T,
    glob = "*uber*seg.txt"
  )) == FALSE) {
    stop(
      "Segment ratio matrix can't be found in the provided directory.
      Please make sure a uber.seg file can be found."
    )
  }
  
  # checking for the existence of more than one uber file
  if (length(fs::dir_ls(
    path = dir,
    recurse = T,
    glob = "*uber*seg.txt"
  )) > 1) {
    stop(
      "More than one uber.seg file can be found at the provided directory.
      Please make sure to only have one sample at that location."
    )
  }
  
  # importing data
  message("Importing segment ratios.")
  dat <- data.table::fread(input = fs::dir_ls(
    path = dir,
    recurse = T,
    glob = "*uber*seg.txt"
  ),
  showProgress = TRUE,
  integer64 = "double") %>%
    as.data.frame()
  
  colnames(dat) <- stringr::str_replace_all(colnames(dat), "\\.", "-")
  
  if (clean_names == TRUE) {
    dat <- janitor::clean_names(dat)
  }
  
  if (remove_Y == TRUE) {
    dat <- dat %>%
      dplyr::filter(chrom != 24)
  }
  
  #saving segment data
  seg_data <- dat %>%
    dplyr::select(-c(chrom,
                     chrompos,
                     abspos))
  
  
  # reading ratios
  message("Importing ratios.")
  dat_rat <- data.table::fread(fs::dir_ls(
    path = dir,
    recurse = T,
    glob = "*uber*ratio.txt"
  ),
  showProgress = TRUE,
  integer64 = "double") 
  
  if (clean_names == TRUE) {
    dat_rat <- janitor::clean_names(dat_rat)
  }
  
  if (remove_Y == TRUE) {
    dat_rat <- dat_rat %>%
      dplyr::filter(chrom != 24)
  }
  
  dat_rat <- dat_rat %>%
    dplyr::select(-c(chrom,
                     chrompos,
                     abspos))
  
  # reading bin counts
  message("Importing bin counts.")
  dat_bin <- data.table::fread(fs::dir_ls(
    path = dir,
    recurse = T,
    glob = "*uber*bin.txt"
  ),
  showProgress = TRUE,
  integer64 = "double") %>%
    as.data.frame()
  
  if (clean_names == TRUE) {
    dat_bin <- janitor::clean_names(dat_bin)
  }
  
  if (remove_Y == TRUE) {
    dat_bin <- dat_bin %>%
      dplyr::filter(chrom != 24)
  }
  
  dat_bin <- dat_bin %>%
    dplyr::select(-c(chrom,
                     chrompos,
                     abspos))
  
  # Fetch the locations (and other informations) of varbins
  rg <- dat %>%
    dplyr::select(c(chrom,
                    chrompos,
                    abspos)) %>%
    as.data.frame()
  
  genome_version <- match.arg(genome_version)
  bin_size <- match.arg(bin_size)
  grlist_varbin <- switch(genome_version,
                          hg19 = varbin_hg19_grangeslist)
  tmp_key <- paste0('res_', bin_size)
  gr_varbin_full <- grlist_varbin[[tmp_key]]
  
  
  rg$chrompos <- rg$chrompos + 1  ## To 1-based
  rg$abspos <- rg$abspos + 1  ## To 1-based
  rg2 <- rg %>% mutate(chrom = ifelse(chrom == "23", "chrX", ifelse(chrom == "24", "chrY", paste0("chr",chrom))))
  
  key_query <- paste0(rg2$chrom, '_', rg2$chrompos)
  key_ref   <- paste0(GenomicRanges::seqnames(gr_varbin_full), '_', IRanges::start(gr_varbin_full))
  idx <- match(key_query, key_ref)
  if (anyNA(idx)) {
    warning('Input ', sum(is.na(idx)),
            'varbins are not recorded in white sheet\n.')
  }
  g <- gr_varbin_full[idx, ]
  g$abspos <- rg2$abspos
  
  # creating scCNA object
  cna_obj <- CopyKit(assay = list(
    segment_ratios = seg_data,
    ratios = dat_rat,
    bin_counts = dat_bin),
    rowRanges = g
  )
  
  #sample name to metadata
  SummarizedExperiment::colData(cna_obj)$sample <- names(seg_data)
  
  # reading metrics
  if (rlang::is_empty(
    fs::dir_ls(path = dir, recurse = T, glob = "*stat_metrics.txt"))){
    warning("No metrics file found. \n
            Metrics are needed if you'd like to run copykit::runMetrics()\n
            Make sure folder metrics with file all_stat_metrics.txt can be found by copykit::runVarbinCNA()")
  } else {
    if (fs::file_exists(
      fs::dir_ls(path = dir, recurse = T, glob = "*stat_metrics.txt"))) {
      message("Importing metrics.")
      dat_metrics <- data.table::fread(
        fs::dir_ls(
          path = dir,
          recurse = T,
          glob = "*stat_metrics.txt"
        ),
        showProgress = TRUE,
        integer64 = "double"
      ) %>%
        janitor::clean_names() %>%
        dplyr::rename(sample = "sample_name") %>%
        as.data.frame()
      
      if (clean_names == TRUE) {
        dat_metrics <- dat_metrics %>%
          dplyr::mutate(sample = janitor::make_clean_names(sample))
      }
      
      # adding metrics to metadata
      # making sure they have the same order
      dat_metrics <- dat_metrics[match(SummarizedExperiment::colData(cna_obj)$sample,
                                       dat_metrics$sample),]
      
      if (identical(dat_metrics$sample,
                    SummarizedExperiment::colData(cna_obj)$sample)) {
        
        SummarizedExperiment::colData(cna_obj)$reads_total <- dat_metrics$total_reads
        SummarizedExperiment::colData(cna_obj)$reads_assigned_bins <- dat_metrics$reads_kept
        SummarizedExperiment::colData(cna_obj)$percentage_duplicates <- round(dat_metrics$dups_removed/dat_metrics$total_reads,2)
        SummarizedExperiment::colData(cna_obj)$median_bin_count <- dat_metrics$median_bin_count
      }
      
    }
  }
  
  if (remove_Y == TRUE) {
    message("Removed ChrY information.")
  }
  
  return(cna_obj)
  
}
filterCells <- function(scCNA, assay = 'segment_ratios',k = 5, resolution = 0.9, BPPARAM = BiocParallel::bpparam()) {
  if (!is.numeric(resolution)) {
    stop("Resolution needs to be a number between 0 and 1")
  }
  if (resolution < 0 || resolution > 1) {
    stop("Resolution needs to be a number between 0 and 1")
  }
  seg <- SummarizedExperiment::assay(scCNA, assay)
  message("Calculating correlation matrix.")
  # correction to avoid correlations calculations with standard deviation zero
  zero_sd_idx <- which(apply(seg, 2, sd) == 0)
  if (length(zero_sd_idx) >= 1) {
    seg[1, zero_sd_idx] <- seg[1, zero_sd_idx] + 1e-3
  }
  
  # calculating correlations
  dst <- parCor(seg, BPPARAM=BPPARAM)
  dst_knn_df <- apply(as.matrix(dst), 1, function(x) {
    mean(sort(x, decreasing = T)[2:(k + 1)])}) %>% tibble::enframe(name = "sample", value = "cor")
  
  dst_knn_df <- dst_knn_df %>% 
    dplyr::mutate(filtered = dplyr::case_when(cor >= resolution ~ "kept", cor < resolution ~ "removed"))
  
  message(
    "Adding information to metadata. Access with colData(scCNA)."
  )
  if (identical(SummarizedExperiment::colData(scCNA)$sample, dst_knn_df$sample)) {
    SummarizedExperiment::colData(scCNA)$filter_corr_value <- round(dst_knn_df$cor, 3)
    SummarizedExperiment::colData(scCNA)$filtered <- dst_knn_df$filtered
  } else
    stop("Sample names do not match metadata sample info. Check colData(scCNA).")
  message("Done.")
  return(scCNA)
}
countBreakpoints <- function(scCNA) {
  # scCNA <- varbin_mtx_tumor2
  # bindings for NSE
  arm <- chrarm <- NULL
  
  rg_chr <- SummarizedExperiment::rowRanges(scCNA) %>%
    as.data.frame() %>%
    dplyr::mutate(chrarm = paste0(seqnames, arm)) %>%
    dplyr::select(chrarm)
  
  dat_seg_cp <- segment_ratios(scCNA)
  # split by chrom
  message("Counting breakpoints.")
  dat_seg_split <- split(dat_seg_cp, dplyr::pull(rg_chr, chrarm))
  
  brkpt_by_chrom <-
    lapply(dat_seg_split, function(x) {
      apply(x, 2, function(i) {
        length(rle(i)$values) - 1
      }) %>%
        unlist()
    })
  
  brkpt_by_chrom_df <- dplyr::bind_rows(brkpt_by_chrom) %>%
    t() %>%
    as.data.frame()
  
  brkpt_count <- rowSums(brkpt_by_chrom_df)
  
  # making sure order is identical
  brkpt_count <- brkpt_count[SummarizedExperiment::colData(scCNA)$sample]
  
  SummarizedExperiment::colData(scCNA)$breakpoint_count <-
    brkpt_count
  
  return(scCNA)
}
l2e.normal.sd <- function(xs){
  # Need at least two values to get a standard deviation
  stopifnot(length(xs) >= 2)
  optim.result <- stats::optimize(
    # L2E loss function
    f=function(sd)
      # "Data part", the sample average of the likelihood
      -2 * mean(stats::dnorm(xs, sd=sd)) +
      # "Theta part", the integral of the squared density
      1/(2*sqrt(pi)*sd),
    # Parameter: standard deviation of the normal distribution fit
    interval = c(0, diff(range(xs))))
  return(optim.result$minimum)
}
overdispersion <- function(v){
  # 3 elements, 2 differences, can find a standard deviation
  stopifnot(length(v) >= 3)
  # Differences between pairs of values
  y <- v[-1]
  x <- v[-length(v)]
  # Normalize the differences using the sum. The result should be around zero,
  # plus or minus square root of the index of dispersion
  vals.unfiltered <- (y-x)/sqrt(y+x)
  # Remove divide by zero cases, and--considering this is supposed to be count
  # data--divide by almost-zero cases
  vals <- vals.unfiltered[y + x  >= 1]
  # Check that there's anything left
  stopifnot(length(vals) >= 2)
  # Assuming most of the normalized differences follow a normal distribution,
  # estimate the standard deviation
  val.sd <- l2e.normal.sd(vals)
  # Square this standard deviation to obtain an estimate of the index of
  # dispersion
  iod <- val.sd^2
  # subtract one to get the overdispersion criteria
  iod.over <- iod -1
  # normalizing by mean bincounts
  iod.norm <- iod.over/mean(v)
  return(iod.norm)
  
}
clonality_log_trinary_neu = function(log_ratio_df = log_ratio_df, lower_cutoff = -0.1, upper_cutoff = 0.1, 
                                     cell_pct = 0.90, neu_pct = 0.90){
  # lower_cutoff <- -0.15
  # upper_cutoff <-  0.15
  # cell_pct <- 0.95
  # neu_pct <- 0.9
  df_mtmap = t(log_ratio_df) %>% data.frame()
  # convert to trinary
  df_mtmap <- apply(df_mtmap, 1, function(x) ifelse(x>=upper_cutoff, 1, ifelse(x<lower_cutoff, -1, 0)))
  df_mtmap <- as.data.frame(t(df_mtmap))
  # count the # events
  df_mtmap$count.amp <- apply(df_mtmap, 1, function(x) length(which(x=="1")))
  df_mtmap$count.neu <- apply(df_mtmap, 1, function(x) length(which(x=="0")))
  df_mtmap$count.del <- apply(df_mtmap, 1, function(x) length(which(x=="-1")))
  nc = ncol(df_mtmap)-3
  # compute the x% cells
  n90 = round(cell_pct*nc, digits = 0)
  n90_neu = round(neu_pct*nc, digits = 0)
  df_mtmap = df_mtmap %>% dplyr::mutate(gene_clonal=if_else(count.amp >= n90 | count.del >= n90, "cCNA", 
                                                            ifelse(count.neu >= n90_neu, "neu","sCNA")))
  return(df_mtmap$gene_clonal)
}
getEventMat <- function(
    scCNA,          # consensus CN matrix of which will be converted to event matrix
    bin_adj = 2,    # number of bins allowed to be adjusted to consider as the same breakpoint
    ploidy_trunc = 8   # maximum integer value, all integer value larger than this will be set to this
){
  
  ## trunc integer matrix
  seg_df = copykit::consensus(scCNA)
  seg_df[seg_df>=ploidy_trunc] = ploidy_trunc
  
  intmat <- SummarizedExperiment::rowRanges(scCNA) %>%
    dplyr::as_tibble() %>%
    dplyr::select(seqnames, start, end) %>%
    cbind(seg_df) %>%
    tibble::remove_rownames()
  
  ## merge segments
  res_int <- as.data.frame(intmat[1,])
  for(i in 2:nrow(intmat)){
    if(identical(as.character(intmat[i,-c(1:3)]), as.character(intmat[i-1,-c(1:3)]))){
      next
    }else{
      res_int <- rbind(res_int, intmat[i,])
    }
  }
  res_int$bin <- as.numeric(rownames(res_int))
  
  ## finding common breakpoints
  res_int_cbp <- as.data.frame(res_int[1,])
  for(i in 2:(nrow(res_int)-1)){
    if(res_int$bin[i+1]-res_int$bin[i]<=bin_adj){
      next
    }else{
      res_int_cbp <- rbind(res_int_cbp, res_int[i,])
    }
  }
  res_int_cbp <- rbind(res_int_cbp, res_int[nrow(res_int),])
  
  res_int_cbp$n.bins=c(res_int_cbp$bin[-1], nrow(scCNA)+1) - res_int_cbp$bin
  
  res_int_cbp$end.pos = intmat$end[as.numeric(res_int_cbp$bin)+res_int_cbp$n.bins-1]
  res_int_cbp$end.chr = intmat$seqnames[as.numeric(res_int_cbp$bin)+res_int_cbp$n.bins-1]
  res_df <- res_int_cbp %>%
    dplyr::rename(start.chr=seqnames, 
                  start.pos=start) %>%
    dplyr::select(starts_with("start"), end.chr, end.pos, bin, n.bins, everything(), -end) %>%
    tibble::remove_rownames()
  
  return(res_df)
}
consensus_genomic_classes <- function(consensus_int) {
  
  # consensus_int <- ploidy_scale(ploidy_VAL, consensus_df)
  
  percent_clonal <- 1
  # percent_extant <- 1 / nrow(consensus_int)
  
  # for every bin
  ps_percents_list <- future.apply::future_apply(consensus_int, 2, function(x) {
    if (any(is.na(x))){
      perc <- 1
    } else {
      perc <- janitor::tabyl(x) %>%
        dplyr::pull(percent)
    }
    
  })
  
  bin_classes <- future.apply::future_lapply(ps_percents_list, function(x) {
    # if (any(x == percent_extant)) {
    #   return("uCNA")
    # } else 
    if (any(x == percent_clonal)) {
      return("cCNA")
    } else
      return("sCNA")
  })
  
  bin_classes <- unlist(unname(bin_classes))
  
  return(bin_classes)
}
link_bam_files <- function(bincounts, path_bam_files, output_path, output_name, target_reads = 500000) {
  
  ###
  # From the names of a uber.sample.bin.txt matrix extracts the cell names,
  # searches for the bam files within a path and outputs a txt with 
  # the path of each located bam file. Only consider files with number of reads
  # higher than `target_reads`
  
  # bincounts: uber.sample.bin matrix
  # path_bam_files: path to the folder where original bam files are located
  # output_path: path to which the output txt will be written
  # output_name: name of the output file
  # target_reads: minimum number of reads a bam file must have to be linked.
  ###
  # bincounts = bin_f
  # path_bam_files = bam_path
  # output_path = output_path
  # output_name = pro_name
  # target_reads = 500000
  
  print(paste("Number of cells in bincounts: ", length(names(bincounts)[-c(1:3)])))
  
  files <- list.files(path_bam_files, full.names = T)
  files <- files[!grepl(".bai", files)] 
  # WARNING, SOME FILES DO NOT HAVE .markdup ON THEIR filenames
  bam_name <- ifelse(stringr::str_detect(files[1], "markdup"), ".sort.markdup.bam", ".sort.bam")
  
  cells <- data.frame(cell = names(bincounts)[-c(1:3)]) %>% mutate(cell = str_replace_all(cell, "\\.", "-")) %>% 
    mutate(cell_bam = paste0(cell, bam_name)) %>% pull(cell_bam)
  
  # https://stackoverflow.com/questions/7597559/grep-using-a-character-vector-with-multiple-patterns
  matches <- unique(grep(paste(cells, collapse = "|"), files, value = T))
  
  if (length(matches) == 0) {
    stop("No matches.")
  }
  
  bam_counts <- BiocParallel::bplapply(matches,function(x) {
    count_bam <- Rsamtools::countBam(x)
    count_vector <- data.frame(file = count_bam$file,reads = count_bam$records)
  })
  
  names(bam_counts) <- matches
  bam_counts_df <- bind_rows(bam_counts, .id = 'filepath') %>% 
    dplyr::filter(reads >= target_reads)
  
  message(paste("Found", length(matches), "bam files"))
  message(paste(length(bam_counts_df$filepath), "had more than", target_reads, "reads."))
  
  write.table(bam_counts_df$filepath, paste0(output_path, output_name), sep = "\t", row.names = F, col.names = F)
}
gini.index <- function(x, n) {
  this.order <- order(x)
  x <- x[this.order]
  n <- n[this.order]
  # Fraction of population with given wealth
  f <- n / sum(n)
  # Cumulative total wealth
  s <- cumsum(x*n)
  # Delayed cumulative total wealth
  ds <- c(0, s[-length(s)])
  # Total wealth
  m <- s[length(s)]
  # Formula from wikipedia
  # https://en.wikipedia.org/w/index.php?title=Gini_coefficient&oldid=968643673#Discrete_probability_distribution
  1 - sum(f * (ds + s)) / m
}
calc_coverage <- function(path) {
  library(readr)
  library(dplyr)
  library(stringr)
  # from a path containing the covhist.txt files for each cell from running the
  # snakemake pipeline returns the breadth of coverage for a sample
  # path <- my_cov
  inpaths <- Sys.glob(paste0(path, "*.covhist.txt"))
  coverage.stats <- tibble(bed_path=inpaths) %>%
    mutate(cellname = stringr::str_extract(basename(bed_path), "^[^.]*")) %>% dplyr::group_by(cellname) %>%
    dplyr::summarize(.groups="keep", readr::read_tsv(bed_path, 
                                              col_names=c("refname", "depth", "count", "refsize", "frac"),
                                              col_types=cols(col_character(), col_double(), col_double(), col_double(), col_double())), ) %>%
    dplyr::filter(refname=="genome") %>%
    dplyr::summarize(breadth = 1 - frac[depth==0], gini_index = gini.index(depth, count), .groups="keep")
  
}
mpd_scTree <- function(df, method = "nj", metric = "manhattan", assay = "ratio", n_threads = parallel::detectCores() / 4) {
  # cores check
  if (n_threads < 1) {n_threads <- 1}
  seg_data <- df
  if (assay == "integer") {
    ## with integers
    message("Using integer data...")
    seg_data <- t(seg_data) %>% as.matrix()
    ## recommend using hamming distance for integer profiles
    distMat <- as.matrix(parallelDist::parDist(seg_data, method= "hamming", diag=T, upper=T,n_threads=n_threads))
    
    if (metric != "hamming") {
      stop("Recommend only using hamming distance for integer profiles")
    }
  } else {
    # with ratios
    message("Using ratio data...")
    seg_data <- t(seg_data) %>% as.data.frame()
    # calculating distance matrix
    message("Calculating distance matrix")
    distMat <- amap::Dist(seg_data, method = metric, nbproc = n_threads)
  }
  
  # ordering cells
  if (method %in% c("nj", "me")) {
    if (method == "nj") {
      message("Creating neighbor-joining tree.")
      tree <- ape::nj(distMat)
    }
    if (method == "me") {
      message("Creating minimum evolution tree.")
      tree <- ape::fastme.bal(distMat)
    }
  } else {
    stop("Currently only nj and me trees are supported.")
  }
  
  n<-length(tree$tip.label)
  ## removing end node
  tree$edge.length[sapply(1:n,function(x,y)   which(y==x),y=tree$edge[,2])] <- 0
  
  return(mean(cophenetic(tree)))  
}
sample_select <- function(df, n, reps = 10, seed = 17){
  if(ncol(df)<n){
    return(NULL)
  }
  
  set.seed(seed)
  l_df <- lapply(1:reps, function(i){df[,sample.int(ncol(df), n, replace = F)]})
  return(l_df)
  
}
countEvents <- function(scCNA, subset, min_bin = 0) {
  # bindings for NSE
  scCNA <- obj
  subset <- "recurrence"
  min_bin <- 1
  arm <- chrarm <- NULL
  rg_chr <- SummarizedExperiment::rowRanges(scCNA) %>% as.data.frame() %>%
    dplyr::mutate(chrarm = paste0(seqnames, arm)) %>% dplyr::select(chrarm)
  meta <- table(scCNA@colData$subclones, scCNA@colData$timepoint)
  dat_seg_cp <- as.data.frame(scCNA@consensus[,names(which(meta[,subset]!=0))])
  
  # split by chrom
  dat_seg_split <- split(dat_seg_cp, dplyr::pull(rg_chr, chrarm))
  brkpt_by_chrom <-
    lapply(dat_seg_split, function(x) {
      apply(x, 2, function(i) {
        sum(rle(i)$lengths>min_bin)
      }) %>% unlist()
    })
  
  brkpt_by_chrom_df <- dplyr::bind_rows(brkpt_by_chrom) %>% t() %>% as.data.frame()
  brkpt_count <- mean(rowSums(brkpt_by_chrom_df))
  return(brkpt_count)
}
countEvents2 <- function(scCNA, subset, min_bin = 0) {
  #---remove split by chromosomal to save events thta ocurred at chromosomal level.
  # bindings for NSE
  # scCNA <- obj
  # subset <- "recurrence"
  # subset <- "primary"
  # min_bin <- 1
  arm <- chrarm <- NULL
  rg_chr <- SummarizedExperiment::rowRanges(scCNA) %>% as.data.frame() %>%
    dplyr::mutate(chrarm = paste0(seqnames, arm)) %>% dplyr::select(chrarm)
  meta <- table(scCNA@colData$subclones, scCNA@colData$timepoint)
  dat_seg_cp <- as.data.frame(scCNA@consensus[,names(which(meta[,subset]!=0))])
  
  # split by chrom
  # dat_seg_split <- split(dat_seg_cp, dplyr::pull(rg_chr, chrarm))
  brkpt_by_chrom <-
    # lapply(dat_seg_split, function(x) {
    apply(dat_seg_cp, 2, function(i) {
      sum(rle(i)$lengths>min_bin)
    }) %>% unlist()
  # })
  
  brkpt_by_chrom_df <- dplyr::bind_rows(brkpt_by_chrom) %>% t() %>% as.data.frame()
  brkpt_count <- mean(rowSums(brkpt_by_chrom_df))
  return(brkpt_count)
}
corr_cells <- function(original,
                       downsampled) {
  
  ###
  # Loops across all cells matching by their names and calculates spearman
  # correlation
  
  
  # original: uber.sample.seg.txt matrix of segment ratios from original files
  # downsampled: uber.sample.seg.txt matrix of segment ratios from downsampled files
  ###
  
  # keeping only cells that were downsampled
  original <- original[,names(downsampled)]
  
  # sanity check
  stopifnot(identical(names(downsampled), names(original)))
  
  # cell names vector
  cell_names <- names(downsampled)[-c(1:3)]
  
  # running correlation
  cor_list <- BiocParallel::bplapply(cell_names, function(x)
    cor(downsampled[,x], original[,x], method = 'spearman')
  )
  
  # binding to a vector
  cor_vector <- do.call(c, cor_list)
  names(cor_vector) <- cell_names
  cor_vector
  
}
create_corr_df <- function(cor_1M,
                           cor_750k,
                           cor_500k,
                           cor_250k,
                           cor_125k,
                           cor_75k,
                           cor_50k) {
  
  ###
  # returns a data frame of all correlations
  ###
  cat(
    length(cor_1M),
    length(cor_750k),
    length(cor_500k),
    length(cor_250k),
    length(cor_125k),
    length(cor_75k),
    length(cor_50k)
  )
  
  data.frame(cell = c(
    names(cor_1M),
    names(cor_750k),
    names(cor_500k),
    names(cor_250k),
    names(cor_125k),
    names(cor_75k),
    names(cor_50k)
  ),
  correlation = c(cor_1M,
                  cor_750k,
                  cor_500k,
                  cor_250k,
                  cor_125k,
                  cor_75k,
                  cor_50k),
  n_reads = c(rep('1M', length(cor_1M)),
              rep('750k', length(cor_750k)),
              rep('500k', length(cor_500k)),
              rep('250k', length(cor_250k)),
              rep('125k', length(cor_125k)),
              rep('75k', length(cor_75k)),
              rep('50k', length(cor_50k))
  )) %>% 
    mutate(n_reads = fct_relevel(n_reads, n_reads_levels))
  
}

dim_res_finder2 <- function(res, dim, seurat_sub, group) {
  #---1 dim -> 1 top10---
  library(cowplot)
  library(ggplot2)
  if(missing(group)){
    group_t <- c("group")
  } else {
    group_t <- group
  }
  j<-1
  l<-1
  plot_s <- list()
  plot_u <- list()
  plot_top <- list()
  plot_rna <- list()
  plot_umi <- list()
  res_t <- res
  dim_t <- dim
  seurat_name <- deparse(substitute(seurat_sub))
  for (i in dim_t) {
    seurat_sub1_1 <- FindNeighbors(seurat_sub, dims = 1:i)
    for (k in res_t) {
      print(paste0(seurat_name, c(" Now is running: dim "),i,c("res "),k))
      seurat_sub1_2 <- FindClusters(seurat_sub1_1, resolution = k)
      seurat_sub1_1_tsne <- RunTSNE(object = seurat_sub1_2, dims = 1:i)#, perplexity = 10)
      seurat_sub1_1_tsu <- RunUMAP(object = seurat_sub1_1_tsne, dims = 1:i)#, perplexity = 10)
      p_tsne<-TSNEPlot(object = seurat_sub1_1_tsne, label = T)+ggtitle(paste0(seurat_name, c("_"), i, c("_singlets_"), k))
      p_umap<-DimPlot(seurat_sub1_1_tsu, reduction = "umap", label = T)+ggtitle(paste0(seurat_name, c("_"), i, c("_singlets_"), k))
      plot_s[[j]]<-p_tsne
      plot_u[[j]]<-p_umap
      j <- j+1
    }
    seurat_sub1_1_tsu.markers <- FindAllMarkers(seurat_sub1_1_tsu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    seurat_sub1_1_tsu.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
    top10 <- seurat_sub1_1_tsu.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
    p_top<-DoHeatmap(seurat_sub1_1_tsu, features = top10$gene)+theme(legend.position="top")+ggtitle(paste0(seurat_name, c("_"), i, c("_singlets_"), k))
    p_rna<-ggplot(seurat_sub1_1_tsu@meta.data, aes(x=seurat_clusters, y=nFeature_RNA))+ geom_boxplot()+
      geom_jitter(aes(color=seurat_clusters), size=0.5)+ggtitle(paste0(seurat_name, c("_"), i, c("_singlets_"), k))
    p_umi<-ggplot(seurat_sub1_1_tsu@meta.data, aes(x=seurat_clusters, y=nCount_RNA))+ geom_boxplot()+
      geom_jitter(aes(color=seurat_clusters), size=0.5)+ggtitle(paste0(seurat_name, c("_"), i, c("_singlets_"), k))
    
    p3<-TSNEPlot(object = seurat_sub1_1_tsne, group.by= group_t)+ggtitle(paste0(seurat_name,c("_L1-L2_singlets_1")))
    p4<-DimPlot(seurat_sub1_1_tsu, reduction = "umap", group.by= group_t)+ggtitle(paste0(seurat_name,c("_L1-L2_singlets_1")))
    plot_top[[l]] <- p_top
    plot_rna[[l]] <- p_rna
    plot_umi[[l]] <- p_umi
    plot_s[[j]]<-p3
    plot_u[[j]]<-p4
    l <- l+1
    j<-j+1
  }
  p_list_tsne<-plot_grid(plotlist = plot_s, ncol = length(res)+1)
  p_list_umap<-plot_grid(plotlist = plot_u, ncol = length(res)+1)
  p_list_top<-plot_grid(plotlist = plot_top, ncol = length(dim_t))
  p_list_rna<-plot_grid(plotlist = plot_rna, ncol = length(dim_t))
  p_list_umi<-plot_grid(plotlist = plot_umi, ncol = length(dim_t))
  plist <- list(p_list_tsne, p_list_umap,p_list_top,p_list_rna,p_list_umi)
  return(plist)
}

runVst2 <- function(scCNA,
                   transformation = c("ft", "log"),
                   assay = 'bincounts') {
  transformation <- match.arg(transformation)
  
  message(paste("Running variance stabilization transformation:",
                transformation))
  
  # recovering assay
  varbin_counts_df <- assay(scCNA, assay)
  
  if (transformation == "ft") {
    counts_df_ft <- as.data.frame(apply(varbin_counts_df,
                                        2,
                                        function(x) sqrt(x) + sqrt(x + 1)))
  }
  
  if (transformation == "log") {
    varbin_counts_df[varbin_counts_df == 0] <- 1e-4
    counts_df_ft <- as.data.frame(apply(varbin_counts_df,
                                        2,
                                        function(x) log(x)))
  }
  
  counts_df_ft <- as.data.frame(counts_df_ft)
  
  S4Vectors::metadata(scCNA)$vst <- transformation
  SummarizedExperiment::assay(scCNA, transformation) <- counts_df_ft
  
  return(scCNA)
}

knnSmooth2 <- function(scCNA,
                      k = 4,
                      BPPARAM = bpparam()) {
  # setup data
  #k = 4
  #scCNA <- varbin_mtx_all_log
  bin <- scCNA@assays@data$bin_counts
  seg <- scCNA@assays@data$segment_ratios
  
  # finding neighbors
  message("Finding neighbors.")
  neighbors <- BiocNeighbors::findKNN(t(seg), k = k)
  
  message(paste("Smoothing cells using k =", k))
  # collect neighbors and sum counts
 # smoothed_bins_list <- BiocParallel::bplapply(seq_along(bin), function(i) {
  #  cells_neighbors_df <- bin[c(i, neighbors$index[i,])]
   # smoothed_cell <- rowSums(cells_neighbors_df)
  #  smoothed_cell
  #})
  smoothed_bins_list <- lapply(seq_along(bin), function(i) {
    cells_neighbors_df <- bin[c(i, neighbors$index[i, ])]
    smoothed_cell <- rowSums(cells_neighbors_df)
    smoothed_cell
  })
  
  # re-adding names
  names(smoothed_bins_list) <- colnames(scCNA)
  
  smoothed_bins_df <- as.data.frame(do.call(cbind,
                                            smoothed_bins_list))
  
  # adding knn smoothed bins to assay and re-running Vst
  assay(scCNA, 'smoothed_bincounts') <- smoothed_bins_df
  
  # re-running vst and segmentation
  scCNA <- runVst2(scCNA, assay = 'smoothed_bincounts')
  scCNA <- runSegmentation(scCNA)
  
  message("Replacing segment_ratios assay.")
  message("Replacing logr assay.")
  
  # re-normalizing
  scCNA <- logNorm(scCNA)
  
  message("Done.")
  
  return(scCNA)
  
}
DoComHeatmap2 <- function(object,genelist,col_palette=NULL,col_title=F,assay='RNA', breaks=c(-2.5,0,2.5)){
  # This version could set breaks
  # prepare matrix
  # make sure the color palette is correct if get weird color annotation on the row
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  dmtx <- as.matrix(GetAssayData(object = object, slot = 'scale.data')[genelist,])
  ident.srt<- sort(object@active.ident)
  dmtx2 <- dmtx[,names(ident.srt)]
  #######----------------------Add Annotations------
  ch_palette <- col_palette %||% yarrr::piratepal("basel")[1:length(unique(ident.srt))]
  names(ch_palette) <- levels(ident.srt)
  
  ha_top=HeatmapAnnotation(foo=anno_block(gp=gpar(fill=ch_palette)), show_annotation_name = F)
  topn <- length(genelist)/length(unique(ident.srt))
  row_s <- rep(unique(ident.srt), each=topn)
  my_palette_h <- c("blue","white","red")
  my_col <- circlize::colorRamp2(breaks = breaks, colors = my_palette_h)
  if(col_title){
    p1 <- Heatmap(dmtx2, name = "Scaled Expression", show_row_dend = F, column_split = ident.srt, row_split = rep(unique(ident.srt), each=topn), 
                  column_title_gp = gpar(col=ch_palette),row_names_gp = gpar(col=ch_palette),column_title_rot = 90,row_title = NULL,top_annotation = ha_top,
                  heatmap_legend_param = list(direction="horizontal", title_position="topcenter"), cluster_columns = FALSE, cluster_rows = FALSE, show_row_names = T, show_column_names = FALSE, use_raster = TRUE, col = my_col) #col = my_col,
    p2 <- draw(p1, heatmap_legend_side="bottom", padding=unit(c(2,2,15,2),"mm"))
  }else{
    p1 <- Heatmap(dmtx2, name = "Scaled Expression", show_row_dend = F, column_split = ident.srt, row_split = rep(unique(ident.srt), each=topn), 
                  column_title_gp = gpar(col=ch_palette),row_names_gp = gpar(col=ch_palette),column_title_rot = 90,row_title = NULL,column_title = NULL,top_annotation = ha_top,
                  heatmap_legend_param = list(direction="horizontal", title_position="topcenter"), cluster_columns = FALSE, cluster_rows = FALSE, show_row_names = T, show_column_names = FALSE, use_raster = TRUE, col = my_col) #col = my_col,
    p2 <- draw(p1, heatmap_legend_side="bottom")
  }
  return(grid::grid.grabExpr(p2))
}

pathways_detector4 <- function(de_gene_rank, genelist_gsea_hall, genelist_gsea_c5_bp, genelist_gsea_kegg, genelist_gsea_rectome, top_n = '25'){
  # de_gene_rank=tcell_reselect_rank
  library(gplots)
  library(gdata)
  
  j<-1
  plot_s <- list()
  plot_top <- list()
  data_s <- list()
  #-------GSEA pathways--
  # genelist_gsea_hall <- geneIds(getGmt("/volumes/USR2/wangkl/db/gsea/human/h.all.v7.5.1.symbols.gmt"))
  # genelist_gsea_c5_bp <- geneIds(getGmt("/volumes/USR2/wangkl/db/gsea/human/c5.go.bp.v7.5.1.symbols.gmt"))
  # genelist_gsea_c5 <- geneIds(getGmt("/volumes/USR2/wangkl/db/gsea/human/c5.all.v7.5.1.symbols.gmt"))
  # #genelist_gsea_all <- geneIds(getGmt("/volumes/USR2/wangkl/db/gsea/human/msigdb.v6.2.symbols.gmt"))
  # genelist_gsea_kegg <- geneIds(getGmt("/volumes/USR2/wangkl/db/gsea/human/c2.cp.kegg.v7.5.1.symbols.gmt"))
  # # genelist_gsea_rectome <- geneIds(getGmt("/volumes/USR2/wangkl/db/gsea/human/c2.cp.reactome.v7.1.symbols.gmt"))
  # genelist_gsea_rectome <- geneIds(getGmt("/volumes/USR2/wangkl/db/gsea/human/c2.cp.reactome.v7.5.1.symbols.gmt"))
  
  genelist_all <- c("genelist_gsea_hall","genelist_gsea_c5_bp","genelist_gsea_kegg","genelist_gsea_rectome")
  obj_name <- deparse(substitute(de_gene_rank))
  #for (i in c("genelist_gsea_hall","genelist_gsea_c2","genelist_gsea_c5","genelist_gsea_all")) {
  for (i in genelist_all) {
    # i="genelist_gsea_c5_bp"
    # i="genelist_gsea_rectome"
    # # #NOTE: higher nperm would give more stable or reliable results, I used to use 10000, but not that great!
    res_fgsea <- fgsea(pathways = get(i), stats = de_gene_rank,minSize=10,maxSize=500)
    # res_fgsea <- fgsea(pathways = get(i), stats = de_gene_rank,minSize=10,maxSize=500)
    
    res_fgsea_srt <- res_fgsea %>% as_tibble() %>% dplyr::arrange(desc(NES))
    #my_palette<-c("#d53e4f","#bababa")
    #scales::show_col(c("#bababa","#d53e4f"))
    top20_PathwaysUp <- head(res_fgsea_srt, n=top_n)
    top20_PathwaysDown <- tail(res_fgsea_srt, n=top_n)
    topPathways <- rbind(top20_PathwaysUp,top20_PathwaysDown)
    for (k in c("res_fgsea_srt","topPathways")) {
      my_path <- get(k)
      if(min(my_path$padj)>0.05){
        my_palette <- c("#bababa")
        assign(paste0("my_palette_",k), my_palette)
      }else if(max(my_path$padj)<0.05){
        my_palette <- c("#d53e4f")
        assign(paste0("my_palette_",k), my_palette)
      }else {
        my_palette<-c("#bababa","#d53e4f")
        assign(paste0("my_palette_",k), my_palette)
      }}
    p1 <- ggplot(res_fgsea_srt, aes(reorder(pathway, NES), NES)) + geom_col(aes(fill=padj<0.05)) + coord_flip() + labs(x="Pathways", y="Normalized GSEA Enrichment Score", title= paste0("GSEA Pathways:\n", i, "\n", obj_name))+scale_fill_manual(values = my_palette_res_fgsea_srt)+theme(axis.text.y = element_text(size=8))
    p2 <- ggplot(topPathways, aes(reorder(pathway, NES), NES)) + geom_col(aes(fill=padj<0.05)) + coord_flip() + labs(x="Pathways", y="Normalized GSEA Enrichment Score", title= paste0("GSEA top ", top_n, " Pathways: \n", i,"\n", obj_name))+scale_fill_manual(values = my_palette_topPathways)+theme(axis.text.y = element_text(size=8))
    plot_s[[j]]<-p1
    plot_top[[j]]<-p2
    data_s[[j]] <- res_fgsea_srt
    j <- j+1
  }
  
  p_list_s <- plot_grid(plotlist = plot_s, ncol = length(genelist_all))
  p_list_top <- plot_grid(plotlist = plot_top, ncol = length(genelist_all))
  
  plist <- list(p_list_s, p_list_top, data_s)
  return(plist)
}
