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


################RNA data preprocessing                                        
Run_SCTV2 <- function(RNA_path1,RNA_path2) {

  rna_path <- RNA_path1 
  seurat_object <- Read10X(data.dir = rna_path)
  seurat_object <- CreateSeuratObject(counts = seurat_object)
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & percent.mt <= 30)
  seurat_object <- seurat_object[rowSums(seurat_object@assays$RNA@counts)>0,]
  
  Chip1 <- seurat_object
  Chip1$chip <- "Chip1"
  
  rna_path <- RNA_path2 
  seurat_object <- Read10X(data.dir = rna_path)
  seurat_object <- CreateSeuratObject(counts = seurat_object)
  
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & percent.mt <= 30)
  seurat_object <- seurat_object[rowSums(seurat_object@assays$RNA@counts)>0,]
  
  Chip2 <- seurat_object
  Chip2$chip <- "Chip2"
  
  object_combined <- merge(Chip1,Chip2,add.cell.ids=c("Chip1","Chip2"))
  object_combined <- DietSeurat(object_combined)
  object_combined <- SCTransform(object_combined, vst.flavor = "v2", verbose = FALSE,vars.to.regress = c("chip")) %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
    FindClusters(resolution = 0.6, verbose = FALSE)
  
  return(object_combined)
}                                        

########gene dosage function
plot_CN_heatmap_clone_state <- function(copykit_obj,assay_name,annotation_variables,color_set,NA_variable_name) {
  #######function to plot integer or segment ratio HT; assay_name: segment_ratios or integer; annotation_variables:rowannotation columns except subclones; annotation variable has to be no NA values; specify the transformed NA value to be set to white color
  color_heat <- structure(
    c(
      "#181C43",
      "#1F77BB",
      "#F0ECEB",
      "#E2C6BE",
      "#D7A290",
      "#CB7D62",
      "#BF583B",
      "#AD2E25",
      "#8D1128",
      "#640D22",
      "#3C0912"
      
    ),
    names = c("0", "1", "2", "3", "4", "5", "6","7","8","9","10"))
  ht_matrix <- t(SummarizedExperiment::assay(copykit_obj, "segment_ratios"))
  ht_matrix <- as.matrix(ht_matrix)
  
  ht_matrix_scale <- scale01(ht_matrix)
  dist_mat <- parallelDist::parallelDist(as.matrix(ht_matrix_scale),method = 'manhattan')
  
  hclust_avg <- fastcluster::hclust(dist_mat, method = 'ward.D2')
  cell_order <- data.frame(hclust_avg$labels,hclust_avg$order)
  colnames(cell_order) <- c("sample","order")
  cell_order <- cell_order[order(cell_order$order),]
  
  copykit_obj@colData$hclust_order <- cell_order$order[match(copykit_obj@colData$sample,cell_order$sample)]
  
  
  mtx_srt <- as.data.frame(copykit_obj@colData)
  
  mtx_srt <- mtx_srt[order(mtx_srt$subclones,mtx_srt$hclust_order),] 
  
  ht_mtx <- t(SummarizedExperiment::assay(copykit_obj, assay_name))
  ht_mtx <- ht_mtx[mtx_srt$sample,]
  if (assay_name=="integer") {
    ht_mtx[ht_mtx>10] <- 10
  } else {  ht_mtx <- log2(ht_mtx + 1e-3)}
  
  
  #----annotation bar--
  anno_mtx <- mtx_srt %>% dplyr::select(c("subclones3",annotation_variables)) 
  rownames(anno_mtx) <- NULL
  
  subclone_col <- setNames(color_set,paste0("c",1:length(color_set)))
  subclone_col <- subclone_col[names(subclone_col) %in% as.character(unique(anno_mtx[,1]))]
  subclones_anno <- rowAnnotation(subclones = anno_mtx$subclones3,
                                  col = list(subclones = subclone_col))
  
  color_list <- list(subclone_col)
  names(color_list) <- "subclones3"
  for (i in 1: length(annotation_variables)){
    annotation_color_i <- setNames(color_set[1:length(unique(anno_mtx[,i+1]))],unique(anno_mtx[,i+1]))
    annotation_color_i[NA_variable_name] <- "white"
    color_list[[i+1]] <- annotation_color_i
    names(color_list)[i+1] <- colnames(anno_mtx)[i+1]
  }
  
  
  cluster_anno <-
    ComplexHeatmap::rowAnnotation(
      df = anno_mtx,
      col = color_list,
      show_annotation_name = T,na_col="white"
    )
  
  ########clone state annotation
  Clone_state_col<- setNames(color_set[1:length(unique(CN_bins_df$clone_state))],unique(CN_bins_df$clone_state))
  CN_DE_annot <- ComplexHeatmap::HeatmapAnnotation(
    clone_state = CN_bins_df$clone_state,
    show_legend = T,
    which = "column",
    col = list(clone_state = Clone_state_col,height = unit(2, "cm"),annotation_name_side="left",border=T)
  )
  
  ###########chr bar annotation 
  ########
  # obtaining data
  seg_data <- t(SummarizedExperiment::assay(copykit_obj, "segment_ratios"))
  
  #checking for duplicated names and making names unique if so
  if (any(duplicated(row.names(seg_data)))) {
    row.names(seg_data) <- make.names(row.names(seg_data), unique = TRUE)
  }
  
  # chromosome bar aesthetic
  chr_ranges <-
    as.data.frame(SummarizedExperiment::rowRanges(copykit_obj))
  chr_lengths <- rle(as.numeric(chr_ranges$seqnames))$lengths
  
  
  chr_binary <- c(rep(c(2, 1), length(chr_lengths) / 2),2)
  
  chr <-
    data.frame(chr = rep.int(x = chr_binary, times = chr_lengths))
  
  # getting lengths for chr numbers annotation
  chr_rl_c <- c(1, cumsum(chr_lengths))
  
  # creating a data frame to calculate rowMeans
  chr_df <-
    data.frame(
      a = chr_rl_c[1:length(chr_rl_c) - 1],
      b = chr_rl_c[2:length(chr_rl_c)]
    )
  chr_l_means <- round(rowMeans(chr_df))
  
  chrom.names <- c(1:22, "X")
  
  # creating the vector for chr number annotations
  v <- vector(length = sum(chr_lengths), mode = "character")
  suppressWarnings(v[chr_l_means] <- chrom.names)
  v[is.na(v)] <- ""
  
  # chr bar with the chr names
  chr_bar_200k <-
    ComplexHeatmap::HeatmapAnnotation(
      chr_text = ComplexHeatmap::anno_text(v[1:ncol(seg_data)],
                                           gp = grid::gpar(fontsize = 14)
      ),
      df = as.character(chr[1:nrow(chr), ]),
      show_legend = FALSE,
      show_annotation_name = FALSE,
      which = "column",
      col = list(df = c("1" = "grey88", "2" = "black"))
    )
  ##############################
  #col = circlize::colorRamp2(breaks = c(-2,0,2), c("dodgerblue3", "white", "firebrick3"))
  if(assay_name=="segment_ratios") {
    p1 <- ComplexHeatmap::Heatmap(as.matrix(ht_mtx), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
                                  row_split = anno_mtx$subclones3, name = "Log2 (Ratio)", show_row_names = F, show_column_names = F, row_title = paste0(nrow(ht_mtx), " single cells"),use_raster = T, raster_quality = 5, col = circlize::colorRamp2(breaks = c(-1,0,1), c("dodgerblue3", "white", "firebrick3")),top_annotation = chr_bar_200k, left_annotation = cluster_anno,bottom_annotation = CN_DE_annot)
  } else {
    
    p1 <- ComplexHeatmap::Heatmap(as.matrix(ht_mtx), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
                                  row_split = anno_mtx$subclones3,name = "Log2 (Ratio)", show_row_names = F, show_column_names = F, row_title = paste0(nrow(ht_mtx), " single cells"),use_raster = T, raster_quality = 5, col = color_heat,top_annotation = chr_bar_200k, left_annotation = cluster_anno,bottom_annotation = CN_DE_annot)
  }
  
  return(p1)
  
}

all_gene_CN_exp_correlation <- function(data,feature_column,assay) {
  ######data: long format gene/bin expression data; assay: seg or integer
  library(dplyr)
  
  gene_symbols <- unique(data[,feature_column])
  n_genes <- length(gene_symbols)
  
  # Pre-allocate vectors for results
  correlations <- numeric(n_genes)
  p_values <- numeric(n_genes)
  mean_CNs <- numeric(n_genes)
  mean_exps <- numeric(n_genes)
  delta_CNs <- numeric(n_genes)
  delta_exps <- numeric(n_genes)
  
  for (i in seq_along(gene_symbols)) {
    gene_symbol <- gene_symbols[i]
    
    gene_integer <- data[data[,feature_column]==gene_symbol, c(assay,"subclone"),drop=F]
    gene_exp <- data[data[,feature_column]==gene_symbol, c("RNA","subclone"),drop=F]
    rownames(gene_exp) <- gene_exp$subclone
    gene_exp <- gene_exp[gene_integer$subclone,]
    
    gene_integer <- gene_integer[,1]
    gene_exp <- gene_exp[,1]
    
    correlation <- cor(gene_integer, gene_exp)
    correlations[i] <- correlation
    p_values[i] <- cor.test(gene_integer, gene_exp)$p.value
    mean_CNs[i] <- mean(gene_integer)
    mean_exps[i] <- mean(gene_exp)
    delta_CNs[i] <- diff(range(gene_integer))
    delta_exps[i] <- diff(range(gene_exp))
    
    if (i %% 50 == 0) print(i)
  }
  
  # Create the resulting data frame
  gene_exp_CN <- data.frame(
    Gene = gene_symbols,
    correlation = correlations,
    p_value = p_values,
    mean_CN = mean_CNs,
    mean_exp = mean_exps,
    delta_CN = delta_CNs,
    delta_exp = delta_exps
  )
  
  return(gene_exp_CN)
}            


calculate_subclone_pair_gene_dosage3 <- function(copykit_obj,seurat_obj,logFC_cutoff,subclone_variable,rna_subclone_variable){
  ######input: copykit tumor object with subclone info (non-matched DNA data is preferred since gives more reliable bin de analysis); seurat_obj; RNA cancer data (matched cell only)
  #######remove global gene dosage part
  seg_df <- assay(copykit_obj,"segment_ratios")
  
  X <- as.matrix(seg_df)
  
  Y <- as.character(copykit_obj@colData[,subclone_variable])
  
  # Define group pairs
  groups <- unique(Y)
  group_pairs <- combn(groups, 2, simplify = FALSE)
  
  # Initialize list to store results
  pairwise_results_dna <- list()
  pairwise_results_rna <- list()
  pairwise_dna_sig <- list()
  pairwise_rna_sig <- list()
  # Loop through group pairs
  for (pair in group_pairs) {
    group1 <- pair[1]
    group2 <- pair[2]
    print(paste0("comparing ",group1, ":", group2))
    # Run Wilcoxon rank-sum test using presto
    result <- wilcoxauc(X, Y, groups_use = c(group1, group2))
    
    # Extract results and store
    result <- result[result$group==group1,]
    
    result$group_pair <- paste(group1, group2, sep = "_")
    result$feature <- str_replace(result$feature,"Feature","")
    pairwise_results_dna[[paste(group1, group2, sep = "_")]] <- result
    
    seurat_obj_subset <- seurat_obj[,seurat_obj@meta.data[,rna_subclone_variable] %in% c(group1,group2)]
    seurat_obj_subset <- seurat_obj_subset[!rowSums(seurat_obj_subset@assays$RNA@counts)==0,]
    mappable_bins <- hg19_mappable_genes$pos[hg19_mappable_genes$gene %in% rownames(seurat_obj_subset)]
    result_rna <- FindMarkers(seurat_obj_subset,ident.1 = group1,ident.2 = group2,group.by = rna_subclone_variable,min.pct=0.2, logfc.threshold = 0)
    result_rna$group_pair <- paste(group1, group2, sep = "_")
    pairwise_results_rna[[paste(group1, group2, sep = "_")]] <- result_rna
    
    dna_up_sig <-  result[result$padj<0.05 & result$logFC >=0.2,]
    if(nrow(dna_up_sig)==0) {dna_up_sig <- NULL} else {
      ########remove unmappable bins 
      dna_up_sig$direction <- "up"
      dna_up_sig <- dna_up_sig[dna_up_sig$feature %in% mappable_bins,]
      
      
      if(nrow(dna_up_sig)==0) {dna_up_sig <- NULL} else {dna_up_sig <- dna_up_sig}
    }
    dna_down_sig <-  result[result$padj<0.05 & result$logFC <= (-0.2),]
    if(nrow(dna_down_sig)==0) {dna_down_sig <- NULL} else {
      ########remove unmappable bins 
      dna_down_sig$direction <- "down"
      dna_down_sig <- dna_down_sig[dna_down_sig$feature %in% mappable_bins,]
      
      if(nrow(dna_down_sig)==0) {dna_down_sig <- NULL} else {dna_down_sig <- dna_down_sig}
    }
    
    
    rna_up_sig <- result_rna[result_rna$p_val_adj<0.05 & result_rna$avg_log2FC >= logFC_cutoff,]
    if(nrow(rna_up_sig)==0) {rna_up_sig <- NULL} else {
      rna_up_sig$gene <- rownames(rna_up_sig)
      rna_up_sig$direction <- "up"
      rna_up_sig$CN_bins <- hg19_mappable_genes$pos[match(rna_up_sig$gene,hg19_mappable_genes$gene)]
      
      
      if(nrow(rna_up_sig)==0) {rna_up_sig <- NULL} else {rna_up_sig <- rna_up_sig}
    }
    
    rna_down_sig <- result_rna[result_rna$p_val_adj<0.05 & result_rna$avg_log2FC <= (logFC_cutoff*(-1)),]
    if(nrow(rna_down_sig)==0) {rna_down_sig <- NULL} else {
      rna_down_sig$gene <- rownames(rna_down_sig)
      rna_down_sig$direction <- "down"
      rna_down_sig$CN_bins <- hg19_mappable_genes$pos[match(rna_down_sig$gene,hg19_mappable_genes$gene)]
      
      if(nrow(rna_down_sig)==0) {rna_down_sig <- NULL} else {rna_down_sig <- rna_down_sig}
    }
    
    dna_up_sig$subclone_gene_dosage <- ifelse(dna_up_sig$feature %in% unique(rna_up_sig$CN_bins),"Positive_dosage",ifelse(dna_up_sig$feature %in% unique(rna_down_sig$CN_bins),"Negative_dosage","No_dosage"))
    dna_down_sig$subclone_gene_dosage <- ifelse(dna_down_sig$feature %in% unique(rna_down_sig$CN_bins),"Positive_dosage",ifelse(dna_down_sig$feature %in% unique(rna_up_sig$CN_bins),"Negative_dosage","No_dosage"))
    
    dna_sig <- rbind(dna_up_sig,dna_down_sig)
    
    rna_up_sig$subclone_gene_dosage <- ifelse(rna_up_sig$CN_bins %in% unique(dna_up_sig$feature),"Positive_dosage",ifelse(rna_up_sig$CN_bins %in% unique(dna_down_sig$feature),"Negative_dosage","No_dosage"))
    rna_down_sig$subclone_gene_dosage <- ifelse(rna_down_sig$CN_bins %in% unique(dna_down_sig$feature),"Positive_dosage",ifelse(rna_down_sig$CN_bins %in% unique(dna_up_sig$feature),"Negative_dosage","No_dosage"))
    
    rna_sig <- rbind(rna_up_sig,rna_down_sig)
    rna_sig <- rna_sig[!is.na(rna_sig$gene),]
    
    pairwise_dna_sig[[paste(group1, group2, sep = "_")]] <- dna_sig
    pairwise_rna_sig[[paste(group1, group2, sep = "_")]] <- rna_sig
  }
  
  # Combine results into a single dataframe
  pairwise_results_dna_df <- do.call(rbind, pairwise_results_dna)
  rownames(pairwise_results_dna_df) <- NULL
  pairwise_results_rna_df <- do.call(rbind, pairwise_results_rna)
  rownames(pairwise_results_rna_df) <- NULL
  pairwise_dna_sig_df <- do.call(rbind, pairwise_dna_sig)
  rownames(pairwise_dna_sig_df) <- NULL
  pairwise_rna_sig_df <- do.call(rbind, pairwise_rna_sig)
  rownames(pairwise_rna_sig_df) <- NULL
  
  final_list <- list("dna_result"=pairwise_results_dna_df,"rna_result"=pairwise_results_rna_df,"dna_sig_result"=pairwise_dna_sig_df,"rna_sig_result"=pairwise_rna_sig_df)
  write_rds(final_list,paste0("subclone_pair_gene_dosage_result_cutoff_",logFC_cutoff,".rds"))
  return(final_list)
}

run_gene_dosage <- function(DNA_data_tumor,RNA_cancer_matched,sample_name,subclone_variable,rna_subclone_variable,rna_cell_variable) 
  {
  #######input: DNA full data (cancer cell only), RNA cancer data (matched data only)
  #sample_name = "BCIS66T"
  #subclone_variable="subclones3"
  #rna_subclone_variable="subclone3"
  #rna_cell_variable <- "RNA_Cell_name"
  
  subclonal_bins_seg <- wilcoxauc(DNA_data_tumor,subclone_variable,assay="segment_ratios")
  
  head(subclonal_bins_seg)
  dim(subclonal_bins_seg)
  length(unique(subclonal_bins_seg$feature))
  subclonal_bins_sig <- subclonal_bins_seg[subclonal_bins_seg$padj<0.05 & abs(subclonal_bins_seg$logFC) >=0.15,]
  nrow(subclonal_bins_sig)
  length(unique(subclonal_bins_sig$feature))
  head(subclonal_bins_sig)
  min(subclonal_bins_sig$logFC)
  
  subclone_columns <- as.numeric(str_replace(unique(subclonal_bins_sig$feature),"Feature",""))
  CN_bins_df <- data.frame(rownames(assay(DNA_data_tumor,"integer")))
  colnames(CN_bins_df) <- "CN_bin"
  CN_bins_df$clone_state <- ifelse(CN_bins_df$CN_bin %in% subclone_columns,"Subclonal","Clonal/Neutral")
  table(CN_bins_df$clone_state)
  
  DNA_data_tumor@colData$RNA_clusters <- as.character(DNA_data_tumor@colData$RNA_clusters)
  DNA_data_tumor@colData$RNA_clusters[is.na(DNA_data_tumor@colData$RNA_clusters)] <- "RNA_missing"
  p1 <- plot_CN_heatmap_clone_state(DNA_data_tumor,"segment_ratios",c("experiment","RNA_clusters"),YZRY_color_set2,"RNA_missing")
  pdf("Tumor_Seg_ratio_CN_bin_clone_state.pdf",height = 10,width = 8,useDingbats = F)
  print(p1)
  dev.off()
  ################step2, subset subclonal bins with mappable genes that are highly expressed
  ########subset RNA data to matched data only
  
  normalized_exp_df <- RNA_cancer_matched@assays$RNA@data
  dim(normalized_exp_df)
  ####furtehr subset to mappable genes
  normalized_exp_df_filtered <- data.frame(normalized_exp_df)
  normalized_exp_df_filtered <- normalized_exp_df_filtered[rownames(normalized_exp_df_filtered) %in% hg19_mappable_genes$gene,]
  
  normalized_exp_df_filtered_t <- data.frame(t(normalized_exp_df_filtered))
  dim(normalized_exp_df_filtered_t)
  normalized_exp_df_filtered_t$subclone <- RNA_cancer_matched$subclone3[match(rownames(normalized_exp_df_filtered_t),colnames(RNA_cancer_matched))]
  table(normalized_exp_df_filtered_t$subclone)
  
  subclone_mean <- normalized_exp_df_filtered_t %>%
    dplyr::group_by(subclone) %>%
    dplyr::summarise(
      across(everything(), 
             list(mean = ~ mean(.x, na.rm = TRUE), 
                  var = ~ var(.x, na.rm = TRUE)), 
             .names = "{.col}_{.fn}")
    )
  subclone_mean <- data.frame(subclone_mean)
  subclone_means <- subclone_mean[,grep("mean",colnames(subclone_mean))]
  dim(subclone_means)
  colnames(subclone_means) <- str_replace(colnames(subclone_means),"_mean","")
  colnames(subclone_means) <- str_replace(colnames(subclone_means),"\\.","-")
  
  #######remove lowly expressed genes
  # calculate max subclonal mean expression for each gene
  row_means <- apply(subclone_means, 2, max)
  
  # Calculate row-wise variances
  row_variances <- apply(normalized_exp_df_filtered, 1, var)
  
  result_exp <- data.frame(row_means, row_variances)
  p2 <- ggplot(result_exp,aes(x=row_means,y=row_variances))+geom_point(size=0.1)+theme_minimal()+xlab("Max average subclonal normalized expression")+ylab("Average gene expression variance")+geom_hline(yintercept = 0.1,color="red",linetype="dashed")+geom_vline(xintercept = 0.1,color="red",linetype="dashed")
  
  pdf("Gene_expression_distribution_plot.pdf",height = 6,width = 6)
  print(p2)
  dev.off()
  
  all_genes <- data.frame(rownames(normalized_exp_df_filtered))
  colnames(all_genes) <- "Gene"
  all_genes$subclonal <- ifelse(all_genes$Gene %in% hg19_mappable_genes$gene[hg19_mappable_genes$pos %in% subclone_columns],"Subclonal","Clonal/neutral")
  table(all_genes$subclonal)
  all_genes$mean_exression <- row_means[match(all_genes$Gene,names(row_means))] 
  min(all_genes$mean_exression)
  max(all_genes$mean_exression)
  all_genes$expression_bins <- ifelse(all_genes$mean_exression<=0.05,"0-0.05",ifelse(all_genes$mean_exression<=0.1 & all_genes$mean_exression>0.05,"0.05-0.1",ifelse(all_genes$mean_exression<=0.5 & all_genes$mean_exression>0.1,"0.1-0.5",ifelse(all_genes$mean_exression<=1 & all_genes$mean_exression>0.5,"0.5-1",ifelse(all_genes$mean_exression<=2 & all_genes$mean_exression>1,"1-2",">2")))))
  table(all_genes$expression_bins)
  
  all_genes$expression_bins <- factor(all_genes$expression_bins,levels = c("0-0.05","0.05-0.1","0.1-0.5","0.5-1","1-2",">2"))
  
  all_genes |> 
    tidyplot(x = expression_bins, color = subclonal) |> 
    add_barstack_absolute() |> adjust_x_axis_title("Expression bins") |> adjust_y_axis_title("Gene count") |> adjust_legend_title("Gene dosage") 
  ###############################identify subclonal gene-dosage sensitive genes
  counts_exp_df <- RNA_cancer_matched@assays$RNA@counts
  counts_exp_df <- counts_exp_df[rownames(counts_exp_df) %in% rownames(normalized_exp_df_filtered),]
  dim(counts_exp_df)
  RNA_exp_obj <- CreateSeuratObject(counts_exp_df)
  subclone_info <- data.frame(DNA_data_tumor@colData[,c(subclone_variable,rna_cell_variable)])
  rownames(subclone_info) <- subclone_info[,2]
  RNA_exp_obj$subclones <- subclone_info[,1][match(colnames(RNA_exp_obj),rownames(subclone_info))]
  table(RNA_exp_obj$subclones)
  RNA_exp_obj$subclones <- droplevels(RNA_exp_obj$subclones)
  
  exp_pseudobulk <- AggregateExpression(RNA_exp_obj,assays = "RNA",slot = "counts",return.seurat = T,group.by = "subclones")
  
  exp_pseudobulk <- data.frame(exp_pseudobulk@assays$RNA@counts)
  meta_bulk <- data.frame(colnames(exp_pseudobulk))
  rownames(meta_bulk) <- meta_bulk$colnames.exp_pseudobulk.
  library(DESeq2)
  exp_dds_bulk <- DESeqDataSetFromMatrix(countData = exp_pseudobulk,
                                         colData = meta_bulk,
                                         design = ~ 1)
  exp_dds_bulk <- vst(exp_dds_bulk, blind=TRUE)
  exp_dds_bulk <- assay(exp_dds_bulk)
  dim(exp_dds_bulk)
  min(exp_dds_bulk)
  max(exp_dds_bulk)
  rownames(exp_dds_bulk)
  exp_dds_bulk <- data.frame(exp_dds_bulk)
  write_rds(exp_dds_bulk,"subclonal_mappable_expr_pseudobulk.rds")
  
  exp_dds_bulk$CN_bins <- hg19_mappable_genes$pos[match(rownames(exp_dds_bulk),hg19_mappable_genes$gene)]
  length(unique(exp_dds_bulk$CN_bins))
  

  
  DNA_data_tumor <- calcConsensus(DNA_data_tumor, consensus_by = subclone_variable, assay = 'segment_ratios',fun = "mean")
  CN_seg_consensus <- data.frame(DNA_data_tumor@consensus)
  rownames(CN_seg_consensus) <- str_replace(rownames(CN_seg_consensus),"V","")
  
  DNA_data_tumor <- calcConsensus(DNA_data_tumor, consensus_by = subclone_variable, assay = 'integer',fun = "mean")
  CN_integer_consensus <- data.frame(DNA_data_tumor@consensus)
  rownames(CN_integer_consensus) <- str_replace(rownames(CN_integer_consensus),"V","")
  
  exp_dds_bulk$gene <- rownames(exp_dds_bulk)
  exp_dds_bulk_long <- gather(exp_dds_bulk,subclones,exp,colnames(exp_dds_bulk)[1]:colnames(exp_dds_bulk)[ncol(exp_dds_bulk)-2],factor_key = T)
  head(exp_dds_bulk_long)
  #########################VERY Important!!!!!!! subclones must be character variable! Otherwise CN assign will be messed up!
  exp_dds_bulk_long$subclones <- as.character(exp_dds_bulk_long$subclones)
  ###################################################
  for (i in 1:nrow(exp_dds_bulk_long)){
    exp_dds_bulk_long$seg[i] <- CN_seg_consensus[exp_dds_bulk_long$CN_bins[i],exp_dds_bulk_long$subclones[i]]
    exp_dds_bulk_long$integer[i] <- CN_integer_consensus[exp_dds_bulk_long$CN_bins[i],exp_dds_bulk_long$subclones[i]]
  }
  head(exp_dds_bulk_long)
  colnames(exp_dds_bulk_long)[3] <- "subclone"
  colnames(exp_dds_bulk_long)[4] <- "RNA"
  
  
  if(length(unique(subclone_info[,1]))>2) {
    gene_seg_cor <- all_gene_CN_exp_correlation(exp_dds_bulk_long,"gene","seg")
    gene_integer_cor <- all_gene_CN_exp_correlation(exp_dds_bulk_long,"gene","integer")
    print(nrow(gene_seg_cor[gene_seg_cor$p_value<=0.05,]))
    print(nrow(gene_integer_cor[gene_integer_cor$p_value<=0.05,]))
    write_rds(gene_seg_cor,"gene_seg_cor.rds")
    write_rds(gene_integer_cor,"gene_integer_cor.rds")
    breast_genes <- unique(c("CKDN2C","AKT3","CKS1B","RUVBL1","C8orf4","FGFR1","BAG4","MTDH","MYC","PAK1","CDK4","MDM2", "FHIT","NEDD9","FOXK2","CCNE2","PLA2G10","GRB7","RPS6KB1","PPM1D","CCNE1","YWHAB","ZNF217","AURKA","PTK6","PPP2R2A","CCND1","NCOA3","TP53","CDKN1A","GADD45A","MDM2","RB1","CDKN2A","CDKN2C","ESR1","SHC1","PGR","PIK3CA","PTEN","BCL2","GATA3","FGFR4","EGFR","ERBB2","EMSY","BRCA1","FOXA1","BRCA2","PIK3CA","TP53","NF1","KMT2C","NCOR1","ATM","MAP3K1","GATA3","PIK3R1","TBX3","ERBB3","ARID1A","MAP2K4","ERBB2","CBFB","CDH1","RUNX1","CHEK2","BRCA1","BRCA2","SMAD4"))
    
    breast_genes <- breast_genes[breast_genes %in% exp_dds_bulk_long$gene]
    breast_genes_df <- data.frame(breast_genes)
    breast_genes_df$cor <- gene_integer_cor$correlation[match(breast_genes_df$breast_genes,gene_integer_cor$Gene)]
    
    breast_genes_df <- breast_genes_df[order(-breast_genes_df$cor),]
    pdf("Breast_genes_integer_gene_dosage.pdf",height = 6,width = 6,useDingbats = F)
    for (i in 1: length(breast_genes_df$breast_genes)){
      p1 <- ggscatter(exp_dds_bulk_long[exp_dds_bulk_long$gene==breast_genes_df$breast_genes[i],], x = "integer", y = "RNA", 
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson",
                      xlab = "Average integer copy number", ylab = "Mean average expression")+theme(aspect.ratio = 1)+geom_point(aes(color = subclone))+scale_color_manual(values =subclone_color,name="DNA subclones")+ggtitle(breast_genes_df$breast_genes[i])
      print(p1)
    }
    dev.off()
    
    breast_genes_df$cor <- gene_seg_cor$correlation[match(breast_genes_df$breast_genes,gene_seg_cor$Gene)]
    
    breast_genes_df <- breast_genes_df[order(-breast_genes_df$cor),]
    
    pdf("Breast_gene_seg_gene_dosage.pdf",height = 6,width = 6,useDingbats = F)
    for (i in 1: length(breast_genes_df$breast_genes)){
      p1 <- ggscatter(exp_dds_bulk_long[exp_dds_bulk_long$gene==breast_genes_df$breast_genes[i],], x = "seg", y = "RNA", 
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson",
                      xlab = "Average CN segment ratios", ylab = "Mean average expression")+theme(aspect.ratio = 1)+geom_point(aes(color = subclone))+scale_color_manual(values =subclone_color,name="DNA subclones")+ggtitle(breast_genes_df$breast_genes[i])
      print(p1)
    }
    dev.off()
    
    ########create gene classification dataframe 
    
    
    ###############################################################################use segment ratios for global gene dosage
    dosage_sensitive_genes <- gene_seg_cor[gene_seg_cor$p_value<0.05 & gene_seg_cor$correlation>0 ,]
    dosage_sensitive_genes <- dosage_sensitive_genes[!is.na(dosage_sensitive_genes$Gene),]
    
    dosage_sensitive_genes_integer <- gene_integer_cor[gene_integer_cor$p_value<0.05 & gene_integer_cor$correlation>0 ,]
    dosage_sensitive_genes_integer <- dosage_sensitive_genes_integer[!is.na(dosage_sensitive_genes_integer$Gene),]
    
    table(dosage_sensitive_genes$Gene %in% dosage_sensitive_genes_integer$Gene)
    
    dosage_sensitive_genes[!dosage_sensitive_genes$Gene %in% dosage_sensitive_genes_integer$Gene,]
    
    dosage_sensitive_genes_integer[!dosage_sensitive_genes_integer$Gene %in% dosage_sensitive_genes$Gene,]
    
    table(dosage_sensitive_genes$Gene %in% all_genes$Gene[all_genes$subclonal=="Subclonal"])
    ###########check genes been detected as DSG when integer CN doesn't change
    nrow(gene_integer_cor[gene_integer_cor$delta_CN==0,])
    table(dosage_sensitive_genes$Gene %in% gene_integer_cor[gene_integer_cor$delta_CN==0,]$Gene)
    table(dosage_sensitive_genes$Gene %in% gene_integer_cor[gene_integer_cor$delta_CN==0,]$Gene, dosage_sensitive_genes$Gene %in% all_genes$Gene[all_genes$subclonal=="Subclonal"])
    
    head(dosage_sensitive_genes)
    dosage_sensitive_genes$CN_integer_dif <- ifelse(dosage_sensitive_genes$Gene %in% gene_integer_cor[gene_integer_cor$delta_CN==0,]$Gene,"No_dif","Dif")
    dosage_sensitive_genes$CN_region <- ifelse(dosage_sensitive_genes$Gene %in% all_genes$Gene[all_genes$subclonal=="Subclonal"],"Subclonal","Clonal/neutral")
    table(dosage_sensitive_genes$CN_integer_dif)
    table(dosage_sensitive_genes$CN_region,dosage_sensitive_genes$CN_integer_dif)
    
    dosage_sensitive_genes$Gene[dosage_sensitive_genes$CN_integer_dif=="No_dif" & dosage_sensitive_genes$CN_region=="Subclonal"]
    
    dosage_sensitive_genes$Gene[!dosage_sensitive_genes$Gene %in% all_genes$Gene[all_genes$subclonal=="Subclonal"]]
    dosage_sensitive_genes <- dosage_sensitive_genes[dosage_sensitive_genes$CN_integer_dif=="Dif" & dosage_sensitive_genes$CN_region== "Subclonal",]
    all_genes$gene_dosage <- ifelse(all_genes$subclonal == "Clonal/neutral" ,"no_dosage","dosage")
    ##########refine dosage genes by moving those w/o integer CN dif to no dosage genes
    all_genes$gene_dosage[all_genes$Gene[all_genes$gene_dosage=="dosage"] %in% gene_integer_cor[gene_integer_cor$delta_CN==0,]$Gene] <- "no_dosage"
    
    all_genes$gene_dosage[all_genes$gene_dosage=="dosage" & (all_genes$Gene %in% dosage_sensitive_genes$Gene)] <- "dosage-sensitive"
    all_genes$gene_dosage[all_genes$gene_dosage=="dosage" & (!all_genes$Gene %in% dosage_sensitive_genes$Gene)] <- "dosage-insensitive"
    
    ##########relabel low expressed genes 
    all_genes$gene_dosage[all_genes$Gene %in% rownames((result_exp[result_exp$row_means<=0.1 & result_exp$row_variances<=0.1,]))] <- "low_expressed"
    
    table(all_genes$gene_dosage)
    
    all_genes$CN_bins <- hg19_mappable_genes$pos[match(all_genes$Gene,hg19_mappable_genes$gene)]
    length(unique(all_genes$CN_bins[all_genes$gene_dosage=="dosage-sensitive"]))
    length(unique(all_genes$CN_bins[all_genes$gene_dosage %in% c("dosage-sensitive","dosage-insensitive")]))
    
    CN_bins_df <- data.frame(rownames(assay(DNA_data_tumor,"integer")))
    colnames(CN_bins_df) <- "CN_bin"
    CN_bins_df$clone_state <- ifelse(CN_bins_df$CN_bin %in% subclone_columns,"dosage_bins","No-dosage bins")
    CN_bins_df$clone_state[CN_bins_df$clone_state=="dosage_bins" & (CN_bins_df$CN_bin %in% unique(all_genes$CN_bins[all_genes$gene_dosage=="dosage-sensitive"]))] <- "Dosage-sensitive bins"
    CN_bins_df$clone_state[CN_bins_df$clone_state=="dosage_bins" & (CN_bins_df$CN_bin %in% unique(all_genes$CN_bins[all_genes$gene_dosage=="dosage-insensitive"]))] <- "Dosage-insensitive bins"
    CN_bins_df$clone_state[CN_bins_df$clone_state=="dosage_bins" & (CN_bins_df$CN_bin %in% unique(all_genes$CN_bins[all_genes$gene_dosage=="low_expressed"]))] <- "Low expression bins"
    CN_bins_df$clone_state[CN_bins_df$clone_state=="dosage_bins"] <- "No expression bins"
    table(CN_bins_df$clone_state)
    
    all_genes$patient <- sample_name
    CN_bins_df$patient <- sample_name
    write_rds(all_genes,"all_genes_seg_classification.rds")
    write_rds(CN_bins_df,"all_bins_seg_classification.rds")
    
    table(all_genes$gene_dosage)
    table(CN_bins_df$clone_state)
    all_genes |> 
      tidyplot(x = patient, color = gene_dosage) |> 
      add_barstack_absolute() |> adjust_x_axis_title("") |> adjust_y_axis_title("Gene count") |> adjust_legend_title("Gene dosage classification") |> add(theme(aspect.ratio = 2.5)) |> save_plot("All_gene_level_dosage_prop_seg.pdf",height = 6,width = 6,units = "in")
    
    CN_bins_df |> 
      tidyplot(x = patient, color = clone_state) |> 
      add_barstack_absolute() |> adjust_x_axis_title("") |> adjust_y_axis_title("CN bin count") |> adjust_legend_title("Bin dosage classification") |> add(theme(aspect.ratio = 2.5)) |>save_plot("All_bin_level_dosage_prop.pdf",height = 6,width = 6,units = "in")
    
    all_genes_subset <- all_genes[all_genes$gene_dosage %in% c("dosage-insensitive","dosage-sensitive"),c("Gene","gene_dosage")] 
    CN_bins_subset <- CN_bins_df[CN_bins_df$clone_state %in% c("Dosage-insensitive bins","Dosage-sensitive bins"),c("CN_bin","clone_state")]
    colnames(all_genes_subset) <- c("feature","gene_dosage")
    colnames(CN_bins_subset) <- c("feature","gene_dosage")
    CN_bins_subset$gene_dosage[CN_bins_subset$gene_dosage=="Dosage-insensitive bins"] <- "dosage-insensitive"
    CN_bins_subset$gene_dosage[CN_bins_subset$gene_dosage=="Dosage-sensitive bins"] <- "dosage-sensitive"
    all_genes_subset$category <- "Gene"
    CN_bins_subset$category <- "CN bin (220kb)"
    
    head(all_genes_subset)
    head(CN_bins_subset)
    
    rbind(all_genes_subset,CN_bins_subset) |> 
      tidyplot(x = category, color = gene_dosage) |> 
      add_barstack_relative(width = 0.4) |> adjust_x_axis_title("") |> adjust_y_axis_title("Proportion in mappable genes and CN bins") |> adjust_legend_title("Dosage effect")   |>save_plot("Sample_mappable_gene_dosage_prop_seg.pdf",height = 6,width = 6,units = "in")
    
    reverse_dosage_genes <- gene_seg_cor[gene_seg_cor$p_value<0.05 & gene_seg_cor$correlation<0 ,]
    reverse_dosage_genes <- reverse_dosage_genes[!is.na(reverse_dosage_genes$Gene),]
    
    write_rds(reverse_dosage_genes,"reverse_dosage_genes.rds")
    
    gene_seg_cor_subset <- gene_seg_cor[gene_seg_cor$Gene %in% all_genes$Gene[all_genes$gene_dosage %in% c("dosage-insensitive","dosage-sensitive")],]
    dim(gene_seg_cor_subset)
    gene_seg_cor_subset$dosage <- all_genes$gene_dosage[match(gene_seg_cor_subset$Gene,all_genes$Gene)]
    p3 <- ggplot(gene_seg_cor_subset, aes(y = delta_CN, x = dosage,color=dosage)) + geom_boxplot(position = position_dodge2(0.9),outlier.shape=NA)+labs(x="",y="Maximum copy number difference", color="Gene dosage")+geom_point(position=position_jitterdodge(),size=0.1)+theme_classic2()+theme(axis.text.x = element_text(size = 12, face = "bold",angle = 45,hjust = 1),axis.text.y = element_text(size = 16,face = "bold"),axis.title=element_text(size=16,face = "bold"),aspect.ratio = 1)+scale_color_manual(values = YZRY_color_set2)+stat_compare_means(aes(group = dosage), label = "p.signif",label.x = 1.5,method = "wilcox.test")
    pdf("CN_diffenrece_DSRvsDIR.pdf",height = 6,width = 6,useDingbats = F)
    print(p3)
    dev.off()
    
    gene_seg_cor_subset$delta_CN_integer <- gene_integer_cor$delta_CN[match(gene_seg_cor_subset$Gene,gene_integer_cor$Gene)]
    gene_seg_cor_subset$delta_CN_class <- gene_seg_cor_subset$delta_CN_integer
    gene_seg_cor_subset$delta_CN_class[gene_seg_cor_subset$delta_CN_integer>=4 & gene_seg_cor_subset$delta_CN_integer <8 ] <- "4-8"
    gene_seg_cor_subset$delta_CN_class[gene_seg_cor_subset$delta_CN_integer>=8] <- "8+"
    table(gene_seg_cor_subset$delta_CN_class)
    
    
    library(tidyplots)
    
    gene_seg_cor_subset |> 
      tidyplot(x = delta_CN_class, color = dosage) |> 
      add_barstack_absolute() |> adjust_x_axis_title("Subclonal copy number difference") |> adjust_y_axis_title("Gene count") |> adjust_legend_title("Gene dosage") |> save_plot("gene_dosage_count_by_CN_dif.pdf")
    
    gene_seg_cor_subset |> 
      tidyplot(x = delta_CN_class, color = dosage) |> 
      add_barstack_relative() |> adjust_x_axis_title("Subclonal copy number difference") |> adjust_y_axis_title("Gene count") |> adjust_legend_title("Gene dosage") |> save_plot("gene_dosage_prop_by_CN_dif.pdf")
    
    gene_seg_cor_subset2 <- gene_seg_cor_subset[gene_seg_cor_subset$dosage=="dosage-sensitive",]
    
    gene_seg_cor_subset2 |> 
      tidyplot(x = delta_CN_class, y = correlation, color = delta_CN_class) |> 
      add_data_points(size=0.1,alpha=0.25) |> 
      add_data_points_jitter(size=0.1,jitter_width = 0.4,alpha=0.2) |>
      add_mean_bar(alpha = 0.5) |> 
      add_sem_errorbar(alpha = 1) |> 
      adjust_x_axis_title("Subclonal copy number difference") |>
      adjust_y_axis_title("Correlation of gene expression and copy number") |>
      adjust_legend_title("") |> save_plot("gene_dosage_effects_by_CN_dif.pdf")
  }
  
  
  write_rds(exp_dds_bulk_long,"exp_dds_bulk_long.rds")
  write_rds(DNA_data_tumor,"DNA_data_tumor.rds")
  write_rds(RNA_cancer_matched,"RNA_cancer_matched.rds")
  

  
  ###########################subclone pair gene-dosage calculation
  gene_dosage_prop <- calculate_subclone_pair_gene_dosage3(DNA_data_tumor,RNA_cancer_matched,0.5,subclone_variable,rna_subclone_variable)
  
  rna_sig_result <- gene_dosage_prop$rna_sig_result
  rna_sig_result <- rna_sig_result[!is.na(rna_sig_result$CN_bins),]
  rna_freq <- data.frame(table(rna_sig_result$group_pair))
  rna_freq$Var1[rna_freq$Freq<5]
  #####remove low group_pairs
  rna_sig_result <- rna_sig_result[!rna_sig_result$group_pair %in% rna_freq$Var1[rna_freq$Freq<5],]
  library(dplyr)
  library(tidyr)
  RNA_subclone_pair_freq <- rna_sig_result %>%
    dplyr::count(group_pair, subclone_gene_dosage) %>%
    tidyr::complete(group_pair, subclone_gene_dosage, fill = list(n = 0))%>%
    dplyr::group_by(group_pair) %>%
    dplyr::mutate(
      frequency = n / sum(n)  # Frequency of each category within the group
    ) %>%
    ungroup() 
  
  
  RNA_subclone_pair_freq$patient <- sample_name
  
  RNA_subclone_pair_freq[RNA_subclone_pair_freq$subclone_gene_dosage=="Positive_dosage",] |> 
    tidyplot(x= patient ,y = frequency,color=patient) |> 
    add_data_points_jitter(size=0.5,jitter_width = 0.2,color="black") |>
    add_mean_bar(alpha = 0.5,width = 0.4) |> 
    add_sem_errorbar(color="black") |> 
    adjust_x_axis_title("") |>
    adjust_y_axis_title("Gene-dosage frequency") |> add(ylim(0,1)) |>
    adjust_legend_title("")|> add(ggtitle(paste0("Cis gene proportion: ",round(mean(RNA_subclone_pair_freq[RNA_subclone_pair_freq$subclone_gene_dosage=="Positive_dosage",]$frequency),2))))  |> save_plot("RNA_subclone_gene_dosage_seg_0.5_logFC.pdf")
  
  write_rds(RNA_subclone_pair_freq,"RNA_subclone_pair_freq.rds")
  
  
  ###############################cis/trans gene-dosage effects (local dosage effect)
  
  rna_sig_result <- gene_dosage_prop$rna_sig_result
  head(rna_sig_result)
  rna_sig_result <- rna_sig_result[!is.na(rna_sig_result$CN_bins),]
  rna_freq <- data.frame(table(rna_sig_result$group_pair))
  rna_freq$Var1[rna_freq$Freq<5]
  #####remove low group_pairs
  rna_sig_result <- rna_sig_result[!rna_sig_result$group_pair %in% rna_freq$Var1[rna_freq$Freq<5],]
  #########summarize count and proportion of genes being identified as positive/negative/no dosage
  #####get count
  summary_df <- rna_sig_result %>%
    group_by(gene, subclone_gene_dosage) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = subclone_gene_dosage, values_from = count, values_fill = 0)
  
  max(summary_df$Positive_dosage)
  
  #####get frequency
  summary_prop_df <- summary_df
  if (ncol(summary_prop_df)==3) {summary_prop_df[,2:3] <- summary_prop_df[,2:3]/length(unique(rna_sig_result$group_pair))} else {
    summary_prop_df[,2:4] <- summary_prop_df[,2:4]/length(unique(rna_sig_result$group_pair))
  }
  
 
  
  
  summary_prop_df <- data.frame(summary_prop_df)
  summary_prop_df$cis_to_reverse <- summary_prop_df$Positive_dosage-summary_prop_df$Negative_dosage
  summary_prop_df$cis_to_trans <- summary_prop_df$Positive_dosage-summary_prop_df$No_dosage
  summary_prop_df$trans_to_reverse <- summary_prop_df$No_dosage-summary_prop_df$Negative_dosage
  summary_prop_df <- summary_prop_df[order(-summary_prop_df$Positive_dosage),]
  head(summary_prop_df)
  
  cis_genes <- summary_prop_df$gene[summary_prop_df$Positive_dosage>0 & summary_prop_df$cis_to_reverse>0 & summary_prop_df$cis_to_trans>0]
  trans_genes <- summary_prop_df$gene[summary_prop_df$No_dosage>0 & summary_prop_df$trans_to_reverse>0 & summary_prop_df$cis_to_trans<0]
  reverse_genes <- summary_prop_df$gene[summary_prop_df$Negative_dosage>0 & summary_prop_df$cis_to_reverse<0 & summary_prop_df$trans_to_reverse<0]
  
  #####show top cis, trans, reverse genes
  summary_prop_df[order(-summary_prop_df$Positive_dosage),][1:20,]$gene
  summary_prop_df[order(-summary_prop_df$No_dosage),][1:20,]
  summary_prop_df[order(-summary_prop_df$Negative_dosage),][1:20,]
  
  summary_prop_df$dosage_class <- ifelse(summary_prop_df$gene %in% cis_genes,"Cis",ifelse(summary_prop_df$gene %in% trans_genes,"Trans",ifelse(summary_prop_df$gene %in% reverse_genes,"Reverse","Uncertain")))
  table(summary_prop_df$dosage_class)
  summary_prop_df$sample <- sample_name
  
  write_rds(summary_prop_df,"cis_trans_reverse_gene.rds")
  
  
} 
                                          
