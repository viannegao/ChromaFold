#library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome)
library(dplyr)
library(Rcpp)
library(base)
library(data.table)
library(dplyr)
library(tidyr)
library(GenomeInfoDb)
library(stats)
library(MASS)
options("scipen"=100, "digits"=4)

get_chr_sizes<-function(gen="Hsapiens",gen_ver="hg19",chrs=NULL){
  genome <- paste("BSgenome.", gen, ".UCSC.", gen_ver, sep = "")
  library(genome, character.only = TRUE)
  if (is.null(chrs)) {
    #get list of chromosomes
    chrs<-get_chrs(gen,gen_ver)
    chrs <- chrs[!chrs %in% c('chrY', 'chrM')]
  }
  sizes<-GenomeInfoDb::seqlengths(get(gen))[chrs]
  return(sizes)
}

get_chrs<-function(gen="Hsapiens",gen_ver="hg19"){
  genome <- paste("BSgenome.", gen, ".UCSC.", gen_ver, sep = "")
  library(genome, character.only = TRUE)
  chrs <- GenomeInfoDb::seqnames(get(gen))[
    vapply(GenomeInfoDb::seqnames(get(gen)), nchar, FUN.VALUE = 1) <= 5]
  chrs <- chrs[!chrs=='chrMT']
  return(chrs)
}

#' hic2hicdc
#'
#' This function converts a .hic file into a HiC-DC readable matrix file for
#' each chromosome
#'@import BSgenome
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@param hic_path path to the .hic file
#'@param bintolen_path path to the relevant feature file
#'@param gen name of the species: e.g., default Hsapiens
#'@param gen_ver genomic assembly version: e.g., default hg19
#'@param binsize binsize in bp, e.g., default 5000
#'@param Dthreshold maximum distance of interactions e.g., default 2000000
#'@param inter Interchromosomal interaction counts added as a covariate or not;
#'default FALSE
#'@param chrs select a subset of chromosomes' e.g., c("chrY","chrM"). Defaults
#'to all chromosomes in the genome.
#'@param output_path the path to the folder you want to place HiCDC-compatible
#'matrices into. Each file will have the .hic file preamble as its name prefix.
#'@param bin_type 'Bins-uniform' if uniformly binned by binsize in
#'bp, or 'Bins-RE-sites' if binned by number of
#'restriction enzyme fragment cutsites!
#'@return paths of a list of matrix files generated for each chromosome
#'readable by HiC-DC
#'@examples hic2hicdc(hic_path="~/HiCDCPlus_files/mm10/cb_ctrl1_map30_raw.hic",
#' bintolen_path=
#' "~/HiCDCPlus_files/mm10_features/features_combo_5kb_bintolen.rds", 
#' straw_path="~/Downloads/straw-R.cpp",
#' gen="Mmusculus",gen_ver="mm10",binsize=5000,
#' inter=FALSE,Dthreshold=2000000,chrs=c("chr1"))
#'@export

hic2hicdc <-
  function(hic_path,
           bintolen_path,
           straw_path,
           gen = "Hsapiens",
           gen_ver = "hg19",
           binsize = 5000,
           Dthreshold = 2000000,
           bin_type = "Bins-uniform",
           inter = FALSE,
           chrs = NULL,
           output_path) {
    options(scipen = 9999)
    #get straw
    suppressWarnings({Rcpp::sourceCpp(path.expand(straw_path))})
    #get list of chromosomes and their sizes
    if (is.null(chrs)) {
      chrom <- get_chrs(gen,gen_ver)
    } else{
      chrom <- chrs
    }
    genome_chromSizes <- data.frame(chr = chrom,
                                    size = get_chr_sizes(gen,gen_ver,chrom))
    #chrs <- gsub('chr', '', chrom) ######################### change chrom format
    chrs <- chrom
    #get bintolen
    if (base::grepl('.txt',bintolen_path,ignore.case = TRUE)){
      bintolen <- data.table::fread(bintolen_path,
                                    sep = "\t",
                                    header = TRUE,
                                    stringsAsFactors = FALSE
      )
    }else if (base::grepl('.rds',bintolen_path,ignore.case = TRUE)){
      bintolen<-readRDS(bintolen_path)
    }
    #prepare bintolen
    bintolen <- bintolen %>% tidyr::separate("bins",
                                             c("chr", "start", "end"),
                                             sep = "-",
                                             remove = FALSE,
                                             convert = TRUE,
                                             extra = "drop",
                                             fill = "warn"
    )
    #modify output path (get string between last / and .)
    hic_preamble <- rev(regmatches(hic_path,regexec('(.*?)/(.*?)\\.',
                                                    hic_path))[[1]])[1]
    output_path <- paste0(output_path, hic_preamble)
    filepaths = NULL
    for (chr in chrs) {
      print(paste0("Processing chromosome: ", chr))
      #process bintolen for the chr
      #filter bintolen to the chr
      bintolen_chr <- bintolen[bintolen$chr == chr, ]
      chrom_size <- genome_chromSizes[
        genome_chromSizes$chr == chr, ]$size
      #get all possible start values for the chr and complete bintolen
      starts <- data.frame(start = seq(1,
                                       genome_chromSizes[genome_chromSizes$chr == chr, ]$size,
                                       binsize))
      bintolen_chr <- dplyr::left_join(starts, bintolen_chr)
      rm(starts)
      bintolen_chr<-bintolen_chr%>%tidyr::replace_na(list(gc=0,map=0,len=0,
                                                          chr=chr))
      bintolen_chr<-bintolen_chr%>%dplyr::mutate(
        start=floor(.data$start/1000)*1000,
        end=pmin(.data$start+ binsize, chrom_size),
        mid=(.data$start + .data$end) / 2,
        bin=paste0(.data$chr,"-",.data$start,"-",.data$end),
        bins=NULL)
      #get all bin combinations within distance threshold
      print("Generating the features matrix")
      numbins<-nrow(bintolen_chr)
      maxbins<-ceiling(Dthreshold/binsize)
      index1<-unlist(sapply(1:numbins,
                            function(x) rep(x,min(maxbins+1,numbins-x+1))))
      index2<-unlist(sapply(1:numbins,
                            function(x) seq(x,min(x+maxbins,numbins),1)))
      bintolen_chr<-dplyr::bind_cols(
        (bintolen_chr%>%dplyr::rename_all(
          function(x) paste0(x,"I")))[index1,],
        (bintolen_chr%>%dplyr::rename_all(
          function(x) paste0(x,"J")))[index2,]
      )
      rm(index1,index2)
      bintolen_chr<-bintolen_chr%>%dplyr::mutate(D=abs(.data$midI-.data$midJ))
      bintolen_chr<-bintolen_chr%>%dplyr::filter(D<=Dthreshold)
      cols_bintolen_chr<-c("binI","binJ","chrI","chrJ","startI","startJ",
                           "endI","endJ","gcI","gcJ","mapI","mapJ","lenI","lenJ","D")
      bintolen_chr<-bintolen_chr[cols_bintolen_chr]
      gc()
      print("Generated. Incorporating counts.")
      count_matrix <- straw(
        norm = 'NONE',
        fname = path.expand(hic_path),
        binsize = binsize,
        chr1loc = chr,
        chr2loc = chr,
        unit = 'BP')
      colnames(count_matrix) <- c("startI", "startJ", "counts")
      count_matrix <-count_matrix %>% dplyr::mutate(
        chrI = chr,
        chrJ = chr,
        startI = .data$startI,
        startJ = .data$startJ
      )
      bintolen_chr <- dplyr::left_join(bintolen_chr,
                                       count_matrix) %>%
        tidyr::replace_na(list(counts = 0))
      rm(count_matrix)
      if (inter) {
        print("Incorporated. Incorporating interchromosomal counts.")
        count_matrix <- NULL
        for (chr2 in chrs) {
          if (chr2 != chr) {
            count_matrix_add <- unique(
              straw(
                norm = 'NONE',
                fname = path.expand(hic_path),
                binsize = binsize,
                chr1loc = chr,
                chr2loc = chr2,
                unit = 'BP'))
            xmax <- max(count_matrix_add[, 1])
            ymax <- max(count_matrix_add[, 2])
            if (abs(chrom_size - xmax) <= abs(chrom_size - ymax)) {
              count_matrix_add <- count_matrix_add[, c(1, 3)]
            } else {
              count_matrix_add <- count_matrix_add[, c(2, 3)]
              
            }
            colnames(count_matrix_add) <- c("start", "inter")
            count_matrix <- dplyr::bind_rows(count_matrix, count_matrix_add)
          }
        }
        rm(count_matrix_add)
        count_matrix <- count_matrix %>% dplyr::group_by(.data$start) %>%
          dplyr::summarize(inter = sum(.data$inter)) %>%
          dplyr::mutate(start = .data$start)
        bintolen_chr <-
          dplyr::left_join(bintolen_chr,
                           count_matrix %>%
                             dplyr::rename("startI" = "start", "interI" = "inter")) %>%
          tidyr::replace_na(list(interI = 0))
        bintolen_chr <-
          dplyr::left_join(bintolen_chr,
                           count_matrix %>%
                             dplyr::rename("startJ" = "start", "interJ" = "inter")) %>%
          tidyr::replace_na(list(interJ = 0))}
      print("Incorporated. Printing to file.")
      filepath <- paste0(output_path, '_', binsize/1000, 'kb_chr', chr, '_matrix.txt')
      data.table::fwrite(bintolen_chr,filepath,
                         quote = FALSE,
                         sep = "\t",
                         row.names = FALSE)
      rm(bintolen_chr)
      gc()
      print(paste0("Chromosome: ", chr, " processed."))
      filepaths <- append(filepaths, paste0(filepath))
    }
    return(filepaths)
  }

hicdcplus_model=function(file="matrix.rds",
                         bin_type = 'Bins-uniform',
                         gen = "Hsapiens",
                         gen_ver = "hg19",
                         covariates = NULL,
                         distance_type = 'spline',
                         output_path = "K562_combined",
                         chrs = "chr1",
                         df = 6,
                         Dmin = 0,
                         Dmax = 2e6,
                         ssize = 0.01,
                         model_file = FALSE){
  
  remove_outliers_nb <- function(dat, mod) {
    mu <- stats::predict(mod, newdata = dat, type = 'response')
    dat <- dat %>% dplyr::mutate(q = stats::qnbinom(0.975,size = mod$theta, mu = mu))
    new.dat <- dat %>% dplyr::filter(.data$counts <= .data$q) %>% dplyr::select(-.data$q)
    return(new.dat)
  }
  
  glm.nb.trycatch<- function(model.formula,data,theta.init=0.0001){
    model<-tryCatch({
      # try fitting distribution as requested by user
      strip_glm(MASS::glm.nb( 
        model.formula, 
        data))
    }, error = function(e) {
      temp.model <- MASS::glm.nb( 
        model.formula, 
        data,
        init.theta = theta.init)
      return(temp.model)
    }
    )
    return(model)
  }
  
  GLM_nb <- function(data,df,bdpts,covariates,distance_type) {
    if (distance_type == 'spline') {
      model.formula<-stats::as.formula(paste0("counts ~",paste(covariates, collapse = " + ")," + splines::bs(D,df=",df,",Boundary.knots = c(",paste(bdpts, collapse = ","),"))"))
    } else {
      model.formula<=stats::as.formula(paste0("counts ~", paste(covariates, collapse = " + "), " + logD"))
    }
    # Remove outliers
    print(2)
    tmp_data <<- data
    mod<-glm.nb.trycatch(model.formula,data)
    new.dat <- remove_outliers_nb(data, mod)
    print(3)
    # Refit the model
    mod<-glm.nb.trycatch(model.formula,new.dat,theta.init = mod$theta)
    print(4)
    return(mod)
  }
  
  Normalize <- function(x) {
    x.new <- ifelse(is.finite(x), (x -
                                     base::mean(x[is.finite(x)])) / stats::sd(x[is.finite(x)]), NA)
    return(x.new)
  }
  
  transform.vec <- function(x, y)
    as.vector(ifelse(x==0 | y==0, 0, Normalize(log(x * y))))
  
  strip_glm = function(m1) {
    m1$data <- NULL
    m1$y <- NULL
    m1$linear.predictors <- NULL
    m1$weights <- NULL
    m1$fitted.values <- NULL
    m1$model <- NULL
    m1$prior.weights <- NULL
    m1$residuals <- NULL
    m1$effects <- NULL
    m1$qr$qr <- NULL
    m1$optim <- NULL
    m1$start <- NULL
    return(m1)
  }
  
  # Main body of algorithm
  options("scipen" = 9999, "digits" = 4)
  
  #set default covariates if need be
  if (is.null(covariates)) {
    covariates <- c("gc", "map", "len")
  }
  if (bin_type=='Bins-RE-sites' & is.null(covariates)) {
    covariates <- c("gc", "map", "width")
  }
  
  #set default chrs if need be
  if (is.null(chrs)) {
    #get list of chromosomes
    chrs<-get_chrs(gen,gen_ver)
    chrs <- chrs[!chrs %in% c('chrY', 'chrM')]
  }
  
  new.x <- data.table::fread(paste0(file)) # reading in data
  new.x <- new.x %>% dplyr::filter(.data$D >= Dmin & .data$D <= Dmax)
  binsize<-min(new.x$D[new.x$D>Dmin&new.x$D%%1000==0])-Dmin
  bdpts <- range(new.x$D)
  ix <- rep(TRUE, nrow(new.x))
  for (cov in covariates) {
    result <- tryCatch({
      x <- as.numeric(unlist(new.x[, paste0(cov, 'I'), with = FALSE]))
      y <- as.numeric(unlist(new.x[, paste0(cov, 'J'), with = FALSE]))
      list(x = x, y = y)}, error = function(e) {
        x <- as.numeric(unlist(new.x[, paste0(cov, 'I')]))
        y <- as.numeric(unlist(new.x[, paste0(cov, 'J')]))
        return(list(x = x, y = y))})
    new.x[, cov] <- transform.vec(result$x, result$y)
    ix <- ix & is.finite(unlist(new.x[, cov, with=FALSE])) #& !is.na(unlist(new.x[, cov, with=FALSE]))
    rm(result)
  }
  dat <- new.x[ix,]
  dat <- dat[dat$gcI>0 & dat$gcJ>0 & dat$mapI>0 & dat$mapJ>0,]
  rm(ix)
  set.seed(1010)
  #stratified sampling for training data
  zeroPairedBins <- dat %>% dplyr::filter(.data$counts == 0)
  countPairedBins <- dat %>% dplyr::filter(.data$counts > 0)
  Intervals <- seq(min(dat$D), to = max(dat$D), by = binsize)
  rm(dat)
  cls.counts <- findInterval(countPairedBins$D, Intervals,
                             rightmost.closed = TRUE)
  cls.zeros <- findInterval(zeroPairedBins$D, Intervals,
                            rightmost.closed = TRUE)
  bins.counts <-split(seq_len(nrow(countPairedBins)), cls.counts)
  bins.zeros <- split(seq_len(nrow(zeroPairedBins)), cls.zeros)
  idx.counts <- unlist(lapply(bins.counts, function(x) {
    sample(x, floor(length(x) * ssize), replace = FALSE)}))
  idx.zeros <- unlist(lapply(bins.zeros, function(x) {
    sample(x, floor(length(x) * ssize), replace = FALSE)}))
  dat <- dplyr::bind_rows(countPairedBins[idx.counts, ],zeroPairedBins[idx.zeros, ])
  rm(cls.counts,cls.zeros,bins.counts,bins.zeros,idx.counts,idx.zeros,Intervals)
  print(0)
  
  stime <- system.time({fit <- GLM_nb(dat, df, bdpts, covariates, distance_type)})[3] / 60
  #normal model
  print(1)
  
  mu <- stats::predict(fit,newdata = new.x,dispersion = fit$theta ** (-1),type = "response")
  sdev <- sqrt(mu + mu ** 2 / fit$theta)
  ix <- !is.finite(mu)
  new.x$pvalue <- 1
  new.x$pvalue[!ix] <- stats::pnbinom(q = new.x$counts[!ix] - 1,size = fit$theta,mu = mu[!ix],lower.tail = FALSE)
  new.x <- new.x %>% dplyr::mutate(mu = mu, sdev = sdev)
  new.x$qvalue <- 1
  new.x$qvalue[!ix] <- stats::p.adjust(new.x$pvalue[!ix],
                                       method = 'fdr') 

  new.x <- new.x %>%
    dplyr::mutate(
      normcounts = .data$counts / .data$mu,
      zvalue = (.data$counts - .data$mu) / .data$sdev ) %>%
    tidyr::replace_na(list(mu = 0,pvalue = 1,qvalue = 1,sdev = 0,zvalue = 0,
                           normcounts = 0))
  default_columns_prefix <- do.call(paste0,
                                    expand.grid(c('bin', 'bin_id', 'chr', 'start', 'end'), c('I', 'J')))
  columns_add <- sort(do.call(paste0,
                              expand.grid(covariates, c('I', 'J'))))
  default_columns_suffix <- c('D','counts','pvalue','qvalue','normcounts',
                              'zvalue',
                              'mu',
                              'sdev')
  column_set<-unique(c(
    default_columns_prefix,
    columns_add,
    default_columns_suffix))
  column_set<-column_set[column_set%in%colnames(new.x)]
  new.x <- new.x[, column_set, with=FALSE]
  if (model_file){
    fitpath<-paste0(output_path, "_Model-Fit_on_", chrs, ".rds")
    saveRDS(fit,  file = fitpath)
  }
  saveRDS(new.x,paste0(output_path, "_Results_on_", chrs, ".rds"))
  return(paste0("Saved results: ", output_path, "_Results_on_", chrs, ".rds"))
}

#' hicdc2hic
#'
#' This function converts HiC-DC Plus significant counts (uniformly binned)
#' for a set of chromosomes back into a .hic file
#'@import BSgenome splines
#'@importFrom dplyr %>%
#'@importFrom rlang .data
#'@param input_path the path to the folder and name prefix where
#'HiCDC-processed matrix text files for each chromosome reside.
#'Should have the same value as \code{output_path}
#'for \code{HiCDCPlus}.
#'@param jarfile path to the Juicebox .jar file. Juicebox command line tools
#'can be downloaded from:
#'\url{https://github.com/aidenlab/juicer/wiki/Download}
#'@param binsize binsize in bp, e.g., default 5000
#'@param gen name of the species: e.g., default \code{'Hsapiens'}
#'@param gen_ver genomic assembly version: e.g., default \code{'hg19'}
#'@param mode What to put to the .hic file
#'as score. Allowable options are: 'pvalue' for -log10 significance p-value,
#''qvalue' for -log10 FDR corrected p-value,
#''normcounts' for raw counts/expected counts, and
#''zvalue' for standardized counts and 'raw' to pass-through
#'raw counts. Defaults to 'zvalue'.
#'@param output_file the path to the .hic file
#'@param chrs select a subset of chromosomes' e.g.,
#'c("chr21","chr22"). Defaults to all chromosomes (except Y and M)
#'in the genome specified.
#'@return path of the.hic file
#'@examples hicdc2hic(input_path="~/Downloads/",
#'jarfile="~/Downloads/juicer_tools_1.13.02.jar",
#'binsize=5000,
#'output_file="~/Downloads/out.hic",chrs=c("chr22"))
#'@export

hicdc2hic<-function(input_path,jarfile,
                    binsize=5000,gen="Hsapiens",gen_ver="hg19",mode='zvalue',
                    output_file="~/Downloads/out.hic",chrs=sprintf("chr%s",c(1:22,"X"))){
  options("scipen" = 9999, "warn" = -1)
  
  #create path to the pre input
  preinputpath<-path.expand(paste0(input_path,'_normcounts.txt'))
  preoutputpath<-path.expand(output_file)
  #get list of chromosomes if not specified
  if (is.null(chrs)) {
    #get list of chromosomes
    chrs<-get_chrs(gen,gen_ver)
    chrs <- chrs[!chrs %in% c('chrY', 'chrM')]
  }
  for (chr in chrs){
    #read the HiC-DC result
    input_path_chr<-paste0(input_path, "_Results_on_", chr, ".rds")
    norm_data<-readRDS(path.expand(input_path_chr))
    #convert to "short with score format" with scores as normalized counts
    if (mode=='normcounts'){
      norm_data<-norm_data%>%dplyr::mutate(score=.data$normcounts)
    }
    if (mode=='pvalue'){
      norm_data<-norm_data%>%dplyr::mutate(score=pmax(-log10(.data$pvalue),0))
    }
    if (mode=='qvalue'){
      norm_data<-norm_data%>%dplyr::mutate(score=pmax(-log10(.data$qvalue),0))
      norm_data$score[!is.finite(norm_data$score)] <- 1
    }
    if (mode=='zvalue'){
      norm_data<-norm_data%>%dplyr::mutate(score=.data$zvalue)
    }
    if(mode=='raw'){
      norm_data<-norm_data%>%dplyr::mutate(score=.data$counts)
    }
    norm_data<-norm_data%>%dplyr::filter(.data$score>-2147400000)%>%
      dplyr::mutate(strI=0,strJ=0,fragI=0,fragJ=1,
                    startI=.data$startI+binsize/2,
                    startJ=.data$startJ+binsize/2, chrI=gsub("chr","",chrI),
                    chrJ=gsub("chr","",chrJ))%>%
      dplyr::select(.data$strI,.data$chrI,.data$startI,.data$fragI,.data$strJ,
                    .data$chrJ,.data$startJ,.data$fragJ,.data$score)
    #dump output to file
    data.table::fwrite(x=norm_data,file=preinputpath,append=TRUE,quote=FALSE,
                       sep="\t",row.names=FALSE,col.names=FALSE)
    print(paste0("Output generated for chr: ", chr))
  }
  #run pre
  if (mode=="zvalue"){
    #make sure negative values get processed
   system2("java", args = c('-Xmx8g', '-jar',
                            path.expand(jarfile), "pre", "-v", "-d", "-n",
                            "-r", binsize,
                            "-m", -2147400000,
                            path.expand(preinputpath), path.expand(preoutputpath),
                            gen_ver))
 } else {
   system2("java", args = c('-Xmx8g', '-jar',
                            path.expand(jarfile), "pre", "-v", "-d",
                            "-r", binsize, path.expand(preinputpath),
                            path.expand(preoutputpath),
                            gen_ver))
 }

  #system2('java',args=c('-Xmx8g', '-jar' ,path.expand(jarfile),
  #'pre','-v','-d',path.expand(preinputpath),path.expand(preoutputpath),
  #gen_ver))
  #remove file
  #system2('rm',args=path.expand(preinputpath))
}









