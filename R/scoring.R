#' Keep one pathway per ligand-receptor pair
#'
#' @param pairs         A ligand-receptor table such as output by \code{pValuesLR} and \code{naiveBayesLR}.
#' @return A table with the same structure as \code{pairs} but reduced to only report one pathway by ligand-receptor pair.
#' Depending on how the pairs were scored, \code{pValuesLR} or \code{naiveBayesLR}, the pathway with
#' smallest P-value or largest log-likelihood ratio is selected.
#'
#' During the execution of \code{pValuesLR} or \code{naiveBayesLR}, ligand-receptor pairs are evaluated
#' in relation with pathways that allow checking receptor downstream correlations. It is thus possible
#' that several pathways are reported for a same pair.
#' @export
#' @examples
#' \dontrun{
#' data(sdc,package="BulkSignalR")
#' sample.types <- rep("tumor",ncol(sdc))
#' sample.types[grep("^N",names(sdc),perl=TRUE)] <- "normal"
#' ds <- prepareDataset(sdc,sample.types)
#' ds <- learnParameters(ds)
#'
#' # since model training takes time and it can be reused, it is advisable to save the model
#' save(ds,file="...",compress="bzip2")
#'
#' ds.LR <- getCorrelatedLR(ds)
#' ds.LR <- checkReceptorSignaling(ds,ds.LR)
#' pairs.p <- pValuesLR(ds.LR,ds$param)
#' pairs.red <- reduceToBestPathway(pairs.p)
#' }
reduceToBestPathway <- function(pairs){

  t <- sum(c("pval","LLR") %in% names(pairs))
  if (t ==2)
    stop("Both P-values and LLR present in the ligand-receptor table")
  if (t == 0)
    stop("P-values or LLR must be available in the ligand-receptor table")

  if ("pval" %in% names(pairs))
    pairs %>% dplyr::group_by(L,R) %>% dplyr::filter(pval==min(pval))
  else
    pairs %>% dplyr::group_by(L,R) %>% dplyr::filter(LLR==max(LLR))

} # reduceToBestPathway


#' Ligand-receptor pair P-value estimation
#'
#' Estimate the P-value of each ligand-receptor pair in the \code{merged.pairs} slot of the
#' ligand-receptor analysis \code{lr}.
#'
#' @param lr         A ligand-receptor analysis such as output by \code{checkReceptorSignaling}.
#' @param param         A list containing the statistical model parameters.
#' @param rank.p        A number between 0 and 1 defining the rank of the last considered target genes.
#' @param signed        A logical indicating whether correlations should be considered with their signs.
#' @param fdr.proc      The procedure for adjusting P-values accoring to \code{\link[multtest]{mt.rawp2adjp}}.
#' @param use.gamma     A logical indicating whther the Gamma distribution should be used instead of the Gaussian.
#' @param best.pw.only  A logical indicating that only the pathway with best P-value should be reported for each ligand-receptor pair.
#' @return A table reporting every ligand-receptor pair with P-values and Q-values (adjusted P-values).
#' @export
#' @examples
#' \dontrun{
#' data(sdc,package="BulkSignalR")
#' sample.types <- rep("tumor",ncol(sdc))
#' sample.types[grep("^N",names(sdc),perl=TRUE)] <- "normal"
#' ds <- prepareDataset(sdc,sample.types)
#' ds <- learnParameters(ds)
#'
#' # since model training takes time and it can be reused, it is advisable to save the model
#' save(ds,file="...",compress="bzip2")
#'
#' ds.LR <- getCorrelatedLR(ds)
#' ds.LR <- checkReceptorSignaling(ds,ds.LR)
#' pairs.p <- pValuesLR(ds.LR,ds$param)
#' }
pValuesLR <- function(lr,param,rank.p=0.75,signed=TRUE,fdr.proc=c("BH","Bonferroni","Holm","Hochberg","SidakSS","SidakSD","BY","ABH","TSBH"),
                      use.gamma=FALSE,best.pw.only=FALSE){

  if (use.gamma && !signed)
    stop("A gamma model with squared correlations is not available")
  if (rank.p<0 || rank.p>1)
    stop("rank.p must lie in [0;1]")
  fdr.proc <- match.arg(fdr.proc)
  pairs <- lr$merged.pairs

  res <- NULL
  for (i in 1:nrow(pairs)){
    pwid <- unlist(strsplit(pairs$pwid[i],split="\\|"))
    pwname <- unlist(strsplit(pairs$pwname[i],split="\\|"))
    tg <- unlist(strsplit(pairs$target.genes[i],split="\\|"))
    spear <- unlist(strsplit(pairs$all.corr[i],split="\\|"))
    len <- as.numeric(unlist(strsplit(pairs$len[i],split="\\|")))

    if (use.gamma)
      p.lr <- 1-stats::pgamma(1+pairs$corr[i],shape=param$LR.0$gamma$k,scale=param$LR.0$gamma$theta)
    else
      p.lr <- 1-stats::pnorm(pairs$corr[i],param$LR.0$norm$mu,param$LR.0$norm$sigma)

    for (k in 1:length(len)){
      spears <- as.numeric(strsplit(spear[k],split=";")[[1]])
      r <- min(max(1,trunc(rank.p*len[k])),len[k])
      if (signed){
        rank.corr <- spears[r]
        if (use.gamma)
          p.rt <- stats::pbinom(r-1,len[k],stats::pgamma(1+rank.corr,shape=param$RT.0$gamma$k,scale=param$RT.0$gamma$theta))
        else
          p.rt <- stats::pbinom(r-1,len[k],stats::pnorm(rank.corr,param$RT.0$norm$mu,param$RT.0$norm$sigma))
      }
      else{
        z.rt <- (spears-param$RT.0$norm$mu)/param$RT.0$norm$sigma
        z.rt.sq <- z.rt**2
        o <- order(z.rt.sq)
        rank.corr <- spears[o][r]
        p.rt <- stats::pbinom(r-1,len[k],stats::pchisq(z.rt.sq[o][r],df=1))
      }
      res <- rbind(res,data.frame(pairs[i,c("L","R")],pw.id=pwid[k],pw.name=pwname[k],pval=p.lr*p.rt,LR.corr=pairs$corr[i],LR.pval=p.lr,
                                  rank=r,len=len[k],rank.corr=rank.corr,RT.pval=p.rt,target.genes=tg[k],all.corr=spear[k],
                                  stringsAsFactors=FALSE))
    }
  }

  rawp <- res$pval
  adj <- multtest::mt.rawp2adjp(rawp,fdr.proc)
  res$qval <- adj$adjp[order(adj$index),fdr.proc]

  if (best.pw.only)
    reduceToBestPathway(res)
  else
    res

} # pValuesLR


#' Ligand-receptor pair LLR estimation
#'
#' Estimate the log-likelihood ratio of every ligand-receptor pair in the \code{merged.pairs} slot of the
#' ligand-receptor analysis \code{lr}.
#'
#' @param lr         A ligand-receptor analysis such as output by \code{checkReceptorSignaling}.
#' @param param         A list containing the statistical model parameters.
#' @param rank.p        A number between 0 and 1 defining the rank of the last considered target genes.
#' @param signed        A logical indicating whether correlations should be considered with their signs.
#' @param use.gamma     A logical indicating whther the Gamma distribution should be used instead of the Gaussian.
#' @param best.pw.only  A logical indicating that only the pathway with best P-value should be reported for each ligand-receptor pair.
#' @return A table reporting every ligand-receptor pair with P-values and Q-values (adjusted P-values).
#' @export
#' @examples
#' \dontrun{
#' data(sdc,package="BulkSignalR")
#' sample.types <- rep("tumor",ncol(sdc))
#' sample.types[grep("^N",names(sdc),perl=TRUE)] <- "normal"
#' ds <- prepareDataset(sdc,sample.types)
#' ds <- learnParameters(ds)
#'
#' # since model training takes time and it can be reused, it is advisable to save the model
#' save(ds,file="...",compress="bzip2")
#'
#' ds.LR <- getCorrelatedLR(ds)
#' ds.LR <- checkReceptorSignaling(ds,ds.LR)
#' pairs.llr <- naiveBayesLR(ds.LR,ds$param)
#' }
naiveBayesLR <- function(lr,param,rank.p=0.75,signed=TRUE,use.gamma=FALSE,best.pw.only=FALSE){

  if (use.gamma && !signed)
    stop("A gamma model with squared correlations is not available")
  if (!is.null(param$learn.altern) && !param$learn.altern)
    stop("No alternative distribution was learnt")
  pairs <- lr$merged.pairs

  eps.m <- sqrt(.Machine$double.eps)
  res <- NULL
  for (i in 1:nrow(pairs)) {
    pwid <- unlist(strsplit(pairs$pwid[i],split="\\|"))
    pwname <- unlist(strsplit(pairs$pwname[i],split="\\|"))
    tg <- unlist(strsplit(pairs$target.genes[i],split="\\|"))
    spear <- unlist(strsplit(pairs$all.corr[i],split="\\|"))
    len <- as.numeric(unlist(strsplit(pairs$len[i],split="\\|")))

    if (use.gamma)
      lr.score <- stats::dgamma(1+pairs$corr[i],shape=param$LR.1$gamma$k,scale=param$LR.1$gamma$theta,log=TRUE)-stats::dgamma(1+pairs$corr[i],shape=param$LR.0$gamma$k,scale=param$LR.0$gamma$theta,log=TRUE)
    else
      lr.score <- stats::dnorm(pairs$corr[i],param$LR.1$norm$mu,param$LR.1$norm$sigma,log=TRUE)-stats::dnorm(pairs$corr[i],param$LR.0$norm$mu,param$LR.0$norm$sigma,log=TRUE)

    for (k in 1:length(len)){
      spears <- as.numeric(strsplit(spear[k],split=";")[[1]])
      r <- min(max(1,trunc(rank.p*len[k])),len[k])
      if (signed){
        rank.corr <- spears[r]
        if (use.gamma){
          p.0 <- stats::pgamma(1+rank.corr,shape=param$RT.0$gamma$k,scale=param$RT.0$gamma$theta)
          p.1 <- stats::pgamma(1+rank.corr,shape=param$RT.1$gamma$k,scale=param$RT.1$gamma$theta)
        }
        else{
          p.0 <- stats::pnorm(rank.corr,param$RT.0$norm$mu,param$RT.0$norm$sigma)
          p.1 <- stats::pnorm(rank.corr,param$RT.1$norm$mu,param$RT.1$norm$sigma)
        }
      }
      else{
        z.rt <- (spears-param$RT.0$norm$mu)/param$RT.0$norm$sigma
        z.rt.sq <- z.rt**2
        o <- order(z.rt.sq)
        rank.corr <- spears[o][r]
        p.0 <- stats::pchisq(z.rt.sq[o][r],df=1)

        z.1 <- (spears-param$RT.1$norm$mu)/param$RT.1$norm$sigma
        z.1.sq <- z.1**2
        p.1 <- stats::pchisq(z.1.sq[o][r],df=1)
      }
      if (p.0+eps.m<=1)
        d.0 <- (stats::pbinom(r-1,len[k],p.0+eps.m)-stats::pbinom(r-1,len[k],p.0))/eps.m
      else
        d.0 <- -(stats::pbinom(r-1,len[k],p.0-eps.m)-stats::pbinom(r-1,len[k],p.0))/eps.m
      if (p.1+eps.m<=1)
        d.1 <- (stats::pbinom(r-1,len[k],p.1+eps.m)-stats::pbinom(r-1,len[k],p.1))/eps.m
      else
        d.1 <- -(stats::pbinom(r-1,len[k],p.1-eps.m)-stats::pbinom(r-1,len[k],p.1))/eps.m
      rt.score <- log(d.1/d.0)
      res <- rbind(res,data.frame(pairs[i,c("L","R")],pw.id=pwid[k],pw.name=pwname[k],LLR=(lr.score+rt.score)/log(10),LR.corr=pairs$corr[i],LR.LLR=lr.score,
                                  rank=r,len=len[k],rank.corr=rank.corr,RT.LLR=rt.score,target.genes=tg[k],all.corr=spear[k],
                                  stringsAsFactors=FALSE))
    }
  }

  if (best.pw.only)
    reduceToBestPathway(res)
  else
    res

} # naiveBayesLR


#' Aggregate the ligands of a same receptor
#'
#' Simplifies a ligand-receptor table to focus on the receptors.
#'
#' @param pairs         A ligand-receptor table such as output by \code{pValuesLR} and \code{naiveBayesLR}.
#' @return A table with the same structure as \code{pairs} but reduced to only report one row per receptor.
#' All the ligands are combined in a semi-colon-separated list surounded by curly braces.
#'
#' To limit the output table to a single line per receptor, this function also selects the pathway
#' with the smallest P-value, respectively the largest log-likelihood ratio. The reported list of target genes
#' and their respective correlations are those of this single selected pathway.
#' @export
#' @examples
#' \dontrun{
#' data(sdc,package="BulkSignalR")
#' sample.types <- rep("tumor",ncol(sdc))
#' sample.types[grep("^N",names(sdc),perl=TRUE)] <- "normal"
#' ds <- prepareDataset(sdc,sample.types)
#' ds <- learnParameters(ds)
#'
#' # since model training takes time and it can be reused, it is advisable to save the model
#' save(ds,file="...",compress="bzip2")
#'
#' ds.LR <- getCorrelatedLR(ds)
#' ds.LR <- checkReceptorSignaling(ds,ds.LR)
#' pairs.p <- pValuesLR(ds.LR,ds$param)
#' pairs.recept <- reduceToReceptor(pairs.p)
#' }
reduceToReceptor <- function(pairs){

  t <- sum(c("pval","LLR") %in% names(pairs))
  if (t ==2)
    stop("Both P-values and LLR present in LR table")
  if (t == 0)
    stop("P-values or LLR must be available in LR table")

  red.mode.L <- (sum(regexpr("^{",pairs$L,perl=TRUE)==1)==nrow(pairs) && sum(regexpr("}$",pairs$L,perl=TRUE)==nchar(pairs$L))==nrow(pairs))
  if (red.mode.L)
    stop("Already reduced to receptor")

  # best LLR or pval per receptor/pathway
  if ("pval" %in% names(pairs))
    best <- pairs %>% dplyr::group_by(R,pw.id) %>% dplyr::filter(pval==min(pval))
  else
    best <- pairs %>% dplyr::group_by(R,pw.id) %>% dplyr::filter(LLR==max(LLR))

  # list of ligands for each receptor/pathway combination
  ligands <- pairs %>% dplyr::group_by(R,pw.id) %>% dplyr::summarize(Ls=paste0('{',paste(L,collapse=';'),'}'))

  # replace ligands by ligand lists
  rpw <- paste(best$R,best$pw.id,sep="|")
  rpwL <- stats::setNames(ligands$Ls,paste(ligands$R,ligands$pw.id,sep="|"))
  best$L <- rpwL[rpw]

  best

} # reduceToReceptor


#' Aggregate the receptors of a same ligand
#'
#' Simplifies a ligand-receptor table to focus on the ligands.
#'
#' @param pairs         A ligand-receptor table such as output by \code{pValuesLR} and \code{naiveBayesLR}.
#' @return A table with the same structure as \code{pairs} but reduced to only report one row per ligand.
#' All the receptors are combined in a semi-colon-separated list surounded by curly braces.
#'
#' To limit the output table to a single line per ligand, this function also selects the pathway
#' with the smallest P-value, respectively the largest log-likelihood ratio. The reported list of target genes
#' and their respective correlations are those of this single selected pathway.
#' @export
#' @examples
#' \dontrun{
#' data(sdc,package="BulkSignalR")
#' sample.types <- rep("tumor",ncol(sdc))
#' sample.types[grep("^N",names(sdc),perl=TRUE)] <- "normal"
#' ds <- prepareDataset(sdc,sample.types)
#' ds <- learnParameters(ds)
#'
#' # since model training takes time and it can be reused, it is advisable to save the model
#' save(ds,file="...",compress="bzip2")
#'
#' ds.LR <- getCorrelatedLR(ds)
#' ds.LR <- checkReceptorSignaling(ds,ds.LR)
#' pairs.p <- pValuesLR(ds.LR,ds$param)
#' pairs.lig <- reduceToLigand(pairs.p)
#' }
reduceToLigand <- function(pairs){

  t <- sum(c("pval","LLR") %in% names(pairs))
  if (t ==2)
    stop("Both P-values and LLR present in LR table")
  if (t == 0)
    stop("P-values or LLR must be available in LR table")

  red.mode.R <- (sum(regexpr("^{",pairs$R,perl=TRUE)==1)==nrow(pairs) && sum(regexpr("}$",pairs$R,perl=TRUE)==nchar(pairs$R))==nrow(pairs))
  if (red.mode.R)
    stop("Already reduced to ligand")

  # best LLR or pval per receptor/pathway
  if ("pval" %in% names(pairs))
    best <- pairs %>% dplyr::group_by(L,pw.id) %>% dplyr::filter(pval==min(pval))
  else
    best <- pairs %>% dplyr::group_by(L,pw.id) %>% dplyr::filter(LLR==max(LLR))

  # list of receptors for each ligand/pathway combination
  receptors <- pairs %>% dplyr::group_by(L,pw.id) %>% dplyr::summarize(Rs=paste0('{',paste(R,collapse=';'),'}'))

  # replace receptor by receptor lists
  lpw <- paste(best$L,best$pw.id,sep="|")
  lpwR <- stats::setNames(receptors$Rs,paste(receptors$L,receptors$pw.id,sep="|"))
  best$R <- lpwR[lpw]

  best

} # reduceToLigand


#' Aggregate ligands and receptors at the pathway level
#'
#' Simplifies a ligand-receptor table to focus on the pathways.
#'
#' @param pairs         A ligand-receptor table such as output by \code{pValuesLR} and \code{naiveBayesLR}.
#' @return A table with the same structure as \code{pairs} but reduced to only report one row per pathway.
#' All the ligands and all the receptors forming pairs related to a certain pathway, are combined in two semi-colon-separated
#' lists surounded by curly braces. The information of which ligand interacted with which receptor is lost.
#'
#' There is no pathway selection, i.e. potentially redundant pathways with different P-values or log-likelihood ratios (LLRs)
#' are all kept in the output. \code{reduceToBestPathway()} can be used to select the one with best P-value or LLR.
#' @export
#' @examples
#' \dontrun{
#' data(sdc,package="BulkSignalR")
#' sample.types <- rep("tumor",ncol(sdc))
#' sample.types[grep("^N",names(sdc),perl=TRUE)] <- "normal"
#' ds <- prepareDataset(sdc,sample.types)
#' ds <- learnParameters(ds)
#'
#' # since model training takes time and it can be reused, it is advisable to save the model
#' save(ds,file="...",compress="bzip2")
#'
#' ds.LR <- getCorrelatedLR(ds)
#' ds.LR <- checkReceptorSignaling(ds,ds.LR)
#' pairs.p <- pValuesLR(ds.LR,ds$param)
#' pairs.pw <- reduceToPathway(pairs.p)
#' pairs.pw.red <- reduceToBestPathway(pairs.pw)
#' }
reduceToPathway <- function(pairs){

  t <- sum(c("pval","LLR") %in% names(pairs))
  if (t ==2)
    stop("Both P-values and LLR present in LR table")
  if (t == 0)
    stop("P-values or LLR must be available in LR table")

  # best LLR or pval per pathway
  if ("pval" %in% names(pairs))
    best <- pairs %>% dplyr::group_by(pw.id) %>% dplyr::filter(pval==min(pval))
  else
    best <- pairs %>% dplyr::group_by(pw.id) %>% dplyr::filter(LLR==max(LLR))

  # list of ligands and receptors for each pathway
  lire <- pairs %>% dplyr::group_by(pw.id) %>% dplyr::summarize(Ls=paste0('{',paste(unique(unlist(strsplit(gsub("}$","",gsub("^{","",L,perl=TRUE)),split=";"))),collapse=';'),'}'),
                                                                Rs=paste0('{',paste(unique(unlist(strsplit(gsub("}$","",gsub("^{","",R,perl=TRUE)),split=";"))),collapse=';'),'}'))

  # replace ligands and receptors by lists
  rpwL <- stats::setNames(lire$Ls,lire$pw.id)
  rpwR <- stats::setNames(lire$Rs,lire$pw.id)
  best$L <- rpwL[best$pw.id]
  best$R <- rpwR[best$pw.id]

  best

} # reduceToPathway


