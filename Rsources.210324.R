removeRows <- function(rowNum, data) {
    newData <- data[-rowNum, , drop = FALSE]
    rownames(newData) <- NULL
    newData
}




getTaxonomy <- function(otus, tax_tab, level, taxRanks,na_str = c( "NA")) {
    ranks <- taxRanks
    sel <- ranks[1:match(level, ranks)]
    inds <- apply(tax_tab[otus,sel], 1, function(x) max(which(!(x %in% na_str))))
    retval <- as.data.frame(tax_tab)[cbind(otus, ranks[inds])]
    retval[inds!=match(level, ranks)] <- paste(na_str[1], retval[inds!=match(level, ranks)], sep="_")
    return(retval)
}


getTaxonomyAll <- function(otus, tax_tab, level, na_str = c( "NA")) {
    rank <- rank
    sel <- rank[1:match(level, rank)]
    inds <- apply(tax_tab[otus,sel], 1, function(x) max(which(!(x %in% na_str))))
    retval <- as.data.frame(tax_tab)[cbind(otus, rank[inds])]
    retval[inds!=match(level, rank)] <- paste(na_str[1], retval[inds!=match(level, rank)], sep="_")
    return(retval)
}



# X is indicator matrix of predictions, Y is indicator matrix of truth
# columns are classes, rows are samples
mcc <- function(preds=NULL, actuals=NULL, x=NULL, y=NULL) {
    # if preds and actuals are provided, x and y will be ignored
    if (!is.null(preds)) {
        nclasses <- length(union(preds, actuals))
        x <- matrix(0, nrow=length(preds), ncol=nclasses)
        y <- matrix(0, nrow=length(actuals), ncol=nclasses)
        x[cbind(1:nrow(x), preds+1)] <- 1
        y[cbind(1:nrow(y), actuals+1)] <- 1
    }
    if (!all(dim(x) == dim(y))) {
        stop("X and Y must have the same dimensions")
    }
    
    cov_biased <- function(x, y) {
        sum(sapply(1:ncol(x), function(k) {
            cov(x[,k], y[,k]) # unbiased estimate with (n-1) denominator as opposed to (n), but cancels out anyways so identical result
        }))
    }
    numerator <- cov_biased(x,y)
    denominator <- sqrt(cov_biased(x,x) * cov_biased(y,y))
    numerator / denominator
}


multiplot <- function(..., plotlist=NULL, file, cols=1, rows=1) {
    require(grid)
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    
    i = 1
    while (i < numPlots) {
        numToPlot <- min(numPlots-i+1, cols*rows)
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(i, i+cols*rows-1), ncol = cols, nrow = rows, byrow=T)
        if (numToPlot==1) {
            print(plots[[i]])
        } else {
            # Set up the page
            grid.newpage()
            pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
            # Make each plot, in the correct location
            for (j in i:(i+numToPlot-1)) {
                # Get the i,j matrix positions of the regions that contain this subplot
                matchidx <- as.data.frame(which(layout == j, arr.ind = TRUE))
                print(plots[[j]], vp = viewport(layout.pos.row = matchidx$row,
                layout.pos.col = matchidx$col))
            }
        }
        i <- i+numToPlot
    }
}

normalizeByRows <- function (df, rsum=1)
{
    while (any(abs((rowSums(df)-rsum))>1e-13)) {
        df <- rsum*(df / rowSums(df))
    }
    return(df)
}
normalizeByCols <- function (df, csum=1, level=NULL, delim="\\|")
{
    if (is.null(level)) {
        while (any(abs((colSums(df)-csum))>1e-13 & colSums(df)!=0, na.rm=T)) {
            missing <- which(colSums(df)==0)
            df <- sweep(df, 2, colSums(df)/csum, "/")
            df[,missing] <- 0
        }
    } else {
        tmp <- df
        tmp$taxa <- rownames(tmp)
        tmp$splitter <- factor(unlist(lapply(rownames(tmp), function(x) unlist(strsplit(x, delim))[level])))
        names <- rownames(tmp)[order(tmp$splitter)]
        tmp <- ddply(tmp, .(splitter), function(x) {
            x <- x[, setdiff(colnames(x), c("taxa", "splitter"))]
            while (any(abs((colSums(x)-csum))>1e-13 & colSums(df)!=0, na.rm=T)) {
                x <- sweep(x, 2, colSums(x)/csum, "/")
            }
            x
        })
        rownames(tmp) <- names
        df <- tmp[, setdiff(colnames(tmp), "splitter")]
    }
    return(df)
}

renameLevelsWithCounts <- function(fvec, originalLevelsAsNames=FALSE) {
    tab <- table(fvec)
    retval <- sprintf("%s (n=%d)", fvec, tab[unlist(lapply(fvec, function(x) match(x, names(tab))))])
    #    newlevels <- sprintf("%s (n=%d)", levels(fvec), tab[levels(fvec)])
    newlevels <- sprintf("%s (n=%d)", levels(fvec), tab[unlist(lapply(names(tab), function(x) which(levels(fvec)==x)))])
    retval <- factor(retval, levels=newlevels)
    if (originalLevelsAsNames) {
        names(retval) <- fvec
    }
    return(retval)
}






# ============================================================
# Tutorial on plotting significant taxa and environmental variables on an NMDS plot using ggplot2
# by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)
# =============================================================

library(vegan)
library(ggplot2)
library(grid)


bv.step <- function(fix.mat, var.mat,
fix.dist.method="bray", var.dist.method="euclidean", correlation.method="spearman",
scale.fix=FALSE, scale.var=TRUE,
max.rho=0.95,
min.delta.rho=0.001,
random.selection=TRUE,
prop.selected.var=0.2,
num.restarts=10,
var.always.include=NULL,
var.exclude=NULL,
output.best=10
){
    
    if(dim(fix.mat)[1] != dim(var.mat)[1]){stop("fixed and variable matrices must have the same number of rows")}
    if(sum(var.always.include %in% var.exclude) > 0){stop("var.always.include and var.exclude share a variable")}
    require(vegan)
    
    if(scale.fix){fix.mat<-scale(fix.mat)}else{fix.mat<-fix.mat}
    if(scale.var){var.mat<-scale(var.mat)}else{var.mat<-var.mat}
    
    fix.dist <- vegdist(as.matrix(fix.mat), method=fix.dist.method)
    
    #an initial removal phase
    var.dist.full <- vegdist(as.matrix(var.mat), method=var.dist.method)
    full.cor <- suppressWarnings(cor.test(fix.dist, var.dist.full, method=correlation.method))$estimate
    var.comb <- combn(1:ncol(var.mat), ncol(var.mat)-1)
    RES <- data.frame(var.excl=rep(NA,ncol(var.comb)), n.var=ncol(var.mat)-1, rho=NA)
    for(i in 1:dim(var.comb)[2]){
        var.dist <- vegdist(as.matrix(var.mat[,var.comb[,i]]), method=var.dist.method)
        temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
        RES$var.excl[i] <- c(1:ncol(var.mat))[-var.comb[,i]]
        RES$rho[i] <- temp$estimate
    }
    delta.rho <- RES$rho - full.cor
    exclude <- sort(unique(c(RES$var.excl[which(abs(delta.rho) < min.delta.rho)], var.exclude)))
    
    if(random.selection){
        num.restarts=num.restarts
        prop.selected.var=prop.selected.var
        prob<-rep(1,ncol(var.mat))
        if(prop.selected.var< 1){
            prob[exclude]<-0
        }
        n.selected.var <- min(sum(prob),prop.selected.var*dim(var.mat)[2])
    } else {
        num.restarts=1
        prop.selected.var=1
        prob<-rep(1,ncol(var.mat))
        n.selected.var <- min(sum(prob),prop.selected.var*dim(var.mat)[2])
    }
    
    RES_TOT <- c()
    for(i in 1:num.restarts){
        step=1
        RES <- data.frame(step=step, step.dir="F", var.incl=NA, n.var=0, rho=0)
        attr(RES$step.dir, "levels") <- c("F","B")
        best.comb <- which.max(RES$rho)
        best.rho <- RES$rho[best.comb]
        delta.rho <- Inf
        selected.var <- sort(unique(c(sample(1:dim(var.mat)[2], n.selected.var, prob=prob), var.always.include)))
        while(best.rho < max.rho & delta.rho > min.delta.rho & RES$n.var[best.comb] < length(selected.var)){
            #forward step
            step.dir="F"
            step=step+1
            var.comb <- combn(selected.var, RES$n.var[best.comb]+1, simplify=FALSE)
            if(RES$n.var[best.comb] == 0){
                var.comb.incl<-1:length(var.comb)
            } else {
                var.keep <- as.numeric(unlist(strsplit(RES$var.incl[best.comb], ",")))
                temp <- NA*1:length(var.comb)
                for(j in 1:length(temp)){
                    temp[j] <- all(var.keep %in% var.comb[[j]])
                }
                var.comb.incl <- which(temp==1)
            }
            
            RES.f <- data.frame(step=rep(step, length(var.comb.incl)), step.dir=step.dir, var.incl=NA, n.var=RES$n.var[best.comb]+1, rho=NA)
            for(f in 1:length(var.comb.incl)){
                var.incl <- var.comb[[var.comb.incl[f]]]
                var.incl <- var.incl[order(var.incl)]
                var.dist <- vegdist(as.matrix(var.mat[,var.incl]), method=var.dist.method)
                temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
                RES.f$var.incl[f] <- paste(var.incl, collapse=",")
                RES.f$rho[f] <- temp$estimate
            }
            
            last.F <- max(which(RES$step.dir=="F"))
            RES <- rbind(RES, RES.f[which.max(RES.f$rho),])
            best.comb <- which.max(RES$rho)
            delta.rho <- RES$rho[best.comb] - best.rho
            best.rho <- RES$rho[best.comb]
            
            if(best.comb == step){
                while(best.comb == step & RES$n.var[best.comb] > 1){
                    #backward step
                    step.dir="B"
                    step <- step+1
                    var.keep <- as.numeric(unlist(strsplit(RES$var.incl[best.comb], ",")))
                    var.comb <- combn(var.keep, RES$n.var[best.comb]-1, simplify=FALSE)
                    RES.b <- data.frame(step=rep(step, length(var.comb)), step.dir=step.dir, var.incl=NA, n.var=RES$n.var[best.comb]-1, rho=NA)
                    for(b in 1:length(var.comb)){
                        var.incl <- var.comb[[b]]
                        var.incl <- var.incl[order(var.incl)]
                        var.dist <- vegdist(as.matrix(var.mat[,var.incl]), method=var.dist.method)
                        temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
                        RES.b$var.incl[b] <- paste(var.incl, collapse=",")
                        RES.b$rho[b] <- temp$estimate
                    }
                    RES <- rbind(RES, RES.b[which.max(RES.b$rho),])
                    best.comb <- which.max(RES$rho)
                    best.rho<- RES$rho[best.comb]
                }
            } else {
                break()
            }
            
        }
        
        RES_TOT <- rbind(RES_TOT, RES[2:dim(RES)[1],])
        print(paste(round((i/num.restarts)*100,3), "% finished"))
    }
    
    RES_TOT <- unique(RES_TOT[,3:5])
    
    
    if(dim(RES_TOT)[1] > output.best){
        order.by.best <- RES_TOT[order(RES_TOT$rho, decreasing=TRUE)[1:output.best],]
    } else {
        order.by.best <-  RES_TOT[order(RES_TOT$rho, decreasing=TRUE), ]
    }
    rownames(order.by.best)<-NULL
    
    order.by.i.comb <- c()
    for(i in 1:length(selected.var)){
        f1 <- which(RES_TOT$n.var==i)
        f2 <- which.max(RES_TOT$rho[f1])
        order.by.i.comb <- rbind(order.by.i.comb, RES_TOT[f1[f2],])
    }
    rownames(order.by.i.comb)<-NULL
    
    if(length(exclude)<1){var.exclude=NULL} else {var.exclude=exclude}
    out <- list(
    order.by.best=order.by.best,
    order.by.i.comb=order.by.i.comb,
    best.model.vars=paste(colnames(var.mat)[as.numeric(unlist(strsplit(order.by.best$var.incl[1], ",")))], collapse=","),
    best.model.rho=order.by.best$rho[1],
    var.always.include=var.always.include,
    var.exclude=var.exclude
    )
    out
    
}

bio.env <- function(fix.mat, var.mat,
fix.dist.method="bray", var.dist.method="euclidean", correlation.method="spearman",
scale.fix=FALSE, scale.var=TRUE,
output.best=10,
var.max=ncol(var.mat)
){
    if(dim(fix.mat)[1] != dim(var.mat)[1]){stop("fixed and variable matrices must have the same number of rows")}
    if(var.max > dim(var.mat)[2]){stop("var.max cannot be larger than the number of variables (columns) in var.mat")}
    
    require(vegan)
    
    combn.sum <- sum(factorial(ncol(var.mat))/(factorial(1:var.max)*factorial(ncol(var.mat)-1:var.max)))
    
    if(scale.fix){fix.mat<-scale(fix.mat)}else{fix.mat<-fix.mat}
    if(scale.var){var.mat<-scale(var.mat)}else{var.mat<-var.mat}
    fix.dist <- vegdist(fix.mat, method=fix.dist.method)
    RES_TOT <- c()
    best.i.comb <- c()
    iter <- 0
    for(i in 1:var.max){
        var.comb <- combn(1:ncol(var.mat), i, simplify=FALSE)
        RES <- data.frame(var.incl=rep(NA, length(var.comb)), n.var=i, rho=0)
        for(f in 1:length(var.comb)){
            iter <- iter+1
            var.dist <- vegdist(as.matrix(var.mat[,var.comb[[f]]]), method=var.dist.method)
            temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
            RES$var.incl[f] <- paste(var.comb[[f]], collapse=",")
            RES$rho[f] <- temp$estimate
            if(iter %% 100 == 0){print(paste(round(iter/combn.sum*100, 3), "% finished"))}
        }
        
        order.rho <- order(RES$rho, decreasing=TRUE)
        best.i.comb <- c(best.i.comb, RES$var.incl[order.rho[1]])
        if(length(order.rho) > output.best){
            RES_TOT <- rbind(RES_TOT, RES[order.rho[1:output.best],])
        } else {
            RES_TOT <- rbind(RES_TOT, RES)
        }
    }
    rownames(RES_TOT)<-NULL
    
    if(dim(RES_TOT)[1] > output.best){
        order.by.best <- order(RES_TOT$rho, decreasing=TRUE)[1:output.best]
    } else {
        order.by.best <- order(RES_TOT$rho, decreasing=TRUE)
    }
    OBB <- RES_TOT[order.by.best,]
    rownames(OBB) <- NULL
    
    order.by.i.comb <- match(best.i.comb, RES_TOT$var.incl)
    OBC <- RES_TOT[order.by.i.comb,]
    rownames(OBC) <- NULL
    
    out <- list(
    order.by.best=OBB,
    order.by.i.comb=OBC,
    best.model.vars=paste(colnames(var.mat)[as.numeric(unlist(strsplit(OBB$var.incl[1], ",")))], collapse=",") ,
    best.model.rho=OBB$rho[1]
    )
    out
}



heatmap.3 = function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
distfun = dist, hclustfun = hclust, dendrogram = c("both",
"row", "column", "none"), reorderfun = function(d, w) reorder(d,
w), symm = FALSE, scale = c("none", "row", "column"),
na.rm = TRUE, revC = identical(Colv, "Rowv"), add.expr, breaks,
symbreaks = any(x < 0, na.rm = TRUE) || scale != "none",
col = "heat.colors", colsep, rowsep, sepcolor = "white",
sepwidth = c(0.05, 0.05), cellnote, notecex = 1, notecol = "cyan",
na.color = par("bg"), trace = c("column", "row", "both",
"none"), tracecol = "cyan", hline = median(breaks), vline = median(breaks),
linecol = tracecol, margins = c(5,5,5,5), ColSideColors, RowSideColors, side.height.fraction=0.3,
cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL,
labCol = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0,
NA), adjCol = c(NA, 0), offsetRow = 0.5, offsetCol = 0.5,
key = TRUE, keysize = 1.5, density.info = c("histogram",
"density", "none"), denscol = tracecol, symkey = any(x <
0, na.rm = TRUE) || symbreaks, densadj = 0.25, key.title = NULL,
key.xlab = NULL, key.ylab = NULL, key.xtickfun = NULL, key.ytickfun = NULL,
key.par = list(), main = NULL, xlab = NULL, ylab = NULL,
lmat = NULL, lhei = NULL, lwid = NULL, ColSideColorsSize = 1, RowSideColorsSize = 1, extrafun = NULL, ...)
{
    library(gtools)
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
    "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
    "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
    else if (all(Colv == "Rowv"))
    Colv <- Rowv
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 4)
    stop("`margins' must be a numeric vector of length 4")
    if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((is.logical(Rowv) && !isTRUE(Rowv)) || (is.null(Rowv))) &&
        (dendrogram %in% c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
            dendrogram <- "column"
            else dendrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
            dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((is.logical(Colv) && !isTRUE(Colv)) || (is.null(Colv))) &&
        (dendrogram %in% c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
            dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
            dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
        if (length(rowInd) > nr || any(rowInd < 1 | rowInd >
        nr))
        stop("Rowv dendrogram doesn't match size of x")
    }
    else if (is.integer(Rowv)) {
        browser()
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorderfun(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorderfun(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
        if (length(colInd) > nc || any(colInd < 1 | colInd >
        nc))
        stop("Colv dendrogram doesn't match size of x")
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorderfun(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorderfun(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
    (1:nr)[rowInd]
    else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
    (1:nc)[colInd]
    else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) <
    1) {
        if (missing(col) || is.function(col))
        breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks)
        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
        length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
    col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
            stop("'ColSideColors' must be a matrix of nrow(x) rows")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
            stop("'RowSideColors' must be a matrix of ncol(x) columns")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        if (!is.matrix(RowSideColors)){
            par(mar = c(margins[1], 0, 0, 0.5))
            image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        } else {
            par(mar = c(margins[1], 0, 0, 0.5))
            rsc = t(RowSideColors[,rowInd, drop=F])
            rsc.colors = matrix()
            rsc.names = names(table(rsc))
            rsc.i = 1
            for (rsc.name in rsc.names) {
                rsc.colors[rsc.i] = rsc.name
                rsc[rsc == rsc.name] = rsc.i
                rsc.i = rsc.i + 1
            }
            rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
            image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
            if (length(rownames(RowSideColors)) > 0) {
                axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
            }
        }
    }
    if (!missing(ColSideColors)) {
        
        if (!is.matrix(ColSideColors)){
            par(mar = c(0.5, 0, 0, margins[4]))
            image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        } else {
            #            par(mar = c(0.5, 0, 0, margins[4]))
            par(mar = c(0.5, margins[2], 0, margins[4]))
            csc = ColSideColors[colInd, , drop=F]
            csc.colors = matrix()
            csc.names = names(table(csc))
            csc.i = 1
            for (csc.name in csc.names) {
                csc.colors[csc.i] = csc.name
                csc[csc == csc.name] = csc.i
                csc.i = csc.i + 1
            }
            csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
            image(csc, col = as.vector(csc.colors), axes = FALSE)
            if (length(colnames(ColSideColors)) > 0) {
                axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
            }
        }
    }
    par(mar = margins)
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
        ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
    c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col,
    breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
    retval$rowDendrogram <- ddr
    if (exists("ddc"))
    retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
        col = na.color, add = TRUE)
    }
    if (is.null(srtCol))
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 +
    offsetCol, tick = 0, cex.axis = cexCol, hadj = adjCol[1],
    padj = adjCol[2])
    else {
        if (is.numeric(srtCol)) {
            if (missing(adjCol) || is.null(adjCol))
            adjCol = c(1, NA)
            xpd.orig <- par("xpd")
            par(xpd = NA)
            xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2,
            tick = 0)
            text(x = xpos, y = par("usr")[3] - (1 + offsetCol) *
            strheight("M"), labels = labCol, adj = adjCol,
            cex = cexCol, srt = srtCol)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtCol ignored.")
    }
    if (is.null(srtRow)) {
        par(mar = c(margins[1L], 0, 0, margins[4L]))
        axis(4, iy, labels = labRow, las = 2, line = -0.5 + offsetRow,
        tick = 0, cex.axis = cexRow, hadj = adjRow[1], padj = adjRow[2])
    }
    else {
        if (is.numeric(srtRow)) {
            xpd.orig <- par("xpd")
            par(xpd = NA)
            ypos <- axis(4, iy, labels = rep("", nr), las = 2,
            line = -0.5, tick = 0)
            text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"),
            y = ypos, labels = labRow, adj = adjRow, cex = cexRow,
            srt = srtRow)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtRow ignored.")
    }
    if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
    if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[4] - 1.25)
    if (!missing(add.expr))
    eval(substitute(add.expr))
    if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0,
    xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) +
    1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) +
    1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) +
    1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1,
    col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i - 0.5 + hline.vals, col = linecol,
                lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
    col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        flag <- try(plot.dendrogram(ddr, horiz = TRUE, axes = FALSE,
        yaxs = "i", leaflab = "none"))
        if ("try-error" %in% class(flag)) {
            cond <- attr(flag, "condition")
            if (!is.null(cond) && conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?")
            stop("Row dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
        }
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[4]))
    if (dendrogram %in% c("both", "column")) {
        flag <- try(plot.dendrogram(ddc, axes = FALSE, xaxs = "i",
        leaflab = "none"))
        if ("try-error" %in% class(flag)) {
            cond <- attr(flag, "condition")
            if (!is.null(cond) && conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?")
            stop("Column dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
        }
    }
    else plot.new()
    if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        mar <- c(5, 4, 2, 1)
        if (!is.null(key.xlab) && is.na(key.xlab))
        mar[1] <- 2
        if (!is.null(key.ylab) && is.na(key.ylab))
        mar[2] <- 2
        if (!is.null(key.title) && is.na(key.title))
        mar[3] <- 1
        par(mar = mar, cex = 0.75, mgp = c(2, 1, 0))
        if (length(key.par) > 0)
        do.call(par, key.par)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
        xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        if (is.null(key.xtickfun)) {
            lv <- pretty(breaks)
            xv <- scale01(as.numeric(lv), min.raw, max.raw)
            xargs <- list(at = xv, labels = lv)
        }
        else {
            xargs <- key.xtickfun()
        }
        xargs$side <- 1
        do.call(axis, xargs)
        if (is.null(key.xlab)) {
            if (scale == "row")
            key.xlab <- "Row Z-Score"
            else if (scale == "column")
            key.xlab <- "Column Z-Score"
            else key.xlab <- "Value"
        }
        if (!is.na(key.xlab)) {
            mtext(side = 1, key.xlab, line = par("mgp")[1], padj = 0.5)
        }
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
            if (is.null(key.ytickfun)) {
                yargs <- list(at = pretty(dens$y)/max(dens$y) *
                0.95, labels = pretty(dens$y))
            }
            else {
                yargs <- key.ytickfun()
            }
            yargs$side <- 2
            do.call(axis, yargs)
            if (is.null(key.title))
            key.title <- "Color Key\nand Density Plot"
            if (!is.na(key.title))
            title(key.title)
            par(cex = 0.5)
            if (is.null(key.ylab))
            key.ylab <- "Density"
            if (!is.na(key.ylab))
            mtext(side = 2, key.ylab, line = par("mgp")[1],
            padj = 0.5)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
            if (is.null(key.ytickfun)) {
                yargs <- list(at = pretty(hy)/max(hy) * 0.95,
                labels = pretty(hy))
            }
            else {
                yargs <- key.ytickfun()
            }
            yargs$side <- 2
            do.call(axis, yargs)
            if (is.null(key.title))
            key.title <- "Color Key\nand Histogram"
            if (!is.na(key.title))
            title(key.title)
            par(cex = 0.5)
            if (is.null(key.ylab))
            key.ylab <- "Count"
            if (!is.na(key.ylab))
            mtext(side = 2, key.ylab, line = par("mgp")[1],
            padj = 0.5)
        }
        else if (is.null(key.title))
        title("Color Key")
        if (trace %in% c("both", "column")) {
            vline.vals <- scale01(vline, min.raw, max.raw)
            if (!is.null(vline)) {
                abline(v = vline.vals, col = linecol, lty = 2)
            }
        }
        if (trace %in% c("both", "row")) {
            hline.vals <- scale01(hline, min.raw, max.raw)
            if (!is.null(hline)) {
                abline(v = hline.vals, col = linecol, lty = 2)
            }
        }
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
    high = retval$breaks[-1], color = retval$col)
    if (!is.null(extrafun))
    extrafun()
    invisible(retval)
}


normalizeByRows <- function (df, rsum=1)
{
    while (any(abs((rowSums(df)-rsum))>1e-13)) {
        df <- rsum*(df / rowSums(df))
    }
    return(df)
}
normalizeByCols <- function (df, csum=1, level=NULL, delim="\\|")
{
    if (is.null(level)) {
        while (any(abs((colSums(df)-csum))>1e-13 & colSums(df)!=0, na.rm=T)) {
            missing <- which(colSums(df)==0)
            df <- sweep(df, 2, colSums(df)/csum, "/")
            df[,missing] <- 0
        }
    } else {
        tmp <- df
        tmp$taxa <- rownames(tmp)
        tmp$splitter <- factor(unlist(lapply(rownames(tmp), function(x) unlist(strsplit(x, delim))[level])))
        names <- rownames(tmp)[order(tmp$splitter)]
        tmp <- ddply(tmp, .(splitter), function(x) {
            x <- x[, setdiff(colnames(x), c("taxa", "splitter"))]
            while (any(abs((colSums(x)-csum))>1e-13 & colSums(df)!=0, na.rm=T)) {
                x <- sweep(x, 2, colSums(x)/csum, "/")
            }
            x
        })
        rownames(tmp) <- names
        df <- tmp[, setdiff(colnames(tmp), "splitter")]
    }
    return(df)
}





data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
}




# Title: Geometric Mean of Pairwise Ratios (GMPR) for Microbiome Sequencing data normalization
# Version: 0.1
# Authors: Jun Chen (chen.jun2@mayo.edu)
# Date: 2017/02/07
# Description: The function calculates the normalizing factors for microbiome sequencing data or, more generally, zeroinflated sequencing data.
# The size factors can be used as offsets in count-based regression models or as divisors to produce normalized data


require(matrixStats)

GMPR <- function (comm, intersect.no = 10, ct.min = 1, trace = TRUE) {
  # Computes the GMPR size factor
  #
  # Args:
  #   comm: a matrix of counts, row - features (OTUs, genes, etc) , column - sample
  #   intersect.no: the minimum number of shared features between sample pair, where the ratio is calculated
  #   ct.min: the minimum number of counts required to calculate ratios
  
  #
  # Returns:
  #   a vector of the size factors with attribute 'NSS'. Samples with distinct sets of features will be output as NA.
  #         NSS:   number of samples with significant sharing (> intersect.no) including itself
  
  # mask counts < ct.min
  comm[comm < ct.min] <- 0
  
  if (is.null(colnames(comm))) {
    colnames(comm) <- paste0('S', 1:ncol(comm))
  }
  
  if (trace) cat('Begin GMPR size factor calculation ...\n')
  
  comm.no <- numeric(ncol(comm))
  gmpr <- sapply(1:ncol(comm),  function(i) {
    if (i %% 50 == 0) {
      cat(i, '\n')
    }
    x <- comm[, i]
    # Compute the pairwise ratio
    pr <- x / comm
    # Handling of the NA, NaN, Inf
    pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
    # Counting the number of non-NA, NaN, Inf
    incl.no <- colSums(!is.na(pr))
    # Calculate the median of PR
    pr.median <- colMedians(pr, na.rm=TRUE)
    # Record the number of samples used for calculating the GMPR
    comm.no[i] <<- sum(incl.no >= intersect.no)
    # Geometric mean of PR median
    if (comm.no[i] > 1) {
      return(exp(mean(log(pr.median[incl.no >= intersect.no]))))
    } else {
      return(NA)
    }
  }
  )
  
  if (sum(is.na(gmpr))) {
    warning(paste0('The following samples\n ', paste(colnames(comm)[is.na(gmpr)], collapse='\n'),
                   '\ndo not share at least ', intersect.no, ' common taxa with the rest samples! ',
                   'For these samples, their size factors are set to be NA! \n',
                   'You may consider removing these samples since they are potentially outliers or negative controls!\n',
                   'You may also consider decreasing the minimum number of intersecting taxa and rerun the procedure!\n'))
  }
  
  if (trace) cat('Completed!\n')
  if (trace) cat('Please watch for the samples with limited sharing with other samples based on NSS! They may be outliers! \n')
  names(gmpr) <- names(comm.no) <- colnames(comm)
  
  attr(gmpr, 'NSS') <- comm.no
  
  return(gmpr)
}
