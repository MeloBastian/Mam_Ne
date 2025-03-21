library(ggplot2)
#data(mtcars)
#corr <- M
#p.mat = PP
#head(corr[, 1:6])
library(ggcorrplot)

###data
#diag sup 144 sp 12psfiltred, diag inf 89 sp 12psfiltred
M <- matrix(nrow=7, ncol=7, data=c(1,-0.28,-0.58,-0.73,-0.65,0.13,-0.18,
                                   -0.61,1,0.31,0.37,0.25,-0.46,0.75,
                                   -0.62,0.44,1,0.59,0.58,-0.19,0.31,
                                   -0.73,0.49,0.59,1,0.61,-0.2,0.32,
                                   -0.7,0.52,0.63,0.68,1,-0.21,0.13,
                                   0.3,-0.49,-0.28,-0.31,-0.2,1,-0.56, 
                                   -0.39,0.7,0.37,0.38,0.36,-0.72,1), byrow=T)
PP=matrix(nrow=7, ncol=7, data=c(1,0.0097,0,0,0,0.85,0.06,
                                 0,1,1,1,0.98,0.00059,1,
                                 0,1,1,1,1,0.018,1,
                                 0,1,1,1,1,0.019,1,
                                 0,1,1,1,1,0.015,0.904,
                                 0.99,0.0004,0.0059,0.0025,0.045,1,0,
                                 0.001,1,1,1,1,0,1), byrow=T)

#diag sup 144sp, 6 sp filtred, diag inf 89 sp 6 sp filtred
M<-matrix(nrow=7, ncol=7, data=c(1,-0.27,-0.58,-0.75,-0.67,0.11,-0.13,
                                 -0.6,1,0.29,0.35,0.24,-0.41,0.73,
                                 -0.62,0.44,1,0.59,0.58,-0.15,0.27,
                                 -0.73,0.48,0.59,1,0.61,-0.09,0.26,
                                 -0.7,0.53,0.63,0.68,1,-0.14,0.1,
                                 0.19,-0.44,-0.23,-0.13,-0.1,1,-0.58,
                                 -0.31,0.67,0.34,0.3,0.29,-0.74,1), byrow = T)
PP<-matrix(nrow=7, ncol=7, data=c(1,0.01,0,0,0,0.8,0.136,
                                 0,1,0.99,1,0.98,0.0009,1,
                                 0,1,1,1,1,0.748,1,
                                 0,1,1,1,1,0.2,0.998,
                                 0,1,1,1,1,0.1,0.832,
                                 0.93,0.0009,0.0148,0.1159,0.172,1,0,
                                 0.007,1,1,0.99,1,0,1
                                 ),byrow=T)

rownames(M)=c("dS", "dN/dS", "mass", "maturity","longevity", "pS", "pN/pS")
colnames(M)=c("dS", "dN/dS", "mass", "maturity","longevity", "pS", "pN/pS")
rownames(PP)=c("dS", "dN/dS", "mass", "maturity","longevity", "pS", "pN/pS")
colnames(PP)=c("dS", "dN/dS", "mass", "maturity","longevity", "pS", "pN/pS")

#diagsup 89 sp low ps, 6 sp filtred, diaginf 144 sp 6 sp filtred !pgls!, pas de dS
M<-matrix(nrow=6, ncol=6, data=c(1,0.41,0.44,0.33,-0.36,0.51,
                                 0.16,1,0.58,0.71,-0.25,0.41,
                                 0.22,0.57,1,0.7,-0.31,0.46,
                                 0.20,0.63,0.63,1,-0.20,0.34,
                                 -0.37,-0.14,-0.19,-0.21,1,-0.72,
                                 0.41,0.27,0.31,0.21,-0.62,1), byrow=T)

matrice <- M #matrice à l'envers...
n <- nrow(matrice)
for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    # Échanger les éléments
    temp <- matrice[i, j]
    matrice[i, j] <- matrice[j, i]
    matrice[j, i] <- temp
  }
}
M<- matrice


PP<-matrix(nrow = 6, ncol=6, data=c(1,0.0003,0.0002,0.007,0.001,2.1e-06,
                                    0.063,1,1.8e-07,1.9e-11,0.034,0.0004,
                                    0.017,1.7e-11,1,1.2e-10,0.014,0.0001,
                                    0.032,2.4e-14,1.5e-13,1,0.103,0.0055,
                                    6.8e-06,0.116,0.039,0.023,1,1.7e-13,
                                    6.6e-07,0.0028,0.00091,0.0242,8.5e-16,1),byrow=T)

matrice <- PP #matrice à l'envers...
n <- nrow(matrice)
for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    # Échanger les éléments
    temp <- matrice[i, j]
    matrice[i, j] <- matrice[j, i]
    matrice[j, i] <- temp
  }
}
PP<- matrice
rownames(M)=c("dN/dS", "mass", "maturity","longevity", "pS", "pN/pS")
colnames(M)=c("dN/dS", "mass", "maturity","longevity", "pS", "pN/pS")
rownames(PP)=c("dN/dS", "mass", "maturity","longevity", "pS", "pN/pS")
colnames(PP)=c("dN/dS", "mass", "maturity","longevity", "pS", "pN/pS")

#diagsup 89 sp , 12 sp filtred, diaginf 144 sp  filtred !pgls!, pas de dS
M<-matrix(nrow = 6, ncol=6, data=c(1,0.16,0.22,0.20,-0.34,0.37,
                                   0.37,1,0.57,0.63,-0.22,0.33,
                                   0.41,0.58,1,0.63,-0.28,0.34,
                                   0.34,0.67,0.69,1,-0.24,0.23,
                                   -0.37,-0.25,-0.42,-0.23,1,-0.61,
                                   0.54,0.36,0.45,0.38,-0.78,1),byrow=T)
PP<-matrix(nrow = 6, ncol=6, data=c(1,0.063,0.017,0.032,5.4e-05,1.4e-05,
                                    0.0007,1,1.8e-11,2.4e-14,0.019,0.0002,
                                    0.0003,7.5e-08,1,1.5e-13,0.0035,0.0004,
                                    0.003,4.9e-11,5.8e-11,1,0.01,0.017,
                                    0.0007,0.031,0.0004,0.059,1,1.3e-14,
                                    2e-07,0.0014,0.0001,0.001,2.2e-16,1), byrow=T)

rownames(M)=c("dN/dS", "mass", "maturity","longevity", "pS", "pN/pS")
colnames(M)=c("dN/dS", "mass", "maturity","longevity", "pS", "pN/pS")
rownames(PP)=c("dN/dS", "mass", "maturity","longevity", "pS", "pN/pS")
colnames(PP)=c("dN/dS", "mass", "maturity","longevity", "pS", "pN/pS")


### La ici la longue fonction qu'on peut bidouiller  et tout en bas un exemple avec le tableau corr et les p_value dans p.mat (la c'est les memes)

ggcorrplot <- function(corr,
                       method = c("square", "circle"),
                       type = c("full", "lower", "upper"),
                       ggtheme = ggplot2::theme_minimal,
                       title = "",
                       show.legend = TRUE,
                       legend.title = "Correlation coefficient",
                       show.diag = NULL,
                       colors = c("blue", "white", "red"),
                       outline.color = "gray",
                       outline.size = 1,
                       hc.order = FALSE,
                       hc.method = "complete",
                       lab = FALSE,
                       lab_col = "black",
                       lab_size = 4,
                       p.mat = NULL,
                       sig.level = 0.05,
                       insig = c("pch", "blank"),
                       pch = 4,
                       pch.col = "black",
                       pch.cex = 5,
                       tl.cex = 12,
                       tl.col = "black",
                       tl.srt = 45,
                       digits = 2,
                       as.is = FALSE) {
  type <- match.arg(type)
  method <- match.arg(method)
  insig <- match.arg(insig)
  if (is.null(show.diag)) {
    if (type == "full") {
      show.diag <- FALSE
    } else {
      show.diag <- FALSE
    }
  }
  if (inherits(corr, "cor_mat")) {
    # cor_mat object from rstatix
    cor.mat <- corr
    corr <- .tibble_to_matrix(cor.mat)
    p.mat <- .tibble_to_matrix(attr(cor.mat, "pvalue"))
  }
  
  if (!is.matrix(corr) & !is.data.frame(corr)) {
    stop("Need a matrix or data frame!")
  }
  corr <- as.matrix(corr)
  
  corr <- base::round(x = corr, digits = digits)
  
  if (hc.order) {
    ord <- .hc_cormat_order(corr, hc.method = hc.method)
    corr <- corr[ord, ord]
    if (!is.null(p.mat)) {
      p.mat <- p.mat[ord, ord]
      p.mat <- base::round(x = p.mat, digits = digits)
    }
  }
  
  if (!show.diag) {
    corr <- .remove_diag(corr)
    p.mat <- .remove_diag(p.mat)
  }
  
  # Get lower or upper triangle
  if (type == "lower") {
    corr <- .get_lower_tri(corr, show.diag)
    p.mat <- .get_lower_tri(p.mat, show.diag)
  } else if (type == "upper") {
    corr <- .get_upper_tri(corr, show.diag)
    p.mat <- .get_upper_tri(p.mat, show.diag)
  }
  
  # Melt corr and pmat
  corr <- reshape2::melt(corr, na.rm = TRUE, as.is = as.is)
  colnames(corr) <- c("Var1", "Var2", "value")
  corr$valuelabel <- corr$value
  corr$pvalue <- rep(NA, nrow(corr))
  corr$signif <- rep(NA, nrow(corr))
  
  if (!is.null(p.mat)) {
    p.mat <- reshape2::melt(p.mat, na.rm = TRUE)
    corr$coef <- corr$value
    corr$pvalue <- p.mat$value
    #corr$signif <- as.numeric(p.mat$value < 0.975 & p.mat$value > 0.025 ) #pp
    corr$signif <- as.numeric(p.mat$value > 0.05 ) #pvalue
    if (insig == "blank") {
      corr$value <- corr$value * (1-corr$signif)
    }
  }
  
  corr$Var1 = factor(corr$Var1,levels=rev(levels(corr$Var1)))
  #corr$Var2 = factor(corr$Var2,levels=rev(levels(corr$Var2)))
  
  corr$abs_corr <- abs(corr$value) * 10
  
  # heatmap
  p <-
    ggplot2::ggplot(
      data = corr,
      mapping = ggplot2::aes_string(x = "Var1", y = "Var2", fill = "value")
    )
  
  # modification based on method
  if (method == "square") {
    p <- p +
      ggplot2::geom_tile(color = outline.color,size=outline.size)
  } else if (method == "circle") {
    p <- p +
      ggplot2::geom_point(
        color = outline.color,
        shape = 21,
        ggplot2::aes_string(size = "abs_corr")
      ) +
      ggplot2::scale_size(range = c(4, 10)) +
      ggplot2::guides(size = "none")
  }
  
  # adding colors
  p <- p + ggplot2::scale_fill_gradient2(
    low = colors[1],
    high = colors[3],
    mid = colors[2],
    midpoint = 0,
    limit = c(-1, 1),
    space = "Lab",
    name = legend.title
  )
  
  # depending on the class of the object, add the specified theme
  if (class(ggtheme)[[1]] == "function") {
    p <- p + ggtheme()
  } else if (class(ggtheme)[[1]] == "theme") {
    p <- p + ggtheme
  }
  
  
  p <- p +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = tl.srt,
        vjust = 1,
        size = tl.cex,
        hjust = 1
      ),
      axis.text.y = ggplot2::element_text(size = tl.cex)
    ) +
    ggplot2::coord_fixed()
  
  label <- round(x = corr[, "valuelabel"], digits = digits)
  labelpvalue <- paste("(", corr[, "pvalue"],")",sep="")
  # if (!is.null(p.mat) & insig == "blank") {
  #   ns <- corr$pvalue > sig.level
  #   if (sum(ns) > 0) label[ns] <- " "
  # }
  
  # matrix cell labels
  if (lab) {
    p <- p +
      ggplot2::geom_text(
        mapping = ggplot2::aes_string(x = "Var1", y = "Var2"),vjust=-0.2,
        label = label,
        color = lab_col,
        size = lab_size
      )+
      ggplot2::geom_text(
        mapping = ggplot2::aes_string(x = "Var1", y = "Var2"),vjust=1,
        label = labelpvalue,
        color = lab_col,
        size = lab_size-1
      )
  }
  
  # matrix cell glyphs
  if (!is.null(p.mat) & insig == "pch") {
    p <- p + ggplot2::geom_point(
      data = p.mat,
      mapping = ggplot2::aes_string(x = "Var1", y = "Var2"),
      shape = pch,
      size = pch.cex,
      color = pch.col
    )
  }
  
  # add titles
  if (title != "") {
    p <- p +
      ggplot2::ggtitle(title)
  }
  
  # removing legend
  if (!show.legend) {
    p <- p +
      ggplot2::theme(legend.position = "none")
  }
  
  # removing panel
  p <- p +
    .no_panel()
  p
}



#' Compute the matrix of correlation p-values
#'
#' @param x numeric matrix or data frame
#' @param ... other arguments to be passed to the function cor.test.
#' @rdname ggcorrplot
#' @export

cor_pmat <- function(x, ...) {
  
  # initializing values
  mat <- as.matrix(x)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  
  # creating the p-value matrix
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- stats::cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  
  # name rows and columns of the p-value matrix
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  
  # return the final matrix
  p.mat
}



#+++++++++++++++++++++++
# Helper Functions
#+++++++++++++++++++++++

# Get lower triangle of the correlation matrix
.get_lower_tri <- function(cormat, show.diag = FALSE) {
  if (is.null(cormat)) {
    return(cormat)
  }
  cormat[upper.tri(cormat)] <- NA
  if (!show.diag) {
    diag(cormat) <- NA
  }
  return(cormat)
}

# Get upper triangle of the correlation matrix
.get_upper_tri <- function(cormat, show.diag = FALSE) {
  if (is.null(cormat)) {
    return(cormat)
  }
  cormat[lower.tri(cormat)] <- NA
  if (!show.diag) {
    diag(cormat) <- NA
  }
  return(cormat)
}

.remove_diag <- function(cormat) {
  if (is.null(cormat)) {
    return(cormat)
  }
  diag(cormat) <- NA
  cormat
}
# hc.order correlation matrix
.hc_cormat_order <- function(cormat, hc.method = "complete") {
  dd <- stats::as.dist((1 - cormat) / 2)
  hc <- stats::hclust(dd, method = hc.method)
  hc$order
}

.no_panel <- function() {
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank()
  )
}


# Convert a tbl to matrix
.tibble_to_matrix <- function(x) {
  x <- as.data.frame(x)
  rownames(x) <- x[, 1]
  x <- x[, -1]
  as.matrix(x)
}

ggcorrplot(M,p.mat = PP,digits=2,lab_size=5,ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"),outline.size=.5,
           insig="blank",
           
           outline.col = "black",lab=T,sig.level=0.05) +
  theme(
    axis.text.x =  element_text(angle = 50,color="black", size=15, family="serif"),
    axis.text.y =  element_text(color="black", size=15, family="serif"),
    title =  element_text(color="black", size=15, family="serif"),
    legend.title = element_text(color="black",vjust=5, size=17, family="serif"),
    legend.text =  element_text(color="black", size=15, family="serif"),
    plot.caption = element_text(hjust = 0.4, face= "italic", size=23)
  )

