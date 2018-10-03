#' @name rdm_IBD
#' @title Calculates Isolation by Distance (IBD) residuals for each pair of individuals in a population
#'
#' @description \code{rdm_IBD}: Reads genotype data for individuals in a population together with the geographic 
#' coordinates of the individuals, models overall IBD trend for the population, and then calculates IBD residuals for each pair of individuals.
#' @param Gen_raw A string specifying the path to a file of genotype data in GENEPOP format, with alleles specified
#'  using three digits. The GENEPOP file must have the extension \emph{.gen}.
#' @param Geo_raw A string specifying the path to a text file showing the name (1st column) and geographical x and y 
#'  coordinates (2nd and 3rd columns) of each individual (each row). The first row needs to specify the column headings, 
#'  and all entries need to be tab-delimited. 
#' @param Dist_method An integer from 1 to 6 specifying the method used to calculate genetic distance among individuals. 
#'  Integers 1 to 6 correspond to the methods incorporated into the following six functions from the
#'  \code{poppr} R package, respectively: \code{diss.dist}, \code{nei.dist}, \code{rogers.dist}, \code{reynolds.dist},
#'  \code{edwards.dist}, and \code{provesti.dist}. The default value is 1.
#' @param IBD_method An integer specifying the method used to calculate the IBD residual for each pair of individuals. 
#'  A value of 1 means that the residuals are calculated by fitting a straight line to all pairs of genetic and 
#'  geographic distances, and measuring the vertical distance from each pair to the fitted line. The fitted straight line 
#'  represents the combinations of genetic and geographic distances that are expected from IBD. A value of 2 means
#'  that a non-linear curve of the form \emph{y = a + b*}(1\emph{-}exp(\emph{-}exp(\emph{c})\emph{*x})) is fitted instead of a straight 
#'  line, where \emph{a}, \emph{b} and \emph{c} are the fitted parameters. The default value is 1. The function 
#'  \code{rdm_IBD} automatically plots the data points together with the fitted line or curve. 
#' @return An object of class \code{dist}, containing a matrix showing the IBD residuals for pairs of individuals.
#' @examples
#' rdm_IBD("files/genotypes.gen", "files/coordinates.txt")
#' @export
#' @importFrom adegenet read.genepop
#' @importFrom poppr diss.dist nei.dist rogers.dist revnolds.dist edwards.dist provesti.dist


rdm_IBD <- function(Gen_raw, Geo_raw, Dist_method = 1, IBD_method = 1){
  
  gen_data <- adegenet::read.genepop(Gen_raw, ncode=3)

  if(Dist_method==1){
    gen_dist <- poppr::diss.dist(gen_data)
  }

  if(Dist_method==2){
    gen_dist <- poppr::nei.dist(gen_data)
  }

  if(Dist_method==3){
    gen_dist <- poppr::rogers.dist(gen_data)
  }

  if(Dist_method==4){
    gen_dist <- poppr::reynolds.dist(gen_data)
  }

  if(Dist_method==5){
    gen_dist <- poppr::edwards.dist(gen_data)
  }

  if(Dist_method==6){
    gen_dist <- poppr::provesti.dist(gen_data)
  }

  sample.points <- read.table(file = Geo_raw, sep = '\t', header = T, row.names = 1, fill = T)

  names(sample.points) <- c('x', 'y')

  geo.dist <- dist(sample.points[,1:2])

  IBD.data <- data.frame(cbind(as.vector(geo.dist), as.vector(gen_dist)))

  names(IBD.data) <- c('geo.dist', 'gen.dist')

  plot(IBD.data)

  if(IBD_method==1){
    IBD.model <- lm(IBD.data$gen.dist~IBD.data$geo.dist)
    IBD.res <- gen_dist - (coef(IBD.model)[1] + coef(IBD.model)[2] * geo.dist)
    lines(sort(IBD.data$geo.dist), (coef(IBD.model)[1] + coef(IBD.model)[2] * sort(IBD.data$geo.dist)), col = "red")
  }

  if(IBD_method==2){
    n_row <- nrow(IBD.data)
    IBD.model.prelim <- NLSstAsymptotic(sortedXyData(x = IBD.data[1:n_row,1], y = IBD.data[1:n_row,2]))
    IBD.model <- nls(formula = gen.dist ~ b0 + b1 * (1 - exp(-exp(lrc) * geo.dist)), data = IBD.data, start = list(b0 = IBD.model.prelim [1], b1 = IBD.model.prelim [2], lrc = IBD.model.prelim [3]))
    IBD.res <- gen_dist - (coef(IBD.model)[1] + coef(IBD.model)[2] * (1 - exp(-exp(coef(IBD.model)[3]) * geo.dist)))
    lines(sort(IBD.data$geo.dist), coef(IBD.model)[1] + coef(IBD.model)[2] * (1 - exp(-exp(coef(IBD.model)[3]) * sort(IBD.data$geo.dist))), col = "red")
  }

  return(IBD.res)

}
