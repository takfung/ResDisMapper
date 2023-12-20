#' @name rdm_resistance
#' @title Produces resistance to dispersal over a landscape, for individuals in a population
#'
#' @description \code{rdm_resistance}: Reads in IBD residuals for pairs of individuals in a population, represented
#'  as line segments connecting the pairs of individuals over a defined landscape, and then applies a novel method 
#'  (Tang et al., Molecular Ecology Resources) to the residuals to calculate resistance to dispersal over the landscape. 
#'  Essentially, the method splits the landscape into a grid of cells and for each cell, calculates the 
#'  resistance as the average of the IBD residuals of all line segments that intersect the cell.
#' @param IBD.res An object of class \code{dist}, containing a matrix with the IBD residuals for pairs of 
#'  individuals in a population.
#' @param Res_SLDF An object of class \code{SpatialLinesDataFrame}, containing the IBD residuals for each pair of 
#'  individuals in a population and the coordinates of line segments connecting the pairs of individuals over the 
#'  defined landscape. This object is produced by the function \code{rdm_residual}.
#' @param nrows Number of grid cells in each row of the raster map of resistance. The default value is 30.
#' @param ncols Number of grid cells in each column of the raster map of resistance. The default value is 30.
#' @param conf_intervals The coverage of the confidence interval generated for each resistance value in a grid cell, 
#'  expressed as a proportion. This interval measures the uncertainty in the (observed) resistance in each cell. 
#'  If a cell has no intersecting line segments or only one, then a confidence interval cannot be calculated and the cell
#'  is no longer considered for further calculation. If the confidence interval does not overlap 0, then the resistance 
#'  is statistically different from 0. The default value is 0.95 (corresponding to 95\% intervals).
#' @param random_rep Number of random resamples of the IBD residuals used to construct the null distribution of 
#'  resistances in each grid cell. For a grid cell with \emph{n} intersecting line segments and \emph{n} corresponding IBD residuals,
#'  a random resample consists of randomly sampling \emph{n} IBD residuals from the set of all IBD residuals over the entire
#'  landscape, without replacement. If the observed resistance is above a threshold percentile of the null
#'  distribution, then the cell is considered to have a high resistance that is statistically significant. If instead the observed 
#'  resistance is below another threshold percentile of the null distribution, then the cell is considered to have a low resistance 
#'  that is statistically significant. The default value is 1,000.
#' @param outputfile Name of the .csv file that is produced, containing resistance information in the output data frame. 
#'  The default is 'resistance_map.csv'.
#' @return A data frame with resistance information for each grid cell in the landscape with more than one 
#'  intersecting line segment. Each row corresponds to one cell, and the eight different columns refer to the 
#'  resistance, the number of intersecting line segments, the x and y coordinates specifying the location of the mid-point of the cell, 
#'  the lower and upper limits of the confidence interval for the resistance, the sign of the product of the lower and upper 
#'  limits of the confidence interval, and the percentile of the null distribution of resistances corresponding 
#'  to the observed resistance. When calculating this data, the function \code{rdm_resistance} prints out messages
#'  showing the current stage of calculation. There are four stages: Calculating (1) the resistances, (2) the
#'  numbers of intersecting line segments, (3) the lower limits of the confidence intervals for the resistances,
#'  and (4) the upper limits of the confidence intervals for the resistances. 
#' @references Tang, Q., Fung, T., Rheindt, F.E. ResDisMapper: An R package for fine-scale mapping of resistance to dispersal. Molecular Ecology Resources.
#' @examples
#' IBD_res <- rdm_IBD("files/genotypes.gen", "files/coordinates.txt")
#' Res_SLDF <- rdm_residuals(IBD_res, "files/coordinates.txt")
#' rdm_resistance(IBD_res, Res_SLDF)
#' @export
#' @importFrom raster raster extent projection rasterize
#' @importFrom Rmisc CI

rdm_resistance <- function(IBD.res, Res_SLDF, nrows = 30, ncols = 30, conf_intervals = 0.95, random_rep = 1000, outputfile = 'resistance_map.csv'){
  
  Res_raster <- raster::raster(raster::extent(Res_SLDF), crs=raster::projection(Res_SLDF), nrows, ncols)
  
  #To calculate the resistance of each grid cell (mean of the IBD residuals corresponding to line segments intersecting the cell).
  
  print("Calculating resistance (1/4)")
  
  Res_mean <- raster::rasterize(Res_SLDF, Res_raster,Res_SLDF@data, fun=mean, progress = 'text')
  
  #To calculate the number of IBD residuals (i.e. intersecting line segments) used to calculate resistance of each grid cell.
  
  print("Counting number of intersecting line segments (2/4)")
  
  Res_count <- raster::rasterize(Res_SLDF, Res_raster, Res_SLDF@data, fun='count', progress = 'text')
  
  #To estimate uncertainty by calculating the confidence interval of resistance in each grid grill.
  
  print("Calculating lower bound of confidence interval (3/4)")
  
  Res_R_ci_L <- raster::rasterize(Res_SLDF, Res_raster, Res_SLDF@data, fun = function(x){L<-Rmisc::CI(x, ci = conf_intervals)
  return(L[3])}, progress = 'text')
  
  print("Calculating upper bound of confidence interval (4/4)")
  
  Res_R_ci_U<-raster::rasterize(Res_SLDF, Res_raster, Res_SLDF@data,fun = function(x){L<-Rmisc::CI(x, ci = conf_intervals)
  return(L[1])}, progress = 'text') 
  
  R.spdf <- as(Res_mean, "SpatialPixelsDataFrame")
  R.df <- as.data.frame(R.spdf)
  
  C.spdf <- as(Res_count, "SpatialPixelsDataFrame")
  C.df <- as.data.frame(C.spdf)
  
  RC.df <- cbind(R.df[,1], C.df)
  RC.df <- RC.df[!(RC.df[,2]==1),]
  
  CI.spdf.l <- as(Res_R_ci_L, "SpatialPixelsDataFrame")
  CI.df.l <- as.data.frame(CI.spdf.l)
  
  CI.spdf.u <- as(Res_R_ci_U, "SpatialPixelsDataFrame")
  CI.df.u <- as.data.frame(CI.spdf.u)
  
  F.df <- cbind(RC.df, CI.df.l[,1], CI.df.u[,1])
  
  F.df [,7] <- sign(F.df[,5]*F.df[,6]) 
  
  Rfunction <- function(x,y){
    iteration <- replicate(random_rep, {
      a<-sample(na.omit(as.numeric(IBD.res)),y)
      mean(a)})
    f <- ecdf(iteration)
    return(f(x))
  }
  
  Perct <- data.frame()
  
  for (i in 1:length(F.df[,1])){
    Perct[i,1] <- Rfunction(F.df[i,1], F.df[i,2])
  }
  
  F.df <- cbind(F.df, Perct)
  
  names(F.df) <- c('resistance', 'intersects', 'x', 'y', 'CI.l', 'CI.u', 'sign.check', 'Prob')
  
  write.csv(F.df, outputfile, row.names=FALSE)
  
  return(F.df)
  
}
