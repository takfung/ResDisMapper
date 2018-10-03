#' @name rdm_mapper
#' @title Visualizes the environmental resistance to dispersal over a landscape, for individuals in a population
#'
#' @description \code{rdm_mapper}: Reads in the environmental resistance to dispersal for individuals of a population, over
#'  a grid of cells defining the landscape containing the individuals, and then creates a raster map of resistance.
#'  Red and green contour lines are drawn on the map, delineating genetic barriers and corridors, respectively. In 
#'  addition, blue contour lines are drawn to delineate areas with resistances that are statistically different
#'  from 0. For further descriptions of the meaning of genetic barrier, genetic corridor, and resistances
#'  that are statistically different from 0, please refer to the description of the function \code{rdm_resistance}.
#' @param F.df A data frame containing eight columns with resistance information for each grid cell in a landscape. 
#'  Each row corresponds to one cell, and the eight different columns refer to the resistance, the number of 
#'  intersecting line segments (lines connecting pairs of individuals in the landscape), the x and y coordinates specifying 
#'  the location of the mid-point of the cell, the lower and upper limits of the confidence interval for the resistance, the 
#'  sign of the product of the upper and lower limits of the confidence interval, and the percentile of the null distribution of resistances 
#'  corresponding to the observed resistance. A .csv file with this information is produced by the function
#'  \code{rdm_resistance}. 
#' @param Geo_raw A string specifying the path to a text file showing the name (1st column) and geographical x and y 
#'  coordinates (2nd and 3rd columns) of each individual (each row). The first row needs to specify the column headings, 
#'  and all entries need to be tab-delimited.
#' @param r_size Specifies the size of each grid cell in the plotted raster map of resistances, which can be adjusted 
#' to match the size of the plotted landscape. Default value of 5. 
#' @param p_signf Specifies the percentiles of the null distribution of resistances used to define genetic barriers and 
#'  corridors. A genetic barrier is a grid cell with a resistance that is above the 100*(\code{1-p_signf})\% percentile, 
#'  whereas a genetic corridor is a grid cell with a resistance that is below the 100*(\code{p_signf})\% percentile.
#'  Default value of 0.05.
#' @param p_size Specifies the size of the sampling points in the plotted raster map of resistances. Default value of 2.
#' @return Function plots a map of resistance, but does not return an output object.
#' @examples
#' F.df <- read.csv("files/resistance_map.csv") 
#' rdm_mapper(F.df, "files/coordinates.txt")
#' @export
#' @import ggplot2

rdm_mapper <- function(F.df, Geo_raw, r_size = 5, p_signf = 0.05, p_size = 2){

  sample.points <- read.table(file = Geo_raw, sep = '\t', header = T, row.names = 1, fill = T)

  p_signf_u = 1 - p_signf

  ggplot2::ggplot(F.df, ggplot2::aes(x = x, y = y, colour=resistance)) + 
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))+
    ggplot2::geom_point(size = r_size, shape = 15) +
    ggplot2::scale_color_gradient2(low="forestgreen", high="red3", mid="white", midpoint=0)+
    ggplot2::stat_contour(ggplot2::aes(z = Prob), bins=4, size=1, col="red", breaks = c(p_signf_u))+
    ggplot2::stat_contour(ggplot2::aes(z = Prob), bins=4, size=1, col="green3", breaks = c(p_signf))+
    ggplot2::stat_contour(ggplot2::aes(z = sign.check), bins=4, size=1, col="blue", breaks=c(0))+
    ggplot2::geom_point(data = sample.points, cex = p_size, shape = 21, color = "black", fill = "yellow",stroke = 2)
  
}
