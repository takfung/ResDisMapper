#' @name rdm_residual
#' @title Creates and plots IBD residuals for pairs of individuals in a population, represented as line segments over a 
#'  defined landscape
#'
#' @description \code{rdm_residual}: Reads the IBD residuals for each pair of individuals in a population (as calculated 
#'  by function \code{rdm_IBD}) together with their geographical coordinates, and then creates a 
#'  3-D graph showing the residuals for each pair of individuals, represented as a line segment over the
#'  landscape considered. The 3-D graph can be rotated to obtain a better view.
#' @param IBD.res An object of class \code{dist}, containing a matrix with the IBD residuals for pairs of 
#'  individuals in a population.
#' @param Geo_raw A string specifying the path to a text file showing the name (1st column) and geographical x and y 
#'  coordinates (2nd and 3rd columns) of each individual (each row). The first row needs to specify the column headings, 
#'  and all entries need to be tab-delimited.
#' @param min.dist IBD residuals are only calculated for pairs of individuals that are separated by a distance greater
#'  than \code{min.dist}. The default value is 1.
#' @param max.dist IBD residuals are only calculated for pairs of individuals that are separated by a distance less
#'  than \code{max.dist}. The default value is Inf.
#' @param n_resolution Specifies the number of cells that the landscape is divided into, for the purposes of creating
#'  the 3-D plot. The landscape is divided into \code{n_resolution} cells along both coordinate axes. Individuals that
#'  appear in the same cell are represented by a single point on the plot. The default value is 50.
#' @param proj The coordinate system that is used for plotting. The default is EPSG 4326.
#' @return An object of class \code{SpatialLinesDataFrame}, containing coordinates of the line segments joining each 
#'  pair of individuals and the corresponding IBD residuals.
#' @examples
#' IBD_res <- rdm_IBD("files/genotypes.gen", "files/coordinates.txt")
#' rdm_residual(IBD_res, "files/coordinates.txt")
#' @export
#' @importFrom rgl open3d points3d segments3d axes3d grid3d title3d
#' @importFrom sp CRS Lines Line SpatialLines SpatialLinesDataFrame

rdm_residual <- function(IBD.res, Geo_raw, min.dist = 1, max.dist = Inf, n_resolution = 50, proj=sp::CRS("+init=epsg:4326")){

  sample.points <- read.table(file = Geo_raw, sep = '\t', header = T, row.names = 1, fill = T)

  names(sample.points) <- c('x', 'y')

  geo.dist <- dist(sample.points[,1:2])

  n <- dim(sample.points)[1]
  x1 <- as.dist(matrix(data = sample.points$x, nrow = n, ncol = n, byrow = F))
  y1 <- as.dist(matrix(data = sample.points$y, nrow = n, ncol = n, byrow = F))
  x2 <- as.dist(matrix(data = sample.points$x, nrow = n, ncol = n, byrow = T))
  y2 <- as.dist(matrix(data = sample.points$y, nrow = n, ncol = n, byrow = T))
  
  Res_segments <- data.frame(as.vector(x1[(geo.dist > min.dist) & (geo.dist < max.dist)]),
    as.vector(y1[(geo.dist > min.dist) & (geo.dist < max.dist)]),
    as.vector(IBD.res[(geo.dist > min.dist) & (geo.dist < max.dist)]),
    as.vector(x2[(geo.dist > min.dist) & (geo.dist < max.dist)]),
    as.vector(y2[(geo.dist > min.dist) & (geo.dist < max.dist)]))

  names(Res_segments) <- c('x1', 'y1', 'z', 'x2', 'y2')

  Res_segments <- Res_segments[complete.cases(Res_segments), ]

  k = (max(Res_segments$z)-min(Res_segments$z))/(max(sample.points$x)-min(sample.points$x))

  rgl::open3d(scale = c(1,1,0.5/k))
  x_vec = as.vector(t(Res_segments[,c(1,4)]))
  y_vec = as.vector(t(Res_segments[,c(2,5)]))
  z_vec = as.vector(t(Res_segments[,c(3,3)]))
  rgl::points3d(x = x_vec, y = y_vec,z = z_vec)
  myColorRamp <- function(colors, values){
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
  }
  z_vec_nonpos_ind = which(z_vec<=0)
  x_vec_low = x_vec[z_vec_nonpos_ind]
  y_vec_low = y_vec[z_vec_nonpos_ind]
  z_vec_low = z_vec[z_vec_nonpos_ind]
  cols_low <- myColorRamp(c("dark green", "white"), z_vec_low) 
  rgl::segments3d(x = x_vec_low, y = y_vec_low, z = z_vec_low, col = cols_low)
  z_vec_pos_ind = which(z_vec>0)
  x_vec_hi = x_vec[z_vec_pos_ind]
  y_vec_hi = y_vec[z_vec_pos_ind]
  z_vec_hi = z_vec[z_vec_pos_ind]
  cols_hi <- myColorRamp(c("white", "red"), z_vec_hi) 
  rgl::segments3d(x = x_vec_hi, y = y_vec_hi, z = z_vec_hi, col = cols_hi)
  rgl::axes3d()
  rgl::grid3d(c("z"), at = NULL, col = "gray", lwd = 1, lty = 1, n = n_resolution)
  rgl::title3d(xlab="X", ylab="Y", zlab="Residuals")

  Res_begin <- cbind(Res_segments$x1, Res_segments$y1)

  Res_end <- cbind(Res_segments$x2, Res_segments$y2)

  l <- vector("list", nrow(Res_begin))

  for (i in seq_along(l)) {
    l[[i]] <- sp::Lines(list(sp::Line(rbind(Res_begin[i, ], Res_end[i,]))), as.character(i))
  }

  Res_lines <- sp::SpatialLines(l, proj4string=proj)

  df <- data.frame(Res_segments$z)

  Res_SLDF <- sp::SpatialLinesDataFrame(Res_lines, data=df)

  return(Res_SLDF)

}
