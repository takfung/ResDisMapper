\name{ResDisMapper-package}
\alias{ResDisMapper-package}
\alias{ResDisMapper}
\docType{package}
\title{
\packageTitle{ResDisMapper}
}
\description{
The Resistance to Dispersal Mapper (ResDisMapper) package provides functions to calculate and map environmental resistance to dispersal over a landscape, for a population of individuals. Firstly, there is a function called \code{rdm_IBD}, which reads in the geographic coordinates and genetic data for the individuals and computes the geographic and genetic distance between each pair of individuals. This function compares the geographic and genetic distance for each pair of individuals to expected values in the case of Isolation by distance (IBD), and calculates the corresponding residuals. Secondly, there is a function \code{rdm_residual}, which takes the residuals and visualizes them using a 3-D plot, and outputs the residuals together with line segments representing pairs of individuals over the landscape. Thirdly, there is a function \code{rdm_resistance}, which uses the residuals and line segments to calculate the environmental resistance to dispersal over the landscape. Lastly, the function \code{rdm_mapper} visualizes the resistance values using a raster map, with contour lines indicating areas of high and low resistance.
}
\author{
\packageAuthor{ResDisMapper}

Maintainer: \packageMaintainer{ResDisMapper}
}
\references{
Tang, Q., Fung, T., Rheindt, F.E. (In review). ResDisMapper: An R package for fine-scale mapping of resistance to dispersal.
}
