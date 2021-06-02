#Clark=read.table("Clark113_20min.dat",head=T)[,c(1,3,2)]
#save(Clark,file="../data/Clark.RData")

#' Covered codend haddock data
#'
#' Covered codend haddock data used in Wileman et al. (1996). They are the total haddock catches from
#' three hauls each of 20 minute duration with 113 mm diamond mesh, originally from Clark
#' (1957. ICES Workshop on selectivity, Lisbon. Paper S25).
#'
#' @docType data
#'
#' @usage data(Clark)
#'
#' @format Dataframe of dimension 37 by 3. Lengthclass (cm) in the first column,
#' and haddock frequencies in the codend and cover in columns two and three, respectively.
#'
#' @keywords datasets
#'
#' @references Wileman et al. (1995) Manual of Methods of Measuring the Selectivity of Towed Fishing Gears.
#' ICES Cooerative Research Report, No 215.

#' @examples
#' data(Clark)
"Clark"

#Pope=read.table("haddock.dat",head=T)[,c(1,3,2)]
#save(Pope,file="../data/Pope.RData")

#' Alternate hauls haddock data
#'
#' Alternate hauls haddock data used in Wileman et al. (1996) and Pope et al. (1975).
#' Codend mesh is 87 mm diamond.
#'
#' @docType data
#'
#' @usage data(Pope)
#'
#' @format Dataframe of dimension 24 by 3. Lengthclass (cm) in the first column,
#' and haddock frequencies in the codend and cover in columns two and three, respectively.
#'
#' @keywords datasets
#'
#' @references Wileman et al. (1995) Manual of Methods of Measuring the Selectivity of Towed Fishing Gears.
#' ICES Cooerative Research Report, No 215.

#' @examples
#' data(Pope)
"Pope"
