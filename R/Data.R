#===============================================================================
#' Gillnet data
#'
#' Counts of trout from a five mesh gang of gillnets
#'
#' @docType data
#'
#' @usage data(Trout)
#'
#' @format list object with components Counts (18 by 6 - first column is lengths) and
#' Meshsize (vector of meshsizes)  respectively.
#'
#' @keywords datasets
#'
#' @references Helser T. E., et al. (1963). Estimating gillnet selectivity
#' using nonlinear response surface regression.
#' Can. J. Fish. Aquat. Sci. 55: 1328-1337.

#' @examples
#' data(Trout)
"Trout"

#===============================================================================
#Clark=read.table("../inst/extdata/Clark113_20min.dat",head=T)
#save(Clark,file="../data/Clark.RData")

#' Covered codend haddock data
#'
#' Covered codend haddock data used in Wileman et al. (1996). They are the total
#' haddock catches from three hauls each of 20 minute duration with 113 mm diamond mesh,
#' originally from Clark (1957. ICES Workshop on selectivity, Lisbon. Paper S25).
#'
#' @docType data
#'
#' @usage data(Clark)
#'
#' @format Dataframe of dimension 37 by 3. Lengthclass (cm) in the first column,
#' and haddock frequencies in the cover and codend in columns two and three, respectively.
#'
#' @keywords datasets
#'
#' @references Wileman et al. (1995) Manual of Methods of Measuring the
#' Selectivity of Towed Fishing Gears.
#' ICES Cooerative Research Report, No 215.

#' @examples
#' data(Clark)
"Clark"

#===============================================================================
#Pope=read.table("../inst/extdata/haddock.dat",head=T)
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
#' and haddock frequencies in the control and experimental codend in
#' columns two and three, respectively.
#'
#' @keywords datasets
#'
#' @references Wileman et al. (1995) Manual of Methods of Measuring the
#' Selectivity of Towed Fishing Gears. ICES Cooerative Research Report, No 215.

#' @examples
#' data(Pope)
"Pope"

#===============================================================================
#Counts=read.table("../inst/extdata/holt.dat",head=F);
#Meshsize=c(13.5,14,14.8,15.4,15.9,16.6,17.8,19)
#Holt=list(Counts=Counts,Meshsize=Meshsize)
#names(Holt$Counts)=c("lgth",paste0("M",Meshs))
#save(Holt,file="../data/Holt.RData")

#' Gillnet data
#'
#' Counts from an eight mesh gang of gillnets
#'
#' @docType data
#'
#' @usage data(Holt)
#'
#' @format list object with components Counts (11 by 9 - first column is lengths) and
#' Meshsize (vector of meshsizes)  respectively.
#'
#' @keywords datasets
#'
#' @references Holt S. J. (1963). A method for determining gear selectivity
#' and its application. ICNAF Special Pubn, 5: 106-115

#' @examples
#' data(Pope)
"Pope"


#===============================================================================
#' Permutation example data
#'
#' Simulated catch data from gears A and B. Used in Appendix 2 of Millar (2024) to
#' demonstrate use of permutation testing.
#'
#' @docType data
#'
#' @usage data(SimCatch)
#'
#' @format dataframe `SimCatch` with 1000 rows and 4 columns.
#'
#' @keywords datasets
#'
#' @references Millar R. B. (2024). Incorrect inference from size-selectivity studies
#' due to widespread misuse of bootstrap intervals.

#' @examples
#' data(SimCatch)
"SimCatch"

