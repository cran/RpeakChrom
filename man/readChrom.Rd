\name{readChrom}
\alias{readChrom}
\title{
Reading chromatographic peak from existing csv or txt files.
}
\description{
Function readChrom allows to read a chromatogram from a csv or txt file with two columns (x, y) and subset the peak of interest if the user gives the limiting retention times. Also, this function can draw the peak to check if the given retention times are correct.
}
\usage{
readChrom(filepath, do.plot=TRUE, t1=0, t2=0)
}
\arguments{
  \item{filepath}{
 string indicating the path to the file. The first column represent the retention times and the second one the intensities.
}
  \item{do.plot}{
 if TRUE, the function prints the chromatogram in the interval selected by t1 and t2. By default, it is TRUE.
  }
  \item{t1}{
 numeric. Filter for measurements with retention time >= t1.
  }
  \item{t2}{
 numeric. Filter for measurements with retention time <= t2.
  }
}
\details{
Setting t1 and t2 allows you to filter your peak in the chromatogram.
To insurance the correct selection of the peak, the argument "do.plot" should be used.
}
\value{
The return value is a data frame containing two columns, the retention time and intensity.
}
\author{
Manuel David Peris, Maria Isabel Alcoriza Balaguer
}
\references{
J. Baeza-Baeza, J. Torres-Lapasio, and M. Garcia-Alvarez-Coque. Approaches to estimate the time and height at the peak maximum in liquid chromatography based on a modified gaussian model. J.Chromatography A, 1218(10):1385-1392, 2011.

R. Caballero, M. Garcia-Alvarez-Coque, and J. Baeza-Baeza. Parabolic-lorentzian modified gaussian model for describing and deconvolving chromatographic peaks. J. Chromatography A, 954:59-76, 2002.

J. Foley and J. Dorsey. Equations for calculation of chromatographic figures of merit for ideal and skewed peaks. Analytical Chemistry, 55:730-737, 1983.

E. Grushka, M. Meyers, and J. Giddings. Moment analysis for the discernment of overlapping
chromatographic peaks. Analytical Chemistry, 42:21-26, 1970.

L. He, S. Wang, and X. Geng. Coating and fusing cell membranes onto a silica surface and their
chromatographic characteristics. Chromatographia, 54:71-76, 2001.

T. Pap and Z. Papai. Application of a new mathematical function for describing chromatographic
peaks. J. Chromatography A, 930:53-60, 2001.

J. van Deemter, F. Zuiderweg, and A. Klinkenberg. Longitudinal diffusion and resistance to mass transfer as causes of nonideality in chromatography. Chemical Engineering Science, 5(6):271-289, 1956.

V.B. Di Marco and G.G. Bombi. Mathematical functions for the representation of chromatographic peaks. Journal of Chromatography A, 931:1-30, 2001.
}
\seealso{
 \code{\link{processPeak}},
 \code{\link{vanDeemterAlternative}},
 \code{\link{vanDeemterClassic}}
}
\examples{
# Substitute the file path argument for a csv or txt file

# To see the whole chromatogram
# peak <- readChrom("example_file.csv", do.plot = TRUE)


# To subset the peak make use of t1 and t2 arguments, for example:
# peak <- readChrom("example_file.csv", do.plot = TRUE, t1 = 2, t2 = 2.5)
}
