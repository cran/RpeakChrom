\name{vanDeemterClassic}
\alias{vanDeemterClassic}
\title{
 Characterization of chromatographic columns using vanDeemter equations.
}
\description{
 Characterization of chromatographic columns using Van Deemter equations.
}
\usage{
vanDeemterClassic(col,  ext, dead, length, A, B, C, Foley=FALSE, GG=FALSE,
    do.plot=TRUE)
}
\arguments{
  \item{col}{
 data frame of the columnar measurements obtained using processPeak function.
}
  \item{ext}{
  data frame of the extracolumnar measurements obtained using processPeak function.
  }
  \item{dead}{
  data frame of the dead marker measurements obtained using processPeak function.
  }
  \item{length}{
  numeric value indicating the column length in mm.
  }
  \item{A}{
  numeric value indicating the initial value of the parameter A from the van Deemter equation.
  }
  \item{B}{
  numeric  value indicating the initial value of the parameter B from the van Deemter equation.
  }
  \item{C}{
  numeric value indicating the initial value of the parameter C from the Van Deemter equation.
  }
  \item{Foley}{
  if TRUE, Foley and Dorsey approach is used to estimate the variance from the half-widths measured at 10\% peak height.
  }
  \item{GG}{
  if TRUE the variance and retention time are calculated by using the Grushka and Giddings approach.
  }
  \item{do.plot}{
  logical
  }
}
\details{
The Van Deemter approach has been widely used in column performance in HPLC from the information obtained in the elution of probe compounds at different flow rates, which relate the column plate height to the linear mobile phase velocity given solute, column and mobile phase composition.
In this function the approaches for obtaining the retention time and variance are based on Grushka and Giddings, or Foley and Dorsey.
The Grushka and Giddings approach make use of the half-widths measured at 60.65\% peak height whereas Foley and Dorsey approach is based on the measurements at 10\% peak height where the peak asymmetry is higher.
The theoretical plate heigh (H) is determined according to the Martin and Synge plate model taking into account the measurement of the extra-column contribution.
}
\value{
A, B and C coefficients of vanDeemter equation already fitted.
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
 \code{\link{readChrom}},
 \code{\link{processPeak}},
 \code{\link{vanDeemterAlternative}}
}
\examples{
ggmetoxi <- vanDeemterClassic(col=parameters_col_metoxi, ext=parameters_ext,
              dead=parameters_dead, length=150, A=6, B=200, C=0.04, GG=TRUE,
              Foley=FALSE, do.plot=TRUE)

foleymetoxi <- vanDeemterClassic(col=parameters_col_metoxi, ext=parameters_ext,
              dead=parameters_dead, length=150, A=6, B=200, C=0.04, GG=FALSE,
              Foley=TRUE, do.plot=TRUE)
}
