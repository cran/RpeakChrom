\name{vanDeemterAlternative}
\alias{vanDeemterAlternative}
\title{
 Characterization of chromatographic columns using a new aproximation to
 vanDeemter equations.
}
\description{
 Characterization of chromatographic columns using a new aproximation to
 vanDeemter equations.
}
\usage{
vanDeemterAlternative(col, ext, dead, length, approachI=FALSE, A, B, C,
approachII=FALSE)
}
\arguments{
  \item{col}{
 data frame of the columnar measurements obtained using processPeak function.
}
  \item{ext}{
  data frame of the extracolumnar measurements obtained using processPeak
  function.
  }
  \item{dead}{
  data frame of the dead marker measurements obtained using processPeak
  function.
  }
  \item{length}{
  numeric value indicating the column length in mm.
  }
  \item{approachI}{
  If TRUE approach I is performed.
  }
  \item{A}{
  numeric value indicating the initial value of the parameter A from the van
  Deemter equation.
  }
  \item{B}{
  numeric  value indicating the initial value of the parameter B from the van
  Deemter equation.
  }
  \item{C}{
  numeric value indicating the initial value of the parameter C from the van
  Deemter equation.
  }
  \item{approachII}{
  If TRUE approachII is performed.
  }
}
\details{
In the ApproachI the parameters A, B and C from the Van Deemter equation are
obtained in two steps. First the variance in volume for a set of compounds
eluted a several flows are linearly fitted versus the retention volume to obtain
the plate height. In the second step, the obtained slopes at several flow rates
are non-linearly correlated with the linear mobile phase velocity. In the
ApproachII in the first step the parabolic behavior for the variance in volume
units for each compound in the set is fitted against the flow rate. In the
second step, the A, B and C coefficients for the different compounds are
linearly correlated with their retention volume. The slopes in the
straight-lines are the model parameters A, B and C in the Van Deemter equation.
}
\value{
For Approach I: list containing 6 items.
    Table: a summary of slope estimated values at several flows.
    Coefficients: A, B and C coefficients already fitted.
    Step I: coefficient of the linear fitting for the first step (R2).
    Correlation: R for the second step of the graphic approach (H vs u,
    non-linear fitting).
    Mean error: for the graphic (H vs u).
    RSE: square root of the estimated variance of the random error for the nls
    graphic.

For Approach II: list containing 8 items
    Table: a summary of slope estimated values at several flows.
    Coefficients: A, B and C coefficients already fitted.
    r2A: R2 for the coefficient A.
    r2B: R2 for the coefficient B.
    r2C: R2 for the coefficient C.
    MREA: mean relative prediction error for the coefficient A.
    MREB: mean relative prediction error for the coefficient B.
    MREC: mean relative prediction error for the coefficient C.
}
\author{
Manuel David Peris, Maria Isabel Alcoriza Balaguer
}
\references{
J. Baeza-Baeza, J. Torres-Lapasio, and M. Garcia-Alvarez-Coque. Approaches to
estimate the time and height at the peak maximum in liquid chromatography based
on a modified gaussian model. J.Chromatography A, 1218(10):1385-1392, 2011.

R. Caballero, M. Garcia-Alvarez-Coque, and J. Baeza-Baeza. Parabolic-lorentzian
modified gaussian model for describing and deconvolving chromatographic peaks.
J. Chromatography A, 954:59-76, 2002.

J. Foley and J. Dorsey. Equations for calculation of chromatographic figures of
merit for ideal and skewed peaks. Analytical Chemistry, 55:730-737, 1983.

E. Grushka, M. Meyers, and J. Giddings. Moment analysis for the discernment of
overlapping chromatographic peaks. Analytical Chemistry, 42:21-26, 1970.

L. He, S. Wang, and X. Geng. Coating and fusing cell membranes onto a silica
surface and their
chromatographic characteristics. Chromatographia, 54:71-76, 2001.

T. Pap and Z. Papai. Application of a new mathematical function for describing
chromatographic peaks. J. Chromatography A, 930:53-60, 2001.

J. van Deemter, F. Zuiderweg, and A. Klinkenberg. Longitudinal diffusion and
resistance to mass transfer as causes of nonideality in chromatography. Chemical
Engineering Science, 5(6):271-289, 1956.

V.B. Di Marco and G.G. Bombi. Mathematical functions for the representation of
chromatographic peaks. Journal of Chromatography A, 931:1-30, 2001.
}
\seealso{
 \code{\link{readChrom}},
 \code{\link{processPeak}},
 \code{\link{vanDeemter}}
}
\examples{
coeff1 <- vanDeemterAlternative(col = col, ext = parameters_ext,
    dead = parameters_dead, length = 150, approachI = TRUE, A = 6, B = 200,
    C = 0.04, approachII = FALSE)

coeff2 <- vanDeemterAlternative(col = col, ext = parameters_ext,
    dead = parameters_dead, length = 150, approachI = FALSE, A = 6, B = 200,
    C = 0.04, approachII = TRUE)
}
