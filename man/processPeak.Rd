\name{processPeak}
\alias{processPeak}
\title{
 Estimating the main parameters of a chromatographic peak.
}
\description{
 Function processPeak uses peak data to calculate the main parameters of the peak: retention time, maximum high, A60, B60, A10, B10 and area.
}
\usage{
processPeak(peak, baseline=FALSE, method, flow=FALSE, compound=FALSE, area=FALSE)
}
\arguments{
  \item{peak}{
  data frame. The input is a peak selected by the readChrom function.
}
  \item{baseline}{
  if TRUE, the function estimates a baseline using asymmetric least squares and subtracts it from the data. By default, FALSE.
  }
  \item{flow}{
  numeric. This value will be written in the output data frame.
  }
  \item{method}{
  string indicating the method used to process the peak.
  "pvmg": the upper part of the peak is filtered and fitted by the PVMG model by non-linear fitting obtaining the time and height at the peak maximum. Then cubic splines is applied obtaining the half-widths measured at 60.65\% and 10\% peak height.
  "splines": interpolation by cubic splines to obtain both half-widths left and right either at 60.65\% and 10\% peak height from the signal of the peak.
  "direct": Interpolation of the signal of the peak to estimate the half-widths.
  }
  \item{compound}{
  string. The name of the compound will be written in the output data frame.
  }
  \item{area}{
  if TRUE, the peak area is estimated using by trapezoidal numerical integration (pracma package). By default, FALSE.
  }
}
\details{
The PVMG model peak (Parabolic Variance Modified Gaussian) is a simplification of the PLMG model (Parabolic Lorentzian Modified Gaussian) removing the Lorentzian function in the PLMG model and using a Gaussian model with a variance showing only a parabolic change with time. This simplified model gives good performance in relatively narrow ranges along the peak elution. It fits accurately the upper region of the peak. The estimation of the time and height at the peak maximum is carried out by non-linear fittings of the chromatographic data to the PVMG equation.

The difference between "pvmg" and "splines" is that the height and time at the peak maximum are calculated either by the PVGM model or by interpolation using the natural signal without being modelled. Then interpolation by cubic splines is applied for obtaining the half-widths at 60.65\% and 10\% of peak height. The method "direct" obtain the half-widths directly from the natural signal by interpolation.

If an error "singular gradient matrix at initial parameter estimates" occurs, the retention times used to subset the peak by the readChrom function need to be fitted again. See examples.
}
\value{
 This function returns a data frame with the following items: compound (if supplied), flow(if supplied), retention time (tr), maximum high (Hmax), A60, B60, A10, B10 and area.
}
\author{
Manuel David Peris, Maria Isabel Alcoriza Balaguer
}
\note{
To get compounds and flows in the output data frame, they need to be witten by the user when using this function.
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
 \code{\link{vanDeemterAlternative}},
 \code{\link{vanDeemterClassic}}
}
\examples{
 ## Reading the file to get the peak
 # peak <- readChrom("file.csv", do.plot = T,
 #   t1=28, t2=29.5)
 # parameters <- processPeak(peak, baseline=FALSE, flow=0.1,
 #   method="pvmg", compound="alanine", area=TRUE)

 ## If when using processPeak it gives an error "singular gradient matrix
 ## at initial parameter estimates", the peak needs to be read again to fit
 ## retention times. For example:

 # peak <- readChrom("file.csv", do.plot = T,
 #    t1=28.5, t2=31)

  parameters <- processPeak(peak, baseline=FALSE, flow=0.1,
    method="pvmg", compound="alanine", area=TRUE)

}
