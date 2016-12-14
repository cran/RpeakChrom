readChrom <- function(filepath, do.plot=TRUE, t1=0, t2=0) {
    if (!file.exists(filepath)) {
        stop("The file doesn\'t exist")
    }
    if (!is.logical(do.plot)){
        stop("do.plot parameter must be logical")
    }
    if(t1 < 0 | t2 < 0 | t1 > t2) {
        stop("No possible values for t1 and t2")
    }

    peak <- read.csv(filepath, header = FALSE)

    if (ncol(peak) != 2) {
        stop("The file must have two columns: x and y")
    }

    if (t2 > t1){
        peak <- peak[which(peak[,1] > t1 & peak[,1] < t2),]
    }

    y <- peak[,2]
    x <- peak[,1]

    if (do.plot == TRUE) {
        # library(ggplot2)
        d <- ggplot(peak,aes(x=x,y=y))+geom_point(color="red",size=0.5,show.legend = TRUE)
        d2 <- d+theme_bw()
        d22 <- d2+ggtitle(filepath)+theme(plot.title=element_text(face="bold"),
                        axis.title = element_text(size = 15)
                        ,axis.title.y=element_text(margin=margin(0,20,0,0)),
                        legend.text = element_text(size = 15),
                        legend.title = element_text(size = 15))+scale_y_continuous(name="Intensity")
        d22 <- d22+xlab("RT")
        d22+theme_minimal()
        d22+scale_x_continuous(expression(t[R] (min)))
        print(d22)
    }

    return(peak)
}


