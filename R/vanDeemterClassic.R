vanDeemterClassic<-function(col,  ext, dead, length, A, B, C, Foley=FALSE, GG=FALSE, do.plot=TRUE) {
    if (!is.numeric(A) | !is.numeric(B) | !is.numeric(C)){
        stop("A, B and C parameters must be numeric")
    }
    if (!is.logical(Foley) | !is.logical(GG)){
        stop("Foley and GG parameters must be logical")
    }
    if (Foley == TRUE & GG == TRUE){
        stop("Just one approach (Foley or GG) at the time is possible")
    }
    if (Foley == FALSE & GG == FALSE){
        stop("Foley or GG parameter must be TRUE")
    }
    if (!is.numeric(length)){
        stop("length parameter must be numeric")
    }
    if (!is.data.frame(col)){
        stop("col parameter must be a data frame")
    }
    if (!is.data.frame(ext)){
        stop("ext parameter must be a data frame")
    }
    if (!is.data.frame(dead)){
        stop("dead parameter must be a data frame")
    }

    col <- col[order(col$flow),]
    ext <- ext[order(ext$flow),]
    dead <- dead[order(dead$flow),]

    t <- table(c(ext$flow, col$flow, dead$flow))
    del <- as.numeric(names(which(t < max(as.vector(t)))))

    ext <- ext[!ext$flow %in% del,]
    col <- col[!col$flow %in% del,]
    dead <- dead[!dead$flow %in% del,]

    Flow <- col[,"flow"]
    tr <- col[,"tr"]
    t0 <- dead[,"tr"]
    tre <- ext[,"tr"]
    velocity2 <- length/t0

    if (Foley == TRUE){
        var <- ((col$A10+col$B10)^2)/(1.764*((col$B10/col$A10)^2)-11.15*(col$B10/col$A10)+28)
        vare <- ((ext$A10+ext$B10)^2)/(1.764*((ext$B10/ext$A10)^2)-11.15*(ext$B10/ext$A10)+28)

        H <- ((var-vare)/(tr*tr-tre*tre))*length*1000

        # library(ggplot2)

        data <- as.data.frame(cbind(velocity2,H))
        d <- ggplot(data,aes(x=data[,1],y=data[,2]))+
            geom_point(color="black",size=2,show.legend = TRUE)+
            geom_smooth(method="nls",
                        formula=y ~ A + (B/x) + C*x,
                        method.args = list(start=c(A=A,B=B,C=C)),
                        se=FALSE)
        d2 <- d+theme_bw()
        d22 <- d2+theme(plot.title=element_text(face="bold"),
                        axis.title = element_text(size = 15)
                        ,axis.title.y=element_text(margin=margin(0,20,0,0)),
                        legend.text = element_text(size = 15),
                        legend.title = element_text(size = 15))+scale_y_continuous(expression(H~(mu~m)))
        d22 <- d22+theme_minimal()
        d22 <- d22+scale_x_continuous(expression(u(mm/min)))
        if (do.plot == TRUE){
            print(d22)
        }

        x <- data[,1]
        y <- data[,2]
        fit <- nls(y ~ A + (B/x) + C*x, start=list(A=A,B=B,C=C))

        return (coefficients(fit))
    }

    if (GG == TRUE){
        col <- col[order(col$flow),]
        ext <- ext[order(ext$flow),]
        dead <- dead[order(dead$flow),]

        t <- table(c(ext$flow, col$flow, dead$flow))
        del <- as.numeric(names(which(t < max(as.vector(t)))))

        ext <- ext[!ext$flow %in% del,]
        col <- col[!col$flow %in% del,]
        dead <- dead[!dead$flow %in% del,]

        var <- (((col$A60^3)+(col$B60^3))/(col$A60+col$B60))-(((2*(col$B60-col$A60)^2))/pi)
        vare <- (((ext$A60^3)+(ext$B60^3))/(ext$A60+ext$B60))-(((2*(ext$B60-ext$A60)^2))/pi)
        tr <- tr+(2*(col$B60-col$A60))/(sqrt(2*pi))
        tre <- tre+(2*(ext$B60-ext$A60))/(sqrt(2*pi))

        H <- ((var-vare)/(tr*tr-tre*tre))*length*1000

        # library(ggplot2)

        data <- as.data.frame(cbind(velocity2,H))
        d <- ggplot(data,aes(x=data[,1],y=data[,2]))+
            geom_point(color="black",size=2,show.legend = TRUE)+
            geom_smooth(method="nls",
                        formula=y ~ A + (B/x) + C*x,
                        method.args = list(start=c(A=A,B=B,C=C)),
                        se=FALSE)
        d2 <- d+theme_bw()
        d22 <- d2+theme(plot.title=element_text(face="bold"),
                        axis.title = element_text(size = 15)
                        ,axis.title.y=element_text(margin=margin(0,20,0,0)),
                        legend.text = element_text(size = 15),
                        legend.title = element_text(size = 15))+scale_y_continuous(expression
                                                                                   (H~(mu~m)))
        d22 <- d22+theme_minimal()
        d22 <- d22+scale_x_continuous(expression(u(mm/min)))
        if (do.plot == TRUE){
            print(d22)
        }

        x <- data[,1]
        y <- data[,2]
        fit <- nls(y ~ A + (B/x) + C*x, start=list(A=A,B=B,C=C))

        return(coefficients(fit))
    }
}
