processPeak <- function(peak, baseline = FALSE, method, flow = FALSE,
    compound = FALSE, area = FALSE){
    if (is.data.frame(peak) == FALSE) {
        stop("peak must be a data.frame with 2 columns")
    }
    if (ncol(peak) != 2) {
        stop("peak must be a data.frame with 2 columns")
    }
    if (baseline != TRUE & baseline != FALSE) {
        stop("baseline parameter must be logical")
    }
    if (method != "direct" & method != "pvmg" & method != "splines") {
        stop("method parameter must be one of these: \"direct\",\"pvmg\",\"splines\"")
    }
    if (is.numeric(flow) == FALSE & flow != FALSE) {
        stop("flow parameter must be numeric or FALSE")
    }
    if (is.character(compound) == FALSE & compound != FALSE) {
        stop("compound parameter must be character or FALSE")
    }
    if (is.logical(area) == FALSE) {
        stop("area parameter must be logical")
    }
    if (baseline == TRUE) {
        peak[, 2] <- baseline.corr(peak[, 2])
        peak <- data.frame(peak)
    }

    y <- peak[, 2]
    x <- peak[, 1]
    xmax <- x[y == max(y)]
    H <- max(y) - min(y)
    x1 <- x[x < xmax][which.min(abs(y[x < xmax] - max(y) * 0.6065))]
    x2 <- x[x > xmax][which.min(abs(y[x > xmax] - max(y) * 0.6065))]
    A60 <- xmax - x1
    B60 <- x2 - xmax
    x3 <- x[x < xmax][which.min(abs(y[x < xmax] - max(y) * 0.1))]
    x4 <- x[x > xmax][which.min(abs(y[x > xmax] - max(y) * 0.1))]
    A10 <- xmax - x3
    B10 <- x4 - xmax
    if (method == "pvmg") {
        izq <- y[x < xmax]
        limizq <- izq[which.min(abs(izq - (abs(max(y) - max(y) * 0.25))))]
        xizq <- x[y == limizq]
        drcha <- y[x > xmax]
        limdrcha <- drcha[which.min(abs(drcha - (abs(max(y) - max(y) * 0.25))))]
        xdrcha <- x[y == limdrcha]
        df <- data.frame(x, y)
        peak20 <- df[df$x < xdrcha & df$x > xizq, ]
        plot(peak20$x, peak20$y, xlab = "RT", ylab = "Intensity")
        title(main = paste("compound :", compound, "/ flow: ", flow),
            col.main = "black")
        x_a <- peak20$x
        y_a <- peak20$y
        s0 <- (A60 + B60)/2
        b <- (B60 - A60) * ((s0^2)/A60 * B60)
        c <- 1 - ((s0^2)/A60 * B60)
        pvmg <- function(x, H, xmax, b, c, s0){H * exp(-0.5 *
                ((x - xmax)^2)/s0^2 + b * (x - xmax) + c * (x - xmax)^2)}
        set.seed(123)
        mp <- nlsLM(y_a ~ pvmg(x_a, H, xmax, b, c, s0), start = list(H = H,
            xmax = xmax, b = b, c = c, s0 = s0), control = list(tol = 0,
                maxiter = 500))
        lines(peak20$x, predict(mp), col = "red")
        cor <- cor(y_a,predict(mp))
        ME <- mean(abs(y_a-predict(mp))/abs(y_a))*100
        summary <- summary(mp)$sigma
        H <- coefficients(mp)[1]
        xmax <- coefficients(mp)[2]
    }
    if (method == "pvmg" | method == "splines"){
        izq2 <- y[x < xmax]
        limizq2 <- izq2[which.min(abs(izq2 - (abs(H - H * 0.9))))]
        xizq2 <- x[y == limizq2]
        drcha2 <- y[x > xmax]
        limdrcha2 <- drcha2[which.min(abs(drcha2 - (abs(H - H * 0.9))))]
        xdrcha2 <- x[y == limdrcha2]
        int <- peak[peak[, 1] < xdrcha2 & peak[, 1] > xizq2,]
        colnames(int) <- c("x", "y")
        eee <- spline(int$x, int$y)
        if (method == "splines") {
            plot(int$x, int$y, xlab = "RT", ylab = "Intensity")
            title(main = paste("compound :", compound, "/ flow: ", flow),
                col.main = "black")
            lines(eee$x, eee$y, col = "red")
            }
        x1 <- eee$x[eee$x < xmax][which.min(abs(eee$y[eee$x < xmax] - H *
                0.6065))]
        x2 <- eee$x[eee$x > xmax][which.min(abs(eee$y[eee$x > xmax] - H *
                0.6065))]
        A60 <- xmax - x1
        B60 <- x2 - xmax
        x3 <- eee$x[eee$x < xmax][which.min(abs(eee$y[eee$x < xmax] - H * 0.1))]
        x4 <- eee$x[eee$x > xmax][which.min(abs(eee$y[eee$x > xmax] - H * 0.1))]
        A10 <- xmax - x3
        B10 <- x4 - xmax
        }
    if (area != FALSE) {
        area <- trapz(x, y)
        }
    if (compound != FALSE & area != FALSE & flow != FALSE) {
        res <- data.frame(compound = as.vector(compound), flow = as.vector(flow),
            tr = as.vector(xmax), Hmax = as.vector(H), A60 = as.vector(A60),
            B60 = as.vector(B60), A10 = as.vector(A10), B10 = as.vector(B10),
            area = as.vector(area), RSE = as.vector(summary), MeanError = ME,
            correlation = cor)
        } else if (compound != FALSE & area == FALSE & flow != FALSE){
            res <- data.frame(compound = as.vector(compound),
                flow = as.vector(flow), tr = as.vector(xmax),
                Hmax = as.vector(H), A60 = as.vector(A60), B60 = as.vector(B60),
                A10 = as.vector(A10), B10 = as.vector(B10),
                RSE = as.vector(summary), MeanError = ME, correlation = cor)
            } else if (compound == FALSE & area == FALSE & flow != FALSE){
                res <- data.frame(flow = as.vector(flow), tr = as.vector(xmax),
                    Hmax = as.vector(H), A60 = as.vector(A60),
                    B60 = as.vector(B60), A10 = as.vector(A10),
                    B10 = as.vector(B10), RSE = as.vector(summary), MeanError=ME,
                    correlation = cor)
                } else if (compound != FALSE & area != FALSE & flow == FALSE){
                    res <- data.frame(compound = as.vector(compound),
                        tr = as.vector(xmax), Hmax = as.vector(H),
                        A60 = as.vector(A60), B60 = as.vector(B60),
                        A10 = as.vector(A10), B10 = as.vector(B10),
                        area = as.vector(area), RSE = as.vector(summary),
                        MeanError = ME, correlation = cor)
                    } else if (compound != FALSE & area == FALSE &
                            flow == FALSE){
                        res <- data.frame(compound = as.vector(compound),
                            tr = as.vector(xmax), Hmax = as.vector(H),
                            A60 = as.vector(A60), B60 = as.vector(B60),
                            A10 = as.vector(A10), B10 = as.vector(B10),
                            RSE = as.vector(summary), MeanError = ME,
                            correlation = cor)
                        } else if (compound == FALSE & area == FALSE &
                                flow == FALSE){
                            res <- data.frame(tr = as.vector(xmax),
                                Hmax = as.vector(H), A60 = as.vector(A60),
                                B60 = as.vector(B60), A10 = as.vector(A10),
                                B10 = as.vector(B10), RSE = as.vector(summary),
                                MeanError = ME, correlation = cor)
                            }
    return(res)
}
