vanDeemterAlternative <- function(col, ext, dead, length, approachI = FALSE, A,
    B, C, approachII = FALSE){
  if (!is.numeric(A) | !is.numeric(B) | !is.numeric(C)) {
    stop("A, B and C parameters must be numeric")
  }
  if (approachI == TRUE & approachII == TRUE) {
    stop("Just one approach at the time is possible")
  }
  if (approachI == FALSE & approachII == FALSE) {
    stop("approachI or approachII parameter must be TRUE")
  }
  if (!is.numeric(length)) {
    stop("length parameter must be numeric")
  }
  if (!is.data.frame(col)) {
    stop("col parameter must be a data frame")
  }
  if (!is.data.frame(ext)) {
    stop("ext parameter must be a data frame")
  }
  if (!is.data.frame(dead)) {
    stop("dead parameter must be a data frame")
  }
  if (approachI == TRUE & approachII == FALSE) {
    col <- col[order(col$flow), ]
    ext <- ext[order(ext$flow), ]
    dead <- dead[order(dead$flow), ]
    t <- table(c(ext$flow, col$flow, dead$flow))
    del <- as.numeric(names(which(t < max(as.vector(t)))))
    ext <- ext[!ext$flow %in% del, ]
    col <- col[!col$flow %in% del, ]
    dead <- dead[!dead$flow %in% del, ]
    flows <- unique(col$flow)

    if (length(flows) == 0) {
      stop("There's no common flows at the given data frames col, ext and dead")
    }

    restep1 <- data.frame()
    summary <- data.frame()
    for (i in 1:length(flows)) {
      ef02 <- subset(ext, ext$flow == flows[i])
      f02 <- subset(col, col$flow == flows[i])
      Flowf02 <- f02[, "flow"]
      trf02 <- f02[, "tr"]
      etrf02 <- ef02[, "tr"]
      trf02 <- trf02 + (2 * (f02$B60 - f02$A60))/(sqrt(2 * pi))
      etrf02 <- etrf02 + (2 * (ef02$B60 - ef02$A60))/(sqrt(2 * pi))
      Vr02 <- trf02 * Flowf02
      eVr02 <- etrf02 * Flowf02
      varf02 <- (((f02$A60^3) + (f02$B60^3))/(f02$A60 + f02$B60)) -
          (((2 * (f02$B60 - f02$A60)^2))/pi)
      y <- varf02 * Flowf02 * Flowf02
      x <- (Vr02 - eVr02)^2
      data <- as.data.frame(cbind(x, y))

      d <- ggplot(data, aes(x = data[, 1], y = data[, 2])) +
        geom_point(color = "black", size = 2, show.legend = TRUE) +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE)
      d2 <- d + theme_bw()
      d22 <- d2 + theme(plot.title = element_text(face = "bold"),
                        axis.title = element_text(size = 15),
          axis.title.y = element_text(margin = margin(0, 20, 0, 0)),
          legend.text = element_text(size = 15),
                        legend.title = element_text(size = 15)) +
          scale_y_continuous(expression(sigma^{2} * ~F^2 ~ (mL^{2})))
      d22 <- d22 + ggtitle(paste("Flow: ", flows[i]))
      d22 <- d22 + theme_minimal()
      d22 <- d22 + scale_x_continuous((V[R] - V[ext])^{2} ~ (mL)^{2})
      print(d22)
      x <- data[, 1]
      y <- data[, 2]
      fit <- lm(y ~ x)
      r2 <- summary(fit)$r.squared
      H <- coefficients(fit)[2] * length * 1000
      restep1d <- cbind(Flow = flows[i],r2=r2)
      restep1dd <- as.data.frame(restep1d)
      restep1 <- rbind(restep1, restep1dd)

      res <- cbind(Flow = flows[i], H = H)
      ress <- as.data.frame((res))
      summary <- rbind(summary, ress)
      rm(H, trf02, varf02, Vr02, x, y, fit, d2, d22, etrf02, eVr02, d, Flowf02)
    }

    appI <- summary
    t0 <- dead[, "tr"]
    velocity2 <- length/t0
    appI <- as.data.frame(appI)
    data <- as.data.frame(cbind(velocity2, appI[, 2]))
    d <- ggplot(data, aes(x = data[, 1], y = data[, 2])) +
      geom_point(color = "black", size = 2, show.legend = TRUE) +
      geom_smooth(method = "nls", formula = y ~ A + (B/x) + C * x,
          method.args = list(start = c(A = A, B = B, C = C)), se = FALSE)
    d2 <- d + theme_bw()
    d22 <- d2 + theme(plot.title = element_text(face = "bold"),
                      axis.title = element_text(size = 15),
        axis.title.y = element_text(margin = margin(0, 20, 0, 0)),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
        scale_y_continuous(expression(H ~ (mu ~ m)))
    d22 <- d22 + theme_minimal()
    d22 <- d22 + scale_x_continuous(expression(u(mm/min)))
    print(d22)
    x <- data[, 1]
    y <- data[, 2]
    fit <- nls(y ~ A + (B/x) + C * x, start = list(A = A, B = B, C = C))

    cor <-cor(y,predict(fit))
    ME <- mean(abs(y-predict(fit))/abs(y))*100
    return(list(table = summary, coefficients = coefficients(fit), step1=restep1,
        correlation=cor, MeanError=ME, RSE=summary(fit)$sigma))
  }

  if (approachI == FALSE & approachII == TRUE) {
    col <- col[order(col$flow), ]
    ext <- ext[order(ext$flow), ]
    dead <- dead[order(dead$flow), ]
    t <- table(c(ext$flow, col$flow, dead$flow))
    del <- as.numeric(names(which(t < max(as.vector(t)))))
    ext <- ext[!ext$flow %in% del, ]
    col <- col[!col$flow %in% del, ]
    dead <- dead[!dead$flow %in% del, ]
    compounds <- unique(col$compound)
    summary <- data.frame()
    for (i in 1:length(compounds)) {
      ecompound <- ext[ext$compound == ext$compound[1],
                       ]
      ccompound <- col[col$compound == compounds[i], ]
      Flow <- ccompound[, "flow"]
      tr <- ccompound[, "tr"]
      etr <- ecompound[, "tr"]
      trr <- tr + (2 * (ccompound$B60 - ccompound$A60))/(sqrt(2 * pi))
      etrr <- etr + (2 * (ecompound$B60 - ecompound$A60))/(sqrt(2 * pi))
      t0 <- dead[, "tr"]
      velocity2 <- length/t0
      varc <- (((ccompound$A60^3) + (ccompound$B60^3))/(ccompound$A60 +
            ccompound$B60)) - (((2 * (ccompound$B60 - ccompound$A60)^2))/pi)
      data <- as.data.frame(cbind(velocity2, varc * length * Flow * Flow * 100))
      d <- ggplot(data, aes(x = data[, 1], y = data[, 2])) +
        geom_point(color = "black", size = 2, show.legend = TRUE) +
        geom_smooth(method = "nls", formula = y ~ A + (B/x) + C * x,
            method.args = list(start = c(A = A, B = B, C = C)), se = FALSE)
      d2 <- d + theme_bw()
      d22 <- d2 + theme(plot.title = element_text(face = "bold"),
          axis.title = element_text(size = 15),
          axis.title.y = element_text(margin = margin(0, 20, 0, 0)),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)) +
          scale_y_continuous(expression(sigma^{2} * ~F^2 ~ L))
      d22 <- d22 + ggtitle(paste("Compound: ", compounds[i]))
      d22 <- d22 + theme_minimal()
      d22 <- d22 + scale_x_continuous(expression(u(mm/min)))
      print(d22)
      x <- data[, 1]
      y <- data[, 2]
      fit <- nls(y ~ A + (B/x) + C * x, start = list(A = A, B = B, C = C))
      cor <-cor(y,predict(fit))
      ME <- mean(abs(y-predict(fit))/abs(y))*100
      RSE <- summary(fit)$sigma
      Av <- coefficients(fit)[1]
      Bv <- coefficients(fit)[2]
      Cv <- coefficients(fit)[3]
      k <- (trr - etrr)/(t0 - etrr)
      res <- cbind(compound = compounds[i], meanK = mean(k), Av, Bv, Cv, cor,
          ME, RSE)
      ress <- as.data.frame((res))
      summary <- rbind(summary, ress)
      rm(Flow, tr, etr, trr, etrr, t0, velocity2, varc, d, d2, d22, x, y, data,
          fit, Av, Bv, Cv, k)
    }
    app2 <- summary
    t0 <- dead[, "tr"]
    Flow0 <- dead[, "flow"]
    data <- as.data.frame(cbind((1 + app2[, 2])^2, app2[, 3:5]))
    d <- ggplot(data) + geom_point(aes(data[, 1], data[,2]), size = 2,
        show.legend = TRUE, color = "blue") +
      geom_smooth(aes(data[, 1], data[, 2]), method = lm, se = FALSE) +
        geom_point(aes(data[, 1], data[, 3]), size = 2, color = "green",
            show.legend = TRUE) + geom_smooth(aes(data[, 1], data[, 3]),
                method = lm, se = FALSE) + geom_point(aes(data[, 1],
                    data[, 4]), size = 2, color = "red", show.legend = TRUE) +
      geom_smooth(aes(data[, 1], data[, 4]), method = lm, se = FALSE)
    d2 <- d + theme_bw()
    d22 <- d2 + theme(plot.title = element_text(face = "bold"),
         axis.title = element_text(size = 15),
        axis.title.y = element_text(margin = margin(0, 20, 0, 0)),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
        scale_y_continuous(expression(Av ~ Bv ~ Cv))
    d22 <- d22 + theme_minimal()
    d22 <- d22 + scale_x_continuous(expression((1 + k)^{2}))
    print(d22)
    x <- data[, 1]
    y1 <- data[, 2]
    y2 <- data[, 3]
    y3 <- data[, 4]
    fit1 <- lm(y1 ~ x)
    fit2 <- lm(y2 ~ x)
    fit3 <- lm(y3 ~ x)
    V0 <- mean(t0 * t0 * Flow0 * Flow0)
    A <- as.vector((coefficients(fit1)[2]/V0) * 10)
    B <- as.vector((coefficients(fit2)[2]/V0) * 10)
    C <- as.vector((coefficients(fit3)[2]/V0) * 10)
    r2.A <- summary(fit1)$r.squared
    r2.B <- summary(fit2)$r.squared
    r2.C <- summary(fit3)$r.squared
    MRE.A <- mean(abs(y1-predict(fit1))/abs(y1))*100
    MRE.B <- mean(abs(y2-predict(fit2))/abs(y2))*100
    MRE.C <- mean(abs(y3-predict(fit3))/abs(y3))*100

    return(list(table = summary, coefficients = c(A = A, B = B, C = C),
        r2A = r2.A, r2B = r2.B, r2C = r2.C, MREA = MRE.A, MREB = MRE.B,
        MREC = MRE.C))
  }
}
