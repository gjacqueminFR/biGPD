

#' Decluster univariate data.
#'
#' @param data A vector of data.
#' @param nbDaysPerYear The number of days considered per year (integer).
#' @param nbYears The number of distinct years.
#' @param threshold The value above which data are supposed to follow a GPD.
#' @param blockSize The size of the blocks used for declustering.
#' @return A list of declustered data: only the maxima of blocks that also exceeds the threshold.
#' @export

UnivariateDeclustering <- function(data, nbDaysPerYear, nbYears, threshold, blockSize) {

  nbBlocksPerYear <- (nbDaysPerYear %/% blockSize) + min(1, nbDaysPerYear %% blockSize)
  dataDeclustered <- c()

  for (year in seq_len(nbYears)) {
    startYear <- 1 + (year - 1) * nbDaysPerYear
    endYear <- year * nbDaysPerYear

    dataYear <- data[startYear: endYear]

    for (block in seq_len(nbBlocksPerYear)) {
      startBlock <- 1 + (block - 1) * blockSize
      endBlock <- min(block * blockSize, nbDaysPerYear)

      maxBlock <- max(dataYear[startBlock: endBlock])

      if (maxBlock >= threshold) {
        dataDeclustered <- append(dataDeclustered, maxBlock)
      }
    }
  }
  return(dataDeclustered)
}


#' Decluster bivariate data.
#'
#' @param data1 A vector of data for the first variable.
#' @param data2 A vector of data for the second variable.
#' @param nbDaysPerYear The number of days considered per year (integer).
#' @param nbYears The number of distinct years.
#' @param thresholds A vector of size 2. The first threshold is for the first variable, and the second value is for the second variable. The threshold is the value above which data are supposed to follow a GPD.
#' @param blockSize The size of the blocks used for declustering.
#' @param logic Either "AND" or "OR". Couples of bloxk maxima are retained if they both exceed their respective threshold in case of "AND", and if at least one of them exceed its respective threshold in case of "OR". Default value to "AND".
#' @return A dataframe with two columns, one for each variable, of declustered data following the logic parameter
#' @export

BivariateDeclustering <- function(data1, data2, nbDaysPerYear, nbYears, thresholds, blockSize, logic) {
  if (missing(logic)) {
    logic <- "AND"
  }

  nbBlocksPerYear <- (nbDaysPerYear %/% blockSize) + min(1, nbDaysPerYear %% blockSize)
  dataDeclustered1 <- c()
  dataDeclustered2 <- c()

  for (year in seq_len(nbYears)) {
    startYear <- 1 + (year - 1) * nbDaysPerYear
    endYear <- year * nbDaysPerYear

    dataYear1 <- data1[startYear: endYear]
    dataYear2 <- data2[startYear: endYear]

    for (block in seq_len(nbBlocksPerYear)) {
      startBlock <- 1 + (block - 1) * blockSize
      endBlock <- min(block * blockSize, nbDaysPerYear)

      maxBlock1 <- max(dataYear1[startBlock: endBlock])
      maxBlock2 <- max(dataYear2[startBlock: endBlock])

      if (logic == "AND") {
        if (maxBlock1 >= thresholds[1] && maxBlock2 >= thresholds[2]) {
          dataDeclustered1 <- append(dataDeclustered1, maxBlock1)
          dataDeclustered2 <- append(dataDeclustered2, maxBlock2)
        }
      } else if (logic == "OR") {
        if (maxBlock1 >= thresholds[1] || maxBlock2 >= thresholds[2]) {
          dataDeclustered1 <- append(dataDeclustered1, maxBlock1)
          dataDeclustered2 <- append(dataDeclustered2, maxBlock2)
        }
      }
    }
  }
  dataBivDecluster <- data.frame("Var1" = dataDeclustered1, "Var2" = dataDeclustered2)
  return(dataBivDecluster)
}


#' Estimate the parameters of the GPD
#'
#' @param data A vector of data.
#' @param threshold The value above which data are supposed to follow a GPD.
#' @importFrom tea gpdFit
#' @return The estimated parameters of the GPD, in a list with the following order: c(threshold, sigma, xi).
#' @export

EstimateGPDParameters <- function(data, threshold) {
  if (is.nan(tea::gpdFit(data, threshold = threshold, method = "mle")$par.ses[["Scale (Intercept)"]])) {
    fit <- tea::gpdFit(data, threshold = threshold, method = "pwm")
  } else {
    fit <- tea::gpdFit(data, threshold = threshold, method = "mle")
  }

  GPDparam <- c(threshold, fit$par.ests[["Scale (Intercept)"]], fit$par.ests[["Shape (Intercept)"]]) # threshold, sigma, xi

  return(GPDparam)
}


#' Select the copula family
#'
#' @param data A dataframe with two columns, one for each variable.
#' @param probaQuantile Data above the threshold of this probability are considered extremes in a peaks-over-threshold approach.
#' @param nbDaysPerYear The number of days considered per year (integer).
#' @param nbYears The number of years considered (integer).
#' @param blockSizes Vector of integers of size 3. The sizes of the blocks considered for the declustering, with the univariates first and the bivariate block size at the end.
#' @param listProbas List of probabilities of quantiles. Have to be high probabilities. Default to c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98).
#' @param listCopulaNumber The list of copula from which the copula is selected. Default to c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 16, 17, 18, 19, 20, 23, 24, 26, 27, 28, 29, 30, 33, 34, 36, 37, 38, 39, 40). See VineCopula documentation.
#' @importFrom VineCopula as.copuladata
#' @importFrom VineCopula BiCopEstList
#' @importFrom VineCopula BiCopName
#' @importFrom tea pgpd
#' @return The copula selected. It's the short name from the VineCopula package.
#' @export

CopulaSelection <- function(data, probaQuantile, nbDaysPerYear, nbYears, blockSizes, listProbas, listCopulaNumber) {
  start("Copula selection")
  if (missing(listProbas)) {
    listProbas <- c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98)
  }
  if (missing(listCopulaNumber)) {
    listCopulaNumber <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 16, 17, 18, 19, 20, 23, 24, 26, 27, 28, 29, 30, 33, 34, 36, 37, 38, 39, 40)
  }

  copulaScore <- data.frame(Families = listCopulaNumber,
                            ScoreSumPlaces = integer(length(listCopulaNumber)),
                            ScoreMeanBIC = numeric(length(listCopulaNumber)),
                            PercentAppearances = numeric(length(listCopulaNumber)),
                            AppearAtSelectedQuantile = logical(length(listCopulaNumber)))
  copulaScore$AppearAtSelectedQuantile <- FALSE

  data1 <- data[, 1]
  data2 <- data[, 2]

  for (q in listProbas) {

    ### first variable
    threshold1 <- quantile(data1, q)
    dataDecluster1 <- UnivariateDeclustering(data1, nbDaysPerYear, nbYears, threshold1, blockSizes[1])
    GPDparam1 <- EstimateGPDParameters(dataDecluster1, threshold1)

    ### second variable
    threshold2 <- quantile(data2, q)
    dataDecluster2 <- UnivariateDeclustering(data2, nbDaysPerYear, nbYears, threshold2, blockSizes[2])
    GPDparam2 <- EstimateGPDParameters(dataDecluster2, threshold2)

    ### Bivariate
    dataBiv <- BivariateDeclustering(data1, data2, nbDaysPerYear, nbYears, c(threshold1, threshold2), blockSizes[3])

    if (length(dataBiv$Var1) > 2) {

      if (GPDparam1[3] < 0) {
        dataBiv[, 1][dataBiv[, 1] >= GPDparam1[1] - (GPDparam1[2] / GPDparam1[3])] <- GPDparam1[1] - (GPDparam1[2] / GPDparam1[3]) - 0.0001
      }
      if (GPDparam2[3] < 0) {
        dataBiv[, 2][dataBiv[, 2] >= GPDparam2[1] - (GPDparam2[2] / GPDparam2[3])] <- GPDparam2[1] - (GPDparam2[2] / GPDparam2[3]) - 0.0001
      }

      dataBiv[, 1] <- tea::pgpd(dataBiv[, 1], GPDparam1[1], GPDparam1[2], GPDparam1[3])
      dataBiv[, 2] <- tea::pgpd(dataBiv[, 2], GPDparam2[1], GPDparam2[2], GPDparam2[3])

      dataCop <- VineCopula::as.copuladata(dataBiv)

      Unif1 <- dataCop[, 1]
      Unif2 <- dataCop[, 2]

      # Copula selection

      if (length(unique(Unif1)) > 2 && length(unique(Unif2)) > 2) {

        listCopula <- VineCopula::BiCopEstList(Unif1, Unif2, familyset = listCopulaNumber)$summary

        orderedCopula <- listCopula[order(listCopula$BIC), ]

        for (i in seq_along(orderedCopula$family)) {

          copulaScore[copulaScore$Families == orderedCopula[i, ]$family, ]$ScoreSumPlaces <- i + copulaScore[copulaScore$Families == orderedCopula[i, ]$family, ]$ScoreSumPlaces

          copulaScore[copulaScore$Families == orderedCopula[i, ]$family, ]$ScoreMeanBIC <- orderedCopula[i, ]$BIC + copulaScore[copulaScore$Families == orderedCopula[i, ]$family, ]$ScoreMeanBIC

          if (as.integer(100 * q) == as.integer(100 * probaQuantile)) {
            copulaScore[copulaScore$Families == orderedCopula[i, ]$family, ]$AppearAtSelectedQuantile <- TRUE
          }

          copulaScore[copulaScore$Families == orderedCopula[i, ]$family, ]$PercentAppearances <- copulaScore[copulaScore$Families == orderedCopula[i, ]$family, ]$PercentAppearances + 1 / length(listProbas)

        }

      }
    }
  }

  if (length(copulaScore$AppearAtSelectedQuantile[copulaScore$AppearAtSelectedQuantile]) > 0) {
    copulaScore <- copulaScore[copulaScore$AppearAtSelectedQuantile, ]
  }
  copulaScore <- copulaScore[copulaScore$PercentAppearances > 0.5, ]
  copulaScore$ScoreMeanBIC <- copulaScore$ScoreMeanBIC / (length(listProbas) * copulaScore$PercentAppearances) # Transform sum of BIC into a mean BIC pondered by the percentage of appearances

  orderedcopulaScore <- copulaScore[order(copulaScore$ScoreMeanBIC), ]
  CandidatesCopula <- orderedcopulaScore[orderedcopulaScore$ScoreMeanBIC < orderedcopulaScore$ScoreMeanBIC[1] + 2 & orderedcopulaScore$ScoreMeanBIC <= 0, ]

  if (length(orderedcopulaScore$Families) < 1) {
    selectedCopula <- "I"
  } else if (CandidatesCopula$ScoreMeanBIC[1] > -0.02) {
    selectedCopula <- "I"
  } else if ("Gumbel" %in% VineCopula::BiCopName(CandidatesCopula$Families, short = FALSE)) {
    if (CandidatesCopula[CandidatesCopula$Families == 4, ]$ScoreMeanBIC < 0) {
      selectedCopula <- "G"
    } else {
      selectedCopula <- VineCopula::BiCopName(CandidatesCopula$Families[1], short = TRUE)
    }
  } else {
    selectedCopula <- VineCopula::BiCopName(CandidatesCopula$Families[1], short = TRUE)
  }

  return(selectedCopula)
}


#' Compute the probability to be under the return level value with a GPD.
#'
#' @param returnLevel The value at which the cdf is computed.
#' @param GPDParam The parameters of the GPD, in a list with the following order: c(threshold, sigma, xi).
#' @importFrom tea pgpd
#' @return The cdf value at the return level value for a GPD.
#' @export

GetGPDproba <- function(returnLevel, GPDParam) {
  if (GPDParam[3] < 0 && returnLevel > GPDParam[1] - GPDParam[2] / GPDParam[3]) {
    GPDproba <- 1
  } else {
    GPDproba <- tea::pgpd(returnLevel, loc = GPDParam[1], scale = GPDParam[2], shape = GPDParam[3])
  }
  return(GPDproba)
}


#' Compute the value corresponding to the cdf probability with a GPD.
#'
#' @param GPDproba The probability at which the inverse cdf is computed.
#' @param GPDParam The parameters of the GPD, in a list with the following order: c(threshold, sigma, xi).
#' @importFrom tea qgpd
#' @return The value at which the cdf equals proba for a GPD.
#' @export

GetGPDvalue <- function(GPDproba, GPDParam) {

  GPDvalue <- tea::qgpd(GPDproba, loc = GPDParam[1], scale = GPDParam[2], shape = GPDParam[3])

  return(GPDvalue)
}


#' Calculate the univariate return period with the GPD approach
#'
#' @param proba The GPD probability corresponding to the return level.
#' @param extremalIndex The extremal index cf UnivariateExtremalIndex.
#' @param nbDaysPerYear The number of days considered per year (integer).
#' @param probaQuantile The probabilityof the value above which data are supposed to follow a GPD. A common value is 0.95.
#' @param h The parameter of non-concurrence (integer).
#' @param probaOccurrence The probability that the values are reached before the return period time. Default value to 0.63.
#' @importFrom tea gpdFit
#' @return The univariate return period with the GPD approach.
#' @export

UnivariateReturnPeriodCopula <- function(proba, extremalIndex, nbDaysPerYear, probaQuantile, h, probaOccurrence) {
  if (missing(probaOccurrence)) {
    probaOccurrence <- 1 - exp(-1)
  }

  ReturnPeriod <- - log(1 - probaOccurrence) * h / (nbDaysPerYear * (1 - proba) * (1 - probaQuantile^(extremalIndex * h)))

  return(ReturnPeriod)
}


#' Calculate the non-concurrent joint excess probability.
#'
#' @param FU1U2 The probability to be below both thresholds (quantiles of probability probaQuantile).
#' @param h The parameter of non-concurrence (integer).
#' @param extremalIndexes A list of size 3, with the extremal index for the first variable in first, the one for the second variable in second and the bivaraite extremal index in third.
#' @param probaQuantile The probability of the value above which data are supposed to follow a GPD. A common value is 0.95.
#' @param nbYears The number of distinct years. Default value is 1.
#' @param Dparam Cf documentation of the dgaps function of the exdex package. Default value is 3.
#' @return The non-concurrent joint excess probability.
#' @export

Hbar <- function(FU1U2, h, extremalIndexes, probaQuantile, nbYears, Dparam) {

  minExtremalIndexBiv <- max((1 - probaQuantile) * extremalIndexes[1] / (1 - FU1U2), (1 - probaQuantile) * extremalIndexes[2] / (1 - FU1U2))
  maxExtremalIndexBiv <- (1 - probaQuantile) * extremalIndexes[1] / (1 - FU1U2) + (1 - probaQuantile) * extremalIndexes[2] / (1 - FU1U2)
  if (extremalIndexes[3] < minExtremalIndexBiv) {
    extremalIndexBiv <- minExtremalIndexBiv
  } else if (extremalIndexes[3] > maxExtremalIndexBiv) {
    extremalIndexBiv <- maxExtremalIndexBiv
  } else {
    extremalIndexBiv <- extremalIndexes[3]
  }

  res <- 1 - probaQuantile^(extremalIndexes[1] * h) - probaQuantile^(extremalIndexes[2] * h) + FU1U2^(extremalIndexBiv * h)

  return(res)
}


#' Calculate the bivariate return period with the copula approach.
#'
#' @param GPDprobas The cdf value at the return level value for a GPD. A vector of size 2, one for each variable.
#' @param h The parameter of non-concurrence (integer).
#' @param extremalIndexes A list of size 3, with the extremal index for the first variable in first, the one for the second variable in second and the bivaraite extremal index in third.
#' @param probaQuantile The probabilityof the value above which data are supposed to follow a GPD. A common value is 0.95.
#' @param copula The copula family. The short name from the VineCopula package is expected.
#' @param FU1U2 The probability to be below both thresholds (quantiles of probability probaQuantile).
#' @param nbDaysPerYear The number of days considered per year (integer).
#' @param nbYears The number of distinct years. Default value is 1.
#' @param Dparam Cf documentation of the dgaps function of the exdex package. Default value is 3.
#' @param probaOccurrence The probability that the values are reached before the return period time. Default value to 0.63.
#' @return The bivariate return period with the copula approach.
#' @export

BivariateReturnPeriodCopula <- function(GPDprobas, h, extremalIndexes, probaQuantile, copula, FU1U2, nbDaysPerYear, nbYears, Dparam, probaOccurrence) {
  if (missing(probaOccurrence)) {
    probaOccurrence <- 1 - exp(-1)
  }
  if (missing(nbYears)) {
    nbYears <- 1
  }
  if (missing(Dparam)) {
    Dparam <- 3
  }

  hbarU1U2 <- Hbar(FU1U2, h, extremalIndexes, probaQuantile, nbYears, Dparam)

  if (is.nan(VineCopula::BiCopCDF(GPDprobas[1], GPDprobas[2], copula))) {
    copValueBiv <- 1
  } else {
    copValueBiv <- VineCopula::BiCopCDF(GPDprobas[1], GPDprobas[2], copula)
  }

  if (is.nan(VineCopula::BiCopCDF(GPDprobas[1], 1, copula))) {
    copValueOne <- 1
  } else {
    copValueOne <- VineCopula::BiCopCDF(GPDprobas[1], 1, copula)
  }

  if (is.nan(VineCopula::BiCopCDF(1, GPDprobas[2], copula))) {
    copValueTwo <- 1
  } else {
    copValueTwo <- VineCopula::BiCopCDF(1, GPDprobas[2], copula)
  }

  probaBiv <- 1 + copValueBiv - copValueOne - copValueTwo
  returnPeriodBivariate <- - log(1 - probaOccurrence) * h / (nbDaysPerYear * hbarU1U2 * probaBiv)

  return(c(returnPeriodBivariate, probaBiv))
}


#' Example of copula approach. Compute univariate and bivariate return periods from 2 time series.
#'
#' @param data A dataframe with two columns, one for each time series. 
#' @param returnLevels A vector of size 2, one for each time series. Return periods correspond to the return levels.
#' @param probaQuantile The probabilityof the value above which data are supposed to follow a GPD. A common value is 0.95.
#' @param nbDaysPerYear The number of days considered per year (integer).
#' @param nbYears The number of distinct years. Default value is 1.
#' @param h The parameter of non-concurrence (integer).
#' @param blockSizes Vector of integers of size 3. The sizes of the blocks considered for the declustering, with the univariates first and the bivariate block size at the end.
#' @param Dparam Cf documentation of the dgaps function of the exdex package. Default value is 3.
#' @param probaOccurrence The probability that the values are reached before the return period time. Default value to 0.63.
#' @param logic Either "AND" or "OR". Couples of bloxk maxima are retained if they both exceed their respective threshold in case of "AND", and if at least one of them exceed its respective threshold in case of "OR". Default value to "AND".
#' @return A list with, in that order: the return period of first variable, the return period of the second variable, the bivariate return period, the non-concurrent joint excess probability, chi and chiBarre.
#' @export

CopulaApproach <- function(data, returnLevels, probaQuantile, nbDaysPerYear, nbYears, h, blockSizes, Dparam, probaOccurrence, logic) {

  data1 <- data[, 1]
  threshold1 <- quantile(data1, probaQuantile)
  data2 <- data[, 2]
  threshold2 <- quantile(data2, probaQuantile)

  extremalIndex1 <- UnivariateExtremalIndex(data1, probaQuantile, nbYears, Dparam)
  extremalIndex2 <- UnivariateExtremalIndex(data2, probaQuantile, nbYears, Dparam)

  if (missing(blockSizes)) {
    blockSize1 <- ceiling(1 / extremalIndex1)
    blockSize2 <- ceiling(1 / extremalIndex2)
  } else {
    if (blockSizes[1] == 0) {
      blockSize1 <- ceiling(1 / extremalIndex1)
    } else {
      blockSize1 <- blockSizes[1]
    }
    if (blockSizes[2] == 0) {
      blockSize2 <- ceiling(1 / extremalIndex2)
    } else {
      blockSize2 <- blockSizes[2]
    }
  }

  dataDecluster1 <- UnivariateDeclustering(data1, nbDaysPerYear, nbYears, threshold1, blockSize1)
  dataDecluster2 <- UnivariateDeclustering(data2, nbDaysPerYear, nbYears, threshold2, blockSize2)

  # Estimate univariate parameters
  GPDparam1 <- EstimateGPDParameters(dataDecluster1, threshold1)
  GPDparam2 <- EstimateGPDParameters(dataDecluster2, threshold2)

  # Calculate univariate return periods
  GPDproba1 <- GetGPDproba(returnLevels[1], GPDparam1)
  GPDproba2 <- GetGPDproba(returnLevels[2], GPDparam2)
  returnPeriod1 <- UnivariateReturnPeriodCopula(GPDproba1, extremalIndex1, nbDaysPerYear, probaQuantile, h, probaOccurrence)
  returnPeriod2 <- UnivariateReturnPeriodCopula(GPDproba2, extremalIndex2, nbDaysPerYear, probaQuantile, h, probaOccurrence)

  # Bivariate

  # Bivariate extremal index
  # Transform data to Frechet margins
  ECDF1 <- ecdf(data1)
  ECDF2 <- ecdf(data2)
  frechet1 <- -1 / log(ECDF1(data1))
  frechet2 <- -1 / log(ECDF2(data2))
  extremalIndexBiv <- BivariateExtremalIndex(frechet1, frechet2, probaQuantile, c(probaQuantile, probaQuantile), nbYears, Dparam)

  # Bivariate declustering
  if (missing(blockSizes)) {
    blockSizeBiv <- ceiling(1 / extremalIndexBiv)
  } else {
    if (blockSizes[3] == 0) {
      blockSizeBiv <- ceiling(1 / extremalIndexBiv)
    } else {
      blockSizeBiv <- blockSizes[3]
    }
  }

  dataBiv <- BivariateDeclustering(data1, data2, nbDaysPerYear, nbYears, c(threshold1, threshold2), blockSizeBiv, logic)

  if (length(unique(dataBiv$Var1)) <= 3 && length(unique(dataBiv$Var2)) <= 3) {
    copula <- VineCopula::BiCop(0, 0) # Independent copula
  } else {

    if (GPDparam1[3] < 0) {
      dataBiv[, 1][dataBiv[, 1] >= GPDparam1[1] - (GPDparam1[2] / GPDparam1[3])] <- GPDparam1[1] - (GPDparam1[2] / GPDparam1[3]) - 0.0001
    }
    if (GPDparam2[3] < 0) {
      dataBiv[, 2][dataBiv[, 2] >= GPDparam2[1] - (GPDparam2[2] / GPDparam2[3])] <- GPDparam2[1] - (GPDparam2[2] / GPDparam2[3]) - 0.0001
    }

    dataBiv[, 1] <- tea::pgpd(dataBiv[, 1], GPDparam1[1], GPDparam1[2], GPDparam1[3])
    dataBiv[, 2] <- tea::pgpd(dataBiv[, 2], GPDparam2[1], GPDparam2[2], GPDparam2[3])

    dataCop <- VineCopula::as.copuladata(dataBiv)

    Unif1 <- dataCop[, 1]
    Unif2 <- dataCop[, 2]

    # Copula estimation
    copulaFamily <- CopulaSelection(data, probaQuantile, nbDaysPerYear, nbYears, c(blockSize1, blockSize2, blockSizeBiv))
    copula <- VineCopula::BiCopEst(Unif1, Unif2, family = VineCopula::BiCopName(copulaFamily), se = TRUE)
  }

  # Calculate bivariate return period
  dataBelow <- data[data[, 1] <= threshold1 & data[, 2] <= threshold2, ]
  FU1U2 <- length(dataBelow[, 1]) / length(data[, 1])

  out <- BivariateReturnPeriodCopula(c(GPDproba1, GPDproba2), h, c(extremalIndex1, extremalIndex2, extremalIndexBiv), probaQuantile, copula, FU1U2, nbDaysPerYear, nbYears, Dparam, probaOccurrence)
  returnPeriodBiv <- out[1]
  probaBiv <- out[2]

  # Transform negative or NA return periods to Inf
  if (is.na(returnPeriodBiv)) {
    returnPeriodBiv <- Inf
    probaBiv <- 0
  } else if (returnPeriodBiv < 0) {
    returnPeriodBiv <- Inf
    probaBiv <- 0
  }

  # Probability of non-concurrent excess
  # probaBiv <- h / (nbDaysPerYear * returnPeriodBiv)

  chi <- max(2 - log(VineCopula::BiCopCDF(0.9999, 0.9999, copula)) / log(0.9999), 0)
  chiBarre <- 2 * log(1 - 0.9999) / log(1 - 2 * 0.9999 + VineCopula::BiCopCDF(0.9999, 0.9999, copula)) - 1

  return(c(returnPeriod1, returnPeriod2, returnPeriodBiv, probaBiv, chi, chiBarre, copulaFamily))
}
