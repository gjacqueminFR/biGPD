
#' Bivariate return period when searching for return levels.
#'
#' @param probaRL A float between 0 and 1. The return levels correspond to the quantiles of this probability.
#' @param probaQuantile A float between 0 and 1. Data above the thresholds of this probability are used to compute the bivariate exceedence probability. The value has a small impact over the probability. It should be close to 1, 0.95 is a classical value. The same value is used for both margins for simplicity.
#' @param FbarU1U2 The probability to be above both thresholds (quantiles of probability probaQuantile).
#' @param extremalIndex1 The univariate extremal indexe of the first variable.
#' @param extremalIndex2 The univariate extremal indexe of the second variable.
#' @param extremalIndexBivOG The bivariate extremal index.
#' @param h The parameter of non-concurrence (integer).
#' @param nbDaysPerYear The number of days considered per year (integer).
#' @param probaOccurrence The probability that the values are reached before the return period time. Default value to 0.63.
#' @return The bivariate return period.
#' @export

ReturnPeriodBiGPDReturnLevels <- function(probaRL, probaQuantile, FbarU1U2, extremalIndex1, extremalIndex2, extremalIndexBivOG, h, nbDaysPerYear, probaOccurrence) {
  FbarX1X2 <- FbarU1U2 * (1 - probaRL) / (1 - probaQuantile)

  # Constraint on the bivariate extremal index
  minExtremalIndexBiv <- max(extremalIndex1 * (1 - probaRL) / (2 - probaRL - probaRL - FbarX1X2), extremalIndex2 * (1 - probaRL) / (2 - probaRL - probaRL - FbarX1X2))
  maxExtremalIndexBiv <- extremalIndex1 * (1 - probaRL) / (2 - probaRL - probaRL - FbarX1X2) + extremalIndex2 * (1 - probaRL) / (2 - probaRL - probaRL - FbarX1X2)
  if (extremalIndexBivOG < minExtremalIndexBiv) {
    extremalIndexBiv <- minExtremalIndexBiv
  } else if (extremalIndexBivOG > maxExtremalIndexBiv) {
    extremalIndexBiv <- maxExtremalIndexBiv
  } else {
    extremalIndexBiv <- extremalIndexBivOG
  }

  returnPeriodBiv <- ReturnPeriodBiGPD(c(probaRL, probaRL, FbarX1X2), c(extremalIndex1, extremalIndex2, extremalIndexBiv), h, nbDaysPerYear, probaOccurrence)

  return(returnPeriodBiv)
}

#' Example of biGPD approach to get return levels from return period. Compute the two return levels with equal univariate probability from a given bivariate return period..
#'
#' @param data A dataframe with two columns, one for each time series.
#' @param returnPeriod An integer, in years.
#' @param EGPDtypes A vector of size 2, one for each time series. Values are 1 or 4. See Naveau et al. (2016) and the package mev for further information.
#' @param initParams1 A vector of parameters to initialize the estimation for the first time series. For type 1, it is (Kappa, Sigma, Xi), for type 4, it is (Prob, Kappa, Delta, Sigma, Xi).
#' @param initParams2 A vector of parameters to initialize the estimation for the second time series. For type 1, it is (Kappa, Sigma, Xi), for type 4, it is (Prob, Kappa, Delta, Sigma, Xi).
#' @param probaQuantile A float between 0 and 1. Data above the thresholds of this probability are used to compute the bivariate exceedence probability. The value has a small impact over the probability. It should be close to 1, 0.95 is a classical value. The same value is used for both margins for simplicity.
#' @param nbDaysPerYear The number of days considered per year (integer).
#' @param nbYears The number of distinct years. Default value is 1.
#' @param h The parameter of non-concurrence (integer).
#' @param Dparam Cf documentation of the dgaps function of the exdex package. Default value is 3.
#' @param probaOccurrence The probability that the values are reached before the return period time. Default value to 0.63.
#' @return A list with, in that order: the return level of first variable, the return level of the second variable, chi and chiBarre.
#' @export

BiGPDApproachReturnLevels <- function(data, returnPeriod, EGPDtypes, initParams1, initParams2, probaQuantile, nbDaysPerYear, nbYears, h, Dparam, probaOccurrence) {

  if (missing(probaOccurrence)) {
    probaOccurrence <- 1 - exp(-1)
  }

  threshold1 <- quantile(data[, 1], probaQuantile)
  threshold2 <- quantile(data[, 2], probaQuantile)

  extremalIndex1 <- UnivariateExtremalIndex(data[, 1], probaQuantile, nbYears, Dparam)
  extremalIndex2 <- UnivariateExtremalIndex(data[, 2], probaQuantile, nbYears, Dparam)

  # Get the probability of no rain, and filter the zeros
  probaZero1 <- length(data[, 1][data[, 1] == 0]) / length(data[, 1])
  dataFiltered1 <- data[, 1][data[, 1] > 0]
  probaZero2 <- length(data[, 2][data[, 2] == 0]) / length(data[, 2])
  dataFiltered2 <- data[, 2][data[, 2] > 0]

  EGPDparam1 <- EstimateEGPDParameters(dataFiltered1, EGPDtypes[1], initParams1)
  EGPDparam2 <- EstimateEGPDParameters(dataFiltered2, EGPDtypes[2], initParams2)

  # Transformation to exponential margins
  dataBiv <- data[data[, 1] > threshold1 | data[, 2] > threshold2, ]

  exp1 <- EGPDtoExp(dataBiv[, 1], EGPDtypes[1], EGPDparam1, probaZero1)
  exp2 <- EGPDtoExp(dataBiv[, 2], EGPDtypes[2], EGPDparam2, probaZero2)

  thresholdExp1 <- EGPDtoExp(threshold1, EGPDtypes[1], EGPDparam1, probaZero1)
  thresholdExp2 <- EGPDtoExp(threshold2, EGPDtypes[2], EGPDparam2, probaZero2)

  # ECDF
  delta <- (exp1 - thresholdExp1) - (exp2 - thresholdExp2)
  ECDFDelta <- ecdf(delta)

  # Estimate bivariate extremal index
  frechet1 <- EGPDtoFrechet(data[, 1], EGPDtypes[1], EGPDparam1, probaZero1)
  frechet2 <- EGPDtoFrechet(data[, 2], EGPDtypes[2], EGPDparam2, probaZero2)

  extremalIndexBivOG <- BivariateExtremalIndex(frechet1, frechet2, probaQuantile, c(0.95, 0.95), nbYears, Dparam) # The value of the univariate probability is not necessary here, the only thing that matters is that proba1 = proba2.

  # Calculate return levels
  dataAbove <- data[data[, 1] >= threshold1 & data[, 2] >= threshold2, ]
  FbarU1U2 <- length(dataAbove[, 1]) / length(data[, 1])

  iteration <- 0
  probaRLmin <- probaQuantile
  probaRLmax <- 1

  if (ReturnPeriodBiGPDReturnLevels(probaRLmin, probaQuantile, FbarU1U2, extremalIndex1, extremalIndex2, extremalIndexBivOG, h, nbDaysPerYear, probaOccurrence) > returnPeriod) {
    probaRL <- probaRLmin
    print("Problem !!")
  } else {
    probaRL <- (probaRLmin + probaRLmax) / 2
    res <- ReturnPeriodBiGPDReturnLevels(probaRL, probaQuantile, FbarU1U2, extremalIndex1, extremalIndex2, extremalIndexBivOG, h, nbDaysPerYear, probaOccurrence) - returnPeriod

    while (abs(res) > 1 && iteration < 1000) {
      if (res > 0) {
        probaRLmax <- probaRL
        probaRL <- (probaRLmin + probaRLmax) / 2
      } else {
        probaRLmin <- probaRL
        probaRL <- (probaRLmin + probaRLmax) / 2
      }
      res <- ReturnPeriodBiGPDReturnLevels(probaRL, probaQuantile, FbarU1U2, extremalIndex1, extremalIndex2, extremalIndexBivOG, h, nbDaysPerYear, probaOccurrence) - returnPeriod
      iteration <- iteration + 1
    }
  }

  returnLevel1 <- ValueEGPD(probaRL, EGPDtypes[1], EGPDparam1, probaZero1)
  returnLevel2 <- ValueEGPD(probaRL, EGPDtypes[2], EGPDparam2, probaZero2)

  # Chi and chi barre values
  chi <- FbarU1U2 / (1 - probaQuantile)
  chiBarre <- (log(1 - 0.9999) - log(chi)) / (log(1 - 0.9999) + log(chi))

  return(c(returnLevel1, returnLevel2, chi, chiBarre))
}


#' Example of copula approach to get return levels from return period. Compute the two return levels with equal univariate probability from a given bivariate return period..
#'
#' @param data A dataframe with two columns, one for each time series.
#' @param returnPeriod An integer, in years.
#' @param probaQuantile The probabilityof the value above which data are supposed to follow a GPD. A common value is 0.95.
#' @param nbDaysPerYear The number of days considered per year (integer).
#' @param nbYears The number of distinct years. Default value is 1.
#' @param h The parameter of non-concurrence (integer).
#' @param blockSizes Vector of integers of size 3. The sizes of the blocks considered for the declustering, with the univariates first and the bivariate block size at the end.
#' @param Dparam Cf documentation of the dgaps function of the exdex package. Default value is 3.
#' @param probaOccurrence The probability that the values are reached before the return period time. Default value to 0.63.
#' @param logic Either "AND" or "OR". Couples of bloxk maxima are retained if they both exceed their respective threshold in case of "AND", and if at least one of them exceed its respective threshold in case of "OR". Default value to "AND".
#' @return A list with, in that order: the return level of first variable, the return level of the second variable, chi and chiBarre.
#' @export

CopulaApproachReturnLevels <- function(data, returnPeriod, probaQuantile, nbDaysPerYear, nbYears, h, blockSizes, Dparam, probaOccurrence, logic) {

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

  # Calculate return levels
  dataBelow <- data[data[, 1] <= threshold1 & data[, 2] <= threshold2, ]
  FU1U2 <- length(dataBelow[, 1]) / length(data[, 1])

  iteration <- 0
  probaRLmin <- probaQuantile
  probaRLmax <- 1
  GPDprobaRL1 <- 1 - (1 - probaRLmin^(h * extremalIndex1)) / (1 - probaQuantile^(h * extremalIndex1))
  GPDprobaRL2 <- 1 - (1 - probaRLmin^(h * extremalIndex2)) / (1 - probaQuantile^(h * extremalIndex2))

  if (BivariateReturnPeriodCopula(c(GPDprobaRL1, GPDprobaRL2), h, c(extremalIndex1, extremalIndex2, extremalIndexBiv), probaQuantile, copula, FU1U2, nbDaysPerYear, nbYears, Dparam, probaOccurrence)[1] > returnPeriod) {
    probaRL <- probaRLmax
  } else {
    probaRL <- (probaRLmin + probaRLmax) / 2
    GPDprobaRL1 <- 1 - (1 - probaRL^(h * extremalIndex1)) / (1 - probaQuantile^(h * extremalIndex1))
    GPDprobaRL2 <- 1 - (1 - probaRL^(h * extremalIndex2)) / (1 - probaQuantile^(h * extremalIndex2))
    res <- BivariateReturnPeriodCopula(c(GPDprobaRL1, GPDprobaRL2), h, c(extremalIndex1, extremalIndex2, extremalIndexBiv), probaQuantile, copula, FU1U2, nbDaysPerYear, nbYears, Dparam, probaOccurrence)[1] - returnPeriod

    while (abs(res) > 1) {
      if (res > 0) {
        probaRLmax <- probaRL
        probaRL <- (probaRLmin + probaRLmax) / 2
      } else {
        probaRLmin <- probaRL
        probaRL <- (probaRLmin + probaRLmax) / 2
      }
      GPDprobaRL1 <- 1 - (1 - probaRL^(h * extremalIndex1)) / (1 - probaQuantile^(h * extremalIndex1))
      GPDprobaRL2 <- 1 - (1 - probaRL^(h * extremalIndex2)) / (1 - probaQuantile^(h * extremalIndex2))
      res <- BivariateReturnPeriodCopula(c(GPDprobaRL1, GPDprobaRL2), h, c(extremalIndex1, extremalIndex2, extremalIndexBiv), probaQuantile, copula, FU1U2, nbDaysPerYear, nbYears, Dparam, probaOccurrence)[1] - returnPeriod
    }
  }

  returnLevel1 <- GetGPDvalue(GPDprobaRL1, GPDparam1)
  returnLevel2 <- GetGPDvalue(GPDprobaRL2, GPDparam2)

  chi <- max(2 - log(VineCopula::BiCopCDF(0.9999, 0.9999, copula)) / log(0.9999), 0)
  chiBarre <- 2 * log(1 - 0.9999) / log(1 - 2 * 0.9999 + VineCopula::BiCopCDF(0.9999, 0.9999, copula)) - 1

  return(c(returnLevel1, returnLevel2, chi, chiBarre))
}


#' Select either the biGPD approach or the copula approach to get return levels from a given bivariate return period.
#'
#' @param data A dataframe with two columns, one for each time series.
#' @param returnPeriod An integer, in years.
#' @param EGPDtypes A vector of size 2, one for each time series. Values are 1 or 4. See Naveau et al. (2016) and the package mev for further information.
#' @param initParams1 A vector of parameters to initialize the estimation for the first time series. For type 1, it is (Kappa, Sigma, Xi), for type 4, it is (Prob, Kappa, Delta, Sigma, Xi).
#' @param initParams2 A vector of parameters to initialize the estimation for the second time series. For type 1, it is (Kappa, Sigma, Xi), for type 4, it is (Prob, Kappa, Delta, Sigma, Xi).
#' @param probaQuantile The probabilityof the value above which data are supposed to follow a GPD. A common value is 0.95.
#' @param nbDaysPerYear The number of days considered per year (integer).
#' @param nbYears The number of distinct years. Default value is 1.
#' @param h The parameter of non-concurrence (integer).
#' @param blockSizes Vector of integers of size 3. The sizes of the blocks considered for the declustering, with the univariates first and the bivariate block size at the end.
#' @param Dparam Cf documentation of the dgaps function of the exdex package. Default value is 3.
#' @param probaOccurrence The probability that the values are reached before the return period time. Default value to 0.63.
#' @param logic Either "AND" or "OR". Couples of bloxk maxima are retained if they both exceed their respective threshold in case of "AND", and if at least one of them exceed its respective threshold in case of "OR". Default value to "AND".
#' @return A list with, in that order: the return level of first variable, the return level of the second variable, chi and chiBarre.
#' @export

GetReturnLevels <- function(data, returnPeriod, EGPDtypes, initParams1, initParams2, probaQuantile, nbDaysPerYear, nbYears, h, blockSizes, Dparam, probaOccurrence, logic) {

  # Estimate chi
  threshold1 <- quantile(data[, 1], probaQuantile)
  threshold2 <- quantile(data[, 2], probaQuantile)
  dataAbove <- data[data[, 1] >= threshold1 & data[, 2] >= threshold2, ]
  FbarU1U2 <- length(dataAbove[, 1]) / length(data[, 1])
  chi <- FbarU1U2 / (1 - probaQuantile)

  if (chi > 0.1) {
    result <- BiGPDApproachReturnLevels(data, returnPeriod, EGPDtypes, initParams1, initParams2, probaQuantile, nbDaysPerYear, nbYears, h, Dparam, probaOccurrence)
  } else {
    result <- CopulaApproachReturnLevels(data, returnPeriod, probaQuantile, nbDaysPerYear, nbYears, h, blockSizes, Dparam, probaOccurrence, logic)
  }

  return(result)
}
