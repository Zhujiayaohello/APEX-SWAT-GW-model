#library(APEXSENSUN)
#----
getExampleFolder <- function(folderName = "Example") {
  if (file.exists(folderName)) {
    stop(paste(folderName)," already exists! Please select a different name for tutorial folder.")
  } else {
    folderLocation <- system.file("extdata", "Exampleee", package = "APEXSENSUN")
    if (file.exists("Exampleee")) {
      stop("Can not copy the tutorial folder!")
    }
    file.copy(from = folderLocation, to = getwd(), overwrite = T, recursive = T)
    file.rename(from = "Exampleee", folderName)
  }
}

#----
inputGen <- function (inputs) {
  
  if (hasArg(inputs)) {
    fileLocation <- system.file("extdata", "Packaged_Test.txt",
                                package = "APEXSENSUN")
    file.copy(from = fileLocation, to = inputs)
  } else {
    globalInputs <- list(sampleSize = 100 ,captionVarSim = "ET", captionVarObs = "ET",
                         startDate = "2002 01 01", endDate = "2003 01 01", labelAPEXExe = "APEX0806",
                         labelWatershedParam = "PARM", labelControlParam = "Apexcont",
                         labelOutputVariableAWP = "APEX001", labelOutputVariableACY = "APEX001",
                         labelOutputVariableDWS = "APEX001", labelObservedVar.txt = "observed.txt",
                         backUpPARM0806.dat = "Example/Back_Up/PARM.dat", backUpAPEXCONT.dat = "Example/Back_Up/Apexcont.dat",
                         folderPathProject = "Example/APEX",folderPathRCodes = getwd(),
                         folderPathObserved = "Example/Observed_Files", folderPathGsaOutputs = "Example/GSA_Outputs",
                         storeFolderPathWatershed = "Example/Generated_Inputs/PARM",
                         storeFolderPathControl = "Example/Generated_Inputs/CONT",
                         calculatedOutputFolderAWP = "Example/Calculated_Outputs/AWP",
                         calculatedOutputFolderACY = "Example/Calculated_Outputs/ACY",
                         calculatedOutputFolderDWS = "Example/Calculated_Outputs/DWS", gsaType = "SOBOL",
                         saParms = list(morrisRFactor = 15, morrisLevels = 8,  # Morris parameters.
                                        sobolOrder = 1,  # Sobol parameter.
                                        ksTestPerf = "NASH", ksTestThreshold = 0.2, ksTestSigmaLevel = 0.05),  #Ks-test.
                         apexPARM = genParmRange(), apexCONT = genContRange())  # Apex model parameters
    return(globalInputs)
  }
  
}

#----
dws2perf <- function(observedFilePath, dwsFolderPath,
                     startDate, endDate,
                     captionDwsVar, captionObsVar,
                     TW="day", perfMatrixFileName=NA) {
  
  perfMatrixFolderPath <- getwd()
  
  sampleSize <- length(list.files(path = dwsFolderPath, pattern = ".SAD"))
  namesList <- list.files(path = dwsFolderPath)
  dwsCaptions <- list.files(path = dwsFolderPath, pattern = "1.SAD")
  shortestCharLength <- .Machine$double.xmax
  for (caption in dwsCaptions) {
    shortestCharLength <- min(shortestCharLength, nchar(caption))
  }
  dwsRootName<- substr(dwsCaptions[1],start = 1 ,stop = (shortestCharLength-5))
  obsTimeseries <- obs2timeseries(observedFilePath, captionObsVar,
                                  startDate, endDate)
  
  perfMat <- data.frame(Index = rep(NaN, sampleSize[1]),
                        RMSE = rep(NaN, sampleSize[1]),
                        NASH = rep(NaN, sampleSize[1]),
                        PBIAS = rep(NaN, sampleSize[1]),
                        MEAN = rep(NaN, sampleSize[1]))
  
  perfMat$Index <- 1:sampleSize[1]
  
  for (i in 1:sampleSize[1]) {
    simTimeseries <- dws2timeseries(paste(dwsFolderPath,
                                          "/", dwsRootName,
                                          toString(i), ".sad",
                                          sep = ""),
                                    captionDwsVar, startDate, endDate)
    
    perfMat['RMSE'][i,1] <- outputRMSEProTW(simTimeseries, obsTimeseries,
                                            startDate, endDate, TW)
    perfMat['NASH'][i,1] <- outputNASHProTW(simTimeseries, obsTimeseries,
                                            startDate, endDate, TW)
    perfMat['PBIAS'][i,1] <- outputPBIASProTW(simTimeseries, obsTimeseries,
                                              startDate, endDate, TW)
    perfMat['MEAN'][i,1] <- outputMEANProTW(simTimeseries,
                                            startDate, endDate, TW)
    
    print(paste("Files Processed...", toString(i), sep = ""))
    
  }
  
  perfMat<-na.omit(perfMat)
  if (!is.na(perfMatrixFileName)) {
    save(perfMat,file = paste(perfMatrixFileName,".RData", sep = ""))
    print("performance measure matrix was written successfully!")
  }
  
  return(perfMat)
}
#----
#
# Helper functions
#
obs2timeseries <- function(observedFilePath, captionObsVar,
                           startDate, endDate) {
  #Creates a dataframe containing observed variable.
  #
  ATemp <- read.table(observedFilePath, header = TRUE)
  
  firstDwsDate <- paste(substr(x = ATemp$Date[1],start = 1,stop = 4),
                        "-", substr(x = ATemp$Date[1],start = 5,stop = 6),
                        "-", substr(x = ATemp$Date[1],start = 7,stop = 8),
                        sep = "")
  
  lastDwsDate <- paste(substr(x= ATemp$Date[nrow(ATemp)], start = 1,stop = 4),
                       "-", substr(x= ATemp$Date[nrow(ATemp)], start = 5, stop = 6),
                       "-", substr(x= ATemp$Date[nrow(ATemp)], start = 7, stop = 8),
                       sep = "")
  
  ATemp$Date <- as.character(seq.Date(as.Date(firstDwsDate),
                                      as.Date(lastDwsDate),
                                      by = "day"))
  dateSeq <- as.character(seq.Date(as.Date(startDate),
                                   as.Date(endDate),
                                   by = "day"))
  idxAllDates <- 1:nrow(ATemp)
  idxStartDate <- idxAllDates[ATemp$Date == startDate]
  idxEndDate <- idxAllDates[ATemp$Date == endDate]
  timeSeriesObs <- array(NaN,
                         dim = c(length(dateSeq),
                                 length(captionObsVar) +1))
  timeSeriesObs[, 1] <- dateSeq
  colNames <- c("Date", captionObsVar)
  colnames(timeSeriesObs) <- colNames
  timeSeriesObs <- as.data.frame(timeSeriesObs)
  
  for (i in (1:(length(colNames) - 1))) {
    
    timeSeriesObs[colNames[i + 1]] <- ATemp[captionObsVar[i]][idxStartDate:idxEndDate,]
  }
  return(timeSeriesObs)
}
#
#
dws2timeseries <- function(fileName, captionDwsVar,
                           startDate, endDate) {
  
  nameLineNumber <-9
  connSource <- file(description = fileName,open = "r+",blocking = TRUE)
  readLineObj <- readLines(con = connSource, n = -1)
  nameLine <- strsplit(x = readLineObj[nameLineNumber], split = " ")
  nameLineVector <- array(data = NaN, dim = c(1, length(nameLine[[1]])))
  
  i=1
  for (elm in nameLine[[1]]) {
    nameLineVector[i] <- nameLine[[1]][i]
    i=i+1
  }
  
  nameLineVector <- nameLineVector[nameLineVector!=""]
  
  ATemp <- as.data.frame(array(data = NaN,
                               dim = c(length(readLineObj)
                                       ,length(nameLineVector)-2)))
  
  colnames(ATemp) <- c("Date",nameLineVector[4:length(nameLineVector)])
  
  firstDwsDate <- paste(toString(strtoi(substr(readLineObj[nameLineNumber+1]
                                               ,start = 1, stop = 5))),
                        "-",
                        toString(strtoi(substr(readLineObj[nameLineNumber+1]
                                               ,start = 6, stop = 9))),
                        "-",
                        toString(strtoi(substr(readLineObj[nameLineNumber+1]
                                               ,start = 10, stop = 13)))
                        ,sep="")
  
  ATemp['Date'] <- as.character(seq.Date(from = as.Date(firstDwsDate),
                                         by = "day",
                                         length.out = length(readLineObj)))
  close(connSource)
  
  try(ATemp <- read.table(fileName, header = TRUE, skip = 8), silent = TRUE)
  ATemp$Date <- as.character(seq.Date(from = as.Date(firstDwsDate),
                                      by = "day",
                                      length.out = nrow(ATemp)))
  
  if (captionDwsVar=="TN") {
    ATemp$TN <- ATemp$QN + ATemp$YN + ATemp$QDRN +
      ATemp$RSFN + ATemp$QRFN + ATemp$SSFN
  }
  
  if(captionDwsVar=="TP") {
    ATemp$TP <- ATemp$QP + ATemp$YP + ATemp$QDRP + ATemp$QRFP
  }
  
  dateSeq <- as.character(seq.Date(as.Date(startDate),
                                   as.Date(endDate),
                                   by = "day"))
  idxAllDates <- 1:nrow(ATemp)
  idxStartDate <- idxAllDates[ATemp$Date == startDate]
  idxEndDate <- idxAllDates[ATemp$Date == endDate]
  timeSeriesDws <- array(NaN,
                         dim = c(length(dateSeq),
                                 length(captionDwsVar) + 1))
  timeSeriesDws[, 1] <- dateSeq
  colNames <- c("Date", captionDwsVar)
  colnames(timeSeriesDws) <- colNames
  timeSeriesDws <- as.data.frame(timeSeriesDws)
  
  for (i in (1:(length(colNames) - 1))) {
    timeSeriesDws[colNames[i + 1]] <- ATemp[captionDwsVar[i]][idxStartDate:idxEndDate,]
  }
  return(timeSeriesDws)
}
#
#
outputNASHProTW <- function (simSeries, obsSeries,
                             startDate, endDate,
                             TW="day") {
  # Calculates Nash-Sutcliffe performance measure.
  simSeriesAgg <- agg4timeseries(timeSeries = simSeries,
                                 startDate = startDate,
                                 endDate = endDate, TW = TW)
  obsSeriesAgg <- agg4timeseries(timeSeries = obsSeries,
                                 startDate = startDate,
                                 endDate = endDate, TW = TW)
  NsNumerator <- sum((simSeriesAgg[[2]] - obsSeriesAgg[[2]])^2)
  NsDenominator <- sum((obsSeriesAgg[[2]] - mean(obsSeriesAgg[[2]], na.rm = TRUE))^2)
  NASH <- 1 - (NsNumerator/NsDenominator)
  
  if (is.na(NASH)) {
    NASH <- -999
  }
  return(NASH)
}
#
#
outputRMSEProTW <- function (simSeries, obsSeries,
                             startDate, endDate, TW="day") {
  # Calculates RMSE performance measure
  simSeriesAgg <- agg4timeseries(timeSeries = simSeries, startDate = startDate,
                                 endDate = endDate, TW = TW)
  obsSeriesAgg <- agg4timeseries(timeSeries = obsSeries, startDate = startDate,
                                 endDate = endDate, TW = TW)
  RMSESeries <- sqrt(mean((simSeriesAgg[[2]] - obsSeriesAgg[[2]])^2, na.rm = TRUE))
  if (is.na(RMSESeries)) {
    RMSESeries <- -999
  }
  return(RMSESeries)
}
#
#
outputMEANProTW <- function (simSeries, startDate, endDate, TW="day") {
  #Calculates mean output (No observed data is used).
  simSeriesAgg <- agg4timeseries(timeSeries = simSeries, startDate = startDate,
                                 endDate = endDate, TW = TW)
  
  MEANSeries <- mean(simSeriesAgg[[2]], na.rm = TRUE)
  if (is.na(MEANSeries)) {
    MEANSeries <- -999
  }
  return(MEANSeries)
}
#
#
outputPBIASProTW <- function (simSeries, obsSeries,
                              startDate, endDate, TW="day") {
  #Calculates PBIAS performance measure
  simSeriesAgg <- agg4timeseries(timeSeries = simSeries, startDate = startDate,
                                 endDate = endDate, TW = TW)
  obsSeriesAgg <- agg4timeseries(timeSeries = obsSeries, startDate = startDate,
                                 endDate = endDate, TW = TW)
  NsNumerator <- sum(obsSeriesAgg[[2]] - simSeriesAgg[[2]])*100
  NsDenominator <- sum(obsSeriesAgg[[2]])
  PBIAS <- NsNumerator/NsDenominator
  if (is.na(PBIAS)) {
    PBIAS <- -999
  }
  return(PBIAS)
}
#
#
agg4timeseries <- function(timeSeries, startDate, endDate, TW="day") {
  #Aggregates daily time series.
  if (TW == "day") {
    return(timeSeries)
  } else {
    aggDateSeq <- as.character(seq.Date(from = as.Date(startDate), to = as.Date(endDate), by = TW))
    idxSeq <- 1:nrow(timeSeries)
    idxSeqStart <- idxSeq[timeSeries$Date %in% aggDateSeq[1:(length(aggDateSeq)-1)]]
    idxSeqEnd <- (idxSeq[timeSeries$Date %in% aggDateSeq][-1]) - 1
    idxSeqEnd[length(idxSeqEnd)] <- idxSeqEnd[length(idxSeqEnd)] + 1
    #Empty data frame for storing aggregated time series...
    aggTimeseries <- as.data.frame(array(data = NaN,dim = c(length(idxSeqStart), 2)))
    names(aggTimeseries) <- names(timeSeries)
    aggTimeseries$Date <- aggDateSeq[1:length(aggDateSeq)-1]
    
    for (i in 1:nrow(aggTimeseries)) {
      aggTimeseries[[2]][i] <- mean(c(timeSeries[[2]][idxSeqStart[i]:idxSeqEnd[i]]), na.rm = TRUE)
    }
    return(aggTimeseries)
    
  }
}

#----
genContRange <- function (controlRangeFile) {
  
  contParamDataframe <- watershedParamDataframe <- as.data.frame(array(data = -1, dim = c(2,37)))
  controlParamNames <- c("RFN", "CO2", "CQN", "PSTX", "YWI", "BTA",
                         "EXPK", "QG", "QCF", "CHSO", "BWD", "FCW",
                         "FPSC", "GWSO", "RFTO", "RFPO", "SATO", "FL",
                         "FW","ANG", "UXP", "DIAM", "ACW", "GZL0",
                         "RTN0", "BXCT", "BYCT", "DTHY", "QTH", "STND",
                         "DRV", "PCO0","RCC0", "CSLT", "BUS1", "BUS2", "BUS3")
  
  names(contParamDataframe) <- controlParamNames
  
  
  if (hasArg(controlRangeFile)) {
    write.table(contParamDataframe, file = paste(controlRangeFile, sep=""), sep = "\t",
                row.names = FALSE)
  }
  contParamDataframe
}

#----
genParmRange <- function(parmRangeFile) {
  watershedParamDataframe <- as.data.frame(array(data = -1, dim = c(2,103)))
  watershedParamDataframeNames <- c("Crop_canopy_PET", "Root_growth_soil", "Water_stress_harvest",
                                    "Water_storage_N", "Soil_water_limit", "Winter_dormancy",
                                    "N_fixation", "Soluble_P_runoff", "Pest_damage_moisture",
                                    "Pest_damage_cover", "Moisture_req_seed_germ", "Soil_evap_coeff",
                                    "Wind_erod_coeff", "Nitrate_leac_ratio", "Runoff_CN_Adj_parm",
                                    "Expand_CN_ret_parm", "Soil_evap_plant_cover", "Sedim_rout_exponent",
                                    "Sedim_rout_coeff", "Runoff_CN_int_abs", "Soluble_C_adsorp_coeff",
                                    "CN_retention_frozen_soil", "Harg_equation_parm", "Pest_leach_ratio",
                                    "Expo_coeff_rainfall", "Matur_frac_spring", "CEC_effect_nitrification",
                                    "N_fixation_limit", "Biological_mix_efficiency", "Soluble_P_exponent",
                                    "Max_depth_bio_mixing", "OrgP_loss_exponent", "MUST_coeff",
                                    "Harg_PET_exponent", "Denit_soil_threshold", "Daily_denit_limit",
                                    "SWAT_delivery_ratio_exponent", "Water_stress_coeff", "Puddling_sat_conduct",
                                    "Groundwater_stor_threshold", "Root_temp_stress_exponent", "SCS_index_coeff",
                                    "Plow_depth", "CN_retention_param", "sediment_rout_travel_coeff",
                                    "RUSLE_c_factor_res", "RUSLE_c_factor_height", "Climate_stress_factor",
                                    "Max_rain_intercept", "Rain_intercept_coeff", "Water_stor_residue_coeff",
                                    "Tillage_residue_decay_rate_coeff", "Microbial_soil_depth_coeff", "N_enrich_coeff",
                                    "N_enrich_rout_exponent", "Fraction_destroyed_burn", "P_enrich_rout_coeff",
                                    "P_enrich_rout_exponent", "P_move_evap_coeff", "Max_days_grazed_rotation",
                                    "Soil_water_up_flow_limit", "Manure_erosion_equation_coeff", "N_enrich_ratio_delivery",
                                    "Dust_distribution_coeff", "RUSLE2_trans_capacity", "RUSLE2_trans_capacity_threshold",
                                    "Dust_distribution_exponent", "Manure_erosion_exponent", "Microbial_top_soil_coeff",
                                    "Microbial_decay_coeff", "Manure_erosion_coeff", "Volt_nitrification_partition_coeff",
                                    "Hydrograph_dev_param", "Partition_N_flow_groundwater", "P_enrich_ratio_deliver_SWAT",
                                    "Stand_dead_fall_rate_coeff", "Runoff_2_delay_pest", "Soil_water_2_delay_tillage",
                                    "Auto_mov_lower_limit", "Nitrification_vol_upper_limit", "Tech_coeff",
                                    "Drainage_lateral_conduct", "P_flux_labile_active_coeff", "P_flux_active_stable_coeff",
                                    "N_salt_evap_coeff", "Water_table_recession_coeff", "Water_table_move_limit",
                                    "Water_table_recession_exponent", "Subsurface_flow_factor", "Flood_evap_limit",
                                    "Runoff_adj_link", "Water_erosion_threshold", "Wind_erosion_threshold",
                                    "Crop_stress_temp_exponent", "Soluble_P_leach_KD", "Unknown1",
                                    "Unknown2", "Unknown3", "Irrigation_cost",
                                    "Lime_cost", "Fuel_cost", "Labor_cost", "Unknown4")
  names(watershedParamDataframe) <- watershedParamDataframeNames
  
  
  
  
  
  if (hasArg(parmRangeFile)) {
    write.table(watershedParamDataframe, file = paste(parmRangeFile, sep=""), sep = "\t",
                row.names = FALSE)
  }
  
  watershedParamDataframe
}

#----
idx2parm <- function (xMat, idxAccept) {
  
  idxParm <- 1:nrow(xMat)
  logicExtract <- idxParm %in% idxAccept
  xMatPost <- xMat[logicExtract, ]
  return(xMatPost)
}

idx2statTimeseries <- function(dwsFolderPath, startDate, endDate,
                               captionDwsVar, TW="day",
                               acceptedIndices="all",
                               probs=c(0.025,0.5,0.975)) {
  
  if (is.character(acceptedIndices)) {
    acceptedIndices=1:length(list.files(path = dwsFolderPath, pattern = ".SAD"))
  }
  
  ensembleMat <- createEnsembleMatrix(dwsFolderPath, startDate, endDate,
                                      captionDwsVar, TW, acceptedIndices)
  
  statTimeseriesMat <- as.data.frame(array(data = NaN,
                                           dim = c(nrow(ensembleMat),
                                                   (length(probs) + 2))))
  quantileNames<-c()
  for (i in (1:length(probs)))  {
    quantileNames[i] <- paste("Quantile_", toString(probs[i]), sep="")
  }
  
  statNames <- c("Date","Mean",quantileNames)
  names(statTimeseriesMat) <- statNames
  statTimeseriesMat['Date'] <- ensembleMat['Date']
  statTimeseriesMat['Mean'] <- apply(X = ensembleMat[,2:ncol(ensembleMat)], 1,
                                     function(x) mean(x = x))
  
  for (i in (1:length(probs))) {
    statTimeseriesMat[quantileNames[i]] <- apply(X = ensembleMat[,2:ncol(ensembleMat)], MARGIN = 1,
                                                 FUN = function (x) stats::quantile(x=x,probs=probs[i], type=1))
  }
  return(statTimeseriesMat)
}
#
#
# Helper Functions
createEnsembleMatrix <- function(dwsFolderPath, startDate, endDate,
                                 captionDwsVar, TW="day", acceptedIndices="all",
                                 timeSeriesFileName="dws_timeseries.txt") {
  #Detecting .DWS root name..
  namesList <- list.files(path = dwsFolderPath,pattern = ".SAD")
  dwsCaptions <- list.files(path = dwsFolderPath,pattern = "1.SAD")
  shortestCharLength <- .Machine$double.xmax
  for (caption in dwsCaptions) {
    shortestCharLength <- min(shortestCharLength, nchar(caption))
  }
  rootName<- substr(dwsCaptions[1], start = 1, stop = (shortestCharLength-5))
  
  if (is.character(acceptedIndices)) {
    acceptedIndices <- 1:length(namesList)
    
  }
  
  originalDws<- dws2timeseries(paste(dwsFolderPath,"/", rootName,toString(1), ".SAD", sep = ""),
                               captionDwsVar, startDate, endDate)
  originalAggDws <- agg4timeseries(originalDws, startDate, endDate,TW)
  ensembleMatrix <- array(data = NaN, dim = c(nrow(originalAggDws), (length(acceptedIndices)+1)))
  colnames(ensembleMatrix) <- c("Date", paste("Sim", acceptedIndices, sep=""))
  ensembleMatrix <- as.data.frame(ensembleMatrix)
  ensembleMatrix['Date'] <- originalAggDws['Date']
  
  for (i in (1:length(acceptedIndices))) {
    tempFileName <- paste(dwsFolderPath, "/", rootName, acceptedIndices[i], ".SAD", sep = "")
    temp_time_series <- dws2timeseries(tempFileName, captionDwsVar, startDate, endDate)
    ensembleMatrix[paste("Sim", acceptedIndices[i], sep="")] <-
      (agg4timeseries(temp_time_series, startDate, endDate, TW))[captionDwsVar]
  }
  return(ensembleMatrix)
}

#----
idx2Timeseries <- function(dwsFolderPath, startDate, endDate,
                           captionDwsVar, TW="day", acceptedIndices="all") {
  #Detecting .DWS root name..
  namesList <- list.files(path = dwsFolderPath,pattern = ".SAD")
  dwsCaptions <- list.files(path = dwsFolderPath,pattern = "1.SAD")
  shortestCharLength <- .Machine$double.xmax
  for (caption in dwsCaptions) {
    shortestCharLength <- min(shortestCharLength,nchar(caption))
  }
  rootName<- substr(dwsCaptions[1], start = 1, stop = (shortestCharLength-5))
  
  if (is.character(acceptedIndices)) {
    acceptedIndices <- 1:length(namesList)
  }
  
  
  originalDws<- dws2timeseries(paste(dwsFolderPath, "/", rootName, toString(1), ".SAD", sep = ""),
                               captionDwsVar, startDate, endDate)
  
  originalAggDws <- agg4timeseries(originalDws, startDate, endDate, TW)
  ensembleMatrix <- array(data = NaN, dim = c(nrow(originalAggDws), (length(acceptedIndices)+1)))
  colnames(ensembleMatrix) <- c("Date", paste("Sim", acceptedIndices, sep=""))
  ensembleMatrix <- as.data.frame(ensembleMatrix)
  ensembleMatrix['Date'] <- originalAggDws['Date']
  
  
  for(i in 1:length(acceptedIndices)) {
    temp_file_name <- paste(dwsFolderPath, "/", rootName, acceptedIndices[i], ".SAD", sep = "")
    temp_time_series <- dws2timeseries(temp_file_name, captionDwsVar, startDate, endDate)
    ensembleMatrix[paste("Sim", acceptedIndices[i], sep="")] <-
      (agg4timeseries(temp_time_series, startDate, endDate,TW))[captionDwsVar]
  }
  return(ensembleMatrix)
}


#----
mc4APEX <- function(input, usingTextFile = FALSE, savingInput4SA = TRUE) {
  if (!isTRUE(usingTextFile)) {
    globalInputs <- input
    sampleSize <- globalInputs$sampleSize
    type <- globalInputs$gsaType
    saParms <- globalInputs$saParms
    
    # Extracting uncertain model parameters from input list
    allBasinParameters <- input$apexPARM
    allControlParameters <- input$apexCONT
    uncertainBasinParameters <- paramExtract(allBasinParameters)
    uncertainControlParameters <- paramExtract(allControlParameters)
    if (ncol(uncertainBasinParameters) !=0 && ncol(uncertainControlParameters) != 0) {
      allUncertainModelParameters <- cbind(uncertainBasinParameters,
                                           uncertainControlParameters)
    }
    else if (ncol(uncertainBasinParameters) !=0 && ncol(uncertainControlParameters) == 0) {
      allUncertainModelParameters <- uncertainBasinParameters
    }
    else if (ncol(uncertainBasinParameters) == 0 && ncol(uncertainControlParameters) != 0) {
      allUncertainModelParameters <- uncertainControlParameters
    } else {
      stop("No uncertain parameter was found!")
    }
    
    # Creating gsaObject...
    gsaObject <- createGsaObject(sampleSize, allUncertainModelParameters,
                                 type, saParms)
    input4SA <- list(gsaObject = gsaObject, basinParms = uncertainBasinParameters,
                     controlParms = uncertainControlParameters)
  } else {
    globalInputs <- text2InputSimulation(input)
    input4SA <- input4MCSimulation(input)
  }
  
  
  if (isTRUE(savingInput4SA)) {
    save(input4SA,file = paste(globalInputs$folderPathGsaOutputs,
                               "/", globalInputs$gsaType, "_",
                               "GSA_WorkSpace", ".RData", sep = ""))
  }
  
  # Running Monte Carlo simulation
  runAPEXModel(input4SA, globalInputs$folderPathProject,
               globalInputs$labelWatershedParam, globalInputs$backUpPARM0806.dat,
               globalInputs$labelControlParam, globalInputs$backUpAPEXCONT.dat,
               globalInputs$labelAPEXExe, globalInputs$folderPathRCodes,
               globalInputs$storeFolderPathWatershed, globalInputs$storeFolderPathControl,
               globalInputs$labelOutputVariableAWP, globalInputs$calculatedOutputFolderAWP,
               globalInputs$labelOutputVariableACY, globalInputs$calculatedOutputFolderACY,
               globalInputs$labelOutputVariableDWS, globalInputs$calculatedOutputFolderDWS)
  
  input4SA
  
}

# Helper Functions
#
#
paramExtract <- function(paramRange) {
  # This function extracts uncertain parameters from a dataframe containing parameter bounds.
  #
  # Args:
  #   paramRange: A dataframe containing parameter bounds.
  #
  #
  # Returns:
  #   A dataframe containing uncertain model parameters only.
  #
  idx_seq <- 1:ncol(paramRange)
  paramRange <- paramRange[idx_seq[as.vector(paramRange[1,] != paramRange[2,])]]
  return(paramRange)
}
#
#
runAPEXModel <- function(input4SA, folderPathProject, labelWatershedParam,
                         backUpPARM0806.dat, labelControlParam, backUpAPEXCONT.dat,
                         labelAPEXExe, folderPathRCodes, storeFolderPathWatershed,
                         storeFolderPathControl, labelOutputVariableAWP, calculatedOutputFolderAWP,
                         labelOutputVariableACY, calculatedOutputFolderACY,
                         labelOutputVariableDWS, calculatedOutputFolderDWS) {
  # Performs Monte Carlo simulation
  #
  # Args:
  #  input4SA: A list contaning a gsaObject and uncertain parameter bounds.
  #  folderPathProject: Folder path to the APEX model folder.
  #  labelWatershedParam: File label containing watershed parameters.
  #  backUpPARM0806.dat: File label containing original watershed parameters.
  #  labelControlParam: File label containing control parameters.
  #  backUpAPEXCONT.dat: File label containing original control parameters.
  #  labelAPEXExe: File label representing APEX executable model.
  #  folderPathRCodes: Working directory
  #  storeFolderPathWatershed: Folder path storing watershed parameters for Monte Carlo simulations.
  #  storeFolderPathControl: Folder path storing control parameters for Monte Carlo simulations.
  #  labelOutputVariableAWP: File label for .AWP output files.
  #  calculatedOutputFolderAWP: Folder path storing .AWP output files for Monte Carlo simulations.
  #  labelOutputVariableACY: File label for .ACY output files.
  #  calculatedOutputFolderACY: Folder path storing .ACY output files for Monte Carlo simulations.
  #  labelOutputVariableDWS: File label for .DWS output files.
  #  calculatedOutputFolderDWS: Folder path storing .DWS output files for Monte Carlo simulations.
  #
  # Returns:
  #  Void
  #
  folderPathRCodes <-getwd()
  sampleSize <- nrow(input4SA$gsaObject$X)
  gsaObject <- input4SA$gsaObject
  uncertainBasinParameters <- input4SA$basinParms
  uncertainControlParameters <- input4SA$controlParms
  
  for (i in 1:sampleSize) {
    tempWatershedFilePath <- paste(folderPathProject,
                                   "/", labelWatershedParam, ".dat", sep = "")
    tempControlFilePath <- paste(folderPathProject,
                                 "/", labelControlParam, ".dat", sep = "")
    if (ncol(uncertainBasinParameters)!=0) {
      basinParm2Write <- as.data.frame(gsaObject$X[i, 1:ncol(uncertainBasinParameters)])
      colnames(basinParm2Write) <- colnames(uncertainBasinParameters)
      writeBasinParameters(basinParm2Write, tempWatershedFilePath, backUpPARM0806.dat)
    }
    
    if (ncol(uncertainControlParameters) != 0) {
      contParm2Write <- as.data.frame(gsaObject$X[i, (ncol(uncertainBasinParameters) +
                                                        1):(ncol(uncertainControlParameters) +
                                                              ncol(uncertainBasinParameters))])
      colnames(contParm2Write) <- colnames(uncertainControlParameters)
      writeControlParameters(contParm2Write, tempControlFilePath, backUpAPEXCONT.dat)
    }
    
    setwd(folderPathProject)
    system(paste(labelAPEXExe, ".exe", sep = ""))
    setwd(folderPathRCodes)
    file.copy(tempWatershedFilePath, paste(storeFolderPathWatershed,
                                           "/", labelWatershedParam, toString(i), ".dat",
                                           sep = ""),overwrite = TRUE)
    file.copy(tempControlFilePath, paste(storeFolderPathControl,
                                         "/", labelControlParam, toString(i),
                                         ".dat", sep = ""),overwrite = TRUE)
    tempOutputFilepathAWP <- paste(folderPathProject,
                                   "/", labelOutputVariableAWP, ".AWP", sep = "")
    file.copy(tempOutputFilepathAWP, paste(calculatedOutputFolderAWP,
                                           "/", labelOutputVariableAWP, toString(i), ".AWP",
                                           sep = ""),overwrite = TRUE)
    tempOutputFilepathACY <- paste(folderPathProject,
                                   "/", labelOutputVariableACY, ".ACY", sep = "")
    file.copy(tempOutputFilepathACY, paste(calculatedOutputFolderACY,
                                           "/", labelOutputVariableACY, toString(i), ".ACY",
                                           sep = ""),overwrite = TRUE)
    tempOutputFilepathDWS <- paste(folderPathProject,
                                   "/", labelOutputVariableDWS, ".SAD", sep = "")
    file.copy(tempOutputFilepathDWS, paste(calculatedOutputFolderDWS,
                                           "/", labelOutputVariableDWS, toString(i), ".SAD",
                                           sep = ""),overwrite = TRUE)
    
  }
}
#
#
createGsaObject <- function (sampleSize, allUncertainModelParameters, type,saParms){
  
  lowerLimits <- as.vector(allUncertainModelParameters[1,], mode = "numeric")
  upperLimits <- as.vector(allUncertainModelParameters[2,], mode = "numeric")
  type = toupper(type)
  if (type == "MORRIS") {
    gsaObject <- sensitivity::morris(model = NULL, factor = names(allUncertainModelParameters),
                                     r = as.numeric(saParms$morrisRFactor) ,
                                     design = list(type = "oat", levels = as.numeric(saParms$morrisLevels),                                                              grid.jump = 1),
                                     binf = lowerLimits, bsup = upperLimits, scale = TRUE)
    gsaObject$X <- as.data.frame(gsaObject$X);
    gsaObject
  }
  else if (type == "SRC") {
    X <- allUncertainModelParameters[-1, ]
    X[1, ] = array(NaN, dim = c(1, length(X)))
    X <- dataFrameAug(X, sampleSize)
    xSize <- dim(X)
    i = 1
    while (i <= xSize[2]) {
      X[i] <- runif(sampleSize, min = lowerLimits[i],
                    max = upperLimits[i])
      i = i + 1
    }
    gsaObject = list(X = X)
  }
  else if (type == "SRRC") {
    X <- allUncertainModelParameters[-1, ]
    X[1, ] = array(NaN, dim = c(1, length(X)))
    X <- dataFrameAug(X, sampleSize)
    xSize <- dim(X)
    i = 1
    while (i <= xSize[2]) {
      X[i] <- runif(sampleSize, min = lowerLimits[i],
                    max = upperLimits[i])
      i = i + 1
    }
    gsaObject = list(X = X)
  }
  else if (type == "SOBOL") {
    X1 <- allUncertainModelParameters[-1, ]
    X2 <- allUncertainModelParameters[-1, ]
    X1[1, ] = array(NaN, dim = c(1, length(X1)))
    X2[1, ] = array(NaN, dim = c(1, length(X2)))
    X1 <- dataFrameAug(X1, sampleSize)
    X2 <- dataFrameAug(X2, sampleSize)
    xSize <- dim(X1)
    i = 1
    while (i <= xSize[2]) {
      X1[i] <- runif(sampleSize, min = lowerLimits[i],
                     max = upperLimits[i])
      X2[i] <- runif(sampleSize, min = lowerLimits[i],
                     max = upperLimits[i])
      i = i + 1
    }
    gsaObject <- sensitivity::sobol(model = NULL, X1 = X1, X2 = X2,
                                    order = saParms$sobolOrder,
                                    nboot = 0, conf = 0.95)
  }
  else if (type == "SOBOL2002") {
    X1 <- allUncertainModelParameters[-1, ]
    X2 <- allUncertainModelParameters[-1, ]
    X1[1, ] = array(NaN, dim = c(1, length(X1)))
    X2[1, ] = array(NaN, dim = c(1, length(X2)))
    X1 <- dataFrameAug(X1, sampleSize)
    X2 <- dataFrameAug(X2, sampleSize)
    xSize <- dim(X1)
    i = 1
    while (i <= xSize[2]) {
      X1[i] <- runif(sampleSize, min = lowerLimits[i],
                     max = upperLimits[i])
      X2[i] <- runif(sampleSize, min = lowerLimits[i],
                     max = upperLimits[i])
      i = i + 1
    }
    gsaObject <- sensitivity::sobol2002(model = NULL, X1 = X1, X2 = X2,
                                        nboot = 0, conf = 0.95)
  }
  else if (type == "SOBOL2007") {
    X1 <- allUncertainModelParameters[-1, ]
    X2 <- allUncertainModelParameters[-1, ]
    X1[1, ] = array(NaN, dim = c(1, length(X1)))
    X2[1, ] = array(NaN, dim = c(1, length(X2)))
    X1 <- dataFrameAug(X1, sampleSize)
    X2 <- dataFrameAug(X2, sampleSize)
    xSize <- dim(X1)
    i = 1
    while (i <= xSize[2]) {
      X1[i] <- runif(sampleSize, min = lowerLimits[i],
                     max = upperLimits[i])
      X2[i] <- runif(sampleSize, min = lowerLimits[i],
                     max = upperLimits[i])
      i = i + 1
    }
    gsaObject <- sensitivity::sobol2007(model = NULL, X1 = X1, X2 = X2,
                                        nboot = 0, conf = 0.95)
  }
  else if (type == "SOBOLEFF") {
    X1 <- allUncertainModelParameters[-1, ]
    X2 <- allUncertainModelParameters[-1, ]
    X1[1, ] = array(NaN, dim = c(1, length(X1)))
    X2[1, ] = array(NaN, dim = c(1, length(X2)))
    X1 <- dataFrameAug(X1, sampleSize)
    X2 <- dataFrameAug(X2, sampleSize)
    xSize <- dim(X1)
    i = 1
    while (i <= xSize[2]) {
      X1[i] <- runif(sampleSize, min = lowerLimits[i],
                     max = upperLimits[i])
      X2[i] <- runif(sampleSize, min = lowerLimits[i],
                     max = upperLimits[i])
      i = i + 1
    }
    gsaObject <- sobolEff(model = NULL, X1 = X1, X2 = X2,
                          order = 1, nboot = 0, conf = 0.95)
  }
  else if (type == "SOBOLJANSEN") {
    X1 <- allUncertainModelParameters[-1, ]
    X2 <- allUncertainModelParameters[-1, ]
    X1[1, ] = array(NaN, dim = c(1, length(X1)))
    X2[1, ] = array(NaN, dim = c(1, length(X2)))
    X1 <- dataFrameAug(X1, sampleSize)
    X2 <- dataFrameAug(X2, sampleSize)
    xSize <- dim(X1)
    i = 1
    while (i <= xSize[2]) {
      X1[i] <- runif(sampleSize, min = lowerLimits[i],
                     max = upperLimits[i])
      X2[i] <- runif(sampleSize, min = lowerLimits[i],
                     max = upperLimits[i])
      i = i + 1
    }
    gsaObject <- sensitivity::soboljansen(model = NULL, X1 = X1, X2 = X2,
                                          nboot = 0, conf = 0.95)
  }
  else if (type == "SOBOLMARA") {
    X <- allUncertainModelParameters[-1, ]
    X[1, ] = array(NaN, dim = c(1, length(X)))
    X <- dataFrameAug(X, sampleSize)
    xSize <- dim(X)
    i = 1
    while (i <= xSize[2]) {
      X[i] <- runif(sampleSize, min = lowerLimits[i],
                    max = upperLimits[i])
      i = i + 1
    }
    gsaObject <- sensitivity::sobolmara(model = NULL, X1 = X)
  }
  else if (type == "SOBOLMARTINEZ") {
    X1 <- allUncertainModelParameters[-1, ]
    X2 <- allUncertainModelParameters[-1, ]
    X1[1, ] = array(NaN, dim = c(1, length(X1)))
    X2[1, ] = array(NaN, dim = c(1, length(X2)))
    X1 <- dataFrameAug(X1, sampleSize)
    X2 <- dataFrameAug(X2, sampleSize)
    xSize <- dim(X1)
    i = 1
    while (i <= xSize[2]) {
      X1[i] <- runif(sampleSize, min = lowerLimits[i],
                     max = upperLimits[i])
      X2[i] <- runif(sampleSize, min = lowerLimits[i],
                     max = upperLimits[i])
      i = i + 1
    }
    gsaObject <- sensitivity::sobolmartinez(model = NULL, X1 = X1, X2 = X2,
                                            nboot = 0, conf = 0.95)
  }
  
  else if (type=="FAST99") {
    
    lowerLimits <- as.vector(allUncertainModelParameters[1,], mode = "numeric")
    upperLimits <- as.vector(allUncertainModelParameters[2,], mode = "numeric")
    limitList <- list()
    i=1
    while (i<=ncol(allUncertainModelParameters)) {
      limitList[[i]] <- list(min=lowerLimits[i],max=upperLimits[i])
      i=i+1
    }
    
    gsaObject <- sensitivity::fast99(model = NULL, factor = names(allUncertainModelParameters),
                                     n = sampleSize, M = 4,
                                     q = rep("qunif",ncol(allUncertainModelParameters)),q.arg = limitList)
    
  }
  else if (type=="KSTEST") {
    
    X <- allUncertainModelParameters[-1, ]
    X[1, ] = array(NaN, dim = c(1, length(X)))
    X <- dataFrameAug(X, sampleSize)
    xSize <- dim(X)
    i = 1
    while (i <= xSize[2]) {
      X[i] <- runif(sampleSize, min = lowerLimits[i],
                    max = upperLimits[i])
      i = i + 1
    }
    gsaObject = list(X = X)
    gsaObject$ksTestPerf <- saParms$ksTestPerf
    gsaObject$ksTestThreshold <- saParms$ksTestThreshold
    gsaObject$ksTestSigmaLevel <- saParms$ksTestSigmaLevel
  } else {
    stop("Please enter a valid method")
  }
  return(gsaObject)
  
}
#
#
text2InputSimulation <- function (inputFile) {
  # Read inputs from a text file.
  globalInputNames <- c("sampleSize","captionVarSim", "captionVarObs",
                        "startDate","endDate", "labelAPEXExe",
                        "labelWatershedParam", "labelControlParam", "labelOutputVariableAWP",
                        "labelOutputVariableACY", "labelOutputVariableDWS", "labelOutputXmat",
                        "labelObservedVar.txt", "labelTemplateBasin.txt", "labelTemplateControl.txt",
                        "backUpPARM0806.dat", "backUpAPEXCONT.dat", "folderPathProject",
                        "folderPathRCodes", "folderPathObserved", "folderPathGsaOutputs",
                        "storeFolderPathWatershed", "storeFolderPathControl", "calculatedOutputFolderAWP",
                        "calculatedOutputFolderACY", "calculatedOutputFolderDWS", "gsaType")
  
  globalInputLocations <- c(10, 15, 17,
                            22, 24, 29,
                            31, 33, 35,
                            37, 39, 41,
                            43, 48, 50,
                            52, 54, 56,
                            58, 60, 62,
                            64, 66, 68,
                            70, 72, 74)
  
  globalInputTypeVector <- c("integer",rep("character",26))
  
  
  saParmsInputNames <- c("morrisRFactor", "morrisLevels", "sobolOrder",
                         "ksTestPerf", "ksTestThreshold", "ksTestSigmaLevel")
  
  saParmsInputLocations <- c(78, 80, 83,
                             86, 88, 90)
  
  saParmsTypeVector <- c(rep("integer",3), "character", rep("double",2))
  
  
  
  # Creating global input list
  globalInputLength <- length(globalInputNames)
  globalInputs <- as.list(x = rep(NA,globalInputLength))
  names(globalInputs) <- globalInputNames
  for (variable in seq(1,globalInputLength)) {
    print(variable)
    globalInputs[variable] <- list2vec(inputFile = inputFile,
                                       lines2Skip = globalInputLocations[variable],
                                       varType = globalInputTypeVector[variable])
  }
  
  # Creating saParms input list
  saParmsInputLength <- length(saParmsInputNames)
  saParms <- as.list(x = rep(NA,saParmsInputLength))
  names(saParms) <- saParmsInputNames
  for (variable in seq(1,saParmsInputLength)) {
    saParms[variable] <- list2vec(inputFile = inputFile,
                                  lines2Skip = saParmsInputLocations[variable],
                                  varType = saParmsTypeVector[variable])
  }
  # Appending saParms to globalInputList
  globalInputs$saParms <- saParms
  
}
#
#
input4MCSimulation <- function (inputFile) {
  
  globalInputs <- text2InputSimulation(inputFile)
  labelTemplateBasin.txt <- globalInputs$labelTemplateBasin.txt
  labelTemplateControl.txt <- globalInputs$labelTemplateControl.txt
  sampleSize =globalInputs$sampleSize
  type=globalInputs$gsaType
  saParms <- globalInputs$saParms
  
  allBasinParameters <- read.table(labelTemplateBasin.txt,
                                   header = TRUE)
  allControlParameters <- read.table(labelTemplateControl.txt,
                                     header = TRUE)
  uncertainBasinParameters <- paramExtract(allBasinParameters)
  uncertainControlParameters <- paramExtract(allControlParameters)
  
  if (ncol(uncertainBasinParameters)!=0 && ncol(uncertainControlParameters)!=0) {
    allUncertainModelParameters <- cbind(uncertainBasinParameters,
                                         uncertainControlParameters)
  }
  else if (ncol(uncertainBasinParameters)!=0 && ncol(uncertainControlParameters)==0) {
    allUncertainModelParameters <- uncertainBasinParameters
    
  }else if (ncol(uncertainBasinParameters)==0 && ncol(uncertainControlParameters)!=0) {
    allUncertainModelParameters <- uncertainControlParameters
    
  }else {
    stop("No uncertain parameter was found!")
  }
  
  gsaObject <- createGsaObject(sampleSize,allUncertainModelParameters,
                               type,saParms)
  
  return(list(gsaObject = gsaObject,
              basinParms = uncertainBasinParameters,
              controlParms = uncertainControlParameters))
}
#
#
dataFrameAug <- function (dataFrame, repNumb) {
  dim_x <- dim(dataFrame)
  i = 1
  dataFrameAugment <- dataFrame
  i = 1
  while (i < repNumb) {
    dataFrameAugment <- rbind(dataFrameAugment, dataFrame,
                              make.row.names = FALSE)
    i = i + 1
  }
  dataFrameAugment
}
#
#
list2vec <- function (input_file, lines2skip, var_type) {
  list_of_vars <- read.table(input_file, skip = lines2skip,
                             colClasses = var_type, nrows = 1)
  vec_of_vars = NULL
  i = 1
  while (i <= length(list_of_vars)) {
    vec_of_vars[i] = list_of_vars[1, i]
    i = i + 1
  }
  return(vec_of_vars)
}


#----
outputIdxCombiner <- function (...) {
  
  listOfArgs <- list(...)
  idxCommon <- listOfArgs[[1]]
  for (elm in listOfArgs) {
    idxCommon <- intersect(idxCommon, elm)
  }
  print(idxCommon)
  idxCommon
}

#----
parmAndContDataframes <- function() {
  
}
#
#
watershed_cont_dataframe <- function() {
  #
  contParamDataframe <- watershedParamDataframe <- as.data.frame(array(data = -1,dim = c(2,37)))
  controlParamNames <- c("RFN", "CO2", "CQN", "PSTX", "YWI", "BTA",
                         "EXPK", "QG", "QCF", "CHSO", "BWD", "FCW",
                         "FPSC", "GWSO", "RFTO", "RFPO", "SATO", "FL",
                         "FW","ANG", "UXP", "DIAM", "ACW", "GZL0",
                         "RTN0", "BXCT", "BYCT", "DTHY", "QTH", "STND",
                         "DRV", "PCO0","RCC0", "CSLT", "BUS1", "BUS2", "BUS3")
  
  names(contParamDataframe) <- controlParamNames
  contParamDataframe
}
#
#
watershed_param_dataframe <- function() {
  #
  watershedParamDataframe <- as.data.frame(array(data = -1,dim = c(2,103)))
  watershedParamDataframeNames <- c("Crop_canopy_PET", "Root_growth_soil", "Water_stress_harvest",
                                    "Water_storage_N", "Soil_water_limit", "Winter_dormancy",
                                    "N_fixation", "Soluble_P_runoff", "Pest_damage_moisture",
                                    "Pest_damage_cover", "Moisture_req_seed_germ", "Soil_evap_coeff",
                                    "Wind_erod_coeff", "Nitrate_leac_ratio", "Runoff_CN_Adj_parm",
                                    "Expand_CN_ret_parm", "Soil_evap_plant_cover", "Sedim_rout_exponent",
                                    "Sedim_rout_coeff", "Runoff_CN_int_abs", "Soluble_C_adsorp_coeff",
                                    "CN_retention_frozen_soil", "Harg_equation_parm", "Pest_leach_ratio",
                                    "Expo_coeff_rainfall", "Matur_frac_spring", "CEC_effect_nitrification",
                                    "N_fixation_limit", "Biological_mix_efficiency", "Soluble_P_exponent",
                                    "Max_depth_bio_mixing", "OrgP_loss_exponent", "MUST_coeff",
                                    "Harg_PET_exponent", "Denit_soil_threshold", "Daily_denit_limit",
                                    "SWAT_delivery_ratio_exponent", "Water_stress_coeff", "Puddling_sat_conduct",
                                    "Groundwater_stor_threshold", "Root_temp_stress_exponent", "SCS_index_coeff",
                                    "Plow_depth", "CN_retention_param", "sediment_rout_travel_coeff",
                                    "RUSLE_c_factor_res", "RUSLE_c_factor_height", "Climate_stress_factor",
                                    "Max_rain_intercept", "Rain_intercept_coeff", "Water_stor_residue_coeff",
                                    "Tillage_residue_decay_rate_coeff", "Microbial_soil_depth_coeff", "N_enrich_coeff",
                                    "N_enrich_rout_exponent", "Fraction_destroyed_burn", "P_enrich_rout_coeff",
                                    "P_enrich_rout_exponent", "P_move_evap_coeff", "Max_days_grazed_rotation",
                                    "Soil_water_up_flow_limit", "Manure_erosion_equation_coeff", "N_enrich_ratio_delivery",
                                    "Dust_distribution_coeff", "RUSLE2_trans_capacity", "RUSLE2_trans_capacity_threshold",
                                    "Dust_distribution_exponent", "Manure_erosion_exponent", "Microbial_top_soil_coeff",
                                    "Microbial_decay_coeff", "Manure_erosion_coeff", "Volt_nitrification_partition_coeff",
                                    "Hydrograph_dev_param", "Partition_N_flow_groundwater", "P_enrich_ratio_deliver_SWAT",
                                    "Stand_dead_fall_rate_coeff", "Runoff_2_delay_pest", "Soil_water_2_delay_tillage",
                                    "Auto_mov_lower_limit", "Nitrification_vol_upper_limit", "Tech_coeff",
                                    "Drainage_lateral_conduct", "P_flux_labile_active_coeff", "P_flux_active_stable_coeff",
                                    "N_salt_evap_coeff", "Water_table_recession_coeff", "Water_table_move_limit",
                                    "Water_table_recession_exponent", "Subsurface_flow_factor", "Flood_evap_limit",
                                    "Runoff_adj_link", "Water_erosion_threshold", "Wind_erosion_threshold",
                                    "Crop_stress_temp_exponent", "Soluble_P_leach_KD", "Unknown1",
                                    "Unknown2", "Unknown3", "Irrigation_cost",
                                    "Lime_cost", "Fuel_cost", "Labor_cost", "Unknown4")
  names(watershedParamDataframe) <- watershedParamDataframeNames
  watershedParamDataframe
}

#----
perf2idx <- function (perfMatrix, lowLimit, upLimit) {
  
  idxAccept <- as.data.frame(array(data = NaN,dim = dim(perfMatrix)))
  names(idxAccept) <- names(perfMatrix)
  i = 1
  j = 1
  while (i <= nrow(perfMatrix)) {
    
    if (lowLimit[1] <= perfMatrix['RMSE'][i,1] && perfMatrix['RMSE'][i,1] <= upLimit[1] &&
        lowLimit[2] <= perfMatrix['NASH'][i,1] && perfMatrix['NASH'][i,1] <= upLimit[2] &&
        lowLimit[3] <= perfMatrix['PBIAS'][i,1] && perfMatrix['PBIAS'][i,1] <= upLimit[3] &&
        lowLimit[4] <= perfMatrix['MEAN'][i,1] && perfMatrix['MEAN'][i,1] <= upLimit[4])  {
      
      idxAccept$Index[j] <- i
      idxAccept$RMSE[j] <- perfMatrix['RMSE'][i,1]
      idxAccept$NASH[j] <- perfMatrix['NASH'][i,1]
      idxAccept$PBIAS[j] <- perfMatrix['PBIAS'][i,1]
      idxAccept$MEAN[j] <- perfMatrix['MEAN'][i,1]
      j = j + 1
      i = i + 1
    } else {
      i = i + 1
    }
  }
  
  return(na.omit(idxAccept))
}

#----
sa4APEX <- function (input,input4SA,usingTextFile=FALSE) {
  if (!isTRUE(usingTextFile)) {
    globalInputs <- input
  }
  
  if(isTRUE(usingTextFile) && missing(input4SA) ) {
    globalInputs <- text2InputSimulation(input)
    load(paste(globalInputs$folderPathGsaOutputs,
               "/", globalInputs$gsaType,"_","GSA_WorkSpace",".RData", sep = ""))
  }
  
  gsaType = globalInputs$gsaType
  gsaObject <- input4SA$gsaObject
  
  if (gsaType=="KSTEST") {
    for (j in 1:length(globalInputs$captionVarSim)) {
      perfFunc <- modelPerfDWS(globalInputs$folderPathObserved,
                               globalInputs$labelObservedVar.txt,
                               globalInputs$startDate[j],
                               globalInputs$endDate[j], 1,
                               nrow(gsaObject$X),
                               globalInputs$calculatedOutputFolderDWS,
                               globalInputs$labelOutputVariableDWS,
                               globalInputs$captionVarSim[j], "Date",
                               globalInputs$captionVarObs[j])
      
      
      subFolderPathGsaOutputs <- paste(globalInputs$folderPathGsaOutputs,
                                       "/", globalInputs$captionVarSim[j], "_",
                                       format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),
                                       sep = "")
      
      dir.create(subFolderPathGsaOutputs)
      write.table(perfFunc,
                  paste(subFolderPathGsaOutputs,
                        "/", "Perform.Fcn", "_",
                        gsaType, "_",
                        format(Sys.time(),"%Y_%B_%d_%H_%M"),
                        ".txt", sep = ""),
                  sep = "\t")
      #
      if (gsaObject$ksTestPerf=="NASH") {
        Y_PF <- perfFunc$NASH
        logicSeq <- (Y_PF > gsaObject$ksTestThreshold[j])
      }
      else if (gsaObject$ksTestPerf=="RMSE") {
        Y_PF <- perfFunc$RMSE
        logicSeq <- (Y_PF < gsaObject$ksTestThreshold[j])
      }
      else if (gsaObject$ksTestPerf=="PBIAS") {
        Y_PF <- perfFunc$PBIAS
        logicSeq <- (abs(Y_PF) < abs(gsaObject$ksTestThreshold[j]))
      } else {
        stop("Please enter a valid performance function:, i.e. NASH, RMSE, or PBIAS")
      }
      
      
      idxSeq <- 1:nrow(gsaObject$X)
      idxAccept <- as.numeric(idxSeq[logicSeq])
      xMat = gsaObject$X
      idxParm <- 1:nrow(xMat)
      logicExtract <- idxParm %in% idxAccept
      xPosterior <- xMat[logicExtract, ]
      dir.create(paste(subFolderPathGsaOutputs, "/",
                       gsaObject$ksTestPerf,
                       sep = ""))
      gsa_analysis(gsaObject,
                   xPosterior,
                   paste(subFolderPathGsaOutputs,
                         "/", gsaObject$ksTestPerf,
                         sep = ""),
                   gsaType)
    }
  } else {
    for (j in 1:length(globalInputs$captionVarSim)) {
      perfFunc <- modelPerfDWS(globalInputs$folderPathObserved,
                               globalInputs$labelObservedVar.txt,
                               globalInputs$startDate[j],
                               globalInputs$endDate[j], 1, nrow(gsaObject$X),
                               globalInputs$calculatedOutputFolderDWS,
                               globalInputs$labelOutputVariableDWS,
                               globalInputs$captionVarSim[j], "Date",
                               globalInputs$captionVarObs[j])
      
      subFolderPathGsaOutputs <- paste(globalInputs$folderPathGsaOutputs, "/",
                                       globalInputs$captionVarSim[j], "_",
                                       format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),
                                       sep = "")
      
      dir.create(subFolderPathGsaOutputs)
      write.table(perfFunc,
                  paste(subFolderPathGsaOutputs, "/",
                        "Perform.Fcn", "_", gsaType,"_",
                        format(Sys.time(), "%Y_%B_%d_%H_%M"),
                        ".txt", sep = ""),
                  sep = "\t")
      
      Y_RMSE <- perfFunc$RMSE
      dir.create(paste(subFolderPathGsaOutputs, "/",
                       "RMSE", sep = ""))
      gsa_analysis(gsaObject,
                   Y_RMSE,
                   paste(subFolderPathGsaOutputs,
                         "/", "RMSE", sep = ""),
                   gsaType)
      
      Y_MEAN <- perfFunc$MEAN
      dir.create(paste(subFolderPathGsaOutputs, "/", "MEAN", sep = ""))
      gsa_analysis(gsaObject,
                   Y_MEAN,
                   paste(subFolderPathGsaOutputs,
                         "/", "MEAN", sep = ""),
                   gsaType)
      
      Y_NASH <- perfFunc$NASH
      dir.create(paste(subFolderPathGsaOutputs, "/", "NASH", sep = ""))
      gsa_analysis(gsaObject,
                   Y_NASH,
                   paste(subFolderPathGsaOutputs,
                         "/", "NASH", sep = ""),
                   gsaType)
      
      Y_PBIAS <- perfFunc$PBIAS
      dir.create(paste(subFolderPathGsaOutputs, "/", "PBIAS", sep = ""))
      gsa_analysis(gsaObject,
                   Y_PBIAS,
                   paste(subFolderPathGsaOutputs,
                         "/", "PBIAS", sep = ""),
                   gsaType)
    }
  }
}
#
#
# Helper Functions...
#
gsa_analysis <- function (GSA_Object, Y, folder_path_GSA_Outputs, GSA_Type) {
  if (GSA_Type == "MORRIS") {
    GSA_Object <- sensitivity::tell(GSA_Object, y = Y)
    Mu_Sigma <- print(GSA_Object)
    
    write.table(Mu_Sigma,
                paste(folder_path_GSA_Outputs,
                      "/", GSA_Type, "_", format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),
                      ".txt", sep = ""),
                sep = "\t")
    
    save(GSA_Object,
         file =paste(folder_path_GSA_Outputs, "/", GSA_Type,
                     "_", "Xmat_and_Yvec", format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),
                     ".RData", sep = ""))
    return(GSA_Object)
  }
  else if (GSA_Type == "SRC") {
    GSA_Object <- sensitivity::src(GSA_Object$X, y = Y, rank = FALSE,
                                   nboot = 0, conf = 0.95)
    SRC_Coeff <- as.data.frame(sensitivity:::print.src(GSA_Object))
    colnames(SRC_Coeff) <- c("SRCs")
    write.table(SRC_Coeff,
                paste(folder_path_GSA_Outputs,
                      "/", GSA_Type, "_",
                      format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),
                      ".txt", sep = ""),
                sep = "\t")
    
    save(GSA_Object,
         file = paste(folder_path_GSA_Outputs, "/", GSA_Type,
                      "_", "Final_GSA_Object", format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),
                      ".RData", sep = ""))
    return(GSA_Object)
  }
  else if (GSA_Type == "SRRC") {
    GSA_Object <- sensitivity::src(GSA_Object$X, y = Y, rank = TRUE,
                                   nboot = 0, conf = 0.95)
    SRC_Coeff <- as.data.frame(sensitivity:::print.src(GSA_Object))
    colnames(SRC_Coeff) <- c("SRRCs")
    write.table(SRC_Coeff,
                paste(folder_path_GSA_Outputs,
                      "/", GSA_Type, "_", format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),
                      ".txt", sep = ""),
                sep = "\t")
    
    save(GSA_Object,
         file = paste(folder_path_GSA_Outputs, "/", GSA_Type,
                      "_", "Final_GSA_Object", format(Sys.time(), "%Y_%B_%d_%H_%M_%S")
                      ,".RData", sep = ""))
    return(GSA_Object)
  }
  else if (GSA_Type == "SOBOL") {
    
    GSA_Object$X1 <- scale(x = GSA_Object$X1,center = TRUE,scale = TRUE)
    GSA_Object$X2 <- scale(x = GSA_Object$X2,center = TRUE,scale = TRUE)
    Y = scale(x = Y,center = TRUE,scale = TRUE)
    GSA_Object <- sensitivity::tell(GSA_Object, y = Y)
    Sobol_Idx <- as.data.frame(sensitivity:::print.sobol(GSA_Object))
    colnames(Sobol_Idx) <- c("Sobol's Indices")
    write.table(Sobol_Idx,
                paste(folder_path_GSA_Outputs,
                      "/", GSA_Type, "_",
                      format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),
                      ".txt", sep = ""),
                sep = "\t")
    
    save(GSA_Object,
         file = paste(folder_path_GSA_Outputs, "/", GSA_Type,
                      "_", "Final_GSA_Object",
                      format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),".RData", sep = ""))
    return(GSA_Object)
  }
  else if (GSA_Type == "SOBOL2002") {
    
    GSA_Object$X1 <- scale(x = GSA_Object$X1,center = TRUE,scale = TRUE)
    GSA_Object$X2 <- scale(x = GSA_Object$X2,center = TRUE,scale = TRUE)
    Y = scale(x = Y,center = TRUE,scale = TRUE)
    GSA_Object <- sensitivity::tell(GSA_Object, y = Y)
    Sobol2002_Idx_S <- GSA_Object$S
    Sobol2002_Idx_T <- GSA_Object$T
    Sobol2002_Idx <- cbind(Sobol2002_Idx_S,Sobol2002_Idx_T)
    colnames(Sobol2002_Idx) <- c("Main Effect","Total Effect")
    
    write.table(Sobol2002_Idx,
                paste(folder_path_GSA_Outputs,
                      "/", GSA_Type, "_",
                      format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),
                      ".txt", sep = ""),
                sep = "\t")
    
    save(GSA_Object,
         file = paste(folder_path_GSA_Outputs,
                      "/", GSA_Type,
                      "_", "Final_GSA_Object",
                      format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),".RData",
                      sep = ""))
  }
  else if (GSA_Type == "SOBOL2007") {
    
    GSA_Object$X1 <- scale(x = GSA_Object$X1,center = TRUE,scale = TRUE)
    GSA_Object$X2 <- scale(x = GSA_Object$X2,center = TRUE,scale = TRUE)
    Y = scale(x = Y,center = TRUE,scale = TRUE)
    GSA_Object <- sensitivity::tell(GSA_Object, y = Y)
    SOBOL2007_Idx_S <- GSA_Object$S
    SOBOL2007_Idx_T <- GSA_Object$T
    SOBOL2007_Idx <- cbind(SOBOL2007_Idx_S,SOBOL2007_Idx_T)
    colnames(SOBOL2007_Idx) <- c("Main Effect","Total Effect")
    
    write.table(SOBOL2007_Idx,
                paste(folder_path_GSA_Outputs,
                      "/", GSA_Type, "_",
                      format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),
                      ".txt", sep = ""),
                sep = "\t")
    
    save(GSA_Object,
         file = paste(folder_path_GSA_Outputs,
                      "/", GSA_Type,
                      "_", "Final_GSA_Object",
                      format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),".RData",
                      sep = ""))
  }
  else if (GSA_Type == "SOBOLEFF") {
    tell(GSA_Object, y = Y)
    Soboleff_Idx <- print(GSA_Object)
    write.table(Soboleff_Idx,
                paste(folder_path_GSA_Outputs,
                      "/", GSA_Type, "_",
                      format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),
                      ".txt", sep = ""),
                sep = "\t")
    
    save(GSA_Object,
         file = paste(folder_path_GSA_Outputs,
                      "/", GSA_Type,
                      "_", "Final_GSA_Object",
                      format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),".RData",
                      sep = ""))
  }
  else if (GSA_Type == "SOBOLJANSEN") {
    
    sensitivity::tell(GSA_Object, y = Y)
    SOBOLJANSEN_Idx_S <- GSA_Object$S
    SOBOLJANSEN_Idx_T <- GSA_Object$T
    SOBOLJANSEN_Idx <- cbind(SOBOLJANSEN_Idx_S,SOBOLJANSEN_Idx_T)
    colnames(SOBOLJANSEN_Idx) <- c("Main Effect","Total Effect")
    
    write.table(SOBOLJANSEN_Idx,
                paste(folder_path_GSA_Outputs,
                      "/", GSA_Type, "_",
                      format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),
                      ".txt", sep = ""),
                sep = "\t")
    
    save(GSA_Object,
         file = paste(folder_path_GSA_Outputs,
                      "/", GSA_Type,
                      "_", "Final_GSA_Object",
                      format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),".RData",
                      sep = ""))
  }
  else if (GSA_Type == "SOBOLMARA") {
    
    GSA_Object$X <- scale(x = GSA_Object$X,center = TRUE,scale = TRUE)
    Y = scale(x = Y,center = TRUE,scale = TRUE)
    sensitivity::tell(GSA_Object, y = Y)
    SOBOLMARA_Idx <- GSA_Object$S
    colnames(SOBOLMARA_Idx) <- c("Main Effect")
    
    write.table(SOBOLMARA_Idx,
                paste(folder_path_GSA_Outputs,
                      "/", GSA_Type, "_",
                      format(Sys.time(), "%Y_%B_%d_%H_%M_S%"),
                      ".txt", sep = ""),
                sep = "\t")
    
    save(GSA_Object,
         file = paste(folder_path_GSA_Outputs,
                      "/", GSA_Type,
                      "_", "Final_GSA_Object",
                      format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),".RData",
                      sep = ""))
  }
  else if (GSA_Type == "SOBOLMARTINEZ") {
    
    sensitivity::tell(GSA_Object, y = Y)
    SOBOLMARTINEZ_Idx_S <-  GSA_Object$S['original']
    SOBOLMARTINEZ_Idx_T <- GSA_Object$T['original']
    SOBOLMARTINEZ_Idx <- cbind(SOBOLMARTINEZ_Idx_S,SOBOLMARTINEZ_Idx_T)
    colnames(SOBOLMARTINEZ_Idx) <- c("Main Effect","Total Effect")
    
    write.table(SOBOLMARTINEZ_Idx,
                paste(folder_path_GSA_Outputs,
                      "/", GSA_Type, "_",
                      format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),
                      ".txt", sep = ""),
                sep = "\t")
    
    save(GSA_Object,
         file = paste(folder_path_GSA_Outputs,
                      "/", GSA_Type,
                      "_", "Final_GSA_Object",
                      format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),".RData",
                      sep = ""))
    
  } else if (GSA_Type=="FAST99") {
    
    
    GSA_Object$X <- scale(x = GSA_Object$X,center = TRUE,scale = TRUE)
    Y = scale(x = Y,center = TRUE,scale = TRUE)
    GSA_Object <- sensitivity::tell(GSA_Object, y = Y)
    FAST99_Idx <- as.data.frame(sensitivity:::print.fast99(GSA_Object))
    colnames(FAST99_Idx) <- c("Main Efffect","Total Effect")
    
    write.table(FAST99_Idx,
                paste(folder_path_GSA_Outputs,
                      "/", GSA_Type, "_",
                      format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),
                      ".txt", sep = ""),
                sep = "\t")
    
    save(GSA_Object,
         file = paste(folder_path_GSA_Outputs,
                      "/", GSA_Type,
                      "_", "Final_GSA_Object",
                      format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),".RData",
                      sep = ""))
    
  } else if (GSA_Type=="KSTEST") {
    
    #Generating an empty vector for storing KS-Test results
    KS_Test <- as.data.frame(array(NaN,dim = c(ncol(GSA_Object$X),3)))
    row.names(KS_Test) <- colnames(GSA_Object$X)
    colnames(KS_Test) <- c("Sensitivity_State","p-value","D-statistic")
    i=1
    X_behave <- Y;
    
    write.table(x = X_behave,
                file = paste(folder_path_GSA_Outputs,
                             "/", GSA_Type,
                             "_", "Behaviour",
                             ".txt", sep = ""),
                sep = "\t",row.names = T)
    
    write.table(x = GSA_Object$X,
                file = paste(folder_path_GSA_Outputs,
                             "/", GSA_Type,
                             "_", "Prior",
                             ".txt", sep = ""),
                sep = "\t",row.names = F)
    
    if (nrow(Y)==0) {
      print("No 'behavior' simulation detected!")
      X_Non_Behaviour = GSA_Object$X
      
    } else {
      idx_seq = 1:nrow(GSA_Object$X)
      idx_non_behavior = idx_seq[!(GSA_Object$X[,1] %in% Y[,1])]
      X_Non_Behaviour <- GSA_Object$X[idx_non_behavior,]
    }
    
    write.table(x = X_Non_Behaviour,
                file = paste(folder_path_GSA_Outputs,
                             "/", GSA_Type,
                             "_", "Non-Behaviour",
                             ".txt", sep = ""),
                sep = "\t",row.names = T)
    
    while (i<=ncol(GSA_Object$X)) {
      ks_test_temp = list(p.value=-999)
      try(ks_test_temp <- ks.test(x=(GSA_Object$X[,i][!(GSA_Object$X[,i]%in%Y[,i])]),
                                  y = Y[,i], alternative = "two.sided"),silent = TRUE)
      
      if((ks_test_temp$p.value <= GSA_Object$KS_TEST_sig_level)&&(ks_test_temp$p.value >= 0) ) {
        KS_Test[i,1] <-  "Sensitive"
        KS_Test[i,2] <-  ks_test_temp$p.value
        KS_Test[i,3] <-  as.numeric(ks_test_temp$statistic)
      }
      else if ( (ks_test_temp$p.value >= GSA_Object$KS_TEST_sig_level)&&(ks_test_temp$p.value <= 1)   ) {
        KS_Test[i,1] <-  "Non-sensitive"
        KS_Test[i,2] <-  ks_test_temp$p.value
        KS_Test[i,3] <-  as.numeric(ks_test_temp$statistic)
      }
      else if (ks_test_temp$p.value < 0){
        KS_Test[i,1] <-  "Inconclusive"
        KS_Test[i,2] <-  NaN
        KS_Test[i,3] <-  NaN
      }
      
      i=i+1
    }
    
    write.table(KS_Test, paste(folder_path_GSA_Outputs,
                               "/", GSA_Type, "_",
                               format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),
                               ".txt", sep = ""),
                sep = "\t")
    
    GSA_Object$Y <- Y
    save(GSA_Object,
         file = paste(folder_path_GSA_Outputs,
                      "/", GSA_Type,
                      "_", "Final_GSA_Object",
                      format(Sys.time(), "%Y_%B_%d_%H_%M_%S"),".RData",
                      sep = ""))
  } else {
    stop("Please Enter a Valid GSA Method")
  }
}
#
#
# Helper Functions




modelPerfDWS <- function (folderPathObserved, labelObservedVar.txt,
                          startDate, endDate, idxFirstFile,
                          idxLastFile, calculatedOutputFolderDWS,
                          labelOutputVariableDWS, captionVarSim,
                          captionDateObs = "Date", captionVarObs) {
  
  
  obsTimeseries <- obsDataframeDaily(paste(folderPathObserved,
                                           "/", labelObservedVar.txt,
                                           sep = ""),
                                     startDate, endDate,
                                     captionDateObs = "Date", captionVarObs)
  
  sampleSize <- (idxLastFile - idxFirstFile) + 1
  perfFunc <- data.frame(RMSE = rep(NaN, sampleSize[1]),
                         MEAN = rep(NaN, sampleSize[1]),
                         NASH = rep(NaN, sampleSize[1]),
                         PBIAS = rep(NaN, sampleSize[1]))
  
  for (i in 1:sampleSize[1]) {
    Sim_timeseries <- simDataframeDWS(paste(calculatedOutputFolderDWS,
                                            "/", labelOutputVariableDWS,
                                            toString(i + idxFirstFile -1),
                                            ".SAD", sep = ""),
                                      startDate, endDate,
                                      captionVarSim)
    
    perfFunc$RMSE[i] <- calculateRMSE(Sim_timeseries, obsTimeseries)
    perfFunc$MEAN[i] <- mean(Sim_timeseries$Signal, na.rm = TRUE)
    perfFunc$NASH[i] <- calculateNASH(Sim_timeseries, obsTimeseries)
    perfFunc$PBIAS[i] <- calculatePBIAS(Sim_timeseries, obsTimeseries)
    print(paste("Files Processed...", toString(i), sep = ""))
    
  }
  
  return(perfFunc)
}
#----
# Helper functions
obsDataframeDaily <- function (obsFile.txt, startDate, endDate,
                               captionDateObs = "Date", captionVarObs) {
  # Creates a dataframe with columns "Date" and "Signal" from a .DWS file.
  #
  #Args:
  #  SITE271.DWS: name of the .DWS file.
  #  startDate: start date.
  #  endDate: end date.
  #  captionVarSim: Variable caption as in the .DWS file.
  #Returns: Simulated series.
  #
  startDate <- gsub(" ", "", startDate)
  endDate <- gsub(" ", "", endDate)
  obsTable <- read.table(obsFile.txt, header = TRUE)
  namesVec <- names(obsTable)
  idxCaptionDateObs <- match(captionDateObs, namesVec)
  idxCaptionVarObs <- match(captionVarObs, namesVec)
  dateString <- as.character(obsTable[[idxCaptionDateObs]])
  signalVector <- obsTable[[idxCaptionVarObs]]
  startIdx <- match(startDate, dateString)
  endIdx <- match(endDate, dateString)
  obsTimeSeriesLength <- (endIdx - startIdx) + 1
  obsTimeSeries <- data.frame(Date = rep(1, times = obsTimeSeriesLength),
                              Signal = rep(1, times = obsTimeSeriesLength))
  obsTimeSeries$Date <- dateString[startIdx:endIdx]
  obsTimeSeries$Signal <- signalVector[startIdx:endIdx]
  obsTimeSeries
}
#----
#
simDataframeDWS <- function(SITE271.SAD, startDate, endDate, captionVarSim) {
  # Creates a dataframe with columns "Date" and "Signal" from a .DWS file.
  #
  #Args:
  #  SITE271.DWS: name of the .DWS file.
  #  startDate: start date.
  #  endDate: end date.
  #  captionVarSim: Variable caption as in the .DWS file.
  #Returns: Simulated series.
  #
  startDate <- gsub(" ", "", startDate)
  endDate <- gsub(" ", "", endDate)
  simTable <- read.table(SITE271.SAD, skip = 9,comment.char = "",header = T,fill=T)
  simTable <- simTable[simTable$ID==1,]            ####select subbasin
  simTable <- simTable[simTable$CPNM=="LCSM",]     ####select crop varieties
  simTable <- simTable[,c("Y","M","D","WYLD","ET")]
  simSize <- dim(simTable)
  dateString <- rep(1, times = simSize[1])
  
  for (i in 1:simSize[1]) {
    yearString <- toString(simTable$Y[i])
    if (simTable$M[i] < 10) {
      monthString <- paste("0", toString(simTable$M[i]), sep = "")
    } else {
      monthString <- toString(simTable$M[i])
    }
    
    if (simTable$D[i] < 10) {
      dayString <- paste("0", toString(simTable$D[i]),sep = "")
      
    } else {
      dayString <- toString(simTable$D[i])
    }
    
    dateString[i] <- paste(yearString, monthString, dayString, sep = "")
  }
  
  startIdx <- match(startDate, dateString)
  endIdx <- match(endDate, dateString)
  simTimeSeriesLength <- (endIdx - startIdx) + 1
  simTimeSeries <- data.frame(Date = rep(1, times = simTimeSeriesLength),
                              Signal = rep(1, times = simTimeSeriesLength))
  
  if (captionVarSim == "WYLD") {
    simTimeSeries$Signal <- simTable$WYLD[startIdx:endIdx]
  }
  else if (captionVarSim == "RFV") {
    simTimeSeries$Signal <- simTable$RFV[startIdx:endIdx]
  }
  else if (captionVarSim == "QDR") {
    simTimeSeries$Signal <- simTable$QDR[startIdx:endIdx]
  }
  else if (captionVarSim == "ET") {
    simTimeSeries$Signal <- simTable$ET[startIdx:endIdx]
  }
  else if (captionVarSim == "Q") {
    simTimeSeries$Signal <- simTable$Q[startIdx:endIdx]
  }
  else if (captionVarSim == "YN") {
    simTimeSeries$Signal <- simTable$YN[startIdx:endIdx]
  }
  else if (captionVarSim == "YP") {
    simTimeSeries$Signal <- simTable$YP[startIdx:endIdx]
  }
  else if (captionVarSim == "QN") {
    simTimeSeries$Signal <- simTable$QN[startIdx:endIdx]
  }
  else if (captionVarSim == "QP") {
    simTimeSeries$Signal <- simTable$QP[startIdx:endIdx]
  }
  else if (captionVarSim == "PRKP") {
    simTimeSeries$Signal <- simTable$PRKP[startIdx:endIdx]
  }
  else if (captionVarSim == "DPRK") {
    simTimeSeries$Signal <- simTable$DPRK[startIdx:endIdx]
  }
  else if (captionVarSim == "QDRN") {
    simTimeSeries$Signal <- simTable$QDRN[startIdx:endIdx]
  }
  else if (captionVarSim == "QDRP") {
    simTimeSeries$Signal <- simTable$QDRP[startIdx:endIdx]
  }
  else if (captionVarSim == "DN") {
    simTimeSeries$Signal <- simTable$DN[startIdx:endIdx]
  }
  else if (captionVarSim == "RSFN") {
    simTimeSeries$Signal <- simTable$RSFN[startIdx:endIdx]
  }
  else if (captionVarSim == "QRFN") {
    simTimeSeries$Signal <- simTable$QRFN[startIdx:endIdx]
  }
  else if (captionVarSim == "QRFP") {
    simTimeSeries$Signal <- simTable$QRFP[startIdx:endIdx]
  }
  else if (captionVarSim == "SSFN") {
    simTimeSeries$Signal <- simTable$SSFN[startIdx:endIdx]
  }
  else if (captionVarSim == "SSF") {
    simTimeSeries$Signal <- simTable$SSF[startIdx:endIdx]
  }
  else if (captionVarSim == "QRF") {
    simTimeSeries$Signal <- simTable$QRF[startIdx:endIdx]
  }
  else if (captionVarSim == "RSSF") {
    simTimeSeries$Signal <- simTable$RSSF[startIdx:endIdx]
  }
  else if (captionVarSim == "TN") {
    simTimeSeries$Signal <- simTable$QN[startIdx:endIdx] +
      simTable$YN[startIdx:endIdx] +
      simTable$QDRN[startIdx:endIdx] +
      simTable$RSFN[startIdx:endIdx] +
      simTable$QRFN[startIdx:endIdx] +
      simTable$SSFN[startIdx:endIdx]
  }
  else if (captionVarSim == "TP") {
    simTimeSeries$Signal <- simTable$QP[startIdx:endIdx] +
      simTable$YP[startIdx:endIdx] +
      simTable$QDRP[startIdx:endIdx] +
      simTable$QRFP[startIdx:endIdx]
  }
  else {
    stop("Please enter a valid variable name, see variables in .SAD files")
  }
  simTimeSeries$Date <- dateString[startIdx:endIdx]
  simTimeSeries
}
#----
#
calculateRMSE <- function (simSeries, obsSeries) {
  # Calculates Root Mean Squar Error (RMSE).
  #
  #Args:
  #  simSeries: Simulated data (i.e., output of a call to simDataframeDWS function).
  #  obsSeries: Observed data (i.e., output of a call to obsDataframeDaily function ).
  #
  #Returns: RMSE
  #
  SE_series <- (simSeries$Signal - obsSeries$Signal)^2
  RMSE <- sqrt(mean(SE_series, na.rm = TRUE))
  return(RMSE)
}
#----
#
calculatePBIAS <- function (simSeries, obsSeries) {
  # Calculates PBIAS performance measure.
  #
  #Args:
  #  simSeries: Simulated data (i.e., output of a call to simDataframeDWS function).
  #  obsSeries: Observed data (i.e., output of a call to obsDataframeDaily function ).
  #
  #Returns: PBIAS
  #
  NS_Numerator <- sum(obsSeries$Signal - simSeries$Signal, na.rm = TRUE) * 100
  NS_Denominator <- sum(obsSeries$Signal, na.rm = TRUE)
  PBIAS <- NS_Numerator/NS_Denominator
  return(PBIAS)
}
#----
#
calculateNASH <- function (simSeries, obsSeries) {
  # Calculates Nash-Sutcliffe (NS) performance measure.
  #
  #Args:
  #  simSeries: Simulated data (i.e., output of a call to simDataframeDWS function).
  #  obsSeries: Observed data (i.e., output of a call to obsDataframeDaily function ).
  #
  #Returns: NS
  #
  NS_Numerator <- sum((simSeries$Signal - obsSeries$Signal)^2,
                      na.rm = TRUE)
  NS_Denominator <- sum((obsSeries$Signal - mean(obsSeries$Signal,
                                                 na.rm = TRUE))^2, na.rm = TRUE)
  NASH <- 1 - (NS_Numerator/NS_Denominator)
  return(NASH)
}

#----
writeBasinParameters <- function (uncertainWatershedParam,paramFile.dat, backUpParamFile.dat) {
  #
  #
  #Loading required functions...
  #source("watershed_param_dataframe.r")
  
  #Setting up connections...
  connSource <- file(description = backUpParamFile.dat,open = "r+",blocking = TRUE)
  connTarget <- file(description = paramFile.dat,open ="r+",blocking = TRUE )
  
  #Reading from back up file into a readLines object...
  readlinesObj <- readLines(con = connSource,n = -1)
  
  
  ####Sub-function
  parm2loc <- function(PARM) {
    switch(names(PARM),
           "Crop_canopy_PET"=return(c(31,1,8,"%8.3f")),
           "Root_growth_soil"=return(c(31,9,16,"%8.3f")),
           "Water_stress_harvest"=return(c(31,17,24,"%8.3f")),
           "Water_storage_N"=return(c(31,25,32,"%8.3f")),
           "Soil_water_limit"=return(c(31,33,40,"%8.3f")),
           "Winter_dormancy"=return(c(31,41,48,"%8.3f")),
           "N_fixation"=return(c(31,49,56,"%8.3f")),
           "Soluble_P_runoff"=return(c(31,57,64,"%8.3f")),
           "Pest_damage_moisture"=return(c(31,65,72,"%8.3f")),
           "Pest_damage_cover"=return(c(31,73,80,"%8.3f")),
           
           #Line 32....
           "Moisture_req_seed_germ"=return(c(32,1,8,"%8.3f")),
           "Soil_evap_coeff"=return(c(32,9,16,"%8.3f")),
           "Wind_erod_coeff"=return(c(32,17,24,"%8.3f")),
           "Nitrate_leac_ratio"=return(c(32,25,32,"%8.3f")),
           "Runoff_CN_Adj_parm"=return(c(32,33,40,"%8.3f")),
           "Expand_CN_ret_parm"=return(c(32,41,48,"%8.3f")),
           "Soil_evap_plant_cover"=return(c(32,49,56,"%8.3f")),
           "Sedim_rout_exponent"=return(c(32,57,64,"%8.3f")),
           "Sedim_rout_coeff"=return(c(32,65,72,"%8.3f")),
           "Runoff_CN_int_abs"=return(c(32,73,80,"%8.3f")),
           #Line 33...
           "Soluble_C_adsorp_coeff"= return(c(33,1,8,"%8.3f")),
           "CN_retention_frozen_soil"=return(c(33,9,16,"%8.3f")),
           "Harg_equation_parm"=return(c(33,17,24,"%8.3f")),
           "Pest_leach_ratio"=return(c(33,25,32,"%8.3f")),
           "Expo_coeff_rainfall"=return(c(33,33,40,"%8.3f")),
           "Matur_frac_spring"=return(c(33,41,48,"%8.3f")),
           "CEC_effect_nitrification"=return(c(33,49,56,"%8.3f")),
           "N_fixation_limit"=return(c(33,57,64,"%8.3f")),
           "Biological_mix_efficiency"=return(c(33,65,72,"%8.3f")),
           "Soluble_P_exponent"=return(c(33,73,80,"%8.3f")),
           #Line 34...
           "Max_depth_bio_mixing"= return(c(34,1,8,"%8.3f")),
           "OrgP_loss_exponent"=return(c(34,9,16,"%8.3f")),
           "MUST_coeff"=return(c(34,17,24,"%8.3f")),
           "Harg_PET_exponent"=return(c(34,25,32,"%8.3f")),
           "Denit_soil_threshold"=return(c(34,33,40,"%8.3f")),
           "Daily_denit_limit"=return(c(34,41,48,"%8.3f")),
           "SWAT_delivery_ratio_exponent"=return(c(34,49,56,"%8.3f")),
           "Water_stress_coeff"=return(c(34,57,64,"%8.3f")),
           "Puddling_sat_conduct"=return(c(34,65,72,"%8.3f")),
           "Groundwater_stor_threshold"=return(c(34,73,80,"%8.3f")),
           #Line 35...
           "Root_temp_stress_exponent" = return(c(35,1,8,"%8.3f")),
           "SCS_index_coeff"=return(c(35,9,16,"%8.3f")),
           "Plow_depth"=return(c(35,17,24,"%8.3f")),
           "CN_retention_param"=return(c(35,25,32,"%8.3f")),
           "sediment_rout_travel_coeff"=return(c(35,33,40,"%8.3f")),
           "RUSLE_c_factor_res"=return(c(35,41,48,"%8.3f")),
           "RUSLE_c_factor_height"=return(c(35,49,56,"%8.3f")),
           "Climate_stress_factor"=return(c(35,57,64,"%8.3f")),
           "Max_rain_intercept"=return(c(35,65,72,"%8.3f")),
           "Rain_intercept_coeff"=return(c(35,73,80,"%8.3f")),
           #Line 36...
           "Water_stor_residue_coeff"=return(c(36,1,8,"%8.3f")),
           "Tillage_residue_decay_rate_coeff"=return(c(36,9,16,"%8.3f")),
           "Microbial_soil_depth_coeff"=return(c(36,17,24,"%8.3f")),
           "N_enrich_coeff"=return(c(36,25,32,"%8.3f")),
           "N_enrich_rout_exponent"=return(c(36,33,40,"%8.3f")),
           "Fraction_destroyed_burn"=return(c(36,41,48,"%8.3f")),
           "P_enrich_rout_coeff"=return(c(36,49,56,"%8.3f")),
           "P_enrich_rout_exponent"=return(c(36,57,64,"%8.3f")),
           "P_move_evap_coeff"=return(c(36,65,72,"%8.3f")),
           "Max_days_grazed_rotation"=return(c(36,73,80,"%8.3f")),
           #Line 37...
           "Soil_water_up_flow_limit"=return(c(37,1,8,"%8.3f")),
           "Manure_erosion_equation_coeff"=return(c(37,9,16,"%8.3f")),
           "N_enrich_ratio_delivery"=return(c(37,17,24,"%8.3f")),
           "Dust_distribution_coeff"=return(c(37,25,32,"%8.3f")),
           "RUSLE2_trans_capacity"=return(c(37,33,40,"%8.3f")),
           "RUSLE2_trans_capacity_threshold"=return(c(37,41,48,"%8.3f")),
           "Dust_distribution_exponent"=return(c(37,49,56,"%8.3f")),
           "Manure_erosion_exponent"=return(c(37,57,64,"%8.3f")),
           "Microbial_top_soil_coeff"=return(c(37,65,72,"%8.3f")),
           'Microbial_decay_coeff'=return(c(37,73,80,"%8.3f")),
           #Line 38...
           "Manure_erosion_coeff"=return(c(38,1,8,"%8.3f")),
           "Volt_nitrification_partition_coeff"=return(c(38,9,16,"%8.3f")),
           "Hydrograph_dev_param"=return(c(38,17,24,"%8.3f")),
           "Partition_N_flow_groundwater"=return(c(38,25,32,"%8.3f")),
           "P_enrich_ratio_deliver_SWAT"=return(c(38,33,40,"%8.3f")),
           "Stand_dead_fall_rate_coeff"=return(c(38,41,48,"%8.3f")),
           "Runoff_2_delay_pest"=return(c(38,49,56,"%8.3f")),
           "Soil_water_2_delay_tillage"=return(c(38,57,64,"%8.3f")),
           "Auto_mov_lower_limit"=return(c(38,65,72,"%8.3f")),
           "Nitrification_vol_upper_limit"=return(c(38,73,80,"%8.3f")),
           #Line 39...
           "Tech_coeff"=return(c(39,1,8,"%8.3f")),
           
           "Drainage_lateral_conduct"=return(c(39,17,24,"%8.3f")),
           "P_flux_labile_active_coeff"=return(c(39,25,32,"%8.3f")),
           "P_flux_active_stable_coeff"=return(c(39,33,40,"%8.3f")),
           "N_salt_evap_coeff"=return(c(39,41,48,"%8.3f")),
           "Water_table_recession_coeff"=return(c(39,49,56,"%8.3f")),
           "Water_table_move_limit"=return(c(39,57,64,"%8.3f")),
           "Water_table_recession_exponent"=return(c(39,65,72,"%8.3f")),
           "Subsurface_flow_factor"=return(c(39,73,80,"%8.3f")),
           #Line 40...
           "Flood_evap_limit"=return(c(40,1,8,"%8.3f")),
           "Runoff_adj_link"=return(c(40,9,16,"%8.3f")),
           "Water_erosion_threshold"=return(c(40,17,24,"%8.3f")),
           "Wind_erosion_threshold"=return(c(40,25,32,"%8.3f")),
           "Crop_stress_temp_exponent"=return(c(40,33,40,"%8.3f")),
           "Soluble_P_leach_KD"=return(c(40,41,48,"%8.3f")),
           
           "Unknown1"=return(c(40,49,56,"%8.3f")),
           "Unknown2"=return(c(40,57,64,"%8.3f")),
           "Unknown3"=return(c(40,65,72,"%8.3f")),
           #Line 41...
           "Irrigation_cost"=return(c(41,1,8,"%8.3f")),
           "Lime_cost"=return(c(41,9,16,"%8.3f")),
           "Fuel_cost"=return(c(41,17,24,"%8.3f")),
           "Labor_cost"=return(c(41,25,32,"%8.3f")),
           "Unknown4"=return(c(41,33,40,"%8.3f"))
    )
    
    
  }
  ####Sub-Function....
  parm2replace <- function(readlinesObj,PARM) {
    loc_vec <- parm2loc(PARM)
    substr(x = readlinesObj[strtoi(loc_vec[1])],start = strtoi(loc_vec[2]),
           stop = strtoi(loc_vec[3])) <- sprintf(fmt = loc_vec[4],PARM[[1]])
    return(readlinesObj)
    
  }
  ################
  i=1
  col_numbers <-ncol(uncertainWatershedParam)
  while(i<=col_numbers) {
    readlinesObj <- parm2replace(readlinesObj,uncertainWatershedParam[i])
    i=i+1
    
  }
  ###Writing the final readLine object...
  writeLines(text = readlinesObj,con = connTarget)
  
  close(connSource)
  close(connTarget)
  
}

#----
writeControlParameters <-function(uncertainContParam, APEXCONT.dat, backUpAPEXCONT.dat) {
  #
  #Loading required functions...
  #source("watershed_cont_dataframe.r")
  
  #Setting up connections...
  connSource <- file(description = backUpAPEXCONT.dat,open = "r+",blocking = TRUE)
  connTarget <- file(description = APEXCONT.dat,open ="r+",blocking = TRUE )
  
  #Reading from back up file into a readLines object...
  readlinesObj <- readLines(con = connSource,n = -1)
  #
  ###Sub-functionn
  cont2loc <- function(APEXCONT) {
    switch(names(APEXCONT),
           
           #Line 3...
           "RFN"= return(c(3,1,8,"%8.2f")),
           "CO2"= return(c(3,9,16,"%8.2f")),
           "CQN"= return(c(3,17,24,"%8.2f")),
           "PSTX"= return(c(3,25,32,"%8.2f")),
           "YWI"= return(c(3,33,40,"%8.2f")),
           "BTA"= return(c(3,41,48,"%8.2f")),
           "EXPK"= return(c(3,49,56,"%8.2f")),
           "QG"= return(c(3,57,64,"%8.2f")),
           "QCF"= return(c(3,65,72,"%8.2f")),
           "CHSO"= return(c(3,73,80,"%8.2f")),
           
           #Line 4...
           "BWD"=return(c(4,1,8,"%8.2f")),
           "FCW"=return(c(4,9,16,"%8.2f")),
           "FPSC"=return(c(4,17,24,"%8.2f")),
           "GWSO"=return(c(4,25,32,"%8.2f")),
           "RFTO"=return(c(4,33,40,"%8.2f")),
           "RFPO"=return(c(4,41,48,"%8.2f")),
           "SATO"=return(c(4,49,56,"%8.2f")),
           "FL"=return(c(4,57,64,"%8.2f")),
           "FW"=return(c(4,65,72,"%8.2f")),
           "ANG"=return(c(4,73,80,"%8.2f")),
           
           #Line 5...
           "UXP"= return(c(5,1,8,"%8.2f")),
           "DIAM"= return(c(5,9,16,"%8.2f")),
           "ACW"= return(c(5,17,24,"%8.2f")),
           "GZL0"= return(c(5,25,32,"%8.2f")),
           "RTN0"=return(c(5,33,40,"%8.2f")),
           "BXCT"= return(c(5,41,48,"%8.2f")),
           "BYCT"=return(c(5,49,56,"%8.2f")),
           "DTHY"=return(c(5,57,64,"%8.2f")),
           "QTH"=return(c(5,65,72,"%8.2f")),
           "STND"=return(c(5,73,80,"%8.2f")),
           
           #Line 6...
           "DRV"=return(c(6,1,8,"%8.2f")),
           "PCO0"=return(c(6,9,16,"%8.2f")),
           "RCC0"=return(c(6,17,24,"%8.2f")),
           "CSLT"=return(c(6,25,32,"%8.2f")),
           "BUS1"=return(c(6,33,40,"%8.2f")),
           "BUS2"=return(c(6,41,48,"%8.2f")),
           "BUS3"=return(c(6,49,56,"%8.2f"))
           
    )
  }
  
  ####Sub-Function....
  cont2replace <- function(readlinesObj,APEXCONT) {
    locVec <- cont2loc(APEXCONT)
    substr(x = readlinesObj[strtoi(locVec[1])],start = strtoi(locVec[2]),
           stop = strtoi(locVec[3])) <- sprintf(fmt = locVec[4],APEXCONT[[1]])
    return(readlinesObj)
    
  }
  ################
  i=1
  colNumbers <-ncol(uncertainContParam)
  while(i<=colNumbers) {
    readlinesObj <- cont2replace(readlinesObj,uncertainContParam[i])
    i=i+1
    
  }
  ###Writing the final readLine object...
  writeLines(text = readlinesObj,con = connTarget)
  
  close(connSource)
  close(connTarget)
}




# Creating a copy of tutorial folder inside the working directory

getExampleFolder()

# 1) Generating a list object with a predefined structure compatible to APEXSENSUN:
globalInput <- APEXSENSUN::inputGen()

# 2) Setting the required inputs (e.g. uncertainty boubds, SA method, sample size, ...)
#
# Setting uncertainty bounds:
globalInput$apexPARM$Root_growth_soil[1] = 0.15
globalInput$apexPARM$Root_growth_soil[2] = 0.2

globalInput$apexPARM$Soil_water_limit[1] = 0
globalInput$apexPARM$Soil_water_limit[2] = 1

globalInput$apexPARM$Soil_evap_coeff[1] = 1.5
globalInput$apexPARM$Soil_evap_coeff[2] = 2.5

globalInput$apexPARM$Soil_evap_plant_cover[1] = 0
globalInput$apexPARM$Soil_evap_plant_cover[2] = 0.5

globalInput$apexPARM$Runoff_CN_int_abs[1] = 0.05
globalInput$apexPARM$Runoff_CN_int_abs[2] = 0.4

globalInput$apexPARM$Max_rain_intercept[1] = 0
globalInput$apexPARM$Max_rain_intercept[2] = 15

globalInput$apexPARM$Rain_intercept_coeff[1] = 0.05
globalInput$apexPARM$Rain_intercept_coeff[2] = 0.3

globalInput$apexPARM$Microbial_top_soil_coeff[1] = 0.1
globalInput$apexPARM$Microbial_top_soil_coeff[2] = 1

globalInput$apexPARM$Microbial_decay_coeff[1] = 0.5
globalInput$apexPARM$Microbial_decay_coeff[2] = 1.5

globalInput$apexPARM$Biological_mix_efficiency[1] = 0.1
globalInput$apexPARM$Biological_mix_efficiency[2] = 0.5

globalInput$apexPARM$Max_depth_bio_mixing[1] = 0.1
globalInput$apexPARM$Max_depth_bio_mixing[2] = 0.3


# SA method and sample size:
#"SRC"、"SRRC", "SOBOL","SOBOL2002","SOBOL2007","SOBOLEFF"、
#"SOBOLJANSEN","SOBOLMARA","SOBOLMARTINEZ","FAST99","KSTEST"
globalInput$gsaType <- "SRRC"
globalInput$sampleSize <- 1000

# 3) Performing Monte Carlo simulation using the setting in globalInput:
input4SA <- mc4APEX(globalInput)

# 4) Calculation of sensitivity indices:
sa4APEX(globalInput,input4SA = input4SA)

# Calculation of performance matrix (containing RMSE, NASH, PBIAS, MEAN) for different Monte Carlo runs:
perfMat <- APEXSENSUN::dws2perf(observedFilePath ="Example/Observed_Files/observed.txt",
                                dwsFolderPath = "Example/Calculated_Outputs/DWS",
                                startDate="2002-01-01", endDate="2003-01-01",
                                captionDwsVar="ET", captionObsVar="ET",TW="week")
#perfMat
# Detecting simulation numbers meeting performance criteria:
acceptedSimulations <- perf2idx(perfMatrix = perfMat,
                                lowLimit = c(0, -25, 0, 0),
                                upLimit = c(10, 100, 25, 5))  
acceptedSimulations
