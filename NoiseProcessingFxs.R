library(gsignal)
library(av)
library(tuneR)
library(hdf5r)
library(rhdf5)
library(ggplot2)
library(lubridate)
library(viridis)
library(scales)


getICESprms<-function(prms){
  
  # Input parameters and flesh out whatever is missing for the file information
  # Metadata group
  # Data group
  
  
  # File group values
  FileGroupVars = c('Email', 'CreationDate', 'StartDate', 'EndDate', 'Institution',
                    'Contact', 'CountryCode', 'StationCode')
  
  
  # Preallocate
  FileGroup <- data.frame(matrix(ncol = length(FileGroupVars), nrow = 1))
  colnames(FileGroup) <- FileGroupVars
  
  # Text for each descriptor
  descriptors <- c(
    Email = "Creator of the HDF5 file/ who holds responsibility for data QA and creation of the submited hdf5 file.:",
    CreationDate = "Date of file creation. UTC DateTime in ISO 8601 format: YYYY-MM-DDThh:mm[:ss] or YYYY-MM-DD hh:mm[:ss] . Seconds are optional.",
    StartDate = "Measurement collection start date. UTC DateTime in ISO 8601 format: YYYY-MM-DDThh:mm[:ss] or YYYY-MM-DD hh:mm[:ss]. Seconds are optional.",
    EndDate = "Measurement collection end date. UTC DateTime in ISO 8601 format: YYYY-MM-DDThh:mm[:ss] or YYYY-MM-DD hh:mm[:ss]. Seconds are optional",
    Institution ='Institution which acquired the data',
    Contact = 'Contact of all future external queries/who submits/holds responsibility for submission ',
    CountryCode = 'Country Code. No context given',
    StationCode ='The station code and its associated coordinates can be found in the ICES station dictionary' 
  )
  
  for (var in FileGroupVars) {
    if (var %in% colnames(prms)) {
      FileGroup[[var]] <- prms[[var]]
    } else {
      value <- readline(descriptors[var])
      FileGroup[[var]] <- value
    }
  }
  
  
  ## Setup the metadata group descriptors
  
  # Text for each descriptor
  descriptors <- c(
    HydrophoneType = 'This field describes the manufacturer and the used hydrophone type/model e.g. Brüell&Kjaer 8106. This field needs to be an array if there are multiple channels (one per channel).',
    HydrophoneSerialNumber = 'e.g. "SN#1234". This field needs to be an array if there are multiple channels (one per channel).',
    RecorderSerialNumber = 'Recorder/data logger type e.g. "Soundtrap"',
    RecorderSerialNumber = 'Recorder serial number e.g. "SN#2345"',
    MeasurementHeight ='Height above the seafloor, in meters',
    MeasuremetnPurpose ='Description of why the continuous underwater noise measurements reported were monitored',
    MeasurementSetup ='Description of deployment. Mandatory in case the purpose is "HELCOM monitoring"',
    RigDesign ='Description of deployment construction. Mandatory in case the purpose is "HELCOM monitoring"',
    FrequencyCount = 'Number of frequency bands',
    FrequencyIndex = 'Third octave band nominal center frequencies. This field is an array of frequencies, with as many columns as the number of frequency bands reported under FrequencyCount',
    FrequencyUnit ='String (10), presume Hz',
    ChannelCount='Number of channels used',
    MeasurementTotalNo = 'Number of measurements',
    MeasurementUnit= 'Unit in which the values are in e.g. dB re 1μPa',
    AveragingTime= 'Averaging time in seconds',
    ProcessingAlgorithm= 'Algorithm used to process the data e.g. computation method for third octave band (fft, filter bank ...)- analysis',
    DataUUID= 'Unique identification number, linking the data submission to the corresponding raw data. It should be used for resubmissions of the same data; matlab function available: uuid = char(java.util.UUID.randomUUID);',
    DatasetVersion= 'Indicates version of the submitted dataset. It should be changed upon resubmission',
    CalibrationProcedure= 'Method used to check the measuring chain. e.g. point calibration with pistonphone, functionality test with microphone and loudspeaker (frequency dependent), or other method used to check the measuring chain. e.g. point calibration with pistonphone, functionality test with microphone and loudspeaker (frequency dependent), or other. Mandatory in case the purpose is "HELCOM monitoring"',
    CalibrationDateTime= 'Date of when the system was last calibrated. Mandatory in case "CalibrationProcedure" is specified UTC DateTime in ISO 8601 format: YYYY-MM-DDThh:mm[:ss] or YYYY-MM-DD hh:mm[:ss]. Seconds are optional.',
    Comments =''
  )
  
  MetaGroupVars = names(descriptors)
  # Preallocate
  MetaGroup <- data.frame(matrix(ncol = length(MetaGroupVars), nrow = 1))
  colnames(MetaGroup) <- MetaGroupVars
  
  
  for (var in MetaGroupVars) {
    if (var %in% colnames(prms)) {
      MetaGroup[[var]] <- prms[[var]]
    } else {
      value <- readline(descriptors[var])
      MetaGroup[[var]] <- value
    }
  }
  
  
  
}

createICESmeta<-function(){
  
}

createAudioDataframe <- function(fileLoc,
                                 nameStringPattern ='\\d{12}', 
                                 lubradiateFormat ="ymdHMS"){
  
  # create df with start and end times of all files
  files <- (list.files(fileLoc, pattern = "\\.wav$"))
  
  audioData =data.frame(files = file.path(fileLoc, files))
  audioData$FileName = files
  audioData$TimeString=lapply(audioData$FileName, str_extract, 
                              pattern =nameStringPattern)
  
  audioData$StartTime =parse_date_time(audioData$TimeString,
                                       lubradiateFormat,
                                       tz='UTC')
  audioData$Duration = as.numeric(sapply(audioData$files, 
                                         av_media_info)["duration", ]) 
  audioData$EndTime = audioData$StartTime+audioData$Duration
  
  
  # Pull sample rate from first file, assume consistnat (required)
  prms = av_media_info(audioData$files[1])$audio
  prms$duration = av_media_info(audioData$files[1])$duration
  colnames(prms)[colnames(prms)=='sample_rate']='Fs' # field standard
  
  dataOut = list(audioData=audioData, prms=prms)
  
  return(dataOut)
}



createWindow<-function(winname, prms){
  
  # get the window functions
  # Input
  # winname - character string either 'none', 'hann', 'Hamming', or 'Blackman'
  # prms - Analysis prameters
  # Returns
  # dataOut - list of the window and the alpha value
  
  
  if (winname== 'none'){
    w = array(1,prms$N);
    alpha = 1                #scaling factor
  }
  
  if (winname == 'Hann'){
    w = (0.5 - 0.5*cos(2*pi*(1:prms$N)/prms$N))
    alpha = 0.5;                #scaling factor
  }
  
  if (winname == 'Hamming'){
    w = (0.54 - 0.46*cos(2*pi*(1:prms$N)/prms$N))
    alpha = 0.55;                #scaling factor
  }
  
  if (winname == 'Blackman'){
    w = (0.42 - 0.5*cos(2*pi*(1:prms$N)/prms$N) + 
           0.08*cos(4*pi*(1:prms$N)/prms$N))
    alpha = 0.42;                #scaling factor
  }
  
  dataOut = list(w, alpha)
  return(dataOut)
}


writeToH5datarH5df<-function(ProjName, instrumentName,
                             dataType='hybridMiliDec', newData, dataStart=1,
                             maxRows=NULL, storagemMode = "double",
                             fillValue = NaN){
  
  # Write the data to a database using the rhdf5 package
  # https://rdrr.io/bioc/rhdf5/man/rhdf5.html
  
  # Input 
  # ProjectName - Initialized H5 files
  # instrumentName - Initialized group for instrument (e.g. soundtrap)
  # dataType - character string of the type of noise level data (e.g. broadband)
  # newData - matrix or vector of the new data to write
  # dataStart-  where in the dataset to start writing the new data
  # maxRows - if initilizing the dataset, total number of expected rows in the 
  # eventual dataset. This can't be updated later. Empty cells filled with 'NULL'
  # Type of storage, e.g. integer, double, character etc.
  # Default fill value- what to populate the dataset with
  
  
  ################################################
  # array or matrix, define datasize
  ############################################
  if(is.null(dim(newData)) || dim(newData)[2] ==1){
    
    # total and chunk dimensiomns
    totalDims=maxRows
    chunkDims = min(10000, length(newData))
    startDims =  dataStart
    countStart = length(newData)
    
  }else{
    
    # total and chunk dimensiomns
    totalDims = c(maxRows, dim(newData)[2])
    chunkDims = dim(newData)
    startDims =  c(dataStart, 1)
    countStart = c(nrow(newData), ncol(newData))
  }
  
  
  ###############################
  # Dataset availability
  #############################
  
  dataDir =  h5dump(ProjName)
  Instrument = dataDir[[instrumentName]]
  

  
  ###############################################
  # determine if dataset exists, if not add it
  ###############################################
  
  if(!dataType %in% names(Instrument)){
    
    # create a dataset using h5df

    # Special case for character values (grrrrr!!!)
    if(storagemMode=='character'){
      
  
      h5write(rep("0000-00-00 00:00:00", totalDims[1]),
           file = ProjName,
           name = paste(instrumentName,dataType,sep="/"))
      
    }else{
      h5createDataset(file = ProjName,
                      dataset = paste(instrumentName,dataType,sep="/"),
                      dims = totalDims,
                      chunk =chunkDims,
                      storage.mode = storagemMode,
                      fillValue= fillValue)
    }
      
  }
  

  # Dataset has been created (or indexed)-populate with new information
  h5write(
    newData,
    file = ProjName,
    name = paste(instrumentName,dataType,sep="/"),
    start = startDims,
    count =countStart)
  
  
}


# Create new or write h5df data
writeToH5data<-function(H5group, dataType='hybridMiliDec',
                        newData, newCols=NULL, newRows=NULL){
  
  # Write the data to a database using the hdf5r package (nicer but ungodly slow)
  
  # Input 
  # H5group - Initialized group for instrument (e.g. soundtrap) within a H5 file
  # dataType - character string of the type of noise level data (e.g. broadband)
  # newData - matrix or vector of the new data to write
  # dataStart-  where in the dataset to start writing the new data
  # maxRows - if initilizing the dataset, total number of expected rows in the 
  # newCols - column names if first run
  
  # if the data isn't in the instrument, add the metric
  if(!dataType %in% H5group$ls()$name){
    
    # not working -trying to intialize dataset to increase speed. Crapps out
    # after 20 6 second files.
    uint2_dt <- h5types$H5T_NATIVE_UINT32$set_size(1)$set_precision(2)$set_sign(h5const$H5T_SGN_NONE)
    space_ds <- H5group$new(dims = c(10, 10), maxdims = c(Inf, 10))
    
    H5group[[dataType]]<-newData
    h5attr(H5group[[dataType]], "colnames") <- newCols
    
    
  }else{
    # add to the existing data
    datDims = H5group[[dataType]]$dims
    
    # some data only one dimension
    if(is.na(ncol(newData)) || is.null(ncol(newData))){
      H5group[[dataType]][datDims+1: datDims+length(newData)]<-newData
    }else{
      # modify the data in the group
      H5group[[dataType]][
        (datDims[1]+1): (datDims[1]+nrow(newData)),1:ncol(newData)]<-newData
    }
  }
  
  # H5group$close()
  
  
}


# Get frequency band limits for hybrid milidecades (usefull)
getHybridBandLimits <-function(prms, f){
  
  # Start parameters
  fCenterStart=455 # above this switch to milidecade below keep 1hz res
  idxNeg =-341
  f0 =1000 # ref frequency 
  
  # frequency centers from min analysis to 455
  fCenters = prms$lcut:fCenterStart   
  Flo= prms$lcut:fCenterStart 
  FhI= prms$lcut:fCenterStart
  
  newCenter=f0
  # Calculate all center frequencies 
  while (newCenter<max(f)){
    
    newCenter= (f0*10^(idxNeg/1000))
    
    # Low and high Figure out the integration bandss
    Flo= c(Flo, (newCenter*10^(-1/2000))) 
    FhI = c(FhI, (newCenter*10^(1/2000)))
    idxNeg = idxNeg+1
    fCenters=c(fCenters, newCenter)
  }
  
  # frequency data for output and validation
  freqLims = data.frame(lowF = (Flo), 
                        center = (fCenters),
                        highF = (FhI))
  
  return(freqLims)
}


# Hybred milidecad bands
calcHybridMiDecade <-function(prms, Psstrimmed, f, w){
  # converts PSD to much smaller milli decade bands
  # exports in decibels
  
  # Input
  # prms - analysis parameters including sample rate (Fs) and window length (N)
  # Psstrimmed - matrix containing spectra of the within the frequency limits 
  # between the user lowcut and nyquist
  # f- vector of frequency values associated with PSS trimmed
  # w- window function
  
  
  # for PSD
  delf = prms$Fs/prms$N  
  B = (1/prms$N)*(sum((w/prms$alpha)^2)) 
  
  
  # Frequency limits
  freqLims<-getHybridBandLimits(prms, f)
  
  
  # below 455, keep frequency bands
  PssKeep = Psstrimmed[,which(f<455)+1]
  
  hybMiliDecade= matrix(0, nrow = dim(PssKeep)[1], 
                        ncol = nrow(freqLims))
  
  hybMiliDecade[1:nrow(PssKeep), 1:ncol(PssKeep)] = PssKeep
  
  start = which(round(f)==round(455))+1
  
  # Step though the upper frequencies and get the mean value
  for(ii in start:nrow(freqLims)-1){
    # Create the sums
    dataSub =  Psstrimmed[,(f>= freqLims$lowF[ii]) & (f<= freqLims$highF[ii])]
    
    if(!is.null(ncol(dataSub))){
      hybMiliDecade[,ii] = rowMeans(dataSub, na.rm=TRUE)
    }else{ hybMiliDecade[,ii] = dataSub}  
    
  }
  
  
  # # delf is now a matrix- should this not change(???)
  # delf = c(1,diff(freqLims$center))
  # 
  # # Repeat the matrix
  # delf = matrix(rep(delf, nrow(hybMiliDecade)),
  #               nrow = nrow(hybMiliDecade), byrow = TRUE ) 
  
  # Third octave level
  hubridBands = 10*log10((1/B)*hybMiliDecade/(prms$pref^2))
  
  # To create the poser spectral density 
  
  dataOut= list(hubridBands,freqLims)
  
}


# Compress in time
welch_compress<-function(prms, Psstrimmed, tt){
  # I'm not 100% sure this is what the welch typically does...
  
  # Apply welch function to the data
  rA = dim(Psstrimmed)[1]; cA=dim(Psstrimmed)[2]  #number of rows in array
  lout = ceiling(rA/prms$welch)                   #length of output array
  AWelch = matrix(0,lout,cA)                        #initialize output array
  tint = ((1-prms$r)*prms$N/prms$Fs)
  
  tcompressed = tt[1:lout]
  
  for (ii in 1:lout){
    stt = tt[1] + (ii-1)*tint*prms$welch            #start time
    ett = stt +dseconds(prms$welch*tint)                  #end time
    
    tidxs = which(tt>= stt & tt<ett)
    
    nowA = Psstrimmed[tidxs,] 
    AWelch[ii,] <- (rowMeans(t(nowA), na.rm = TRUE))	#convert to dB
    tcompressed[ii] <- stt+dseconds(prms$welch*tint/2)		#assign time index
  }
  
  
  dataOut = list(AWelch, tcompressed)
  return(dataOut)
  
}


# Function to create microphone/hydrophone frequency response
calibrationMatrix<- function(prms, dims){
  
  # Produce a matrix of end-to-end frequency response values. 
  # Input
  # prms - analysis parameters
  # dims - dimensions of fft grid
  
  
  freqResp = prms$freqCal[[1]]
   
    
  if(is.data.frame(freqResp)){
     # frequencies out 
    freqOut =round(floor(prms$Fs/2)*seq(1/(prms$N/2),1, length.out =dims[1]))
    
    # Calculate the frequency response over the output frequencies from the PSS 
    
    #fit linear regression model using data frame
    model <- lm(E2E ~ f, data = freqResp)
    
    #interpolate y value based on x value of 13
    y_new = approx(freqResp$f, freqResp$E2E, xout=freqOut)$y
    
    freqCal =  matrix(rep(y_new, dims[2]), nrow = dims[2], byrow = TRUE)
    
   
  }else{
    
    # single value
    freqCal = matrix(freqResp, dims[1], dims[2])
    
  }
    return(freqCal)
    }



# Return the PSD fom which all band levels are calculated
calcPSS <- function(audioData, ii=1, prms, w, freqResp =NULL) {
  
  # audiotdata - dataframe of audio files with full directory in the 'files' 
  # column and the UTC start time of each file in the Start time
  # ii - index (which audiofile to read)
  # prms- dataframe of parameters
  # w - windo function
  # freqResp - matrix for the e2e frequency response of the hydrophone 
  
  delf = prms$Fs/prms$N  
  B = (1/prms$N)*(sum((w/prms$alpha)^2)) 
  
  # Create a tuner object to read the data
  yyObj = tuneR::readWave(audioData$files[ii])
  
  if(length(yyObj@left)>1){
    yy <- yyObj@left / 2^(yyObj@bit -1)
  }else if(length(yyObj@left)>1){
    yy <- yyObj@right / 2^(yyObj@bit -1)
  }else{
    print('No data')
    stop()
  }
  
  
  # remove dc offset
  if(prms$rmDC == TRUE){ 
    yy = yy-mean(yy)}
  
  # Break the segment into frames
  xgrid = gsignal::buffer(yy,
                          prms$N,
                          ceiling(prms$N*prms$r), 'nodelay')
  
  #remove final segment if not full
  if (any(xgrid[,dim(xgrid)[2]]==0)) { # trim xgrid and the window grid
    xgrid = xgrid[,1:dim(xgrid)[2]-1]}
  
  xgrid=t(xgrid)
  
  M = dim(xgrid)[1] 
  
  # windowing grid, should only have to figure out once
  windowGrid = matrix(rep(w/prms$alpha, M), 
                      nrow = M, byrow = TRUE)
  
  # apply the window function and transpose
  xgrid = xgrid*windowGrid
  
  
  # FFT calculation and calibration
  X = abs(apply(xgrid, 1, stats::fft))#calculate DFT of each data segment
  
  # if the calibration values are single point or a frequency response
  # Get the frequency calibration matrix
  calMat = (calibrationMatrix(prms = prms, dims=dim(X)))
  
  
  # Calibration values 
  X = X/(10^(calMat/20)) # apply the calibration here
  P = t((X/prms$N)^2)  
  
  
  # power spectral density
  Pss = 2*P[,2:(floor(prms$N/2)+1)]
  
  # frequency bins
  f = floor(prms$Fs/2)*seq(1/(prms$N/2),1, length.out =dim(Pss)[2])
  
  # index of user defined frequency bounds
  fIdx = which(f>= prms$lcut & f<= prms$hcut)
  
  
  intDur =(1-prms$r)*prms$N/prms$Fs # time between each spec bin
  dur = M*intDur-intDur # total duration of the PSS in seconds
  tt = seq(0, dur, by = intDur) 
  # timestamp for each PSS section
  tt = tt+ audioData$StartTime[ii]+lubridate::seconds(prms$secSkip)
  
  # Trim the PSS to the frequency limits 
  f = f[fIdx];                   #frequency bins in user-defined range
  
  Psstrimmed = Pss[, fIdx] # trimmed PSS
  
  # welch compress
  welchOut = welch_compress(prms, Psstrimmed, tt)
  compressedPss=welchOut[[1]]
  ttout=welchOut[[2]]
  
  
  dataOut = list(compressedPss, ttout, round(f))
  
  return(dataOut)
}


# Calculate broadband 
calcBroadband <- function(prms,Psstrimmed) {
  
  abroad = 10*log10(rowSums(Psstrimmed, na.rm = TRUE)/(prms$pref^2))
  
  
  return(abroad)
}

# Calculate decade bands
calcDecadeBands = function(prms,Psstrimmed, f, w){
  # Psstrimmed= dataOut[[1]]
  # tt= dataOut[[2]]
  # f = dataOut[[3]]
  # delf = dataOut[[4]]
  # B = dataOut[[5]]
  
  
  
  delf = prms$Fs/prms$N  
  B = (1/prms$N)*(sum((w/prms$alpha)^2)) 
  M = dim(Psstrimmed)[1] 
  
  
  #### Calculate Decade Bands ####
  # Number decade bands
  nZeros = c(floor(log10(max(c(prms$lcut, min(f)))+1))+1:
               floor(log10(min(c(prms$hcut, min(f)))+1))+1)
  
  
  nZeros = seq(floor(log10(max(c(prms$lcut, min(f))+1))) + 1,
               floor(log10(min(c(prms$hcut, max(f))+1))) + 1, by=1);
  
  # Uper and lower frequencies
  fUpper  = 10^nZeros
  fLower = fUpper/10;
  
  # Make sure the lowest value is greater than the
  # highcut
  fLower[1] = max(floor(prms$lcut), fLower[1])
  
  # NUmber of decade bands
  nDecadeBands = length(nZeros)
  
  # Preallocate decade band indicies
  fidxlow = rep(0, length(fLower))
  fidxhigh= fidxlow
  
  # Get the lower and upper index of each frequency band
  # this uses the
  for (ii in 1:nDecadeBands)
  {
    fidxs = range(which(f< fUpper[ii] &
                          f>= fLower[ii]))
    fidxhigh[ii]= fidxs[1]
    fidxlow[ii] = fidxs[2] }
  
  
  # Preallocate decade bands then populate
  decadeBands = matrix(0, nrow =dim(Psstrimmed)[1], ncol = nDecadeBands);
  
  for (ii in 1:(nDecadeBands))
  {
    decadeBands[,ii] = 10*log10(rowSums(
      Psstrimmed[,fidxlow[ii]:fidxhigh[ii]], na.rm = TRUE)/(prms$pref^2))}
  
  
  
  freqBands = data.frame(fLow =fLower, fHigh= fUpper)
  outData = list(decadeBands, freqBands)
  return(outData)
  
}

# Calcualte custom bandwidth
calcCustomBand = function(prms,Psstrimmed, f, flow,fhigh){
  # Return the SPL of the user-specified band
  # prms- input parameters
  # f - frequencies of the PSS
  # Psstrimmed -  the trimmed spectral analysis
  # flow - minimum frequency over which to summ (hz)
  # fhigh - max frequency over which to summ (hz)
  
  
  #### Calculate Decade Bands ####
  floIdx = which(f>-flow)[1]
  fhiIdx = which(f>fhigh)[1]
    customBand = 10*log10(rowSums(
      Psstrimmed[,floIdx:fhiIdx], na.rm = TRUE)/(prms$pref^2))

  outData = list(customBand )
  return(outData)
  
}


# Calcualte PSD
calcPsd<- function(prms, Psstrimmed, w){
  
  # Calcualte the power spectral density
  # Psstrimmed -  the trimmed spectral analysis
  # w - array of length N with the window values
  
  B <- (1/prms$N)*(sum((w/prms$alpha)^2))	
  delf <- prms$Fs/prms$N;					
  a <- 10*log10((1/(delf*B))*Psstrimmed[prms$lcut:prms$hcut,]/(prms$pref^2))-prms$freqCal
	
}

# Calculate third octave bands
calcThirdOctBands = function(prms, Psstrimmed,f, w){
  # Calculate third octave bands
  # prms- input parameters
  # Psstrimmed -  the trimmed spectral analysis
  # f - array of frequencies (Hz) coinciding with the PSStrimmed columns
  
  delf = prms$Fs/prms$N  
  B = (1/prms$N)*(sum((w/prms$alpha)^2)) 
  M = dim(Psstrimmed)[1] 
  
  low13BAND = ifelse(prms$lcut < 25, 25, prms$lcut)
  
  lobandf = floor(log10(low13BAND))   #lowest power of 10 frequency for 1/3
  
  # octave band computation
  hibandf = ceiling(log10(prms$hcut))    #highest ""
  nband = 10*(hibandf-lobandf)+1 #number of 1/3-octave bands
  fc = matrix(data =0, nrow=1, ncol=nband)            #initialise 1/3-octave frequency vector
  fc[1] = 10^lobandf             #lowest frequency = lowest power of 10
  
  
  
  # Calculate centre frequencies (corresponds to EQUATION 13) #calculate 1/3 octave centre
  for (i in 2:nband){
    fc[i] = fc[i-1]*10^0.1}     # frequencies to (at least) precision
  # of ANSI standard
  
  #crop frequency vector to frequency
  fc = fc[(fc >= low13BAND & fc <= prms$hcut)]
  nfc = length(fc)               #number of 1/3 octave bands
  
  # Calculate boundary frequencies of each band (EQUATIONS 14-15)
  fb = fc*10^-0.05;               #lower bounds of 1/3 octave bands
  fb[nfc+1] = fc[nfc]*10^0.05;    #upper bound of highest band (upper
  #   bounds of previous bands are lower
  #   bounds of next band up in freq.)
  if (max(fb) > prms$hcut){
    nfc = nfc-1
    fc = fc[1:nfc]}             #   remove highest band
  
  
  # Calculate 1/3-octave band levels (corresponds to EQUATION 16)
  P13 = matrix(data = 0, 
               nrow = dim(Psstrimmed)[1], 
               ncol = nfc)             #initialise TOL array
  
  for (ii in 1:nfc ){                  
    
    # frequency range
    frange = which(f >= fb[ii] & f < fb[ii+1])
    for (jj in 1:M){
      fcl = sum(Psstrimmed[jj,frange],na.rm = TRUE)
      P13[jj,ii] = fcl}
  }
  
  # Third octave level
  TOL = 10*log10((1/B)*P13/(prms$pref^2))
  
  outData = list(TOL, round(fc))
  return(outData)
  
}

# create LTSAs from the h5database
makeLTSA <- function(instrument_group, averagingPeriod='60 min'){
  
  # input, 
  # instrument_group- h5df instrument
  # hybrid decades are stored
  # averagingPeriod - lubridate styled character string of averaging duration

  frequencyData <- instrument_group[['hybridDecFreqHz']]
  timeStamps <- instrument_group[['DateTime']]
  
  # fill values, figure out wher ethe data end
  maxDataLen <- min(which(timeStamps[1:timeStamps$maxdims] ==  
                            "0000-00-00 00:00:00"))-1
  
  timeStamps <- timeStamps[1:maxDataLen]
  timeStamps<- lubridate::ymd_hms(timeStamps)
  
  # the hybrid milidecades and the center frequencies
  data<- instrument_group[['hybridMiliDecLevels']]
  
  # Spacing- for longer datasets 60 min common
  medianSpacing <- lubridate::duration(averagingPeriod)
  
  # Create the dataframe from which plotting will happen
  tstart<- (timeStamps[1])
  timeChunks <-seq(tstart, (timeStamps[maxDataLen]+duration('1 minute')), 
                   by = medianSpacing)
  
  # pre-allocate the data and format for ggplot
  ggData <- expand.grid(freq = frequencyData[],
                        time = tt )
  ggData$Median <- NaN
  
  
  for (ii in 1:((length(timeChunks)-1))){
    
    chunIdx <- which(timeStamps>= timeChunks[ii] & timeStamps< timeChunks[ii+1])
    
    # Pull the data and get the medin for each frequency
    dataChunk <- data[chunIdx, ]
    medVals <- apply(dataChunk, 2, median, na.rm=T)
    
    # fill in the dataframe
    idxStart <- ((ii-1)*frequencyData$dims)+1
    idxStop <- idxStart+frequencyData$dims
    ggData$Median[idxStart:idxStop]<- medVals
    
  }
  
  
  # We need to use geom_rect to define fredquency widths
  hybridFreqDiff <- c(1, diff(frequencyData[]))
  ggData$fmin <- ggData$freq-(hybridFreqDiff*0.5)
  ggData$fmax <- ggData$freq+(hybridFreqDiff*0.5)
  
  ggData<-ggData[!is.nan(ggData$Median),]
  ggData=ggData[ggData$Median>-2 & ggData$Median<2,]
  
  # Create the plot
  p<- ggplot(ggData)+
    geom_rect(aes(ymin=fmin,
                  ymax=fmax,
                  xmin =time+medianSpacing*.5, 
                  xmax= time-medianSpacing*.5,
                  fill = Median))+
    scale_fill_gradientn(colors = viridis::viridis(10), 
                         na.value = "transparent",
                         name = 'Median Level')+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))+
    annotation_logticks(sides = 'l')+  
    ylab('Frequency Hz')+
    xlab('Date')+
    theme_bw()  
  

  return(p)
} 

# create hourly figure
makeHourly<-function(instrument_group){

  # Create  figure with hours on the y axis and frequency on the x
  # input, 
  # instrument_group- h5df instrument
  # hybrid decades are stored
  
  frequencyData <- instrument_group[['hybridDecFreqHz']]
  timeStamps <- instrument_group[['DateTime']]
  
  
  # fill values, figure out wher ethe data end
  maxDataLen <- min(which(timeStamps[1:timeStamps$maxdims] ==  
                            "0000-00-00 00:00:00"))-1
  
  timeStamps <- timeStamps[1:maxDataLen]
  timeStamps<- ymd_hms(timeStamps)
  
  # Initialize the data frame for the intended output style
  ggdataDaily <- expand.grid(freqs = round(frequencyData[]),
                             hours=0:24)
  ggdataDaily$Median <- NaN
  
  # We need to use geom_rect to define fredquency widths
  hybridFreqDiff <- c(1, diff(frequencyData[1:frequencyData$dims]))
  ggdataDaily$fmin <- ggdataDaily$freqs-(hybridFreqDiff*0.5)
  ggdataDaily$fmax <- ggdataDaily$freqs+(hybridFreqDiff*0.5)
  
  # Pull out a bit of the time data and figure out which indices are going to 
  # represent which days
  
  for(ii in 1:24){
    
    # get the index of the hour to load
    hrIdx <- which(hour(timeStamps)==(ii-1))
    
    
    # Load the hybrid levels-one chunk at a time
    dataChunk <- data[hrIdx,]
    dataChunk[is.infinite(dataChunk)]<-NA
    
    # median values for each frequency level
    medVals <- apply(dataChunk, 2, median, na.rm=T)
    ggdataDaily$Median[ggdataDaily$hours==(ii-1)]<-medVals
    
  }
  
  
  # Spiff it up like prom night
  p<-ggplot(ggdataDaily)+
    geom_rect(aes(xmin=fmin,
                  xmax=fmax,
                  ymin =hours, 
                  ymax= hours+1,
                  fill = Median))+
    scale_fill_gradientn(colors = viridis::viridis(10), 
                         na.value = "transparent",
                         name = 'Median Level')+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))+
    annotation_logticks(sides = 'b')+  
    xlab('Frequency Hz')+
    ylab('Hour of Day')+
    theme_bw()  
  
  return(p)
}

# create the PSD distribution plot
makePSDDist<-function(instrument_group){
  
  timeStamps <- instrument_group[['DateTime']]
  
  # Find end of usable data (note we made the bigger than it needed to be)
  maxDataLen <- min(min(which(timeStamps[1:timeStamps$maxdims] == 
                            "0000-00-00 00:00:00"))-1,
                    timeStamps$maxdims)
  
  frequencyData <- instrument_group[['hybridDecFreqHz']]
  timeStamps<- ymd_hms(timeStamps[1:maxDataLen],tz = 'UTC')
  
  
  
  #################################################
  ## PSD plot
  #################################################
  
  
  # the hybrid milidecades and the center frequencies
  data<- instrument_group[['hybridMiliDecLevels']]
  
  
  
  # Dimensions of the dataset
  M<- maxDataLen#data$dims[1]
  N <-data$dims[2]
  
  # Decibel spacing- for now pick sensible ranges. Ideally pull from databse
  dBint <-1
  dbMin<- 0
  dbMax<- 200
  dbVals <- seq(dbMin, dbMax, by = dBint)
  
  
  # Create a dataframe for the matrix psd data (will populate with counts from the
  # database)
  ggData <- data.frame(
    freqs= sort(rep(frequencyData[1:N], length(dbVals)-1)),
    counts = rep(NaN, data$dims[2]*(length(dbVals)-1)),
    levelbins =  rep(dbVals[-1], data$dims[2]))
  
  # Step through the data in chunks and count the values in each 
  # frequency/ sound level bin
  for(ii in 1:N){
    
    # Read a chunk of the dataset
    dataChunk <- data[,ii]
    dataChunk <- dataChunk[!is.infinite(dataChunk)]
    
    # Create the histogram counts
    histDatat <- hist(dataChunk,breaks = dbVals,
                      plot = FALSE,na.rm=TRUE)$counts
    
    idxStart<- ((ii-1)*length(histDatat))+1
    idxStop <- idxStart+length(histDatat)-1
    ggData$counts[idxStart:idxStop]<-histDatat
    
  }
  
  # set NAN
  ggData$counts[ggData$counts==0]<-NA
  
  # Ccreate dummy frequency for plottinng
  ggData$freqHack <- sort(rep(1:data$dims[2],
                              length(dbVals)-1))
  
  # Normalize for density
  ggData$vals<- ggData$counts/(dBint*(M))
  
  # We need to use geom_rect to define fredquency widths
  hybridFreqDiff <- c(1, diff(frequencyData[1:frequencyData$dims]))
  ggData$fmin <- ggData$freqs-(hybridFreqDiff*0.5)
  ggData$fax <- ggData$freqs+(hybridFreqDiff*0.5)
  ggData$levelMin <- ggData$levelbins-(0.5*dBint)
  ggData$levelMax <- ggData$levelbins+(0.5*dBint)
  
  # Remove data for which there were no data
  ggdataCleaned <- ggData[!is.na(ggData$counts),]
  
  # Spiff it up like prom night
  p<-ggplot(ggData)+
    geom_rect(aes(xmin=fmin,xmax=fax,
                  ymin =levelMin, ymax= levelMax,
                  fill = vals))+
    scale_fill_gradientn(colors = viridis::viridis(10), 
                         na.value = "transparent",
                         name = 'Emp. Prob. Den')+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))+
    annotation_logticks(sides = 'b')+  
    xlab('Frequency Hz')+
    ylab('PSD dB re 1uPa^2/Hz')+
    theme_bw()
  
  
  return(p)
}

calcMetrics<-function(prms, Psstrimmed, f, w){
  # Driver function to return a matrix of 
  # metrics - list of metrics to calculate including hybrid, broadband, decadeband
  # and thirdoct
  # prms- dataframe of parameters
  # Psstrimmed - power spectrum
  # f - frequencies associated with the PSS
  
  metrics<-prms$metrics
  
  
  outList = list()
  if('hybrid' %in% metrics){
    # Hybrid milidecade from PSS
    hybridMilidecade = calcHybridMiDecade(prms, Psstrimmed, f, w)
    
    # add key value pair
    outList[['hybLevels']]<-hybridMilidecade[[1]]
    outList[['hybFreqs']]<-hybridMilidecade[[2]]
  }
  
  if('thirdoct' %in% metrics){
    # Third Ocatave Bands (checked)
    thirdOctBands= calcThirdOctBands(prms, Psstrimmed,f, w)
    
    # add key value pair
    outList[['thridOctLevels']]<-thirdOctBands[[1]]
    outList[['thirdOctF']]<-thirdOctBands[[2]]
  }
  
  if('decadeband' %in% metrics){
    # Decade bands
    DecadeBands =calcDecadeBands(prms, Psstrimmed,f, w)
    
    # add key value pair
    outList[['DecadeBandsLevels']]<-DecadeBands[[1]]
    outList[['DecadeBandsF']]<-DecadeBands[[2]]
  }
  
  if('broadband' %in% metrics){
    # Broadband (checked)
    BroadBandLevels =calcBroadband(prms,Psstrimmed)
    
    # add key value pair
    outList[['BroadBandLevels']]<-BroadBandLevels}
  
  return(outList)
  
}

  
# # Function for processing the data within the GUI
# analyzeTheData <- function(audioData,  prms, w, ProjName, instrumentName, input){
#   
#   # if file doesn't already exist, create it
#   if(!file.exists(ProjName)){
#     h5createFile(ProjName)
#   }
#   
#   # create group for location 1 (here we are going with just the hydrophone serial number)
#   h5createGroup(ProjName, instrumentName)
#   
#   # add the parameter and the files
#   h5write(
#     prms,
#     file =ProjName,
#     paste(instrumentName,"Parms",sep="/"))
#   
#   # Write all the audio information to the database for idiot checking in
#   # the future
#   h5write(
#     audioData$FileName,
#     file = ProjName,
#     paste(instrumentName,"Files",sep="/"))
#   
#   # Estimate the number of rows to pre-allocate the datasets
#   dataSetLenght = sum(ceiling(audioData$Duration/60))
#   
#   
#   idStart =1
#     for(ii in 1:nrow(audioData)){
#       
#       # Calculate the PSS within the user defined range, time stamps, and frequency
#       # vector.
#       dataOut = calcPSS(audioData, ii, prms, w)
#       
#       
#       Psstrimmed= dataOut[[1]]
#       tt= dataOut[[2]]
#       f = dataOut[[3]]
#       avPSD= 10*log10(Psstrimmed)
#       
#       # Hybrid milidecade from PSS
#       if("Hybrid Milidecade" %in% input$analysis_type){
#         hybridMilidecade = calcHybridMiDecade(prms, Psstrimmed, f, w)
#         hybLevels = hybridMilidecade[[1]]
#         hybFreqs = hybridMilidecade[[2]]}
#       
#       # Third Ocatave Bands (checked)
#       if("OneThird Octave Bands" %in% input$analysis_type){
#         thirdOctBands= calcThirdOctBands(prms, Psstrimmed,f, w)
#         thridOctLevels = thirdOctBands[[1]]
#         thirdOctF = thirdOctBands[[2]]}
#       
#       # Decade bands
#       if("Decade Band" %in% input$analysis_type){
#         DecadeBands =calcDecadeBands(prms,Psstrimmed, f, w)
#         DecadeBandsLevels =DecadeBands[[1]]
#         DecadeBandsF = DecadeBands[[2]]}
#       
#       # Broadband (checked)
#       if("Broadband" %in% input$analysis_type){
#         BroadBandLevels =calcBroadband(prms,Psstrimmed)}
#       
#       # Custom Bands
#       if("Custom Band" %in% input$analysis_type){
#         customLevels = calcCustomBand(prms,Psstrimmed, f,
#                                       input$custLoF,input$custHiF)[[1]]}
#       
#       # Duration of the output analysis
#       countLen = length(tt)
#       
#       # First run, add add the frequency information
#       if(ii==1){
#         
#         
#         # Write the hybrid frequencies
#         if("Hybrid Milidecade" %in% input$analysis_type){
#           writeToH5datarH5df(ProjName, instrumentName,
#                              dataType='hybridDecFreqHz',
#                              newData = hybFreqs$center,
#                              dataStart=1,
#                              maxRows=nrow(hybFreqs),
#                              storagemMode='integer')}
#         
#         # Write the hybrid frequencies
#         if("Decade Band" %in% input$analysis_type){
#           writeToH5datarH5df(ProjName, instrumentName,
#                              dataType='decadeFreqHz',
#                              newData = DecadeBandsF$fLow,
#                              dataStart=1,
#                              maxRows=nrow(DecadeBandsF),
#                              storagemMode='integer')}
#         
#         # Write the third-octave frequencies
#         if("OneThird Octave Bands" %in% input$analysis_type){
#           writeToH5datarH5df(ProjName, instrumentName,
#                              dataType='thirdOctFreqHz',
#                              newData = thirdOctF,
#                              dataStart=1,
#                              maxRows= length(thirdOctF),
#                              storagemMode='integer')}
#         
#         # Custom Bands
#         if("Custom Band" %in% input$analysis_type){
#           writeToH5datarH5df(ProjName, instrumentName,
#                              dataType='customBand',
#                              newData = c(input$custLoF,input$custHiF),
#                              dataStart=1,
#                              maxRows= length(thirdOctF),
#                              storagemMode='integer')}
#         
#         ###################################################
#         # Add new data- will automatically create dataset on first row
#         ###################################################
#         
#         
#         # Write the timestamps
#         writeToH5datarH5df(ProjName, instrumentName,
#                            dataType='DateTime',
#                            newData = as.matrix(as.character(tt)),
#                            dataStart=idStart,
#                            maxRows<-dataSetLenght,
#                            storagemMode<- 'character')
#         
#         # write the hybrid milidecade levels
#         if("Hybrid Milidecade" %in% input$analysis_type){
#           writeToH5datarH5df(ProjName, instrumentName,
#                              dataType<-'hybridMiliDecLevels',
#                              newData <- hybLevels,
#                              dataStart<-idStart,
#                              maxRows<-dataSetLenght)}
#         
#         # write the third octave band levels
#         if("OneThird Octave Bands" %in% input$analysis_type){
#           writeToH5datarH5df(ProjName, instrumentName,
#                              dataType<-'thirdOctLevels',
#                              newData <- thridOctLevels,
#                              dataStart<-idStart,
#                              maxRows<-dataSetLenght)}
#         
#         # write the decade band levels
#         if("Decade Band" %in% input$analysis_type){
#           writeToH5datarH5df(ProjName, instrumentName,
#                              dataType<-'decadeLevels',
#                              newData <- DecadeBandsLevels,
#                              dataStart<-idStart,
#                              maxRows<-dataSetLenght)}
#         
#         
#         # write the broadband levels
#         if("broadband" %in% input$analysis_type){
#           writeToH5datarH5df(ProjName, instrumentName,
#                              dataType<-'thirdOctLevels',
#                              newData <- thridOctLevels,
#                              dataStart<-idStart,
#                              maxRows<-dataSetLenght)}
#         
#         # write the broadband levels
#         if("Custom Band" %in% input$analysis_type){
#           writeToH5datarH5df(ProjName, instrumentName,
#                              dataType<-'customLevels',
#                              newData <- customLevels,
#                              dataStart<-idStart,
#                              maxRows<-dataSetLenght)}
#         
#         
#         idStart =idStart+countLen
#         print(ii)
#         #incProgress(1/nrow(audioData), detail = paste("File", ii, 'of ', nrow(audioData)))
#       }
#       
#     }
# }
#   
# 
