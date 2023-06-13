

rm(list = ls())
library(here)
library(lubridate)
library(stringr)
library(rhdf5)
library(tuneR)


#Define path to functions that we will use
source(here('NoiseProcessingFxs.r'))

# Scratch directory to set up noise analysis metrics following merchant 2012 (I think)


################################################################
# Load list of sound files and start times, assumes all audio files are the 
# same duration
################################################################

# Assumes 
# 1) file names in UTC
# 2) All sample rates are the same

fileLoc = 'D:\\RECORDINGS\\ADRIFT_001_CENSOR_12kHz'
fileLoc ='D:\\RECORDINGS\\ADRIFT_002_CENSOR_12kHz'
#fileLoc = 'C:\\Users\\kaitlin.palmer\\Downloads\\mee312330-sup-0002-datas1\\PAMGuide'

fileLoc ='D:\\Recordings\\ADRIFT_024_CENSOR\\ST6664'
fileLoc = 'C:\\Users\\kaitlin.palmer\\Desktop/ADRIFTLF/ADRIFT_024'

# Jasco Pattern
nameStringPattern = '\\d{8}T\\d{6}.\\d{3}'
nameStringFormat =" %Y%m%dT%H%M%OS"

# Soundtrap Pattern
nameStringPattern = '\\d{12}'
nameStringFormat ='%y%m%d%H%M%OS'
lubradiateFormat ="ymdHMS"

# Set up the audiodataframe and iniitalize parameters
dataOut = createAudioDataframe(fileLoc,nameStringPattern = nameStringPattern,
                                 lubradiateFormat = lubradiateFormat)
prms=dataOut$prms
audioData= dataOut$audioData

######################################################################
# Set up the analysis parameters
###################################################################


##
# User defined parameters
##

prms$N = prms$Fs # fft lenght in samples
prms$r = 0.5 # overlap in percent
prms$aveSec = 60 # averaging window (seconds)
prms$secSkip = 0 # number of seconds to skip when loading the file (soundtraps...)

# ratio of new to original window lengths in Welch method (avg over aveSec secs)
prms$welch = prms$aveSec*(prms$Fs/prms$N)/(1-prms$r)

# Analysis frequencies (hz). 
prms$lcut =10
prms$hcut = prms$Fs/2

# Reference pressure, set to 1 for water 20 for air
prms$pref =1

# remove DC offset by subtracting mean(yy) from yy
prms$rmDC = TRUE

# Create window functions as defined by NM 2014
dataOut = createWindow('Hann', prms)
w =dataOut[[1]]
prms$alpha = dataOut[[2]]

# Which metrics to calculate, case sensitive
prms$metrics = t(list('hybrid', 'broadband', 'decadeband', 'thirdoct'))


###################################
# Calibration, either end to end or frequency response
###########################################

# Frequency response example, make sure start/end frequencies cover
# low frequency analysis end and nyquest 
freqResp = data.frame(f = round(10:prms$Fs/2))
freqResp$E2E = seq(-185, -165, length.out = nrow(freqResp))
prms$freqCal= list(freqResp)

# End to end example
prms$freqCal= -161.1
prms$StartDate=format(min(audioData$StartTime), "%Y-%m-%dT%H:%M:%SZ")
prms$EndDate = format(max(audioData$EndTime), "%Y-%m-%dT%H:%M:%SZ")


################################################
# Other metrics for ICES/completeness
#####################################################

# Date the analysis was run
prms$DateRun <- now(tzone = "UTC")

# as ICES 
prms$Comments <- 'This is a test of the emergency noise analysis system'

# Measurement depth in meters
prms$MeasurementHeight = 5


prms$HydrophoneType = 'HTI'
prms$HydrophoneSerialNumber ='856163'
prms$RecorderType = 'Soundtrap'
prms$RecorderSerialNumber = "ST6664"
prms$MeasurementPurpose = 'OTH' # Example data not part of the ICES
prms$MeasurementSetup = 'AUT'
prms$RigDesign='VMC'
prms$FrequencyCount
prms$FrequencyIndex
prms$FrequencyUnit = 'Hz'
prms$ChannelCount = 1
prms$MeasurementTotalNo = NA
prms$MeasurementUnit = 'dB re 1μPa'
prms$AveragingTime = prms$aveSec
prms$ProcessingAlgorithm = 'MER'
prms$CalibrationProcedure ='OTH'



################################################################
# Set up the H5DF file and group
################################################################

# Project name with current date (latter mostly for debugging)
ProjName = paste0("Adrift_2022",format(Sys.time(), "%Y-%m-%d%H%M%S"), '.h5')


# Instrument name, here I'm using soundtrap IDs
instrumentName = "ST5987"
instrumentName = "ST6664"

# Create the  hdf5 file and fill out meta and audiofiles
h5createFile(ProjName)

# create group for location 1 (here we are going with just the hydrophone serial number)
h5createGroup(ProjName, instrumentName)

# add the parameter and the files
h5write(
  prms,
  file =ProjName,
  paste(instrumentName,"Parms",sep="/"))

h5write(
  audioData$FileName,
  file = ProjName,
  paste(instrumentName,"Files",sep="/"))

# Estimate the number of rows to pre-allocate the datasets
dataSetLenght = sum(ceiling(audioData$Duration/60))




# ######################################################################
# # parallel version to make the h5df
# #######################################################################
# 
# audioDatasub = audioData[1:5,]
# 
# processAudiotoH5DF<-function(prms, audioData, w){
#   # Calculate the PSS within the user defined range, time stamps, and frequency
#   # vector.
#   dataOut = calcPSS(audioData, ii, prms, w)
#   
#   
#   Psstrimmed= dataOut[[1]]
#   tt= dataOut[[2]]
#   f = dataOut[[3]]
#   avPSD= 10*log10(Psstrimmed)
#   
#   # Hybrid milidecade from PSS
#   hybridMilidecade = calcHybridMiDecade(prms, Psstrimmed, f, w)
#   hybLevels = hybridMilidecade[[1]]
#   hybFreqs = hybridMilidecade[[2]]
#   
#   # Third Ocatave Bands (checked)
#   thirdOctBands= calcThirdOctBands(prms, Psstrimmed, tt,f)
#   thridOctLevels = thirdOctBands[[1]]
#   thirdOctF = thirdOctBands[[2]]
#   
#   # Decade bands
#   DecadeBands =calcDecadeBands(prms,dataOut[[1]], dataOut[[2]],dataOut[[3]])
#   DecadeBandsLevels =DecadeBands[[1]]
#   DecadeBandsF = DecadeBands[[2]]
#   
#   # Broadband (checked)
#   BroadBandLevels =calcBroadband(prms,dataOut[[1]])
#   
#   countLen = length(tt)
# }
# 
# num_cores <- detectCores()
# cl <- makeCluster(num_cores)
# parLapply(cl, file_paths, calc_bblvl_write_csv, csv_file_path)
# stopCluster(cl)
# 
# #
# tstart = Sys.time()
# results <- parLapply(cl = parallelCluster,
#                      dataLocs1,
#                      processAudiotoH5DF)


######################################################################
# Step through the soundfiles, create spectrogram and write to hdf5 file
#######################################################################
idStart =1


for(ii in 1:nrow(audioData)){
  
  # Calculate the PSS within the user defined range, time stamps, and frequency
  # vector.
  dataOut = calcPSS(audioData, ii, prms, w)
  
  
  Psstrimmed= dataOut[[1]]
  tt= dataOut[[2]]
  f = dataOut[[3]]
  avPSD= 10*log10(Psstrimmed)
  
  
  allMetrics<-calcMetrics(prms, Psstrimmed, f, w)

  # # Hybrid milidecade from PSS
  # hybridMilidecade = calcHybridMiDecade(prms, Psstrimmed, f, w)
  # hybLevels = hybridMilidecade[[1]]
  # hybFreqs = hybridMilidecade[[2]]
  # 
  # # Third Ocatave Bands (checked)
  # thirdOctBands= calcThirdOctBands(prms, Psstrimmed,f, w)
  # thridOctLevels = thirdOctBands[[1]]
  # thirdOctF = thirdOctBands[[2]]
  # 
  # # Decade bands
  # DecadeBands =calcDecadeBands(prms, Psstrimmed,f, w)
  # DecadeBandsLevels =DecadeBands[[1]]
  # DecadeBandsF = DecadeBands[[2]]
  # 
  # # Broadband (checked)
  # BroadBandLevels =calcBroadband(prms,Psstrimmed)
  
  countLen = length(tt)

  # First run, add add the frequency information
  if(ii==1){

    if('hybrid' %in% prms$metrics)

    # Write the hybrid frequencies
    writeToH5datarH5df(ProjName, instrumentName,
                       dataType='hybridDecFreqHz',
                       newData = round(hybFreqs$center),
                       dataStart=1,
                       maxRows=nrow(hybFreqs),
                       storagemMode='double')

    # Write the hybrid frequencies
    writeToH5datarH5df(ProjName, instrumentName,
                       dataType='decadeFreqHz',
                       newData = DecadeBandsF$fLow,
                       dataStart=1,
                       maxRows=nrow(DecadeBandsF),
                       storagemMode='integer')

    # Write the third-octave frequencies
    writeToH5datarH5df(ProjName, instrumentName,
                       dataType='thirdOctFreqHz',
                       newData = thirdOctF,
                       dataStart=1,
                       maxRows= length(thirdOctF),
                       storagemMode='integer')}
  ###################################################
  # Add new data- will automatically create dataset on first row
  ###################################################


  # Write the timestamps
  writeToH5datarH5df(ProjName, instrumentName,
                     dataType='DateTime',
                     newData = as.matrix(as.character(tt)),
                     dataStart=idStart,
                     maxRows<-dataSetLenght,
                     storagemMode<- 'character')

  # write the hybrid milidecade levels
  writeToH5datarH5df(ProjName, instrumentName,
                     dataType<-'hybridMiliDecLevels',
                     newData <- hybLevels,
                     dataStart<-idStart,
                     maxRows<-dataSetLenght)

  # write the third octave band levels
  writeToH5datarH5df(ProjName, instrumentName,
                     dataType<-'thirdOctLevels',
                     newData <- thridOctLevels,
                     dataStart<-idStart,
                     maxRows<-dataSetLenght)
  
  # write the decade band levels
  writeToH5datarH5df(ProjName, instrumentName,
                     dataType<-'decadeLevels',
                     newData <- DecadeBandsLevels,
                     dataStart<-idStart,
                     maxRows<-dataSetLenght)

  idStart =idStart+countLen
  
  print(ii)
}

H5Fclose(ProjName)

#######################################################################
# Opening projects and plotting
#######################################################################
h5closeAll()
rm(list = ls())
library(ggplot2)
library(scales)
library(hdf5r)
library(lubridate)
source('NoiseProcessingFxs.R')
ProjName <-  "Adrift_2022.h5"
instrumentName <- "ST6664"

#open the files
my_file <- h5file(ProjName, "r")
instrument_group <- my_file[[instrumentName]]


#################################################
## PSD plot
#################################################

p<-makePSDDist(instrument_group)
print(p)



#######################################################################
# Create Hrly fig 
#######################################################################

#open the files
my_file <- h5file(ProjName, "r")
instrument_group <- my_file[[instrumentName]]
p<- makeHourly(instrument_group)

print(p)

#######################################################################
# Create LTSA - median hourly
#######################################################################
# Make an LTSA by stepping through the code and creating the median hourly levels

#open the files
my_file <- h5file(ProjName, "r")
instrument_group <- my_file[[instrumentName]]

p<- makeLTSA(instrument_group, averagingPeriod = '90 min')

print(p)



############################################################################
# Update ICES metadata group
############################################################################
my_file <- h5file('C:/Users/kaitlin.palmer/Downloads/ICES-Continuous-Underwater-Noise-format (1)/Sample.h5', "r")
metaGroup<-my_file[['Metadata']]

for(ii in 1:length(metaGroup$names)){
  metaName = metaGroup$names[ii]
  
  ds<-metaGroup[[metaName]]
}



