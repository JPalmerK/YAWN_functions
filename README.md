# SoundAnalytics: Yet Another Way to implement Noise analyis (YAWN)

![Rplot01](https://user-images.githubusercontent.com/28478110/221385856-2d31cb30-b8af-4dc7-81f8-66f9c8141849.png)

____________
## TL;DR
An R-based methodology for storing noise metrics in [HDF5 (Hierarchical Data Format version 5)](https://www.hdfgroup.org/solutions/hdf5/)  files. The code is modularized so a user can pick and chose which metrics to store. This repository is licensed as open source software to help the community deal with data irregularities. 

Pre-written metrics
1) Hybrid-millidecade band levels
2) Third-octave band levels
3) Decade band levels
4) Broadband level

Pre-written plotting functions. These returns standard plots in ggplot format.
1) LTSA (Long-term Spectral Average)
2) Probability distribution
3) Hourly noise level distributions
_______________________
## Introduction

Measuring ambient noise levels is an important part of many ecological studies, particularly in the marine environment where noise levels are both an emerging conservation concern for many signaling or listening species and a hindrance to passive acoustic monitoring for signals from soniferous species. Anthropogenic noise has been implicated in the stranding of deep diving marine species such as beaked whales and as a stressor in other marine mammals including endangered killer whales and right whales. With ongoing acoustic monitoring efforts for many vocally active species, quantitative assessment of background noise levels is also key in understanding changes in detection range which could bias monitoring or real-time conservation efforts.

There are plethora of free and paid services used to calculate noise metrics including PAMGuard (www.Pamguard.org), [Triton](https://www.cetus.ucsd.edu/technologies_triton.html), and [PAMGuide](https://sourceforge.net/projects/pamguide/). Each of these systems have their benefits and limitations and bioacousicians often find themselves in need of modification for their own specific data needs. For instance, Triton and its associate Remoras are also free and there is a compiled version that does not require MATLAB. However, any customization does require a MATLAB license which is frequently cost prohibitive, especially for researchers in developing countries. PAMGuard is both free and an industry standard but can be difficult to work with and JAVA is not commonly known among biologists. This makes it challenging for researchers to troubleshoot without reaching out to a small but dedicated team of maintainers. Finally, PAMGuide was written in both R and MATLAB and includes a Matlab-GUI, allowing for user-friendly interface. 

The [PAMGuide paper (Merchant et. al 2014)](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12330) was well received, with over 200 citations, in no small part because of the well-documented and published code provided by the author. However, it too has limitations have required extensive modifications for many long-term noise projects and, after 9 years, the code was due for some updates. Principle among the limitations the speed of analysis and the storage options initially provided (either mat files or .csv files). These become ungainly at best and untenable at worst when working with large, multi-instrument, or multi-year arrays. 

An ideal storage solution would allow storage for large, multi-level datasets with descriptors, metadata, efficient access speed and accessibility across multiple platforms. Presently HDF5 databases meet these criteria. Such formats are platform independent, self-described, and open source (sharing is caring).

Additionally [Martin et al., 2021](https://static1.squarespace.com/static/52aa2773e4b0f29916f46675/t/6033d0181ce4934ad7c3d913/1614008346204/Martin_et_al_2021_Hybrid+millidecade+spectra_practical+format+for+ambient+data+exchange.pdf) provided a efficient methodology for storing and sharing large database of sound metrics. These metrics provide an efficient methodology for sharing large-scale and long-term datasets. 

## YAWN goals, features, & status

The principal goal of YAWN is to produce a reliable and flexible system for recording noise metrics. In achieving this goal I required that the system be built in a well-established language within the biological field ([R](https://www.r-project.org/)) which is free and open source software under the GPL license. These characteristics make the language as accessable to as many researchers as possible. I also required that the data storage should also allow for multi-level organization, storage of large files, be accessable to multiple software packages/languages, and ultimately itegrate with the sturdy metadata managment system, [Tethys](https://tethys.sdsu.edu/). 

In this repository I've used the R-based FFT and level metric calculations by Merchant et. al (2014), but updated the functions and wrappers to allow for procesing of large numbers of sound files and storing in a HDF5 database.   

This repository contains a modidified version of [Merchant et. al's 2014](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12330) noise analysis tools. The major modifications between this version and the origional incude:

1) Saving data to HDF5 files rather than CSV (comma-separated variable) files. This allows for much larger data sets to be stored efficiently, including multiple deployments within a single study;
2) Saving large PSD (power spectral density) results as hybrid milidecades (e.g. [Martin et al. 2021](https://static1.squarespace.com/static/52aa2773e4b0f29916f46675/t/6033d0181ce4934ad7c3d913/1614008346204/Martin_et_al_2021_Hybrid+millidecade+spectra_practical+format+for+ambient+data+exchange.pdf));
3) Allowing user-specified anlysis parameters to be saved as part of database;
4) Modularizing as many aspects of the code as possible in order to integrate with SHINY apps (in the future) and allow for readability;
5) Vectorizing where possible to increase speed;
6) Ability to exclude the first N seconds of the file (this is useful for Soundtraps);
7) Example code for constructing popular noise level plots from the HDF5 files.

Computational results have been validated with test data against the most recent version of PAMGuide, so the code should provide consistant results.

This codebase is still in developmnet (half-baked) and is not ready for release. Use at your own risk.

___

## Limitations

The package used to load audio only handles WAV files (:sob:). The AV package was initially attempted because of the variety of sound files it incorporates, but it was not possible to create validated levels. At least not within the author's available time/patienence! Suggestions on more flexible approaches are welcome. 

This system does not work well with high frequency audio data (i.e. high sample rates). R has not been optimized for this and the approach is memory intensive. I do not recommend for data sampled at rates greater than 48 kHz. 

I have not implemented all of the features of the origional PAMGuide as they do not pertain to much of the my work.

It is not possible to add to HDF5 datasets after the initial creation as the dimensions are set at the start of the run. This means all sound files for which you intend to be in the same dataset must be in one folder. The implementation created here pulls from a scratch directory...

Unlike the original version, this code is not set up to read chunks of sound into memory, instead it must read the whole file. As such, very large and/or high frequency data are likely to gum up the works, as it were. This may or may not be something that gets updated in the future.

At present, only end-to-end calibration values are accepted.

A code review had not been preformed. Volunteers welcome. Snacks provided. 

___
## To do (in no particular order)
1) ~~Implement figure code with hybrid mili-decade structure and HDF5 data (this is the first priortiy actually, not immediately clear the best approach for doing this)~~
2) Incorporate into SHINY app or similar to allow easy access for new/non-acousticians- *Not happening. Too slow. Better establish csv file with suggested parameters.*
3) ~~low users to load frequency response of hydropone/microphone system rather than single value~~ fx structure updated, simulated data there
4) ~~Modularize write to HDF5 section~~
5) ~~Match window options with NM version~~
6) ~~Ensure that delta frequency calculations are correct for decade, third octave, and hybrid bands as the frequency bandwidth varries. Ensure that this does not result in incorrest PS**D** values~~
7) Ensure consistency with [ICES](https://www.ices.dk/data/data-portals/Pages/Continuous-Noise.aspx) database 
8) Implement parallel processing
 
___
## Notes from the trenches

**Window Functions**

The hann window from the R gsignal library varies considerably from the R function outlined in [Merchant et. al (2014)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12330). This is also true in the MATLAB version of PAMGuide, i.e. hann(5) ≠ (0.5 - 0.5*cos(2*pi*(1:5)/5)). In toying around with this, the results have varied by several dB. I have used the published values rather than gsignal for both precision with other modules using this code and reliance on as few packages as possible.

**Welch Compress**

The Welch compression of the power spectra using the published PAMGuide version of the code varies from system welch compressions. I have opted for using the published version for consistency. This could use more teasing out.

**End-to-end Calibration**


Note that **calibration values must be RMS.** Thus, implementing the SoundTrap manufacturer's peak-to-peak calibration values will result in noise metrics that are a few dB off.  


____________________
# Tutorial


## Brief Tutorial

This tutorial provides concrete example of setting up and using a series of functions to process sound data. Ultimately, this should be a GUI, but for the meantime, some R knowledge will be required.

## Setting up the sound files and the project

This code is structured such that data are exported to HDF5 files. These large databases are capable of handling multiple types of data and grouping datasets. For our purposes we often evaluate noise level data from multiple instruments across years and locations. Thus, the example code sets up the data assuming a higher project level and hydrophones as the lower level. Within the hydrophone names we have the metrics of interest (e.g. broadband, hybrid-milidecade, and third octave sound levels). These are arranged as datasets. Additionally, the analysis parameters, time stamps, center frequencies of the third octave and hybrid-milidecade bands are saved as independent datasets within each hydrophone.

If you don't want to muck around navigating to the functions directory, open the project and make sure you have the `here` library installed. The rest should just work.

```{r}
rm(list = ls())
library(here)
library(lubridate)
library(stringr)
library(rhdf5)


#Define path to functions that we will use
source(here('NoiseProcessingFxs.r'))
```


## Audio file directory

At present, this must be wav files. Which is a bummer.

This section of code creates a dataframe with the audio file names and locations and the audio start time (assuming UTC). The data locally are typically in one of a few formats, so two examples are provided here. This identifying names will be put into a function.


```{r}
################################################################
# Load list of sound files and start times, assumes all audio files are the 
# same duration
################################################################

# Assumes 
# 1) file names in UTC
# 2) All sample rates are the same

fileLoc = 'D:\\RECORDINGS\\SoundtrapST5987'

# Jasco Pattern
nameStringPattern = '\\d{8}T\\d{6}.\\d{3}'
nameStringFormat =" %Y%m%dT%H%M%OS"

# Soundtrap Pattern
nameStringPattern = '\\d{12}'
nameStringFormat ='%y%m%d%H%M%OS'
lubradiateFormat ="ymdHMS"

# Set up the audiodataframe 
audioData = createAudioDataframe(fileLoc,nameStringPattern = nameStringPattern,
                                 lubradiateFormat = lubradiateFormat)


```


## Analysis Parameters

The following section of code sets up the analysis parameters. These are standard metrics in acoustics such as sample rate (fs), FFT length (samples, sorry Marie), and the averaging duration (welch) in seconds. Additional parameters are the option of removing the DC offset of the file once the sound is loaded by subtracting the mean. This was in the author's original version, but here is set to an option (recommended to leave as the default: TRUE). In the parameters the 'metrics' value is a list of which of the included calculations should be returned. For all metric included, set value to 'all'. Othwise options are 'hybrid', 'broadbabd', 'decadeband', and 'thirdoct', and 'psd'. 

```{r}
#Which metrics to calculate, case sensitive but order agnostic
prms$metrics = t(list('hybrid', 'broadband', 'decadeband', 'thirdoct', 'psd'))
```

The calibration value can either be end-to-end in decibels or a frequency response. Here I've simulated a frequency response calibration. This could otherwise be uploaded from a CSV or similar. Note that if the calibration values are less than the Nyquist frequency, the user needs to add a value at fs/2 to their calibration dataframe. The last section creates windowing values (vector) and defines alpha in accordance with the functions published in Merchant et al. (2014).

```{r}
################################################################
# Set up FFT metrics
################################################################
# Set up needed parameters

fileLoc = 'D:\\RECORDINGS\\ADRIFT_001_CENSOR_12kHz'
files <- (list.files(fileLoc[1], pattern = "\\.wav$"))

# Pull sample rate from first file, assume consistent (required)
prms = av_media_info(file.path(fileLoc,files[1]))$audio
prms$duration = av_media_info(file.path(fileLoc,files[1]))$duration
colnames(prms)[colnames(prms)=='sample_rate']='Fs' # field standard
prms$N = prms$Fs # fft length in samples
prms$r = 0.5 # overlap in percent
prms$aveSec = 60 # averaging window (seconds)
prms$secSkip = 0 # number of seconds to skip when loading the file (soundtraps...)

# ratio of new to original window lengths in Welch method (avg over aveSec secs)
prms$welch = prms$aveSec*(prms$Fs/prms$N)/(1-prms$r)

# Analysis frequencies (hz). I do not reccomend setting the hi-cut to less than fs/2
prms$lcut =10
prms$hcut = prms$Fs

# Reference pressure, set to 1 for water
prms$pref =1

# remove DC offset by subtracting mean(yy) from yy where yy is the audio signal
prms$rmDC = TRUE

# calibration, either end-to-end or frequency response
prms$freqCal= -175.5

# simulate frequency response (load from CSV or similar)
freqResp = data.frame(f = round(10:prms$Fs/2))
freqResp$E2E = seq(-185, -165, length.out = nrow(freqResp))
prms$freqCal= list(freqResp)

# Date the analysis was run
prms$DateRun = now(tzone = "UTC")


########################################################
# Create window functions as defined by NM 2014
####################################################

w = windowFunctions('hann', prms)[[1]]
prms$alpha = windowFunctions('hann', prms)[[1]]

# Which metrics to calculate, case sensitive, order does not matter
prms$metrics = t(list('hybrid', 'broadband', 'decadeband', 'thirdoct', 'psd'))
```

## Database Initialization

The following section of code sets up the HDF5 database including the study name and writes some of the data that we've already established (parameters and the audio data frame). These could be useful in troubleshooting the data after it has been processed.

The other important thing that happens in this bit of code is estimating the total number of rows that will ultimately be in the dataset. If you are planning on adding to the dataset at a later date (e.g. an instrument comes ashore, is refurbished, & then returned to the field) then you need to ensure that the length is greater or equal to your ultimate length. In this case I know that I will not be adding to this dataset so I'm using the length of the time in minutes (**tt**) later on.

```{r}
# Database name and instrument, adding current computer time for debugging
ProjName = paste0("AcousticStudy",format(Sys.time(), "%Y-%m-%d%H%M%S"), '.h5')
instrumentName = "ST5987"

# Create the  hdf5 file and fill out meta and audiofiles
h5createFile(ProjName)

# create group for location 1 (here we are going with just the hydrophone serial number). 
h5createGroup(ProjName, instrumentName)

# Write the parameters database to the instrument as well as the dataframe containing all the files used
h5write(
  prms,
  file =ProjName,
  paste(instrumentName,"Parms",sep="/"))

h5write(
  audioData$FileName,
  file = ProjName,
  paste(instrumentName,"Files",sep="/"))

##########################################################################
# Guesstimate total duration of the database, this should be revisited but works
# for now. Note this *must* be as long or longer than your intended data. HDF5 is preallocated.
##########################################################################
# figure out maximum dimensions of the output data
timesAll = seq(floor_date(min(audioData$StartTime), 'minute'),
               ceiling_date(max(audioData$EndTime), 'minute'),
               by = paste(prms$welch/2, 'sec'))

```

## Create the Noise Metrics

The system runs by loading each audio file, calculating the PSS using the calcPss function then summing over different bands to calculate third octave, decade, hybrid-milidecade bands etc. Existing functions for these metrics and expected names for the plotting functions are:

- **calcHybridMiDecade** - hybrid-milidecade bands - datatype-'hybridMiliDecLevels'
- **calcBroadband**- broadband levels - datatype-'broadbandLevels'
- **calcDecadeBands** - decade band levels- datatype-'decadeLevels'
- **calcCustomBand** - custom band level (define low and high frequency limits) - datatype- not defined suggest lowToHighHzLevels:)
- **calcThirdOctBands** -third octave band levels- datatype- 'thirdOctLevels'
- **calcPsd** - Power Spectral Density-  datatype- 'psdLevels'
- **calcMetrics** Driver function for the above functions. 

The pre-defined names only pertain to hybrid milidecade and theird octave levels because I've written convienience functions to produce figures from these datasets. In reality, you can call them whatever you want and pull the relevent portions of the make figures code and adapt to your specific dataset. 

A custom function has been created to write the data to the database. The dataType is the name of the dataset, newData is the result of the analysis so will be the octave band levels, time, third octave band levels etc. Datastart is the index of where, in the dataset, the new data will be written. Thus if the analysis results in in a 6x5000 matrix for 6 minutes of the anlaysis the dataStart will be 1 on the first iteration followed by 7, 13, 19 and so forth. The data are organized such that columns represent frequencies and rows represent time. Max rows is the total number of rows that we expect in the dataset (as above). In my understanding its better to overestimate this, as I've done here. Storage mode is the type of data represented in the dataset. The default is double.

Within the noise analysis loop, this function is used only to write the timestamps to the HDF5 file. Otherwise it has been replaced by the driver functions **writeMetricPrms** and **writeAllMetrics** which take in the metrics calculated for each file and automatically assigns them to the correct sheet in the database.

```{r}
  # Write the timestamps
  writeDataToHDF5(ProjName, instrumentName,
                     dataType='timeUTC', 
                     newData = as.matrix(as.character(tt)), 
                     dataStart=((ii-1)*length(tt))+1,
                     maxRows=length(timesAll), 
                     storagemMode= 'character')
```

Putting the above functions in a loop, we now step through the the \*.wav file in the audiodata dataframe, calculates the PSS then passes that matrix to the measurement functions. On the first iteration of the loop, it calculates the center and frequency limits of the analysis bands given the user-defined frequency limits and writes that to the dataset. After the first iteration, new data are added to the same dataset.

```{r}
######################################################################
# Step through the soundfiles, create spectrogram and write to hdf5 file
#######################################################################
idStart =1
for(ii in 1:nrow(audioData)){

  # Calculate the PSS within the user defined range, time stamps, and frequency
  # vector. Returns list [1] power spectrum over the defined frequency range [2] time
  # stamp for each calcuation(PSS rows) [3] vector of frequencies for each calculation (PSS columns)  
  # Calculate the PSS within the user defined range, time stamps, and frequency
  # vector.
  dataOut = calcPSS(audioData, ii, prms, w)
  
  
  Psstrimmed= dataOut[[1]]
  tt= dataOut[[2]]
  f = dataOut[[3]]
  avPSD= 10*log10(Psstrimmed)
  
  # Calculate the metrics defined earlier from the Power spectrum
  allMetrics<-calcMetrics(prms, Psstrimmed, f, w)
  
  # Figure out how long the resulting output is 
  countLen = length(tt)

  ######################################
  # 1) Write the the frequency information for each of the user defined metics.
  # The driver function (writeMetricPrms) is used to call write writeDataToHDF5 for each of the noise metrics
  ###############################################

   if(isTRUE(writePrmsflag)){
    # write the frequency cetners etc for the calculated metrics
    writeMetricPrms(prms, allMetrics, ProjName, instrumentName)
    
    # only on first iteration
    writePrmsflag=FALSE
  }

  ######################################################
  # 2) Write the timestamps
  #  For each of the audiofiles, record the time of each noise measure. In UTC. Or should be if diligince on user end is implemented
 #############################################################

  writeDataToHDF5(ProjName, instrumentName,
                     dataType='DateTime',
                     newData = as.matrix(as.character(tt)),
                     dataStart=idStart,
                     maxRows<-dataSetLenght,
                     storagemMode<- 'character')
  

 #######################################################################
 # 3) Write the metrics to the HDF5 dataframe
 # Last driver function. Write all of the calculated metricies to the HDF5 file. Users can alteratively use
 # writeDataToHDF5() to write specific metrics.
 ######################################################################
  # Driver function to clean things up. Writes all calculated metrics to approperiate location 
  writeAllMetrics(prms, allMetrics, ProjName, 
                     instrumentName, dataStart=idStart, 
                     dataSetLenght= length(allMetrics$hybFreqs$lowF))

  idStart =idStart+countLen
  
  print(ii)
}

H5Fclose(ProjName)



```

## Create Figures from the HDF5 Database

Now that we have a database with our instrument, how do we get the data out? The following section of code shows one way to do that and create a PSD plot.

```{r}
#######################################################################
# Opening projects and plotting
#######################################################################
rm(list = ls())
library(hdf5r) # I have found this package more intuitive but it cannot write to the databse efficiently

ProjName =  "AcousticStudy.h5"
instrumentName = "ST5987"

#open the files
hdfFile <- h5file(ProjName, "r")
instrument_group <- hdfFile[[instrumentName]]


#################################################
## PSD plot
#################################################

p<-makePSDDist(instrument_group)
print(p)



```

![Rplot01](https://user-images.githubusercontent.com/28478110/221385886-ed93fc74-1203-490a-9415-cbaa2bc776b6.png)


And now to make the hourly plot. The first step is, again, to create the dataframe where the data we use to make the plot is going to be stored.

```{r}

#######################################################################
# Create Hrly fig 
#######################################################################

#open the files
my_file <- h5file(ProjName, "r")
instrument_group <- my_file[[instrumentName]]
p<- makeHourly(instrument_group)

print(p)

```


![heatMap](https://user-images.githubusercontent.com/28478110/221384820-00d52ffd-0762-4109-bad4-20afb404cd9e.png)


## Create and LTSA 

In the last example, we create an LTSA from the hybrid-milidecad bands by taking the median houly sound levels.

```{r}
#######################################################################
# Create LTSA - median hourly
#######################################################################

# This has been turned into a function for convienience 

#open the files
my_file <- h5file(ProjName, "r")
instrument_group <- my_file[[instrumentName]]

# Return a ggplot object for any user required adjustments and printing
p<- makeLTSA(instrument_group, averagingPeriod = '90 min')
print(p)


```
![LTSA](https://user-images.githubusercontent.com/28478110/221433499-e689701f-9480-4014-bff1-fd7a3d2d25a9.png)


