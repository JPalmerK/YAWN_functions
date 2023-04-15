


# This script contains relevent references for working with h5df files that will
# be used to create the noise functions. Handy reference for the futre

library(hdf5r)

# Database name and instrument
ProjName = paste0("Adrift",format(Sys.time(), "%Y-%m-%d%H%M%S"), '.h5')
ProjName ='Adrift2023-02-19122754.h5'
instrumentName = "ST5987"
dataType = 'hybridMiliDec' #which measurement metrics

# Create a new h5 file
Adrift.h5 <- H5File$new(ProjName, mode = "w")
Adrift.h5

# Create the instrument group
instrument.grp <- Adrift.h5$create_group(instrumentName)
Adrift.h5$ls()

# Add the data group and attributes
instrument.grp[['hybridMiliDec']]<-hybLevels
h5attr(instrument.grp[['hybridMiliDec']], "colnames") <- as.character(
  round(hybFreqs$center))

# Extract tdata and column names
hybridMiliDec_ds <- instrument.grp[["hybridMiliDec"]]
aa= hybridMiliDec_ds[1:5,]
colnames(aa)=h5attr(instrument.grp[['hybridMiliDec']], "colnames") 

# Add to data
hybridMiliDec_ds$dims
hybridMiliDec_ds$maxdims
instrument.grp[['hybridMiliDec']][7:12,]<-hybLevels
hybridMiliDec_ds$dims

# Determine if group and data exists
"ST5987" %in% names(Adrift.h5)
'hybridMiliDec' %in% instrument.grp$ls()$name

# close the file
Adrift.h5$close_all()
h5close(Adrift.h5)

##########################################################################
# For existing file
##########################################################################


#open the files
my_file <- h5file(ProjName, "r")

# open the group
instrument_group = my_file[[instrumentName]]

# get the dimensions
datDims = instrument_group[[dataType]]$dims


# modify the data in the group
instrument.grp[[dataType]][datDims+1: datDims+nrow(hybLevels),]<-hybLevels



# Adrift.h5$close_all()
h5close(Adrift.h5)

Adrift.h5$ls()

################################################
# hdf5r is struggling to write, try rhdf5 again
################################################

CurrentDate = fileDates[1]
prms$DateRun = now(tzone = "UTC")


# Database name and instrument
ProjName = paste0("Adrift",format(Sys.time(), "%Y-%m-%d%H%M%S"), '.h5')
instrumentName = "ST5987"


# Create the  hdf5 file and fill out meta and audiofiles
h5createFile(ProjName)

# create group for location 1
h5createGroup(ProjName, instrumentName)

# list groups
h5ls(ProjName)$name
