library(hdf5r)
library(rhdf5)
ProjName = 'C:/Users/kaitlin.palmer/Downloads/ICES-Continuous-Underwater-Noise-format (1)/Sample.h5'


#Group name
H5open(ProjName)

groupName<-'Meta2'

# create group for location 1 (here we are going with just the hydrophone serial number)
h5createGroup(ProjName, groupName)

##########################################
# Use prms to populate the metadata
############################################

## AveragingTime
h5createDataset(file = ProjName,
                dataset = paste(groupName,'AveragingTime',sep="/"),
                dims = 1,
                storage.mode = 'integer')
# Dataset has been created (or indexed)-populate with new information
h5write(
  prms$aveSec,
  file = ProjName,
  name = paste(groupName,'AveragingTime',sep="/"),
  start = 1, count =1)

## CalibrationDateTime
h5createDataset(file = ProjName,
                dataset = paste(groupName,'CalibrationDateTime',sep="/"),
                dims = 1,
                storage.mode = 'character')

# Dataset has been created (or indexed)-populate with new information
h5write(
  "0000-00-00 00:00:00",
  file = ProjName,
  name = paste(groupName,'CalibrationDateTime',sep="/"),
  start = 1, count =1)

## CalibrationProcedure
h5createDataset(file = ProjName,
                dataset = paste(groupName,'CalibrationProcedure',sep="/"),
                dims = 1, storage.mode = 'character')

# Dataset has been created (or indexed)-populate with new information
h5write(
  "Soundtrap Calibrations",
  file = ProjName,
  name = paste(groupName,'CalibrationProcedure',sep="/"),
  start = 1, count =1)

## Channel Count
h5createDataset(file = ProjName,
                dataset = paste(groupName,'ChannelCount',sep="/"),
                dims = 1, storage.mode = 'integer')

# Dataset has been created (or indexed)-populate with new information
h5write(
  prms$channels,
  file = ProjName,
  name = paste(groupName,'ChannelCount',sep="/"),
  start = 1, count =1)

## Comments
h5createDataset(file = ProjName,
                dataset = paste(groupName,'Comments',sep="/"),
                dims = 1, storage.mode = 'character')

# Dataset has been created (or indexed)-populate with new information
h5write(
  prms$Comments,
  file = ProjName,
  name = paste(groupName,'Comments',sep="/"),
  start = 1, count =1)

## DataUUID
h5createDataset(file = ProjName,
                dataset = paste(groupName,'DataUUID',sep="/"),
                dims = 1, storage.mode = 'character')

# Dataset has been created (or indexed)-populate with new information
h5write(
  'Matlab Fx available',
  file = ProjName,
  name = paste(groupName,'DataUUID',sep="/"),
  start = 1, count =1)

## DatasetVersion
h5createDataset(file = ProjName,
                dataset = paste(groupName,'DatasetVersion',sep="/"),
                dims = 1, storage.mode = 'character')

# Dataset has been created (or indexed)-populate with new information
h5write(
  '1.0',
  file = ProjName,
  name = paste(groupName,'DatasetVersion',sep="/"),
  start = 1, count =1)

## Frequency Count
h5createDataset(file = ProjName,
                dataset = paste(groupName,'FrequencyCount',sep="/"),
                dims = 1, storage.mode = 'integer')

# Dataset has been created (or indexed)-populate with new information
h5write(
  length(thirdOctF),
  file = ProjName,
  name = paste(groupName,'FrequencyCount',sep="/"),
  start = 1, count =1)

## Frequency Index
h5createDataset(file = ProjName,
                dataset = paste(groupName,'FrequencyIndex',sep="/"),
                dims = 1, storage.mode = 'double')

# Dataset has been created (or indexed)-populate with new information
h5write(
  length(thirdOctF),
  file = ProjName,
  name = paste(groupName,'FrequencyCount',sep="/"),
  start = 1, count =1)

## Frequency Unit
h5createDataset(file = ProjName,
                dataset = paste(groupName,'FrequencyUnit',sep="/"),
                dims = 1, storage.mode = 'character')

# Dataset has been created (or indexed)-populate with new information
h5write(
  'Hz',
  file = ProjName,
  name = paste(groupName,'FrequencyUnit',sep="/"),
  start = 1, count =1)


## HydrophoneSerialNumber
h5createDataset(file = ProjName,
                dataset = paste(groupName,'HydrophoneSerialNumber',sep="/"),
                dims = prms$channels, storage.mode = 'character')

# Dataset has been created (or indexed)-populate with new information
h5write(
  'BNKBlabla',
  file = ProjName,
  name = paste(groupName,'HydrophoneSerialNumber',sep="/"),
  start = 1, count =1)


## HydrophoneType
h5createDataset(file = ProjName,
                dataset = paste(groupName,'HydrophoneType',sep="/"),
                dims = prms$channels, storage.mode = 'character')

# Dataset has been created (or indexed)-populate with new information
h5write(
  'HTI',
  file = ProjName,
  name = paste(groupName,'HydrophoneType',sep="/"),
  start = 1, count =1)

## MeasurementHeight
h5createDataset(file = ProjName,
                dataset = paste(groupName,'MeasurementHeight',sep="/"),
                dims = prms$channels, storage.mode = 'double')

# Dataset has been created (or indexed)-populate with new information
h5write(
  prms$MeasurementHeight,
  file = ProjName,
  name = paste(groupName,'MeasurementHeight',sep="/"),
  start = 1, count = prms$channels)


## MeasurementPurpose
h5createDataset(file = ProjName,
                dataset = paste(groupName,'MeasurementPurpose',sep="/"),
                dims = 1, storage.mode = 'character')

# Dataset has been created (or indexed)-populate with new information
h5write(
  'Measuring Stuff!',
  file = ProjName,
  name = paste(groupName,'MeasurementPurpose',sep="/"),
  start = 1, count =1)



## MeasurementSetup
h5createDataset(file = ProjName,
                dataset = paste(groupName,'MeasurementSetup',sep="/"),
                dims = 1, storage.mode = 'character')

# Dataset has been created (or indexed)-populate with new information
h5write(
  'Dipped Hydrophone',
  file = ProjName,
  name = paste(groupName,'MeasurementSetup',sep="/"),
  start = 1, count =1)

## MeasurementTotalNo
h5createDataset(file = ProjName,
                dataset = paste(groupName,'MeasurementTotalNo',sep="/"),
                dims = 1, storage.mode = 'integer')

# Dataset has been created (or indexed)-populate with new information
h5write(
  9999,
  file = ProjName,
  name = paste(groupName,'MeasurementTotalNo',sep="/"),
  start = 1, count =1)

## MeasurementUnit
h5createDataset(file = ProjName,
                dataset = paste(groupName,'MeasurementUnit',sep="/"),
                dims = 1, storage.mode = 'character')

# Dataset has been created (or indexed)-populate with new information
h5write(
  'SPL',
  file = ProjName,
  name = paste(groupName,'MeasurementUnit',sep="/"),
  start = 1, count =1)


## ProcessingAlgorithm
h5createDataset(file = ProjName,
                dataset = paste(groupName,'ProcessingAlgorithm',sep="/"),
                dims = 1, storage.mode = 'character')

# Dataset has been created (or indexed)-populate with new information
h5write(
  'BIAS',
  file = ProjName,
  name = paste(groupName,'ProcessingAlgorithm',sep="/"),
  start = 1, count =1)


## RecorderSerialNumber
h5createDataset(file = ProjName,
                dataset = paste(groupName,'RecorderSerialNumber',sep="/"),
                dims = 1, storage.mode = 'character')

# Dataset has been created (or indexed)-populate with new information
h5write(
  'Soundtrap',
  file = ProjName,
  name = paste(groupName,'RecorderSerialNumber',sep="/"),
  start = 1, count =1)

## RigDesign
h5createDataset(file = ProjName,
                dataset = paste(groupName,'RigDesign',sep="/"),
                dims = 1, storage.mode = 'character')

# Dataset has been created (or indexed)-populate with new information
h5write(
  'Soundtrap',
  file = ProjName,
  name = paste(groupName,'RigDesign',sep="/"),
  start = 1, count =1)


