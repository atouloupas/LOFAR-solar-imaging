#!/bin/bash

# Set number of threads
export OPENBLAS_NUM_THREADS=1

ANT="HBA"
DATE="20230310"

# Define the folder containing the MS files and the source database
DATA_FOLDER="/project/lofarsw/Data/EVENT_LBA/${ANT}_files/${DATE}"
MS_FILES=$(ls -d ${DATA_FOLDER}/*.MS)

# # Copy MS files to convert from raw to proper MS format (otherwise taql commands don't work)
# for MS_FILE in $MS_FILES; do
#     DP3 msin=$MS_FILE \
#         msout=$MS_FILE.CPY \
#         steps=[]
# done

# MS_FILES_CPY=$(ls -d ${DATA_FOLDER}/*.CPY)

# for MS_FILE in $MS_FILES_CPY; do
#     if [[ -e "$MS_FILE" ]]; then
#         # Update columns
#         taql update $MS_FILE set FLAG_ROW=false
#         taql update $MS_FILE set WEIGHT_SPECTRUM=1
#         taql update $MS_FILE set FLAG=false
#         taql update $MS_FILE set FLAG=true where ANTENNA1=ANTENNA2

#         # Step 1: Average
#         echo "Running averager step..."
#         DP3 msin=$MS_FILE \
#             steps=["averager"] \
#             msout=$MS_FILE.AVG \
#             averager.freqresolution=500000 \
#             averager.timeresolution=60
#     fi
# done

# Merge the averaged Sun files
DP3 msin=[$DATA_FOLDER/L883050_SAP000*.AVG] \
    steps=[] \
    msout="${ANT}_joined_sun_${DATE}.MS"

# Merge the averaged Cas A files
# DP3 msin=[$DATA_FOLDER/*SAP001*.AVG] \
#     steps=[] \
#     msout="${ANT}_joined_cas_${DATE}.MS"

echo "Merging completed."
