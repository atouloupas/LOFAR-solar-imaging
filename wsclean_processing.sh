#!/bin/bash

# Set number of threads
export OPENBLAS_NUM_THREADS=1

# Parameters
FILE="HBA_files/20230311/L883060_SAP000_SB230_uv.MS.CPY.AVG"
# FILE="${ANT}_joined_cas_${DATE}.MS"
TYPE="calibrated"
INTERVALS_OUT=15
TIME_SLOT=0
# CHAN_MIN=19
# CHAN_MAX=20

BASE_NAME="${FILE##*/}"
BASE_NAME="${BASE_NAME%%.*}"
OUTPUT_NAME="images_fits/"$BASE_NAME

# Assign data column based on image type
if [ "$TYPE" = "model" ]; then
    DATA_COLUMN=MODEL_DATA
elif [ "$TYPE" = "calibrated" ]; then
    DATA_COLUMN=CALIBRATED_DATA
else
    DATA_COLUMN=DATA
fi

FILE_NAME="${BASE_NAME}-${TYPE}"
mkdir -p "images_fits/${FILE_NAME}/"

OUTPUT_NAME="images_fits/${FILE_NAME}/${FILE_NAME}"

# Predicted model data
wsclean -name $OUTPUT_NAME \
    -mem 50 \
    -weight briggs 0.2 \
    -scale 30asec \
    -size 512 512 \
    -pol I \
    -multiscale \
    -data-column $DATA_COLUMN \
    -niter 5000 \
    -no-reorder \
    -no-update-model-required \
    -mgain 0.8 \
    -maxuvw-m 3000 \
    -auto-mask 6 \
    -auto-threshold 2 \
    -fit-beam \
    -make-psf \
    -intervals-out $INTERVALS_OUT \
    $FILE \
    # -channel-range $CHAN_MIN $CHAN_MAX

echo "WSClean processing completed."

if [ $INTERVALS_OUT = 1 ]; then
    echo "Created FITS file: ${OUTPUT_NAME}-image.fits"
    python3 helioprojective_plot.py $OUTPUT_NAME"-image.fits"
else
    echo "Created folder of FITS files: "
    python3 reproject.py $OUTPUT_NAME-t*-image.fits images_png/$FILE_NAME/
fi