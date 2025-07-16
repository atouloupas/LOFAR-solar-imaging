#!/bin/bash

# Set number of threads
export OPENBLAS_NUM_THREADS=1

# Define the antenna and date
ANT="HBA"
DATE="20230311"

# Define the folder containing the MS files and the source database
SOURCE_DB="Cas-sources.txt"
SOLUTIONS_FOLDER="./solutions_${ANT}"
mkdir -p $SOLUTIONS_FOLDER

# Check if the source database file exists
if [[ ! -f "$SOURCE_DB" ]]; then
    echo "Source database file '$SOURCE_DB' not found. Exiting."
    exit 1
fi

# CAL_FILE="${ANT}_joined_cas_${DATE}.MS"
# SUN_FILE="${ANT}_joined_sun_${DATE}.MS"
CAL_FILE="${ANT}_files/${DATE}/L883062_SAP000_SB190_uv.MS.CPY.AVG"
SUN_FILE="${ANT}_files/${DATE}/L883060_SAP000_SB190_uv.MS.CPY.AVG"
# Extract the filename without the path and extension for creating parmdb file
BASE_NAME="${SUN_FILE##*/}"
BASE_NAME="${BASE_NAME%%.*}"
PARMDB_FILE="solutions_${BASE_NAME}.h5"

# List of tables
tables=("$CAL_FILE" "$SUN_FILE")

# Loop over each table and apply TAQL commands
for table in "${tables[@]}"; do
    taql "update $table set FLAG_ROW=false"
    taql "update $table set WEIGHT_SPECTRUM=1"
    taql "update $table set FLAG=false"
    taql "update $table set FLAG=true where ANTENNA1=ANTENNA2"
done

# Check if the file exists (to handle cases where the glob finds no matches)
if [[ -e "$CAL_FILE" ]]; then
    echo "Processing file: $CAL_FILE"

    # Step 1: Predict
    echo "Running predict step..."
    DP3 steps=["predict"] \
        msin="$CAL_FILE"/ \
        msout=. \
        msout.datacolumn=MODEL_DATA \
        predict.sourcedb="$SOURCE_DB" \
        predict.usebeammodel=true

    if [[ $? -eq 0 ]]; then
        echo "Predict step successful for: $CAL_FILE"
    else
        echo "Predict step failed for: $CAL_FILE"
        break
    fi

    
    # Produce and read outlier station list
    python3 station_outliers.py $SUN_FILE
    mapfile -t stations < stations/stations_${BASE_NAME}.txt

    # Build baseline flag list
    baseline_flags=""
    for station in "${stations[@]}"; do
        baseline_flags+="[$station,*],[*,$station],"
    done
    baseline_flags="${baseline_flags%,}"

    # Step 2: Gain Calibration
    # flag_times.timeofday=[13:13:00..13:15:25]
    echo "Running gaincal step..."
    DP3 steps=[preflagger,flag_times,gaincal] \
        preflagger.baseline=[$baseline_flags,[RS*,RS*]] \
        flag_times.type=preflagger \
        flag_times.timeofday=[] \
        gaincal.caltype=diagonal \
        gaincal.onebeamperpatch=true \
        gaincal.parmdb="$SOLUTIONS_FOLDER/$PARMDB_FILE" \
        msin="$CAL_FILE"/ \
        gaincal.sourcedb="$SOURCE_DB" \
        msout=. \
        msout.datacolumn=CALIBRATED_DATA \
        gaincal.applysolution=true \
        gaincal.usebeammodel=true \
        gaincal.uvlambdamin=120

    if [[ $? -eq 0 ]]; then
        echo "Gaincal step successful for: $CAL_FILE"
    else
        echo "Gaincal step failed for: $CAL_FILE"
    fi


    # Step 3: Applying the Solution
    echo "Applying solution..."
    python3 interpolate_solutions.py "$SOLUTIONS_FOLDER/$PARMDB_FILE" --interpolate --save

    DP3 steps=[applyphase,applyamp,applybeam] \
        applybeam.updateweights=true \
        applyphase.type=applycal \
        applyphase.interpolation=nearest \
        applyphase.updateweights=true \
        applyphase.parmdb="$SOLUTIONS_FOLDER/$PARMDB_FILE" \
        applyphase.correction=phase000 \
        applyphase.soltab=sol000 \
        \
        applyamp.type=applycal \
        applyamp.interpolation=nearest \
        applyamp.updateweights=true \
        applyamp.parmdb="$SOLUTIONS_FOLDER/$PARMDB_FILE" \
        applyamp.correction=amplitude000 \
        applyamp.soltab=sol000 \
        \
        msin="$SUN_FILE" \
        msout=. \
        msout.datacolumn=CALIBRATED_DATA

    if [[ $? -eq 0 ]]; then
        echo "Applying solution step successful for: $SUN_FILE"
    else
        echo "Applying solution step failed for: $SUN_FILE"
    fi

else
    echo "No matching files found."
fi

echo "Calibration process completed."
