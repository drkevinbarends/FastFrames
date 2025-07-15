#!/bin/bash

runs=(
    "tWZOFSR"
    "tWZSFSR"
    "tWZLooseSR"
    "ttZCR"
    "ZZbCR"
)

for param in "${runs[@]}"; do
    export region="$param"

    echo "Running with region: $region"

    rm -rf /eos/user/k/kebarend/tWZ/ntuples/FastFrames/FinalSelection/FourLeptonChannel/$region/*.root

    # Create a temporary config file
    temp_config="FastFrames/temp3LepConfig.yml"
    envsubst < FastFrames/tWZ4LepConfig.yml > "$temp_config"

    # Run Python with the temporary config file
    python3 FastFrames/python/FastFrames.py -c "$temp_config" --step n

    echo "Completed run with PARAMETER: $region"
done