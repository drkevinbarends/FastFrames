#!/bin/bash

runs=(
    "Preselection"
    "LooseLeptons"
    "TightLeptons"
    "LeptonPt"
    "LeptonCharge"
    "OSSF"
    "ZCandidates"
    "Jets"
    "Bjets"
)

for param in "${runs[@]}"; do
    export region="$param"

    echo "Running with region: $region"

    # Remove existing root files
    rm -f /eos/user/k/kebarend/tWZ/ntuples/FastFrames/Preselection/$region/*.root

    # Create a temporary config file
    temp_config="FastFrames/tempPreselection_config.yml"
    envsubst < FastFrames/preselection_config.yml > "$temp_config"

    # Run Python with the temporary config file
    python3 FastFrames/python/FastFrames.py -c "$temp_config" --step n

    echo "Completed run with PARAMETER: $region"
done