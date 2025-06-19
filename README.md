# FastFrames

FastFrames is a lightweight and efficient framework designed for high-performance data processing and analysis. It provides tools to streamline workflows and optimize computational tasks.

## Features

- **High Performance**: Optimized for speed and efficiency.
- **Modular Design**: Easily extendable and customizable.
- **User-Friendly**: Simple APIs for quick integration.

For more information, please visit the FastFrames documentation [link](https://atlas-project-topreconstruction.web.cern.ch/fastframesdocumentation/tutorial/)

## Getting Started

To start using the package, source the `initialSetup.sh` script:

```bash
source initialSetup.sh
```

If you have already done the initial setup and are in a new session or restarting your session, run:

```bash
source setup.sh
```

This will set up the necessary environment variables and dependencies. Then either move your input ntuples into the input folder or adjust the input directory in the following two files:

removeFoldersWithNoNtuples.py
removeNtuplesWithNoRecoTree.py

The reason you need to these files is because there might be files with no "reco" TTree from the TopCPToolkit output. This looks in each file and removes them if there is no "reco" TTree found, then removes the folder if there are no more files. FastFrames cannot handle files with no "reco" TTree. Then run:

```bash
source generate_file_metadata.sh
```

To run the pipeline, source the `runFastFramesPipeline.sh` script:

```bash
source runFastFramesPipeline.sh
```

This will apply the final selection cuts for the tWZ analysis over the various Signal and Control regions and store the output in an ntuple.

**Note:** Before running the scripts, ensure that the file paths inside `setup.sh` and `runFastFramesPipeline.sh` are updated to reflect their correct locations on your machine.

## Creating another branch

To create a new branch, first ensure you are on the branch you want to base it on. Then run:

```bash
git checkout -b <new-branch-name>
```

Replace `<new-branch-name>` with the desired name for your branch. This will create and switch to the new branch.

## Compiling Changes

If you make any changes to the `tWZClass.h` or `tWZClass.cc` files, you must run the `setup.sh` script in the main directory to compile and link your changes:

```bash
source setup.sh
```