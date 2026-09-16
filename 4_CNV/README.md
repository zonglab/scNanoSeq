
# Ginkgo CNV Calling Guide

Welcome to the standalone version of Ginkgo for Copy Number Variation (CNV) calling. For more detailed information,   
please visit the [Ginkgo GitHub page](https://github.com/robertaboukhalil/ginkgo/blob/master/README.standalone).

## Getting Started

### 1. Installation

First, you need to install Ginkgo. Simply follow the instructions provided in the README.md file included with Ginkgo. Important: Ensure you update the value of `main` in both `scripts/process.R` and `scripts/reclust.R`.

### 2. Prepare Your Workspace

Create a new directory within the `uploads` directory in your Ginkgo installation:

```bash
mkdir ~/src/ginkgo/uploads/new_run_123
```

Note: The directory name `new_run_123` represents the `$GINKGO_USER_ID`. While the main server typically uses a 20-digit random number, you can use any identifier you prefer.

### 3. Manage Your Data Files

Copy your compressed BED files (output from `bamTobed`) into the new directory:

```bash
cp *.bed.gz ~/src/ginkgo/uploads/new_run_123
```

### 4. List the Cells

Create a file named "list" that includes all the BED files:

```bash
ls ~/src/ginkgo/uploads/new_run_123 | grep .bed.gz$ > ~/src/ginkgo/uploads/new_run_123/list
```

### 5. Configuration File

Set up a configuration file with all necessary options for running Ginkgo. This step can be challenging due to limited documentation, but you can start with an example file:

```bash
cp ~/src/ginkgo/config.example ~/src/ginkgo/uploads/new_run_123/
```

Edit this file with your preferred text editor (e.g., vi):

```bash
vi ~/src/ginkgo/uploads/new_run_123/config.example
```

### 6. Execute Ginkgo

Finally, to run Ginkgo, use the following commands:

```bash
cd ~/src/ginkgo/
./scripts/analyze.sh new_run_123
```

## Additional Information

We've also updated the critical CN calling script, `process.R`, which is available as `modified_process.R` for enhanced performance.
Please note that we have removed the part for clustering and modified the output names to be more ordered and meaningful.
However, this will not be compatible with the downstream analysis of Ginkgo pipeline.
