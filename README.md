# EmbryonicCollectiveExpressionCore

A Julia package providing core analysis tools for Ascidian single-cell RNA-seq data. This library includes clustering routines, empirical statistic computations, and maximum-entropy model fitting. The main package, demonstrating in detail how to use these tools is at [https://github.com/Dom-Skinner/EmbryonicCollectiveExpression](https://github.com/Dom-Skinner/EmbryonicCollectiveExpression)

---

## Prerequisites

1. **Julia** (tested on 1.10.4)
2. **Python** 3.x (for SciKit‑Learn integration)
3. **Sanity** executable in `bin/` (Linux: `Sanity`; macOS: `Sanity_macOS`)

---

## Installation

### 1. Clone the repository

```bash
# Create the dev directory if it doesn't exist
mkdir -p ~/.julia/dev/

# Clone into Julia's development folder
git clone https://github.com/Dom-Skinner/EmbryonicCollectiveExpressionCore.git ~/.julia/dev/EmbryonicCollectiveExpressionCore
```

### 2. Register the package in Julia

Launch the Julia REPL and run:

```julia
using Pkg
Pkg.activate(homedir()*"/.julia/dev/EmbryonicCollectiveExpressionCore")
Pkg.develop(path=homedir()*"/.julia/dev/EmbryonicCollectiveExpressionCore")
```

### 3. Install Python dependencies

Ensure you have Python 3 installed, then:

```bash
python -m pip install --upgrade scikit-learn
```

#### Configure PyCall to use your Python

Before building PyCall, set the `PYTHON` environment variable to your Python executable:

```bash
export PYTHON=/full/path/to/python
```

In the Julia REPL:

```julia
using Pkg
Pkg.build("PyCall")
using PyCall
println(PyCall.python)  # Should print the path you set
pyimport("sklearn.svm")  # Verify that SciKit-Learn loads correctly
```

If you encounter issues, refer to the [PyCall documentation](https://github.com/JuliaPy/PyCall.jl#readme).

### 4. Install Sanity executable

Download (or compile your own) `Sanity` binaries from:

```
https://github.com/jmbreda/Sanity/tree/master/bin
```

Place the appropriate executable in this repo’s `bin/` directory:

* `Sanity` for Linux
* `Sanity_macOS` for macOS (this is the non-linux default)

---

## Usage

After installation, you can use in your Julia scripts or REPL:

```julia
using EmbryonicCollectiveExpressionCore

# Example: run PCA
GE = load("path_to_data.h5")
nPC = 3
PC_GE  = PCA_GeneExpression(GE,nPC,true,false,true) 
```

Refer to main package for detailed examples.

---

This package was last tested on Julia **1.10.4**.

---