# Autodock Vina-GAG

## Overview

**Autodock Vina-GAG** is an enhanced version of **Autodock Vina-Carb** specifically tailored for the molecular docking of **Glycosaminoglycans (GAGs)**.

GAGs are complex, highly charged biopolymers, and accurately modeling their interactions with proteins is critical in structural biology. This tool significantly improves the performance and accuracy of GAG docking simulations by explicitly introducing **charged interactions** into the scoring function, a feature vital for realistically modeling the behavior of these complex, polyanionic molecules.

---

## Core Features

* **Enhanced Scoring Function:** Modification of the AutoDock Vina scoring function to include explicit terms for charged interactions (electrostatic potential), which are essential for GAG-protein binding thermodynamics.
* **Glycosaminoglycan Support:** Optimized specifically for large, flexible carbohydrate ligands like GAGs.
* **Performance:** Maintains the speed and efficiency of the original AutoDock Vina engine while providing higher accuracy for GAG-related systems.
* **C++ Implementation:** Built on a robust, high-performance C++ codebase.

---

## Installation

As **Autodock Vina-GAG** is a C++-based application, it requires compilation from source. The following steps outline the general process, which is standard for Vina-based projects.

### Prerequisites

You will need the following dependencies installed on your system:

1.  A **C++ Compiler** (e.g., GCC/g++ that supports C++11 or later).
2.  **Make** (or equivalent build system).
3.  **Boost Libraries:** Specifically, the Boost Program Options, Boost Threads, and Boost System libraries.

### Building from Source (Inferred)

1.  **Clone the Repository:**
    ```bash
    git clone [https://github.com/amikasood/VC_GAG.git](https://github.com/amikasood/VC_GAG.git)
    cd VC_GAG
    ```

2.  **Compile the Source:**
    Given the structure, the compilation process likely uses a `Makefile` or similar script within the `build` directory. You will typically run a command such as:
    ```bash
    # This command is based on standard Vina compilation practices.
    make release
    
    # Or navigate into the build directory and use the build script
    # cd build
    # ./build.sh 
    ```
    Upon successful compilation, the executable file, typically named `vina_gag`, should be generated in the main directory or the `bin/` subdirectory.

---

## Usage

**Autodock Vina-GAG** uses a command-line interface and is designed to be compatible with standard PDBQT file formats for molecular docking.

### Input Files

1.  **Receptor File:** The protein structure (e.g., in PDBQT format).
2.  **Ligand File:** The GAG molecule (e.g., in PDBQT format).
3.  **Configuration File** (Recommended): A text file specifying docking parameters.

### Basic Command Example

To perform a docking run, you will typically use the executable along with a configuration file (`config.txt`):

```bash
./vina_gag --config config.txt
