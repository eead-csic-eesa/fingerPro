# Evolution of the fingerPro Package

Comparative analysis of the `fingerPro` package versions back to the original standalone C program. Each phase documents the state, features, and development metrics of the package, presented in reverse chronological order.

---

## fingerPro 2.1
**API Refinement and CRAN Preparation**

This version focuses on standardizing the package for CRAN, refining the public API, and improving data ingestion workflows.

### Key Features
* **Parameter Standardization**: The `seed` parameter is renamed to `rng_init` across the entire package (`unmix`, `CR`, `CTS_explore`, `CTS_select`, `CI`, etc.) to clearly indicate its role in initializing the random generator.
* **Semantic Renaming**: 
    * `CTS_seeds()` is now `CTS_explore()` to better describe the exploration of tracer combinations.
    * `CTS_error()` is now `CTS_select()`, featuring an improved interface that returns a filtered dataset.
* **New Utility Functions**: 
    * `read_database()`: A centralized entry point that handles file reading, validation, mixture selection, and zero-value detection.
    * `validate_results()`: Diagnostic function for post-unmixing analysis.
* **Function Visibility**: Internal validation functions are now hidden from the public API using `@keywords internal`.

---

## fingerPro 2.0
**Refinement and Variability Propagation**

A major development cycle aimed at polishing the package for CRAN. While the initial objective was to include the Conservative Balance (`CB`) method, the project evolved into a deep revision of almost all existing functions and documentation.

### Key Features
* **Advanced Unmixing API**: Supports `variability` (SD vs SEM), Linear Variability Propagation (`lvp`), and Latin Hypercube `resolution`.
* **C++ Backend Expansion**: Core C++ code includes `unmix_c_lvp` and variability correction helpers (`correct_sem` and `correct_sd`).
* **Technical Metrics**: 
    * **900 lines of code** added or modified.
    * **500 lines** of new documentation and help files.
* **Workflow**: Transitioned to a Pull-Request (PR) based development model to facilitate collaborative updates.

---

## fingerPro 1.4
**Restructuring and Standardization

This phase modernized the codebase by gathering dispersed functions and adopting consistent naming conventions.

### Key Features
* **Naming Convention**: Transitioned all files and functions to `snake_case` (e.g., `box_plot.R`, `DFA_test.R`).
* **Method Generalization**: `CR` and `CTS` logic rewritten to support any number of sources (previously limited to 3 or 4).
* **New Methods**: Introduced Conservative Balance (`CB_method`) for isotopic tracers (initial implementation).
* **Technical Metrics**: 
    * **2,400 lines of code** added or modified.
    * **560 lines** of new documentation.

---

## fingerPro 1.3
**Major Expansion: Selection and Logic**

A significant growth phase that introduced the core tracer selection framework.

### Key Features
* **Codebase Size**: Total of **3,200 lines of R code** (800 of which are documentation) and **1,620 lines of C code**.
* **Tracer Selection Suite**: Introduced Consensus Ranking (`CR`), Consistency Error (`CTS`), and the Conservativeness Index (`CI`).
* **Analytic Solving**: Added C++ support for `least_squares_c` and "triangles" methods.
* **Pre-calculated Data**: Added support for pre-averaged data (mean/SD/n) via the `Means` parameter.

---

## fingerPro 1.2 / 1.1
**Maintenance and Metadata**
* **1.2**: Updated `README.md` to include CRAN/Zenodo badges and formal DOI citations.
* **1.1**: Maintenance release; identical to the 1.0 code base.

---

## fingerPro 1.0
**The Full R Package Transition**

This version marked the transformation from a single-function bridge into a comprehensive statistical toolset.
* **Algorithm**: New optimization using Latin Hypercube sampling with iterative variability corrections.
* **Statistical Testing**: Built-in functions for Kruskal-Wallis (KW) and Discriminant Function Analysis (DFA).
* **Graphic Suite**: Implemented a wide range of plots using `ggplot2`, including PCA biplots, LDA plots, and contribution density plots.

---

## fingerPro 0.9
**The Rcpp Bridge**

The first iteration as an R package, originally named `fingerprinting`. It provided a single Rcpp wrapper around the original C algorithm, implementing a Goodness of Fit (GOF).

---

## fingerPro-C 0.8
**Standalone C Program**

The original implementation of the unmixing logic.
* **Interface**: Command-line execution using `input.txt` and `solution.txt`.
* **Dependencies**: GNU Scientific Library (GSL).
* **Algorithm**: Random weight generation (sorted uniforms) with quadratic GOF optimization.
