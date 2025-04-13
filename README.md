# CC Singlecell

A project focused on using single-cell RNA sequencing data (scRNA-seq) and pseudo time to improve colon cancer diagnosis and outcomes.
If you find our work helpful please use the citation below:
>Andrew Willems, Nicholas Panchy, and Tian Hong. Using Single-Cell RNA Sequencing and MicroRNA Targeting Data to Improve Colorectal Cancer Survival Prediction. (2023) Cells 12(2):228

## Project Organization

```bash
CCsinglecell/
├── R/                        # Refactored, optimized R package code
│   ├── cell_dataset_builder.R
│   ├── cox_model.R
│   ├── mad_calculator.R
│   ├── etc...
├── man/                      # Documentation (auto-generated)
├── vignettes/                # Usage tutorials
├── tests/                    # Unit tests
├── inst/
│   ├── bin/                  # Command-line tools
│   ├── extdata/              # Example data
│   └── legacy/               # Original code for reproducibility
│       ├── Code/             # Original code files
│       ├── Data/             # Original data files
│       └── paper.Rmd         # Original analysis
├── data/                     # Package datasets
├── DESCRIPTION               # Package metadata
├── NAMESPACE                 # Export declarations
└── README.md                 # Package documentation
```
