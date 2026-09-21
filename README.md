# multiomics-platform

Python package for programmatic access to the Human Cancer-Associated TMT Proteome Atlas — a curated cross-repository registry of high-multiplexing TMT/TMTpro proteomics studies in human samples.

## About

The Human Cancer-Associated TMT Proteome Atlas (v1.0) integrates 126 datasets from PRIDE, PDC/CPTAC, MassIVE, iProX, GTEx-Proteome and CCLE Proteomics, covering 15,958 annotated samples. This package provides a Python API to search, filter and export curated metadata and quantitative tables.

**Associated publication:**  
Gordeeva A.I. et al. Human Cancer-Associated TMT Proteome Atlas. *Molecular & Cellular Proteomics* (submitted).

**Atlas repository:** https://github.com/arinaatom-cyber/tmt-projects  
**Web interface:** https://human-cancser-tmt-proteome-atlas.streamlit.app/

## Installation

```bash
pip install git+https://github.com/arinaatom-cyber/multiomics-platform.git
```

## Quick start

```python
from proteomics_explorer import ProteomicsExplorer

# Initialize and load project registry
explorer = ProteomicsExplorer()

# List available projects
projects = explorer.list_projects(limit=10)
print(projects)

# Search by gene or protein
result = explorer.search("TP53")
print(result)

# Filter by cancer type
luad = explorer.list_projects(tumor_type="Lung adenocarcinoma")
print(f"LUAD datasets: {len(luad)}")

# Export results
explorer.df.to_csv("results.csv", index=False)
explorer.df.to_excel("results.xlsx", index=False)
```

## Data access

Each dataset in the atlas includes:
- **summary** — project metadata (repository, PMID, TMT scheme, platform)
- **annotation** — channel-level sample annotation (up to 29 fields per TMT channel)
- **result** — quantitative protein matrix with harmonized HGNC identifiers

Project files are available at:  
https://github.com/arinaatom-cyber/tmt-projects/tree/main/Projects

## Requirements

- Python ≥ 3.10
- pandas ≥ 2.0
- openpyxl ≥ 3.1
- requests ≥ 2.31

## Citation

If you use this package or the atlas in your research, please cite:

> Gordeeva A.I. et al. Human Cancer-Associated TMT Proteome Atlas. *Molecular & Cellular Proteomics* (submitted). github.com/arinaatom-cyber/multiomics-platform

## License

MIT
