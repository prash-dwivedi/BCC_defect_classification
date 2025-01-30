# BCC Defect Classification Script

## Overview
This script is designed for use with **OVITO (Open Visualization Tool)** to classify and visualize **defects in a BCC crystal structure**. It automatically detects and labels **bulk atoms, vacancies, dislocations, twins, planar faults, and surfaces** based on atomic coordination, centrosymmetry, and common neighbor analysis.

By automating defect classification, this script enhances efficiency in **dislocation analysis, defect evolution studies, and mechanical behavior investigations** in materials such as tungsten.

## Citation
If you use this script in your research, please cite the following articles:

- DOI: [10.1016/j.jnucmat.2024.155289](https://doi.org/10.1016/j.jnucmat.2024.155289)  
- DOI: [10.1016/j.jnucmat.2024.155042](https://doi.org/10.1016/j.jnucmat.2024.155042)  

Proper citation ensures recognition of the research that contributed to this work.

## Usage

1. **Open OVITO** and load your atomic structure file.
2. **Attach this script** as a Python modifier in the pipeline.
3. **The script will classify defects** and assign colors accordingly.
4. **Defects are visualized** with the following color mapping:
   - **Gray** - Bulk atoms
   - **Blue** - Surface atoms
   - **Red** - Vacancies
   - **Green** - Dislocations
   - **Yellow** - Twins
   - **Orange** - Planar Faults
   - **Dark Gray** - Unidentified defects

## Requirements

- **OVITO (Open Visualization Tool)**  
- **Python** (for scripting in OVITO)  
- **NumPy** library  

## How It Works

1. The script reads atomic properties such as **coordination number, centrosymmetry, and common neighbor analysis (CNA)** from the current data frame.
2. It applies **defect classification logic** to identify different defects based on their atomic environment.
3. It assigns a **"Defect Type" property** to each atom and applies **color mapping** for visualization.
4. The classified data is displayed **directly in OVITO** for further analysis.

## License

This project is **open-source** and available under the **MIT License**. You are free to use, modify, and distribute this script while providing proper attribution.

---
