# BCC_defect_classification.py
#
# Developed by: Prashant Dwivedi
"""
This script is designed for use with OVITO (Open Visualization Tool) to classify 
and visualize defects in a BCC crystal structure. The script applies various 
crystallographic analysis methods to identify different defect types such as 
dislocations, vacancies, twins, planar faults, and surfaces.

Why this script is important:
--------------------------------
Defect characterization in materials is critical for understanding their 
mechanical and structural properties. Traditional defect analysis can be 
computationally expensive and may require extensive post-processing.

This script streamlines the classification of common defects in BCC structures 
(such as tungsten) by utilizing precomputed atomic properties from OVITO’s 
modifiers (CNA, centrosymmetry, and coordination number). Instead of manually 
analyzing simulation data, this automated process accelerates defect detection, 
providing immediate insights into material behavior and deformation mechanisms.

Key Features:
--------------------------------
- Identifies bulk, surface atoms, vacancies, dislocations, twins, and planar faults.
- Uses coordination number, centrosymmetry, and common neighbor analysis (CNA) 
  to classify defects.
- Assigns a defect type to each atom based on its local atomic environment.
- Applies a color map to visually differentiate defect types in OVITO.

Usage:
--------------------------------
1. Load your atomic configuration file into OVITO.
2. Attach this script to the pipeline as a Python modifier.
3. The script will classify defects and assign colors accordingly.
4. The classified atoms can be further analyzed in OVITO.

Dependencies:
--------------------------------
This script requires the OVITO Python API and NumPy for computations.
"""

from ovito.data import *
from ovito.modifiers import *
import numpy as np

def modify(frame, data):
    # Define defect types
    blk = 0  # Bulk
    srf = 1  # Surface
    vcn = 2  # Vacancy
    dsl = 3  # Dislocation
    twn = 4  # Twin
    plf = 5  # Planar Fault
    els = 6  # Unidentified
    
    # Define cutoff radii and lattice parameters
    alat = 3.16  # Example lattice parameter for BCC Tungsten
    nn2_cutoff = (np.sqrt(2) + 1)/2 * alat
    
    # Step 1: Compute initial properties
    csp = CentroSymmetryModifier(num_neighbors=8)
    data.apply(csp)
    cna = CommonNeighborAnalysisModifier(mode=CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff)
    data.apply(cna)
    coord = CoordinationNumberModifier(cutoff=nn2_cutoff)
    data.apply(coord)
    
    # Prepare neighbor finding (exclude perfect BCC atoms)
    neighbor_finder = CutoffNeighborFinder(nn2_cutoff, data)
    
    # Fetch required properties
    atom_coord = data.particles_['Coordination']
    atom_cna = data.particles_['Structure Type']
    atom_csp = data.particles_['Centrosymmetry']
    
    # Create new defect property
    defect_property = data.particles_.create_property("Defect Type", dtype=int)
    defect_property[:] = blk  # Initialize all particles as bulk
    
    # Helper functions
    def get_neighbors(index):
        return [neigh.index for neigh in neighbor_finder.find(index) 
                if atom_cna[neigh.index] != 3 or atom_coord[neigh.index] != 14]
    
    # 1. Surface Identification
    def is_surface(i):
        if atom_cna[i] != 3:
            if atom_coord[i] <= 11:
                defect_property[i] = srf
                return True
            # Check if neighbors are surface atoms
            neighbors = get_neighbors(i)
            count = 0
            for n in neighbors:
                if defect_property[n] == srf or (atom_cna[n] != 3 and atom_coord[n] <= 11):
                    count += 1
            if count >= 4:
                defect_property[i] = srf
                return True
        return False
    
    # 2. Dislocation Identification
    def is_dislo(i):
        coord = atom_coord[i]
        if coord >= 12 and coord != 14:
            neighbors = get_neighbors(i)
            nr_14 = sum(atom_coord[n] == 14 for n in neighbors)
            nr_non14 = len(neighbors) - nr_14
            if nr_non14 > nr_14:
                defect_property[i] = dsl
                return True
        elif coord == 14:
            neighbors = get_neighbors(i)
            nr_non14 = sum(atom_coord[n] >= 12 and atom_coord[n] != 14 for n in neighbors)
            nr_14 = sum(atom_coord[n] == 14 for n in neighbors)
            if nr_non14 >= 4 and nr_14 <= 6:
                defect_property[i] = dsl
                return True
        return False
    
    # 3. Vacancy Identification
    def is_vac(i):
        coord = atom_coord[i]
        csp = atom_csp[i]
        # Mono-vacancy checks
        if coord == 13 and csp < 1:
            neighbors = get_neighbors(i)
            nr_13 = sum(atom_coord[n] == 13 and atom_csp[n] > 4 for n in neighbors)
            nr_12 = sum(atom_coord[n] == 12 and atom_csp[n] > 4 for n in neighbors)
            if (nr_13 == 4 or (nr_13 == 2 and nr_12 == 2)):
                defect_property[i] = vcn
                return True
        elif coord == 13 and csp > 4:
            neighbors = get_neighbors(i)
            nr_13_4 = sum(atom_coord[n] == 13 and atom_csp[n] > 4 for n in neighbors)
            if nr_13_4 >= 4:
                defect_property[i] = vcn
                return True
        # Di-vacancy checks
        elif coord == 12 and csp > 4:
            neighbors = get_neighbors(i)
            nr_12_4 = sum(atom_coord[n] == 12 and atom_csp[n] > 4 for n in neighbors)
            nr_13_4 = sum(atom_coord[n] == 13 and atom_csp[n] > 4 for n in neighbors)
            if nr_12_4 == 1 and nr_13_4 == 4:
                defect_property[i] = vcn
                return True
        return False
    
    # 4. Twin Identification
    def is_twin(i):
        coord = atom_coord[i]
        csp = atom_csp[i]
        if coord == 13 and csp > 4.5:
            neighbors = get_neighbors(i)
            nr_13 = sum(atom_coord[n] == 13 for n in neighbors)
            nr_14 = sum(atom_coord[n] == 14 for n in neighbors)
            if nr_13 == 5 and nr_14 == 2:
                defect_property[i] = twn
                return True
        elif coord == 14 and csp > 8:
            defect_property[i] = twn
            return True
        elif coord == 14:
            neighbors = get_neighbors(i)
            nr_13 = sum(atom_coord[n] == 13 for n in neighbors)
            nr_14 = sum(atom_coord[n] == 14 for n in neighbors)
            if nr_13 >= 4 and nr_14 >= 2:
                defect_property[i] = twn
                return True
        return False
    
    # 5. Planar Fault Identification
    def is_planarfault(i):
        coord = atom_coord[i]
        if atom_cna[i] != 3:
            neighbors = get_neighbors(i)
            nr_12 = sum(atom_coord[n] == 12 for n in neighbors)
            nr_13 = sum(atom_coord[n] == 13 for n in neighbors)
            if coord == 12:
                if (nr_12 >= 3 and nr_12 <= 6) and (nr_13 >= 7 and nr_13 <= 9):
                    defect_property[i] = plf
                    return True
            elif coord == 13:
                if nr_12 + nr_13 == 9 and nr_13 >= 7:
                    defect_property[i] = plf
                    return True
        return False
    
    # Step 2: Identify defects
    for i in range(len(atom_coord)):
        if atom_cna[i] == 3 and atom_coord[i] == 14:
            defect_property[i] = blk  # Perfect BCC
        else:
            if (is_surface(i) or is_dislo(i) or is_vac(i) or is_twin(i) or is_planarfault(i)):
                continue
            else:
                defect_property[i] = els  # Unidentified
    
    # Step 3: Visualization and Selection
    color_map = {
        blk: (0.8, 0.8, 0.8),  # Gray - Bulk
        srf: (0.0, 0.0, 1.0),  # Blue - Surface
        vcn: (1.0, 0.0, 0.0),  # Red - Vacancy
        dsl: (0.0, 1.0, 0.0),  # Green - Dislocation
        twn: (1.0, 1.0, 0.0),  # Yellow - Twin
        plf: (1.0, 0.5, 0.0),  # Orange - Planar Fault
        els: (0.5, 0.5, 0.5)   # Dark Gray - Unidentified
    }
    
    # Assign colors based on defect type
    color_property = np.array([color_map[defect] for defect in defect_property])
    data.particles_.create_property('Color', data=color_property)
    print("Defect classification completed!")
