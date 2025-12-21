# Estimate Rashba-like spin splitting parameters within the parabolic-band approximation according to VASP outputs.


## Introduction
The script estimates Rashba-like spin splitting parameters by the following two steps. First, it extracts the k-point coordinates and band energies at the valence band maximum (VBM) and conduction band minimum (CBM), as well as those at the nearest high-symmetry k-points to the VBM and CBM. From these data, the Rashba-like spin splitting energy $E_R$ and the momentum offset $k_R$ can be determined. Then it calculates spin splitting coefficient $\lambda_R$ by $\lambda_R=2E_R/k_R$. Note that this script is only applicable to two-band models. Empty results will be returned in the following cases:  

(i) The system is metallic.  

(ii) The VBM(CBM) is degenerate.  

(iii) The degeneracy at the nearest high-symmetry k-point to the VBM(CBM) is not 2.

(iv) More than two bands are involved along the relevant k-path (from the VBM(CBM) k-point to its nearest high-symmetry k-point), in which case the two-band-model may not be suitable.


## Requirements 
1. Pymatgen  
   
2. PyProcar  
   
3. Numpy


## Input args
vasp_output_dir: The path that saves the VASP outputs (at least vasprun.xml, PROCAR, EIGENVAL, KPOINTS files), the default value is the current directory (i.e., "./").

neighbor_bands_factor : The factor that multiplies the Rashba-like spin-splitting energy to determine the number of bands involved along the relevant k-path. The default value is 0.5.


## Output
A json file named "Rashba_parameters_results.json" will be saved which contains the following results:

1. band_gap: band_gap of the system  

2. VBM_kpoint: k-point coordinates (fractional) at the VBM  

3. VBM_near_high_symm_kpoint: name of the nearest high-symmetry k-point to the VBM  

4. VBM_Rashba_energy (eV): Rashba-like spin splitting energy $E_R$  

5. VBM_Rashba_momentum_offset (1/Angst): momentum offset $k_R$  

6. VBM_Rashba_coefficient (eV Angst): spin splitting coefficient calculated by $\lambda_R=2E_R/k_R$  

7. VBM_spin: spin polarization at the VBM  

8. nearest_VB_spin: spin polarization of the valence band immediately below the VBM at the VBM_kpoint

9. Similar 2-8 outputs for the CBM 


## Attention
1. This script is applicable only to two-band systems. Therefore, the electronic band structure should be inspected in advance to verify that the two-band model is appropriate for the system under consideration.

2. An energy tolerance of 0.001 eV is adopted for identifying band degeneracy.
   
3. The "nearest high-symmetry k-point to the VBM(CBM)" refers to the high-symmetry point along the k-point path specified in the "KPOINTS" file that is closest, in terms of k-point index, to the k-point at which the VBM(CBM) occurs.


## Example 
An example is provided in the "example" folder. To run it, execute the following commands in the terminal:
```
   cd example; cd input
   python ../../Rashba_like_spin_splitting_parameters_estimate.py --vasp_output_dir ./
```
The results will be saved to the file "Rashba_parameters_results.json".

## Additional files
All 168,000 generated structures mentioned in our work titled "High-Throughput Inverse Design of Two-Dimensional Ferroelectric Semiconductors with Giant Rashba-like Splitting Energy" are provided in "generated_structures" folder.