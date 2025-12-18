# Estimate Rashba-like spin splitting parameters within the parabolic-band approximation according to VASP outputs.

## introduction
The script can extract the k-point coordinates and band energies at the valence band maximum (VBM) and conduction band minimum (CBM), as well as those at the nearest high-symmetry k-points to the VBM and CBM. This step determines the Rashba like spin splitting energy $E_R$ and the momentum offset $k_R$. Then the spin splitting coefficient $\lambda_R$ will be calculated by  $\lambda_R=2E_R/k_R$. Note that the empty results will be shown in the following situations:  

(i) The material is metallic.  

(ii) The VBM(CBM) is degenerate.  

(ii) More than two bands are involved at the nearest high-symmetry k-point to the VBM(CBM), in which case the two-band-model may not be suitable.

## requirements 
1. Pymatgen  
   
2. PyProcar  
   
3. numpy

## input args
vasp_output_dir: The path that saves the VASP outputs (at least vasprun.xml, PROCAR, EIGENVAL, KPOINTS files), the default value is the current directory (i.e., "./").

## ouput
A json file will be saved which contains the estimated Rashba-like spin splitting parameters.


band_gap: band_gap of the material  

VBM_kpoint: k-point coordinates (fractional) at the VBM  

VBM_near_high_symm_kpoint: name of the nearest high-symmetry k-point to the VBM  

VBM_Rashba_energy (eV): Rashba-like spin splitting energy $E_R$  

VBM_Rashba_momentum_offset (1/Angst): momentum offset $k_R$  

VBM_Rashba_coefficient (eV Angst): spin splitting coefficient calculated by  $\lambda_R=2E_R/k_R$  

VBM_spin: spin polarization at the VBM  

nearest_VB_spin: spin polarization of the valence band immediately below the VBM at the VBM_kpoint

CBM outputs are similar to VBM.

## An example 
cd example; python ../Rashba_like_spin_splitting_parameters_estimate.py --vasp_output_dir ./
