import numpy as np
import os
import json
import argparse

from pymatgen.io.vasp import Vasprun
from pymatgen.io.vasp.outputs import Eigenval
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.outputs import Vasprun

try:
    from pyprocar.io.procarparser import ProcarParser
except ImportError:
    raise ImportError(
        "pyprocar is not installed. Please install it using 'pip install pyprocar'."
    )


def get_neighbor_bands_number_at_kpoint(
    eigenvalfile, kpoint_index, to_check_energy, energy_tolerance=0.01
):
    """Get the number of bands neighboring a specific energy at a given k-point."""
    temp_neighbor_bands_number = 0
    kpoint_eigen_energy_list = eigenvalfile.eigenvalues[Spin(1)][kpoint_index, :, 0]
    for x in kpoint_eigen_energy_list:
        if abs(to_check_energy - x) < energy_tolerance:
            temp_neighbor_bands_number += 1
    return temp_neighbor_bands_number


def check_vbm_cbm_high_symmetry_kpoints(band_vasprun):
    """
    Check the VBM and CBM of a semiconductor at high symmetry k-points.
    Args:
        rootpath (str): The path to the directory containing the VASP output files.
    Returns:
        dict: A dictionary containing the VBM and CBM information.
    """

    bs_data = band_vasprun.get_band_structure(line_mode=True)
    lattice = band_vasprun.final_structure.lattice.matrix
    vbm_info = bs_data.get_vbm()
    cbm_info = bs_data.get_cbm()
    band_gap = (bs_data.get_band_gap())["energy"]

    if band_gap < 1e-6:
        print("conductor")
        return None
    else:
        print("semiconductor")

    high_symm_k_names = []
    high_symm_k_idx = []
    for x in bs_data.branches:
        high_symm_k_names.append(x["name"].split("-")[0])
        high_symm_k_idx.append(x["start_index"])
        high_symm_k_names.append(x["name"].split("-")[1])
        high_symm_k_idx.append(x["end_index"])

    temp_vbm_idx = min(
        range(len(high_symm_k_idx)),
        key=lambda i: abs(high_symm_k_idx[i] - vbm_info["kpoint_index"][0]),
    )
    temp_cbm_idx = min(
        range(len(high_symm_k_idx)),
        key=lambda i: abs(high_symm_k_idx[i] - cbm_info["kpoint_index"][0]),
    )
    high_symm_index_vbm = high_symm_k_idx[temp_vbm_idx]
    vbm_high_symm_kpoint = high_symm_k_names[temp_vbm_idx]

    high_symm_index_cbm = high_symm_k_idx[temp_cbm_idx]
    cbm_high_symm_kpoint = high_symm_k_names[temp_cbm_idx]

    vbm_high_symm_kpoint_frac_coord = bs_data.kpoints[high_symm_index_vbm].frac_coords
    cbm_high_symm_kpoint_frac_coord = bs_data.kpoints[high_symm_index_cbm].frac_coords

    vbm_cbm_high_symm_info = {
        "VBM": {
            "kpoint": vbm_info["kpoint"].frac_coords,
            "energy": vbm_info["energy"],
            "high_symm_kpoint": vbm_high_symm_kpoint,
            "high_symm_index": high_symm_index_vbm,
            "high_symm_kpoint_frac_coord": vbm_high_symm_kpoint_frac_coord,
        },
        "CBM": {
            "kpoint": cbm_info["kpoint"].frac_coords,
            "energy": cbm_info["energy"],
            "high_symm_kpoint": cbm_high_symm_kpoint,
            "high_symm_index": high_symm_index_cbm,
            "high_symm_kpoint_frac_coord": cbm_high_symm_kpoint_frac_coord,
        },
    }
    return vbm_cbm_high_symm_info


def check_vbm_cbm_split_semiconductor(rootpath):
    """_summary_

    Args:
        rootpath (str): The path to the directory containing the VASP output files.

    Returns:
        dict: A dictionary containing the VBM and CBM Rashba parameters information.
    """
    band_vasprun = Vasprun(
        os.path.join(rootpath, "vasprun.xml"), parse_projected_eigen=True
    )
    parser_base_scf = ProcarParser()
    parser_base_scf.readFile(os.path.join(rootpath, "PROCAR"), permissive=True)
    eigenvalfile_scf = Eigenval(
        os.path.join(rootpath, "EIGENVAL"), separate_spins=False
    )
    spd = parser_base_scf.spd
    bs_data = band_vasprun.get_band_structure(line_mode=True)
    vbm_info = bs_data.get_vbm()
    cbm_info = bs_data.get_cbm()
    band_gap = (bs_data.get_band_gap())["energy"]
    if band_gap < 1e-6:
        print("The material is conductor, empty dictionary will be return")
        vbm_dict = {}
        cbm_dict = {}
        return {"VBM": vbm_dict, "CBM": cbm_dict}

    vbm_cbm_high_symm_info = check_vbm_cbm_high_symmetry_kpoints(band_vasprun)
    vbm_high_symm_kpoint_index = vbm_cbm_high_symm_info["VBM"]["high_symm_index"]
    cbm_high_symm_kpoint_index = vbm_cbm_high_symm_info["CBM"]["high_symm_index"]

    if len(vbm_info["band_index"][Spin(1)]) > 1:
        print(
            "The degeneracy at VBM is larger than 1, an empty dictionary will be shown"
        )
        vbm_dict = {}

    else:
        vbm_energy = (
            eigenvalfile_scf.eigenvalues[Spin(1)][
                vbm_info["kpoint_index"][0], vbm_info["band_index"][Spin(1)][0]
            ]
        )[0]

        high_symm_energy_vbm = (
            eigenvalfile_scf.eigenvalues[Spin(1)][
                vbm_high_symm_kpoint_index, (vbm_info["band_index"][Spin(1)][0])
            ]
        )[0]

        vbm_neighbor_bands_number = get_neighbor_bands_number_at_kpoint(
            eigenvalfile_scf,
            vbm_high_symm_kpoint_index,
            high_symm_energy_vbm,
            abs(vbm_energy - high_symm_energy_vbm) / 2,
        )

        if vbm_neighbor_bands_number > 2:
            print(
                "Neighbor bands number at VBM nearest high symmetry point is larger than 2, two-band-model is not suitable, an empty dictionary will be shown"
            )
            vbm_dict = {}
        else:
            vbm_S = spd[
                vbm_info["kpoint_index"][0],
                (vbm_info["band_index"][Spin(1)][0]),
                1:,
                -1,
                -1,
            ]
            vbm_S_1 = spd[
                vbm_info["kpoint_index"][0],
                (vbm_info["band_index"][Spin(1)][0] - 1),
                1:,
                -1,
                -1,
            ]

            delta_energy_vbm = abs(float(vbm_energy - high_symm_energy_vbm))

            vbm_kpoint = vbm_info["kpoint"].frac_coords
            delta_k = np.linalg.norm(
                vbm_info["kpoint"].cart_coords
                - bs_data.kpoints[vbm_high_symm_kpoint_index].cart_coords
            )
            Rashba_alpha_vbm = 2 * (delta_energy_vbm) / delta_k
            vbm_dict = {
                "VBM_kpoint": f"{vbm_kpoint}",
                "VBM_near_high_symm_kpoint": vbm_cbm_high_symm_info["VBM"][
                    "high_symm_kpoint"
                ],
                "VBM_Rashba_energy (eV)": delta_energy_vbm,
                "VBM_Rashba_momentum_offset (1/Angst)": delta_k,
                "VBM_Rashba_coefficient (eV Angst)": Rashba_alpha_vbm,
                "VBM_spin": vbm_S.tolist(),
                "nearest_VB_spin": vbm_S_1.tolist(),
            }

    if len(cbm_info["band_index"][Spin(1)]) > 1:
        print(
            "The degeneracy at CBM is larger than 1, an empty dictionary will be shown."
        )
        cbm_dict = {}
    else:
        cbm_energy = (
            eigenvalfile_scf.eigenvalues[Spin(1)][
                cbm_info["kpoint_index"][0], cbm_info["band_index"][Spin(1)][0]
            ]
        )[0]

        high_symm_energy_cbm = (
            eigenvalfile_scf.eigenvalues[Spin(1)][
                cbm_high_symm_kpoint_index, (cbm_info["band_index"][Spin(1)][0])
            ]
        )[0]

        cbm_neighbor_bands_number = get_neighbor_bands_number_at_kpoint(
            eigenvalfile_scf,
            cbm_high_symm_kpoint_index,
            high_symm_energy_cbm,
            abs(cbm_energy - high_symm_energy_cbm) / 2,
        )

        if cbm_neighbor_bands_number > 2:
            print(
                "Neighbor bands number at CBM nearest high symmetry point is larger than 2, two-band-model is not suitable, an empty dictionary will be shown"
            )
            cbm_dict = {}
        else:
            cbm_S = spd[
                cbm_info["kpoint_index"][0],
                (cbm_info["band_index"][Spin(1)][0]),
                1:,
                -1,
                -1,
            ]
            cbm_S_1 = spd[
                cbm_info["kpoint_index"][0],
                (cbm_info["band_index"][Spin(1)][0] + 1),
                1:,
                -1,
                -1,
            ]

            delta_energy_cbm = abs(float(-cbm_energy + high_symm_energy_cbm))

            cbm_kpoint = cbm_info["kpoint"].frac_coords
            delta_k = np.linalg.norm(
                cbm_info["kpoint"].cart_coords
                - bs_data.kpoints[cbm_high_symm_kpoint_index].cart_coords
            )
            Rashba_alpha_cbm = 2 * (delta_energy_cbm) / delta_k
            cbm_dict = {
                "CBM_kpoint": f"{cbm_kpoint}",
                "CBM_near_high_symm_kpoint": vbm_cbm_high_symm_info["CBM"][
                    "high_symm_kpoint"
                ],
                "CBM_Rashba_energy (eV)": delta_energy_cbm,
                "CBM_Rashba_momentum_offset (1/Angst)": delta_k,
                "CBM_Rashba_coefficient (eV Angst)": Rashba_alpha_cbm,
                "CBM_spin": cbm_S.tolist(),
                "nearest_CB_spin": cbm_S_1.tolist(),
            }

    return {"band_gap": band_gap, "VBM": vbm_dict, "CBM": cbm_dict}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Estimate Rashba spin splitting parameters within the parabolic-band approximation according to VASP outputs."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--vasp_output_dir",
        default="./",
        type=lambda x: x.rstrip("/\\"),
        help="where VASP outputs exist (default: current directory)",
    )

    args = parser.parse_args()

    rootpath = args.vasp_output_dir
    result = check_vbm_cbm_split_semiconductor(rootpath)
    print("=======================================================")
    print("Rashba-like spin splitting parameters estimation results:")
    print("band gap:", result["band_gap"])
    print("VBM results:", result["VBM"])
    print("CBM results:", result["CBM"])
    print("\n")
    print("The results will be saved to Rashba_parameters_results.json")
    print("=======================================================")

    to_output_dict = {
        "band_gap": result["band_gap"],
        "VBM": result["VBM"],
        "CBM": result["CBM"],
    }
    with open("Rashba_parameters_results.json", "w") as f:
        json.dump(to_output_dict, f)
