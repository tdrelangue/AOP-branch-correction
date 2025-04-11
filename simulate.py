from concurrent.futures import ThreadPoolExecutor
from matplotlib import pyplot as plt
from branch_correction import run_goat_dose_response, fill_method,find_critical_path,run_dose_response_on_partial_AOP,clean_up_AOP
from create_AOP import add_AOP_variable_by_keid
import os
import numpy as np
# Démo de la dose-response (avec et sans correction)
# En supposant que run_goat_dose_response est importé et fonctionnel

# Dummy run_goat_dose_response wrapper
def simulate_dose_response(doses, corrected=True):
    probs = []
    for d in doses:
        _, prob = run_goat_dose_response(
            AOP_id="17",
            dose=d,
            probability_values=None,
            AC50_values=None,
            calculated_node="AO0",
            method=fill_method.DUMMY,
        )
        probs.append(prob)
    return probs



def compute_proba_from_id(real_ac50_values, doses, aop_id=17, crit_path=None):
    results = []
    for d in doses:
        _, p = run_goat_dose_response(aop_id, d, AC50_values=real_ac50_values, selected_nodes=crit_path)

        results.append(p)
    return results

def compute_proba_from_aop(doses, aop, crit_path=None):
    results = []
    for d in doses:
        aop = clean_up_AOP(aop)
        p = run_dose_response_on_partial_AOP(aop, d, selected_nodes=crit_path)
        results.append(p)
    return results

def run_simulation_from_id(ac50_values, aop_id, thread =True):
    output_dir = "AOP_branch_corrector/Graphs"
    os.makedirs(output_dir, exist_ok=True)
    if isinstance(ac50_values, list):
        for index, real_ac50_values in enumerate(ac50_values):
            doses = np.arange(0, 3020, 20)
            proba_corrected = []
            proba_uncorrected = []
            # Je garde le crit path de cote il changera pas
            AOP, _ = run_goat_dose_response(aop_id, 0, method=fill_method.DUMMY)
            crit_path = find_critical_path(AOP)
            if thread:
                from concurrent.futures import ThreadPoolExecutor
                # Exécuter les deux en parallèle et récupérer les résultats
                with ThreadPoolExecutor(max_workers=2) as executor:
                    future1 = executor.submit(compute_proba_from_id, real_ac50_values, doses, aop_id)
                    future2 = executor.submit(compute_proba_from_id, real_ac50_values, doses, aop_id, crit_path)
                    proba_corrected = future1.result()
                    proba_uncorrected = future2.result()
            else:
                proba_corrected = compute_proba_from_id( real_ac50_values, doses, aop_id)
                proba_uncorrected = compute_proba_from_id(real_ac50_values, doses, aop_id, crit_path)

            plt.figure(figsize=(10,6))
            plt.plot(doses, proba_corrected, label="With branch correction", linewidth=2)
            plt.plot(doses, proba_uncorrected, label="Without branch correction", linestyle="--", linewidth=2)
            plt.xlabel("Concentration ($\\mu$M)")
            plt.ylabel("AO activation probability")
            plt.title(f"Dose-response for AOP {aop_id}")
            plt.legend()
            plt.grid(True)
            plt.tight_layout()
            if index == 0:
                plt.savefig(os.path.join(output_dir, f"aop{aop_id}_dose_response.png"), dpi=300)
            else:
                plt.savefig(os.path.join(output_dir, f"aop{aop_id}_dose_response({index}).png"), dpi=300)
            plt.show()
    elif isinstance(ac50_values, dict):
        run_simulation_from_id([ac50_values], aop_id)

    
def run_simulation_from_aop(ac50_values, AOP, aop_id, thread =True):
    output_dir = "AOP_branch_corrector/Graphs"
    os.makedirs(output_dir, exist_ok=True)
    if isinstance(ac50_values, list):
        for index, real_ac50_values in enumerate(ac50_values):
            new_AOP = {node: AOP[node] for node in AOP }
            add_AOP_variable_by_keid(AOP,real_ac50_values)

            doses = np.arange(0, 3020, 20)
            proba_corrected = []
            proba_uncorrected = []
            # Je garde le crit path de cote il changera pas
            crit_path = find_critical_path(AOP)
            print(crit_path)

            if thread:
                from concurrent.futures import ThreadPoolExecutor
                # Exécuter les deux en parallèle et récupérer les résultats
                with ThreadPoolExecutor(max_workers=2) as executor:
                    future1 = executor.submit(compute_proba_from_aop, doses, clean_up_AOP(new_AOP))
                    future2 = executor.submit(compute_proba_from_aop, doses, clean_up_AOP(new_AOP), crit_path)
                    proba_corrected = future1.result()
                    proba_uncorrected = future2.result()
            else:
                proba_corrected = compute_proba_from_aop(doses, clean_up_AOP(new_AOP))
                proba_uncorrected = compute_proba_from_aop(doses, clean_up_AOP(new_AOP), crit_path)

            plt.figure(figsize=(10,6))
            plt.plot(doses, proba_corrected, label="With branch correction", linewidth=2)
            plt.plot(doses, proba_uncorrected, label="Without branch correction", linestyle="--", linewidth=2)
            plt.xlabel("Concentration ($\\mu$M)")
            plt.ylabel("AO activation probability")
            plt.title(f"Dose-response for AOP {aop_id}")
            plt.legend()
            plt.grid(True)
            plt.tight_layout()
            if index == 0:
                plt.savefig(os.path.join(output_dir, f"aop{aop_id}_dose_response.png"), dpi=300)
            else:
                plt.savefig(os.path.join(output_dir, f"aop{aop_id}_dose_response({index}).png"), dpi=300)
            plt.show()
    elif isinstance(ac50_values, dict):
        run_simulation_from_aop([ac50_values],AOP, aop_id)

if __name__ == "__main__":
    # ac50_values = [{ #enrichment
    #     "369": 50,
    #     "850": 50,
    #     "844": 21.81779,
    #     "845": 50,
    #     "846": 9.92126,
    #     "18": 57.4041
    # },
    # { #aop-wiki alone
    #     "369": 50,
    #     "850": 50,
    #     "844": 50,
    #     "845": 50,
    #     "846": 50,
    #     "18": 50
    # }]
    # aop_id = 131

    # AOP = {
    #     "AO0": {
    #         "KE_id": "369",
    #         "connections": [],
    #         "genes": [
    #             "OCA2",
    #             "UROD"
    #         ],
    #         "name": "Uroporphyria"
    #     },
    #     "KE1": {
    #         "AC50": 50,
    #         "KE_id": "850",
    #         "connections": [
    #             "KE2"
    #         ],
    #         "genes": [
    #             "BBS9",
    #             "CYP1A2",
    #             "MS4A1"
    #         ],
    #         "name": "Induction, CYP1A2/CYP1A5"
    #     },
    #     "KE2": {
    #         "KE_id": "844",
    #         "connections": [
    #             "KE3"
    #         ],
    #         "genes": [
    #             "UROD"
    #         ],
    #         "name": "Oxidation, Uroporphyrinogen"
    #     },
    #     "KE3": {
    #         "KE_id": "845",
    #         "connections": [
    #             "KE4"
    #         ],
    #         "genes": [
    #             "CYP1A2",
    #             "ROS1",
    #             "UROD"
    #         ],
    #         "name": "Inhibition, UROD"
    #     },
    #     "KE4": {
    #         "KE_id": "846",
    #         "connections": [
    #             "AO0"
    #         ],
    #         "genes": [
    #             ""
    #         ],
    #         "name": "Accumulation, Highly carboxylated porphyrins"
    #     },
    #     "MIE0": {
    #         "KE_id": "18",
    #         "connections": [
    #             "KE1"
    #         ],
    #         "name": "Activation, AhR"
    #     }
    # }
    
    # add_AOP_variable_by_keid(AOP, ac50_values[0])
    # p = run_simulation_from_aop(ac50_values, AOP, aop_id, thread = False)
    # print(AOP, p)
    # real_ac50_values = {
    # "MIE0": 50.0,          # KE 1487
    # "KE1": 70.45238499999999,  # KE 1538
    # "KE2": 50.0,           # KE 1392
    # "KE3": 50.0,           # KE 1488
    # "KE4": 35.0,           # KE 55
    # "KE5": 26.081485,      # KE 188
    # "KE6": 50.0,           # KE 1492
    # "KE7": 41.44697,       # KE 1493
    # "KE8": 37.5,           # KE 386
    # "AO0": 50.0            # AO 341
    # }
    # aop_id=17

    ac50_values = { #enrichement
        "1982": 0.15,
        "1262": 0.15,
        "1190": 0.15,
        "1196": 0.15,
        "1241": 0.15,
        "149": 0.15,
        "1971": 0.15,
        "1376": 0.15,
        "18": 0.00026,
    }

    aop_id=439

    AOP = {
        "AO0": {
            "KE_id": "1982",
            "connections": [],
            "genes": [
                ""
            ],
            "name": "metastatic breast cancer"
        },
        "KE1": {
            "KE_id": "1262",
            "connections": [
                "KE6"
            ],
            "genes": [
                "BCL2",
                "CCDC6",
                "CDKN1A",
                "FGFR1",
                "H3F3AP6",
                "H4C6",
                "HDAC1",
                "NFKB1",
                "PRDX2",
                "TCEAL1",
                "TP53"
            ],
            "name": "Apoptosis"
        },
        "KE2": {
            "KE_id": "1190",
            "connections": [
                "KE7"
            ],
            "genes": [
                "VEGFA"
            ],
            "name": "Increased, Migration (Endothelial Cells)"
        },
        "KE3": {
            "KE_id": "1196",
            "connections": [
                "AO0"
            ],
            "genes": [
                "MMRN1",
                "VEGFA"
            ],
            "name": "Increased, Invasion"
        },
        "KE4": {
            "KE_id": "1241",
            "connections": [
                "KE3"
            ],
            "genes": [
                "ITK",
                "MMRN1",
                "SLC22A3"
            ],
            "name": "Increased, Motility"
        },
        "KE5": {
            "KE_id": "149",
            "connections": [
                "KE3",
                "KE7"
            ],
            "genes": [
                "CXCL8",
                "IL18",
                "IL6",
                "ITK",
                "MMRN1",
                "ROS1",
                "SLC22A3",
                "TNF",
                "VEGFA",
                "ALOX5"
            ],
            "name": "Increase, Inflammation"
        },
        "KE6": {
            "KE_id": "1971",
            "connections": [
                "AO0"
            ],
            "genes": [
                "AKT1",
                "BCL2",
                "BRCA1",
                "BRCA2",
                "CDKN2A",
                "EGF",
                "EGFR",
                "ERBB2",
                "EREG",
                "ITK",
                "MMRN1",
                "MYC",
                "PIK3CA",
                "SLC22A3",
                "TP53",
                "VEGFA"
            ],
            "name": "Increased, tumor growth"
        },
        "KE7": {
            "KE_id": "1376",
            "connections": [
                "AO0"
            ],
            "genes": [
                "VEGFA"
            ],
            "name": "Increase, angiogenesis"
        },
        "MIE0": {
            "KE_id": "18",
            "connections": [
                "KE1",
                "KE2",
                "KE3",
                "KE4",
                "KE5"
            ],
            "genes": [
                "AHR",
                "AIP",
                "ARNT",
                "AURKAIP1",
                "BLOC1S6",
                "CYP2A13",
                "CYP2A6",
                "CYP2A7",
                "EPB42",
                "GCNT2",
                "NCOA3",
                "NRIP1",
                "PSMD9",
                "PTGES3",
                "SRC",
                "TMED10",
                "ZNRD2"
            ],
            "name": "Activation, AhR"
        }
    }

    p = run_simulation_from_aop(ac50_values, AOP, aop_id)


