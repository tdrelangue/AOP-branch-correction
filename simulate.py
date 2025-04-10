from matplotlib import pyplot as plt
from branch_correction import run_goat_dose_response, fill_method,find_critical_path
import os
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



def compute_proba(real_ac50_values, doses, aop_id=17, crit_path=None):
    results = []
    for d in doses:
        _, p = run_goat_dose_response(aop_id, d, AC50_values=real_ac50_values, selected_nodes=crit_path)
        results.append(p)
    return results

def run_simulation(ac50_values, aop_id):
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

            # Exécuter les deux en parallèle et récupérer les résultats
            with ThreadPoolExecutor(max_workers=2) as executor:
                future1 = executor.submit(compute_proba, real_ac50_values, doses, aop_id)
                future2 = executor.submit(compute_proba, real_ac50_values, doses, aop_id, crit_path)
                proba_corrected = future1.result()
                proba_uncorrected = future2.result()

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
        run_simulation([ac50_values], aop_id)

if __name__ == "__main__":
    import numpy as np
    from concurrent.futures import ThreadPoolExecutor
    ac50_values = [{ #enrichment
        "AO0": 50,
        "KE1": 50,
        "KE2": 21.81779,
        "KE3": 50,
        "KE4": 9.92126,
        "MIE0": 57.4041
    },
    { #aop-wiki alone
        "AO0": 50,
        "KE1": 50,
        "KE2": 50,
        "KE3": 50,
        "KE4": 50,
        "MIE0": 50
    }]
    aop_id = 131

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

    # real_ac50_values = { #enrichement
    #     "AO0": 0.15,
    #     "KE1": 0.15,
    #     "KE2": 0.15,
    #     "KE3": 0.15,
    #     "KE4": 0.15,
    #     "KE5": 0.15,
    #     "KE6": 0.15,
    #     "KE7": 0.15,
    #     "MIE0": 0.00026,
    # }
    # aop_id=439

    run_simulation(ac50_values, aop_id)

