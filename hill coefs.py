from branch_correction import hill_activation

import numpy as np
from scipy.stats import linregress

def infer_hill_from_ac_values(ac_dict):
    """
    Infer Hill coefficient from given ACx values where x is percent response (e.g., 5, 10, 20, 50).
    
    Parameters:
        ac_dict (dict): Mapping of response percent (int) to ACx dose (float), e.g.
                        {5: ac5, 10: ac10, 20: ac20, 50: ac50}
    
    Returns:
        hill_coefficient (float): Inferred Hill coefficient
        r_value (float): Correlation coefficient (goodness of fit)
    """
    # Extract values and sort by response level
    responses = sorted(ac_dict.keys())
    doses = [ac_dict[r] for r in responses]

    # Convert response % to fraction
    E = np.array(responses) / 100
    D = np.array(doses)
    
    ac50 = ac_dict[50]
    
    # Log-transform for linearisation
    X = np.log(D / ac50)
    Y = np.log(E / (1 - E))

    slope, intercept, r_value, p_value, std_err = linregress(X, Y)

    return slope, r_value
ac_values = {
    5: 0.5,
    10: 1.0,
    20: 2.0,
    50: 5.0
}

hill_coefficient, r = infer_hill_from_ac_values(ac_values)
print(f"Inferred Hill coefficient: {hill_coefficient:.3f} (R={r:.3f})")
