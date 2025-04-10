import numpy as np
import json

from pluggy import Result

def hill_model(x, tp, g_a, g_w):
    """
    Computes the Hill model response for a given x.
    Parameters:
    x  : float or array-like, concentration or dose values
    tp : float, top response parameter
    g_a: float, activation parameter
    g_w: float, slope parameter

    Returns:
    float or array-like: Hill model response
    """
    return tp / (1 + 10 ** ((g_a - x) * g_w))

def gain_loss_model(x, tp, g_a, g_w, l_a, l_w):
    """
    Computes the Gain-Loss model response for a given x.
    Parameters:
    x  : float or array-like, concentration or dose values
    tp : float, top response parameter
    g_a: float, gain activation parameter
    g_w: float, gain slope parameter
    l_a: float, loss activation parameter
    l_w: float, loss slope parameter

    Returns:
    float or array-like: Gain-Loss model response
    """
    gain_component = hill_model(x, tp, g_a, g_w)
    loss_component = 1 / (1 + 10 ** ((x - l_a) * l_w))
    return tp * gain_component * loss_component

def poly1_model(x, a):
    """
    Computes the response for the Poly1 (Linear Polynomial) model.
    
    Parameters:
    x  : float or array-like, concentration or dose values
    a  : float, linear coefficient

    Returns:
    float or array-like: Poly1 model response
    """
    return a * x

def poly2_model(x, a, b):
    """
    Computes the response for the Poly2 (Quadratic Polynomial) model.
    
    Parameters:
    x  : float or array-like, concentration or dose values
    a  : float, linear coefficient
    b  : float, quadratic coefficient

    Returns:
    float or array-like: Poly2 model response
    """
    return a * x + b * x**2

def exp1_model(x, a):
    return np.exp(a * x)

def exp2_model(x, a, b):
    return a * np.exp(b * x)

def exp3_model(x, a, b, p):
    return a * np.exp(b * x**p)

def exp4_model(x, a, b, c, d):
    """
    Computes the response for the Exp4 (Four-Parameter Exponential) model.
    
    Parameters:
    x  : float or array-like, concentration or dose values
    a  : float, minimum asymptote (response at zero concentration)
    b  : float, Hill slope (steepness of the curve)
    c  : float, inflection point (EC50, the concentration at half-maximum response)
    d  : float, maximum asymptote (response at infinite concentration)

    Returns:
    float or array-like: Exp4 model response
    """
    return d + ((a - d) / (1 + (x / c) ** b))

def exp5_model(x, tp, a, ga, er, p):
    return tp + (a - tp) * np.exp(-np.exp(ga * (x - er)**p))

def find_aic_values_from_assay(entry):
    models = ["hill", "gnls", "poly1", "poly2", "exp2", "exp3", "exp4", "exp5"]
    aic_values = {model: entry["mc4Param"].get(f"{model}_aic", float("inf")) for model in models}
    return aic_values

def find_best_model_from_assay(entry):
    aic_values = find_aic_values_from_assay(entry)
    best_model = min(aic_values, key=aic_values.get)  # Select the model with the lowest AIC
    return best_model, aic_values

def calculate_value_from_model(best_assay, best_model, dose):
    # Define required parameters for each model type
    model_params = {
        "gnls": ["tp", "ga", "er", "la", "lw"],
        "hill": ["tp", "ga", "er"],
        "exp1": ["a"],
        "exp2": ["a", "b"],
        "exp3": ["a", "b", "c"],
        "exp4": ["a", "b", "c", "d"],
        "exp5": ["a", "b", "c", "d", "e"],
        "poly1": ["a"],
        "poly2": ["a", "b"],
    }

    params = {param: best_assay["mc4Param"].get(f"{best_model}_{param}") for param in model_params.get(best_model, [])}

    # Check for missing parameters
    if any(value is None for value in params.values()):
        raise ValueError(f"Missing parameters for {best_model} model: {params}")

    # Compute model response based on best model
    if best_model == "gnls":
        model_responses = gain_loss_model(dose, params["tp"], params["ga"], params["er"], params["la"], params["lw"])
    elif best_model == "hill":
        model_responses = hill_model(dose, params["tp"], params["ga"], params["er"])
    elif best_model == "exp1":
        model_responses = exp1_model(dose, params["a"])
    elif best_model == "exp2":
        model_responses = exp2_model(dose, params["a"], params["b"])
    elif best_model == "exp3":
        model_responses = exp3_model(dose, params["a"], params["b"], params["c"])
    elif best_model == "exp4":
        model_responses = exp4_model(dose, params["a"], params["b"], params["c"], params["d"])
    elif best_model == "exp5":
        model_responses = exp5_model(dose, params["a"], params["b"], params["c"], params["d"], params["e"])
    elif best_model == "poly1":
        model_responses = poly1_model(dose, params["a"])
    elif best_model == "poly2":
        model_responses = poly2_model(dose, params["a"], params["b"])
    else:
        raise ValueError(f"Model {best_model} is not implemented.")
    return model_responses

def find_best_model_for_key_event(assays: list):
    results = []
    for assay in assays:
        aic_values = find_aic_values_from_assay(assay)
        results.append(aic_values)
    print(results)
    # Extract best models and sort by AIC
    best_models = []

    for idx, assay in enumerate(results):

        # Remove None values before comparison
        valid_aic_values = {model: aic for model, aic in assay.items() if aic is not None}

        if not valid_aic_values:  # If all values were None, skip this assay
            print(f"Skipping Assay {idx} (All AIC values are None)")
            continue

        best_model = min(valid_aic_values, key=valid_aic_values.get)  # Get the model with the lowest AIC
        best_aic = valid_aic_values[best_model]  # Get the lowest AIC value
        best_models.append((idx, best_model, best_aic))  # Store as tuple (id, model, aic)

    # Sort the list by AIC values
    best_models_sorted = sorted(best_models, key=lambda x: x[2])

    # Convert to a list of tuples (id, model)
    ordered_models = [(idx, model) for idx, model, _ in best_models_sorted]

    return ordered_models

def calculate_best_value(assays, ordered_models: list, dose: int):
    for model in ordered_models:
        try:
            print("", model)
            num = calculate_value_from_model(assays[model[0]], model[1], dose)
            return num #if num >0 and num < 1 else 0 if num < 0 else 1
        except:
            pass 

def load_json(file_path: str):
    with open(file_path, "r") as file:
        data = json.load(file)
    return data

def reduceData(data, chnm: str="Acrylamide",dtsxid: str=None):
    reduced = []
    for assay in data:
        if assay["chnm"] == chnm:
            reduced.append(assay)
    return reduced

data = load_json("C:/Users/tdrelangue/OneDrive/Bureau/ASSAYSDATA/63.json")
data = reduceData(data)
print(data)
ordered = find_best_model_for_key_event(data)
print(calculate_best_value(data, ordered, -50000000))