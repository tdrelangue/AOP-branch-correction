import unittest
from unittest.mock import patch
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from branch_correction import (
    find_useless_nodes, 
    perform_branch_correction,
    run_dose_response_on_partial_AOP,
    calculate_node_activation_probability,
    complete_ac50_values,
    find_critical_path
)
import numpy as np


class TestAOPFunctions(unittest.TestCase):

    def setUp(self):
        """
        Initialize data for AOP graph and AC50 values.
        This is done before every test.
        """
        self.AOP = {
        "MIE": {"name": "Binding, Thiol\\seleno-proteins involved in protection against oxidative stress","connections": ["KE1"], "genes": ["OCA2"]},
        "KE1": {"name": "Decreased protection against oxidative stress", "connections": ["KE2"], "genes": ["OCA2"]},
        "KE2": {"name": "Oxidative stress", "connections": ["KE3", "KE4"], "genes": []},
        "KE3": {"name": "Glutamate homeostasis", "connections": ["KE4"], "genes": []},
        "KE4": {"name": "Increase cell injury/death", "connections": ["KE5", "KE6", "KE8"], "genes": []},
        "KE5": {"name": "Neuroinflammation", "connections": ["KE4"], "genes": []},
        "KE6": {"name": "Tissue resident cell activation", "connections": [], "genes": []},
        "KE7": {"name": "increased proinflammatory mediators", "connections": ["KE4"], "genes": []},
        "KE8": {"name": "Decreased neural network function", "connections": ["AO","KE3"], "genes": []},
        "AO": {"name": "Impairment learning and memory", "connections": [], "genes": ["AOX1"]},
    }
        AC50 = {node: 20 for node in self.AOP}
        complete_ac50_values(self.AOP)
    
    def test_find_useless_nodes(self):
        """
        Test the identification of useless nodes that form loops.
        """
        useless_nodes = find_useless_nodes(self.AOP)

        # Expect KE2 to be part of a useless loop
        self.assertIn("KE5", useless_nodes)
    
    def test_perform_branch_correction(self):
        """
        Test the correct branch correction in the AOP graph.
        """
        def prod(numbers):
            prod=1
            for item in numbers:
                prod*=item
            return prod

        # Setup mock branch data
        critical_path = find_critical_path(self.AOP)
        branch_path = ["KE2", "KE3"]
        endpoint = "KE4"


        # Calculate initial probabilities
        calculate_node_activation_probability(AOP=self.AOP, dose=1000)
        original_proba = prod([self.AOP[node]["P(prior|event)"] for node in self.AOP])

        # Perform branch correction
        perform_branch_correction(self.AOP, critical_path, branch_path, endpoint)

        # Check if the probability has changed (i.e., branch correction applied)
        updated_proba = prod([self.AOP[node]["P(prior|event)"] for node in self.AOP])
        
        # Assert the probability has been updated
        self.assertNotEqual(original_proba, updated_proba)

    def test_calculate_node_activation_probability(self):
        """
        Test node activation probability calculation.
        """
        dose = 0
        # Calculate node activation probabilities
        calculate_node_activation_probability(self.AOP, dose)

        # Check that the calculated probabilities are within the expected range (0 to 1)
        for node, details in self.AOP.items():
            prob = details.get("P(prior|event)")
            self.assertGreaterEqual(prob, 0)
            self.assertLessEqual(prob, 1)
        
        dose = 10000000
        # Calculate node activation probabilities
        calculate_node_activation_probability(self.AOP, dose)

        # Check that the calculated probabilities are within the expected range (0 to 1)
        for node, details in self.AOP.items():
            prob = details.get("P(prior|event)")
            self.assertGreaterEqual(prob, 0)
            self.assertLessEqual(prob, 1)

    def test_run_dose_response_on_partial_AOP(self):
        """
        Test the response calculation on a partial AOP graph.
        """
        selected_nodes = ["MIE", "KE1", "KE2", "KE4", "KE5","KE8", "AO"]
        dose = 1000
        # Run dose response on selected nodes
        result = run_dose_response_on_partial_AOP(self.AOP, dose, selected_nodes)

        # Verify that the result is a valid probability (between 0 and 1)
        self.assertGreaterEqual(result, 0)
        self.assertLessEqual(result, 1)

if __name__ == '__main__':
    unittest.main()
