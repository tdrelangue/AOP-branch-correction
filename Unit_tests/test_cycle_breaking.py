import unittest
import sys
import os
# Add the root directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from branch_correction import find_useless_connections, remove_useless_connections, build_graph

class TestCycleBreaking(unittest.TestCase):
    def setUp(self):
        # AOP examples for testing
        self.aop1 = {
            "KE1": {"name": "Step 1", "connections": ["KE2"], "genes": []},
            "KE2": {"name": "Step 2", "connections": ["KE3"], "genes": []},
            "KE3": {"name": "Step 3", "connections": ["KE4"], "genes": []},
            "KE4": {"name": "Step 4", "connections": ["KE5", "KE6", "KE3"], "genes": []},  # Cycle back to KE3
            "KE5": {"name": "Step 5", "connections": ["KE4"], "genes": []},  # Cycle back to KE4
            "KE6": {"name": "Step 6", "connections": ["AO"], "genes": []},
            "AO": {"name": "Adverse Outcome", "connections": [], "genes": []},
        }

        self.aop2 = {
            "KE1": {"name": "Start", "connections": ["KE2"], "genes": []},
            "KE2": {"name": "Middle", "connections": ["KE3"], "genes": []},
            "KE3": {"name": "End", "connections": [], "genes": []},
        }

        self.aop_empty = {}

    def test_find_useless_connections_simple(self):
        connections_to_break = find_useless_connections(self.aop1)
        expected_connections = [('KE5', 'KE4'), ('KE4', 'KE3')]
        self.assertCountEqual(connections_to_break, expected_connections)

    def test_find_useless_connections_no_cycles(self):
        connections_to_break = find_useless_connections(self.aop2)
        self.assertEqual(connections_to_break, [])

    def test_remove_useless_connections(self):
        # Expected result after removing useless connections
        expected_aop = {
            "KE1": {"name": "Step 1", "connections": ["KE2"], "genes": []},
            "KE2": {"name": "Step 2", "connections": ["KE3"], "genes": []},
            "KE3": {"name": "Step 3", "connections": ["KE4"], "genes": []},
            "KE4": {"name": "Step 4", "connections": ["KE5", "KE6"], "genes": []},
            "KE5": {"name": "Step 5", "connections": [], "genes": []},
            "KE6": {"name": "Step 6", "connections": ["AO"], "genes": []},
            "AO": {"name": "Adverse Outcome", "connections": [], "genes": []},
        }

        # Remove connections and verify the updated AOP
        remove_useless_connections(self.aop1)
        self.assertEqual(self.aop1, expected_aop)

    def test_empty_aop(self):
        connections_to_break = find_useless_connections(self.aop_empty)
        self.assertEqual(connections_to_break, [])
        remove_useless_connections(self.aop_empty)
        self.assertEqual(self.aop_empty, {})

if __name__ == "__main__":
    unittest.main()
