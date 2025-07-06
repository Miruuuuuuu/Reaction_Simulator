from rdkit import Chem
from rdkit.Chem import AllChem
import sys
import os
from typing import Dict, List, Tuple

class ChemicalReactionBalancer:
    def __init__(self):
        self.rdkit_failed = False
        
    def parse_compound(self, formula: str) -> Dict[str, int]:
        """Parse chemical formula into element counts without RDKit"""
        element_counts = {}
        i = 0
        while i < len(formula):
            # Find element symbol (1 or 2 characters)
            if i + 1 < len(formula) and formula[i+1].islower():
                element = formula[i:i+2]
                i += 2
            else:
                element = formula[i]
                i += 1
            
            # Find number
            num = ""
            while i < len(formula) and formula[i].isdigit():
                num += formula[i]
                i += 1
            
            count = int(num) if num else 1
            element_counts[element] = element_counts.get(element, 0) + count
            
        return element_counts

    def read_reaction_data(self) -> dict:
        """Read reaction data from text file"""
        with open("reaction_data.txt", "r") as f:
            # Read symbols
            symbols = [f.readline().strip() for _ in range(4)]
            
            # Read atom counts
            atom_counts = [int(f.readline().strip()) for _ in range(4)]
            
            # Read ion charges
            num_charges = int(f.readline().strip())
            ion_charges = [int(f.readline().strip()) for _ in range(num_charges)]
            
            # Read positive ions
            num_positive = int(f.readline().strip())
            positive_ions = [f.readline().strip() for _ in range(num_positive)]
            
            # Read negative ions
            num_negative = int(f.readline().strip())
            negative_ions = [f.readline().strip() for _ in range(num_negative)]
            
        return {
            "symbols": symbols,
            "atomCounts": atom_counts,
            "ionCharges": ion_charges,
            "positive": positive_ions,
            "negative": negative_ions
        }

    def write_results(self, results: dict):
        """Write results to text file"""
        with open("reaction_results.txt", "w") as f:
            # Write status
            f.write(f"{results['status']}\n")
            
            if results['status'] == 'error':
                f.write(f"{results['message']}\n")
                return
                
            # Write balance status
            f.write(f"{str(results['is_balanced']).lower()}\n")
            
            if not results['is_balanced']:
                # Write suggestions
                suggestions = results['corrections']['suggestions']
                f.write(f"{len(suggestions)}\n")
                for suggestion in suggestions:
                    f.write(f"{suggestion}\n")
                    
                # Write corrected atom counts
                for count in results['corrections']['atom_counts']:
                    f.write(f"{count}\n")

    # Keep other methods the same (verify_with_rdkit, calculate_total_atoms, etc.)

    def process_reaction(self):
        """Main processing function"""
        try:
            # Read input data
            data = self.read_reaction_data()
            
            # Process the reaction
            results = self.balance_reaction(data)
            
            # Write results
            self.write_results(results)
            return True
            
        except Exception as e:
            # Write error results
            error_results = {
                "status": "error",
                "message": str(e)
            }
            self.write_results(error_results)
            return False

def main():
    balancer = ChemicalReactionBalancer()
    success = balancer.process_reaction()
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()