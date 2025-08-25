from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

# User inputs
Ngen = 3  # Number of generations
reaction_library = pd.DataFrame({
    'Scheme_name': ['Scheme1', 'Scheme2'],
    'Reaction_expression': ['[O:1]C(=O)[C:2].[O:3][C:4]>>[O:1]C(=O)[C:2][O:3][C:4]', '[C:1][O:2]>>[C:1]=[O:2]'],
    'Rank': [1, 2],
    'Reactivity_rule': [True, True],  # Simplified for demonstration
    'Selectivity_rule': [False, False]  # Simplified for demonstration
})

# Initialize parent molecule
mol = ['c1ccccc1OCCOC(=O)CC']  # Parent SMILES
gen = [0]
ParentID = [0]
F_deg = [0.0]
production = [1.0]
accumulation = [0.0]
N_prod = [1]
mol_index = 0
prod_index = 0

# Initialize reaction objects
reactions = [AllChem.ReactionFromSmarts(exp) for exp in reaction_library['Reaction_expression']]

# Main loop for generations
from itertools import combinations
for i in range(1, Ngen + 1):
    m = prod_index
    current_gen_indices = range(mol_index, mol_index + N_prod[i - 1])
    reactant_mols = [Chem.MolFromSmiles(mol[j]) for j in current_gen_indices]
    for p in range(N_prod[i - 1]):
        j = mol_index + p
        Np = 0
        # Ensure F_deg is long enough
        if j >= len(F_deg):
            F_deg.append(0.0)
        else:
            F_deg[j] = 0.0

        for k in range(len(reactions)):
            reaction = reactions[k]
            if reaction_library['Reactivity_rule'][k]:  # Simplified rule check
                F = 7 ** reaction_library['Rank'][k]
                num_reactants = reaction.GetNumReactantTemplates()
                # Generate all possible reactant combinations for this reaction from current generation
                for reactant_tuple_indices in combinations(current_gen_indices, num_reactants):
                    reactant_tuple = tuple(Chem.MolFromSmiles(mol[idx]) for idx in reactant_tuple_indices)
                    try:
                        products = reaction.RunReactants(reactant_tuple)
                    except Exception:
                        continue
                    for product_set in products:
                        for product in product_set:
                            product_smiles = Chem.MolToSmiles(product)
                            Np += 1
                            # Use the first reactant as parent for bookkeeping
                            parent_idx = reactant_tuple_indices[0]
                            # Extend lists for new molecule
                            mol.append(product_smiles)
                            gen.append(i)
                            ParentID.append(parent_idx)
                            # Ensure F_deg is long enough for parent_idx
                            while parent_idx >= len(F_deg):
                                F_deg.append(0.0)
                            F_deg[parent_idx] += F
                            production.append(production[parent_idx] * F / F_deg[parent_idx])
                            accumulation.append(0.0)  # Placeholder for accumulation calculation
                            N_prod.append(Np)
    mol_index += N_prod[i - 1]
    prod_index = mol_index + N_prod[i]

# Calculate accumulation values
for j in range(1, mol_index):
    accumulation[j] = production[j] - production[ParentID[j]] * F_deg[j] / F_deg[ParentID[j]]

# Output results
for idx, smiles in enumerate(mol):
    print(f"Molecule {idx}: SMILES={smiles}, Generation={gen[idx]}, Production={production[idx]}, Accumulation={accumulation[idx]}")