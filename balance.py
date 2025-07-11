import re
import sympy as sp
from typing import List, Dict

def get_molecule_atom_counts_and_element_names(molecules: List[str]) -> Dict[str, int]:
    """
    Takes a list of molecules and returns a list of dictionaries
    where each dictionary maps between chemical symbol and its count
    for a given molecule in the input list.
    """
    molecule_atom_counts = []
    element_names = set()
    for molecule in molecules:
        molecule_atom_count = {}
        separated_molecule = re.findall(r'[A-Z][a-z]*\d*', molecule)
        for element_name_and_count in separated_molecule:
            element_name = re.search(r'[A-Z][a-z]*', element_name_and_count).group()
            # Add the element name to the set of elements.
            element_names.add(element_name)
            if element_name not in molecule_atom_count:
                molecule_atom_count[element_name] = 0
            element_count = re.search(r'\d+', element_name_and_count)
            # If the element has a subscript, add it to the list. Otherwise set it to 1.
            if element_count:
                molecule_atom_count[element_name] += int(element_count.group())
            else:
                molecule_atom_count[element_name] += 1
        molecule_atom_counts.append(molecule_atom_count)
    return molecule_atom_counts, element_names
    

def balance(eq):
    """
    Takes a string of a chemical equation and returns a string of the balanced equation.
    The equation must have no leading coefficients for it to work properly.
    """
    n = len(re.findall(r'\+|\=', eq)) + 1 # The number of molecules in the chemical equation.
    
    # Sides is a list containing strings of the two sides of the equation.
    sides = re.findall(r'[^\=]+', eq)

    # Separate the left and right hand sides into a list of compound names.
    left_molecules = re.findall(r'[^\s\+]+', sides[0])
    right_molecules = re.findall(r'[^\s\+]+', sides[1])
    
    # Separate the compound names into a list of dicts containing each molecules element counts.
    # Also get the element names from each side of the equation and combine them to get the complete list.
    left_element_counts, left_element_names = get_molecule_atom_counts_and_element_names(left_molecules)
    right_element_counts, right_element_names = get_molecule_atom_counts_and_element_names(right_molecules)
    element_names = list(left_element_names.union(right_element_names))
    
    # Create a matrix to solve the chemical equation.
    # Entries from the left side of the equation are positive, and
    # entries from the right side of the equation are negative in the matrix.
    left_length = len(left_element_counts)
    matrix = [[0 for _ in range(n+1)] for _ in range(len(element_names))]
    for i in range(len(element_names)):
        element_name = element_names[i]
        for j in range(len(left_element_counts)):
            if element_name in left_element_counts[j]:
                matrix[i][j] += left_element_counts[j][element_name]
        for j in range(len(right_element_counts)):
            if element_name in right_element_counts[j]:
                matrix[i][j + left_length] -= right_element_counts[j][element_name]
    
    # Solve matrix.
    sympy_matrix = sp.Matrix(matrix)
    syms = sp.symbols([f"x{num}" for num in range(n)])
    sols = sp.solve_linear_system(sympy_matrix, *syms)
    
    # Substitute in 1 for any free variable.
    for item in sols:
        for j in range(len(syms)):
            sols[item] = sols[item].subs(syms[j], 1)
    
    # If a variable is missing from the solutions dict, add 1 as its value.
    for sym in syms:
        if sym not in sols:
            sols[sym] = sp.Rational(1, 1)

    # Create a list of denominators to find the lcm.
    denominators = [x.q for x in sols.values()]
    lcm = sp.lcm(denominators)

    # Multiply each solution by the lcm.
    for sym in sols:
        sols[sym] *= lcm
    solutions_list = list(sols.values())

    # Create string to return.
    left_side = []
    for i in range(left_length):
        count = str(solutions_list[i]) if solutions_list[i] > 1 else ""
        left_side.append("".join([count, left_molecules[i]]))
    left_side = "+".join(left_side)
    right_side = []
    for i in range(len(right_molecules)):
        count = str(solutions_list[i + left_length]) if solutions_list[i + left_length] > 1 else ""
        right_side.append("".join([count, right_molecules[i]]))
    right_side = "+".join(right_side)
    return_string = "=".join([left_side, right_side])
    return return_string

assert balance("CH4   + O2  =  CO2 +   H2O") == "CH4+2O2=CO2+2H2O"
assert balance("H2+O2=H2O") == "2H2+O2=2H2O"
assert balance("CH3OH=CH3OCH3+H2O") == "2CH3OH=CH3OCH3+H2O"