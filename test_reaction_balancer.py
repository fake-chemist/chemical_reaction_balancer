import pytest
import numpy as np
from reaction_balancer import ChemicalEquation

def test_parse_equation():
    eq = ChemicalEquation('H2 + O2 -> H2O')
    assert eq.reactants == ['H2', 'O2']
    assert eq.products == ['H2O']

def test_get_elements():
    eq = ChemicalEquation('H2 + O2 -> H2O')
    elements = eq._get_elements()
    
    assert elements == ['H', 'O']

def test_get_element_count():
    eq = ChemicalEquation('H2 + O2 -> H2O')
    assert eq._get_element_count('H2') == {'H': 2}
    assert eq._get_element_count('O2') == {'O': 2}
    assert eq._get_element_count('H2O') == {'H': 2, 'O': 1}


def test_get_element_count_2():
    eq = ChemicalEquation('C2H5OH + O2 -> CO2 + H2O')
    assert eq._get_element_count('C2H5OH') == {'C': 2, 'H': 6, 'O': 1}


def test_get_matrix():
    ce = ChemicalEquation("H2 + O2 -> H2O")
    matrix = ce._get_matrix()
    expected_matrix = np.array([[2, 0, -2], [0, 2, -1]])
    assert np.array_equal(matrix, expected_matrix)


def test_get_matrix_2():
    ce = ChemicalEquation("C2H5OH + O2 -> CO2 + H2O")
    matrix = ce._get_matrix()
    # C
    # H
    # O 
    expected_matrix = np.array([
        [2, 0, -1, 0],
        [6, 0, 0, -2],
        [1, 2, -2, -1],
        ])
    assert np.array_equal(matrix, expected_matrix)


def test_balance_matrix_2x1():
    eq = ChemicalEquation('H2 + O2 -> H2O')
    balanced_matrix = eq.solve_system()
    assert balanced_matrix == [2,1,2]
    assert balanced_matrix[0] == 2

def test_balance_matrix_2x2():
    eq = ChemicalEquation('C2H5OH + O2 -> CO2 + H2O')
    balanced_matrix = eq.solve_system()
    assert balanced_matrix == [1,3,2,3]
    assert balanced_matrix[0] == 1
    assert balanced_matrix[1] == 3


def test_balnace_reaction():
    eq = ChemicalEquation('C2H5OH + O2 -> CO2 + H2O')
    balanced_eq = eq.balance_rxn()
    # 2C2H5OH + 3O2 -> 4CO2 + H2O
    # C - 4 , 4
    # H = 12, 2
    # O = 3, 9 
    assert balanced_eq == 'C2H5OH + 3O2 -> 2CO2 + 3H2O'

def test_balnace_reaction_2():
    eq = ChemicalEquation('H2 + O2 -> H2O')
    balanced_eq = eq.balance_rxn()
    assert balanced_eq == '2H2 + O2 -> 2H2O'

