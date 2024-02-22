import re 
import numpy as np
import pulp
import argparse

class ChemicalEquation: 
    def __init__(self, equation):
        """
        Initializes a ChemicalEquation object.

        Parameters:
        - equation (str): The chemical equation to be balanced.

        Attributes:
        - equation (str): The chemical equation.
        - reactants (list): List of reactants in the equation.
        - products (list): List of products in the equation.
        - elements (list): List of unique elements present in the equation.
        - matrix (numpy.ndarray): Matrix representation of the equation.

        """
        self.equation = equation
        self.reactants, self.products = self._parse_equation(equation)
        self.elements = self._get_elements()
        self.matrix = self._get_matrix()

    def _parse_equation(self, equation):
        """
        Parse a chemical equation into reactants and products.

        This method takes a string representing a chemical equation, splits it into reactants and products, 
        and returns two lists containing the reactants and products, respectively. Each reactant and product 
        is stripped of leading and trailing whitespace.

        Parameters:
        equation (str): A string representing a chemical equation. Reactants and products should be separated 
        by '->', and individual reactants/products should be separated by '+'.

        Returns:
        tuple: A tuple containing two lists. The first list contains the reactants, and the second list 
        contains the products.
        """
        reactants, products = equation.split('->')
        reactants = reactants.split('+')
        reactants = [reactant.strip() for reactant in reactants]
        products = products.split('+')
        products = [product.strip() for product in products]
        return reactants, products
    

    def _get_elements(self):
        elements = []
        for compound in self.reactants + self.products:
            for element in re.findall('[A-Z][a-z]*', compound):
                if element not in elements:
                    elements.append(element)
        return elements
    
    def _get_element_count(self, compound):
        """
        Get the count of each element in a chemical compound.

        This method takes a string representing a chemical compound and returns a dictionary where the keys 
        are elements and the values are the counts of those elements in the compound.

        Parameters:
        compound (str): A string representing a chemical compound. Elements should be represented by their 
        chemical symbols, and counts should be represented by integers following the symbols. If no count is 
        provided for an element, it is assumed to be 1.

        Returns:
        dict: A dictionary where the keys are elements and the values are the counts of those elements in 
        the compound.
        """
        element_count = {}
        for element in re.findall(r'[A-Z][a-z]*\d*', compound):
            element_name, count = re.match(r'([A-Z][a-z]*)(\d*)', element).groups()
            count = int(count) if count else 1
            if element_name in element_count:
                element_count[element_name] += count
            else:
                element_count[element_name] = count
        return element_count



    def _get_matrix(self):

        """
        Generate a matrix representation of the chemical equation.

        This method creates a matrix where each row corresponds to an element in the equation, and each column 
        corresponds to a reactant or product. The value in each cell is the count of the element in the 
        reactant/product, with reactants being negative and products being positive.

        Returns:
        np.array: A 2D numpy array representing the chemical equation. The number of rows is equal to the 
        number of unique elements in the equation, and the number of columns is equal to the total number of 
        reactants and products. Each cell contains the count of an element in a reactant/product, with 
        reactants being negative and products being positive.
        """

        matrix = np.zeros((len(self.elements), len(self.reactants) + len(self.products)), dtype=int)
        for i, reactant in enumerate(self.reactants):
            for element, count in self._get_element_count(reactant).items():
                matrix[self.elements.index(element), i] = int(count)
        for i, product in enumerate(self.products):
            for element, count in self._get_element_count(product).items():
                #Negative sign added to the count because we're moving the prodcuts to the left side 
                matrix[self.elements.index(element), i + len(self.reactants)] = -int(count)
        return matrix



    def solve_system(self):
        """
        A, B = self.matrix[:, :len(self.reactants)], self.matrix[:, len(self.reactants):]
        A_inv = np.linalg.pinv(A)
        coefficients = np.round(A_inv @ B).astype(int)
        
        An example of where context matters. The above code will solve the matrix for 
        the coefficients of the reactants and products. However, the coefficients are 
        not bound by areal world constraints. This can be seeing by this producing 
        negative coefficients for this equation: 'H2 + O2 -> H2O'.

         '-1H2 + 0O2 -> -1H2O'
        """
        # Create the 'problem'
        matrix = self._get_matrix()
        problem = pulp.LpProblem("Problem", pulp.LpMinimize)

        # Variables
        variables = [pulp.LpVariable(f'var{i}', lowBound=1) for i in range(matrix.shape[1])]

        # Objective function
        problem += 0  # No real objective function, we just want to satisfy constraints

        # Constraints
        for row in matrix:
            problem += pulp.lpSum(row[i]*variables[i] for i in range(len(row))) == 0

        # Solve the problem
        status = problem.solve(pulp.PULP_CBC_CMD(msg=False))

        return variables

    def balance_rxn(self):
        coefficients = self.solve_system()
        coefficients = [int(var.varValue) for var in coefficients]

        # Construct the balanced equation
        balanced_equation = ""
        for i, reactant in enumerate(self.reactants):
            if i > 0:
                balanced_equation += " + "
            # If the coefficient is 1, don't print it
            if coefficients[i] != 1:
                balanced_equation += str(coefficients[i])
            balanced_equation += reactant.strip()

        balanced_equation += " -> "

        for i, product in enumerate(self.products):
            if i > 0:
                balanced_equation += " + "
            # If the coefficient is 1, don't print it
            if coefficients[i + len(self.reactants)] != 1:
                balanced_equation += str(coefficients[i + len(self.reactants)])
            balanced_equation += product.strip()

        return balanced_equation



def main():
    parser = argparse.ArgumentParser(description='Balance a chemical equation.')
    parser.add_argument('equation', type=str, help='The chemical equation to balance.')
    args = parser.parse_args()

    eq = ChemicalEquation(args.equation)
    eq_balanced = eq.balance_rxn()
    print(args.equation)
    print("_______________________________")
    print(eq_balanced)

if __name__ == '__main__':
    main()