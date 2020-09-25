from main import PendRob

if __name__ == '__main__':

    # Calculate Z for a single competent
    compounds = {'oxygen': 1.0}
    calculator = PendRob(compounds)
    calculator.calculate_Z(pressure=1, temperature=100)
    print(compounds)
    print(calculator.z_mix)
    print()

    # Calculate molar volume of mixture
    compounds = {'ethane': 0.5, 'n-butane': 0.5}
    calculator = PendRob(compounds)
    calculator.calculate_V(pressure=2, temperature=100)
    print(compounds)
    print(calculator.V)
    print()

    # Calculate partial vapor fugacities of compounds in mixture
    compounds = {'methane': 0.65, 'ethane': 0.20, 'propane': 0.15}
    calculator = PendRob(compounds)
    calculator.calculate_fv(pressure=2, temperature=100)
    print(compounds)
    print(calculator.fv)
    print()
