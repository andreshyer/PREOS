from pendRob import PendRob

if __name__ == '__main__':

    # # Calculate Z for a single competent
    # compounds = {'oxygen': 1.0}
    # calculator = PendRob(compounds)
    # calculator.calculate_Z(pressure=1, temperature=100)
    # print(compounds)
    # print(calculator.z_mix)
    # print()
    #
    # # Calculate Z for a mixture
    # compounds = {'ethane': 0.5, 'n-butane': 0.5}
    # calculator = PendRob(compounds)
    # calculator.calculate_Z(pressure=1, temperature=100)
    # print(compounds)
    # print(calculator.z_mix)
    # print()
    #
    # # Calculate molar volume of mixture
    # compounds = {'ethane': 0.5, 'n-butane': 0.5}
    # calculator = PendRob(compounds)
    # calculator.calculate_V(pressure=2, temperature=100)
    # print(compounds)
    # print(calculator.V)
    # print()
    #
    # # Calculate partial vapor fugacities of compounds in mixture
    # compounds = {'methane': 0.65, 'ethane': 0.20, 'propane': 0.15}
    # calculator = PendRob(compounds)
    # calculator.calculate_fv(pressure=2, temperature=100)
    # print(compounds)
    # print(calculator.fv)
    # print()

    # Generate df for range of fugacities
    compounds = {'methane': 0.65, 'ethane': 0.20, 'propane': 0.15}
    calculator = PendRob(compounds)
    df = calculator.generate_table_for_fv(pressures=[1, 5, 25], temperature=100, file='output.csv')
    print(compounds)
    print(df)
