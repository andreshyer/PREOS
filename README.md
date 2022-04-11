## Peng-Robinson Equation of State Solver

This script is meant to be able to solve primarly the partial vapor fugacities of
the compounds in a mixture using the PREOS and equations found in "Chemical, Biochemical,
and Engineering Thermodynamics" by Stanley I. Sandler. 

This solver can also solve the Z and molar volumes of a mixture both in vapor 
and liquid phase. 

## Solving for Z of a single component
    compounds = {'oxygen': 1.0}
    calculator = PendRob(compounds)
    calculator.calculate_Z(pressure=1, temperature=100)
    print(calculator.z_mix)
    >>> {'single state': 0.9997068078533051}
   
## Solving for Z of a mixture
    compounds = {'ethane': 0.5, 'n-butane': 0.5}
    calculator = PendRob(compounds)
    calculator.calculate_Z(pressure=1, temperature=100)
    print(calculator.z_mix)
    >>> {'single state': 0.9913218856617976}

## Solving for molar volume of a mixture
    compounds = {'ethane': 0.5, 'n-butane': 0.5}
    calculator = PendRob(compounds)
    calculator.calculate_V(pressure=2, temperature=100)
    print(calculator.V)
    >>> {'single state': 0.015241545952714833}
   
## Solving for partial vapor fugacities of compounds in mixture
    compounds = {'methane': 0.65, 'ethane': 0.20, 'propane': 0.15}
    calculator = PendRob(compounds)
    calculator.calculate_fv(pressure=2, temperature=100)
    print(calculator.fv)
    >>> {'methane': 1.2980602636892002, 'ethane': 0.3967781774811893, 'propane': 0.29602876757692975}
    
## Solving for partial vapor fugacities of compounds in mixture in range of pressures
    compounds = {'methane': 0.65, 'ethane': 0.20, 'propane': 0.15}
    calculator = PendRob(compounds)
    df = calculator.generate_table_for_fv(pressures=[2, 7, 20], temperature=100, file='output.csv')
    print(df)
    >>> Pressure (bar)  ...  Vapor Fugacity of propane (bar)
    0               2  ...                           0.2960
    1               7  ...                           1.0022
    2              20  ...                           2.6266
