## Peng-Robinson Equation of State Solver

This script is meant to be able to solve primarly the partial vapor fugacities of
the compounds in a mixture using the PREOS and equations found in "Chemical, Biochemical,
and Engineering Thermodynamics" by Stanley I. Sandler. 

This solver can also solve the Z and molar volumes of a mixture both in vapor 
and liquid phase. 

This script was made for a project in CLSE 305 at VCU, if someone in the future finds
this script, please do not copy it. This was put on GitHub so anyone wanting to use
the PREOS can easily do so.

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
    >>> {'vapor': 0.9912607006581391, 'liquid': 0.08645540842653841}

## Solving for molar volume of a mixture
    compounds = {'ethane': 0.5, 'n-butane': 0.5}
    calculator = PendRob(compounds)
    calculator.calculate_V(pressure=2, temperature=100)
    print(calculator.V)
    >>> {'vapor': 0.015237604385571256, 'liquid': 0.001938170603008456}
   
## Solving for partial vapor fugacities of compounds in mixture
    compounds = {'methane': 0.65, 'ethane': 0.20, 'propane': 0.15}
    calculator = PendRob(compounds)
    calculator.calculate_fv(pressure=2, temperature=100)
    print(calculator.fv)
    >>> {'methane': 1.3077096361812484, 'ethane': 0.4032211641428493, 'propane': 0.3027229210121678}
    
## Solving for partial vapor fugacities of compounds in mixture in range of pressures
    compounds = {'methane': 0.65, 'ethane': 0.20, 'propane': 0.15}
    calculator = PendRob(compounds)
    df = calculator.generate_table_for_fv(pressures=[1, 5, 15], temperature=100, file='output.csv')
    print(df)
    >>>  Pressure (bar)       Z  ...  Vapor Fugacity of ethane  Vapor Fugacity of propane
    0               1  0.9977  ...                    0.2008                     0.1507
    1               5  0.9885  ...                    1.0204                     0.7673
    2              15  0.9657  ...                    3.1923                     2.4136