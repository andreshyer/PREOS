## Peng-Robinson Equation of State Solver

This script is meant to be able to solve primarly the partial vapor fugacities of
compound in a mixture using the PREOS and equations found in "Chemical, Biochemical,
and Engineering Thermodynamics" by Stanley I. Sandler. 

This solver can also solve the Z and molar volumes of a mixture both in vapor 
and liquid phase. 

This script was made for a project in CLSE 305 at VCU, if someone in the future finds
this script, please do not copy it. This was put on GitHub for anyone wanting to use
the PREOS easily can do so.

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