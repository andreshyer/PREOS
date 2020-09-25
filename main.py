import pandas as pd
import numpy as np
from math import log, sqrt, exp


class PendRob:

    def __init__(self, compounds):

        """

        This the init function, main purpose of this function is to gather data in csv files and define
        variables to be used in other parts of the script.

        As a little note, and function that is formatted __something__(), is meant to only be used by the class,
        and not the user. Anything that is formatted something() is okay to use.

        Also, all equations and tables referenced are from "Chemical, Biochemical, and Engineering Thermodynamics"
        by Stanly I. Sandler

        :param compounds: A dict of compounds, where keys are compound names found in preso_ab.csv and values
        are the vapor mole fraction
        """

        self.compounds = compounds

        self.R = 8.314  # In Pa * m^3 * K^-1 * mol^-1
        self.raw_Tc_Pc_w_data = pd.read_csv('preso_ab.csv')
        self.raw_Tc_Pc_w_data = self.raw_Tc_Pc_w_data.rename(columns={'Tc (K)': 'Tc', 'Pc (MPa)': 'Pc'})
        self.raw_binary_data = pd.read_csv('binary_interaction.csv')

        y_sum = 0
        for compound, yi in self.compounds.items():
            compound_row = self.raw_Tc_Pc_w_data.loc[self.raw_Tc_Pc_w_data['compound'] == compound]
            if len(compound_row) == 0:
                raise TypeError(f'Compound {compound} not found')
            y_sum += yi
        y_sum = round(y_sum, 2)
        if y_sum != 1:
            raise TypeError(f'Mole fractions do not add up to 1, they add up to {y_sum}')

        # Define holding variables
        self.a_data = dict(zip(self.compounds.keys(), [np.nan] * len(self.compounds)))
        self.b_data = dict(zip(self.compounds.keys(), [np.nan] * len(self.compounds)))
        self.a_mix = np.nan
        self.b_mix = np.nan
        self.A_mix = np.nan
        self.B_mix = np.nan
        self.z_mix = np.nan
        self.pressure = np.nan
        self.temperature = np.nan
        self.fv = dict(zip(self.compounds.keys(), [np.nan] * len(self.compounds)))

    def __calculate_a_b__(self, compound):

        """
        This equation used for this can be found in pg.250, equations 6.7-1 - 6.7-4

        :param compound: Individual compound to calculate a_ii and b_ii for
        :return: a_ii and b_ii
        """

        # Pull out information from csv files about the compound
        compound_data = self.raw_Tc_Pc_w_data[self.raw_Tc_Pc_w_data['compound'] == compound]
        compound_data = compound_data.to_dict('records')[0]
        Tc = compound_data['Tc']
        Pc = compound_data['Pc'] * 10**6  # Convert MPa to Pa
        w = compound_data['w']

        # calculate a and b
        k = 0.37464 + 1.54226*w - 0.26992*(w**2)
        alpha = (1 + k * (1 - sqrt(self.temperature/Tc)))**2

        a = 0.45724 * alpha * (self.R*Tc)**2/Pc
        b = 0.07780 * self.R*Tc/Pc

        return a, b

    def __calculate_a_b_mix__(self):
        """
        Goal is to calculate the a_mix and b_mix.
        This uses equation

        :return:
        """
        a_mix = 0
        b_mix = 0
        for compound, y1 in self.compounds.items():
            b_mix += y1 * self.b_data[compound]
            for sub_compound, y2 in self.compounds.items():
                k = self.__fetch_kij__(compound, sub_compound)
                a1 = self.a_data[compound]
                a2 = self.a_data[sub_compound]
                a_mix += (y1*y2) * (1-k) * sqrt(a1*a2)
        return a_mix, b_mix

    def __fetch_kij__(self, compound, sub_compound):
        k = 0
        if compound in list(self.raw_binary_data.columns):
            k = self.raw_binary_data[['compound', compound]]
            k = k.loc[k['compound'] == sub_compound]
            if len(k) > 0:
                k = k.to_dict('records')[0]
                k = k[compound]
                if np.isnan(k):
                    k = 0
        return k

    def calculate_Z(self, pressure, temperature, compound=None):

        """
        Goal is to solve for Z in PREOS of mixtures

        where:
        Z^3 + uZ^2 + wZ + v = 0

        u = -1 + B
        w = A - 3B^2 - 2B
        v = -AB + B^2 + B^3

        A = aP / (RT)^2
        B = bP / RT

        Z = PV / RT

        where a, b, A, B, and Z are all mixture values

        This is the main part of the script that does most of the heavy lifting

        :param pressure: absolute pressure of mixture in bar
        :param temperature: temperture of mixture in C
        :param compound: if not None, calculate Z only for a specific compound in mixture and return Z
        :return: z constant for mixture at conditions given
        """

        self.pressure, self.temperature = self.__convert_variables__(pressure, temperature)
        if self.pressure < 0:
            raise ValueError("Absolute pressure can not less than 0")
        if self.temperature < 0:
            raise ValueError("Temperature cannot be less than 0K")

        for sub_compound in self.compounds:
            self.a_data[sub_compound], self.b_data[sub_compound] = self.__calculate_a_b__(sub_compound)

        if compound is None:
            self.a_mix, self.b_mix = self.__calculate_a_b_mix__()
            self.A_mix, self.B_mix = self.__calculate_A_B__(self.a_mix, self.b_mix)
            u, w, v = self.__calculate_u_w_v__(self.A_mix, self.B_mix)
            self.z_mix = self.__calculate_z_roots__(u, w, v)

        else:
            a, b = self.a_data[compound], self.b_data[compound]
            A, B = self.__calculate_A_B__(a, b)
            u, w, v = self.__calculate_u_w_v__(A, B)
            z = self.__calculate_z_roots__(u, w, v)
            return z

    @staticmethod
    def __convert_variables__(pressure, temperature):
        pressure = pressure * 10**5  # Convert bar to Pa
        temperature = temperature + 273.15  # Convert C to K
        return pressure, temperature

    def __calculate_A_B__(self, a, b):
        A = (a * self.pressure) / ((self.R * self.temperature) ** 2)
        B = (b * self.pressure) / (self.R * self.temperature)
        return A, B

    @staticmethod
    def __calculate_u_w_v__(A, B):
        u = -1 + B
        w = A - 3 * (B ** 2) - 2 * B
        v = -(A * B) + B ** 2 + B ** 3
        return u, w, v

    @staticmethod
    def __calculate_z_roots__(u, w, v):
        roots = np.roots([1, u, v, w])
        roots = roots[np.isreal(roots)]
        roots = np.real(roots)
        if len(roots) == 3:
            roots = roots[roots > 0]
            roots = {'vapor': max(roots), 'liquid': min(roots)}
        elif len(roots) == 1:
            roots = {'single state': roots[0]}
        elif len(roots) == 0:
            raise ValueError("No roots for Z found at given conditions")
        else:
            raise Exception("Something broke, there should not be only two real roots")
        return roots

    def calculate_V(self, pressure, temperature):

        """
        :param pressure: absolute pressure of mixture in bar
        :param temperature: temperture of mixture in C
        :return: molar volume in m^3/mol for mixture at conditions given
        """

        all_z = self.calculate_Z(pressure, temperature)
        V = {}
        for state, z in all_z.items():
            V[state] = z * self.R * self.temperature / self.pressure
        return V

    def calculate_fv(self, pressure, temperature):

        self.calculate_Z(pressure, temperature)
        z_mix = self.__pull_z_value__(self.z_mix, 'vapor')

        for compound, yi in self.compounds.items():
            ai = self.a_data[compound]
            bi = self.b_data[compound]
            Ai, Bi = self.__calculate_Ai_Bi__(ai, bi)
            sum_yj_Aij = self.__calculate_sum_yj_Aij__(compound)

            try:
                term1 = (Bi * (z_mix - 1)/self.B_mix)
                term2 = log(z_mix - self.B_mix)
                term3 = self.A_mix/(2 * sqrt(2) * self.B_mix)
                term4 = (2*sum_yj_Aij/self.A_mix) - Bi/self.B_mix

                term5_1 = z_mix + ((1 + sqrt(2)) * self.B_mix)
                term5_2 = z_mix + ((1 - sqrt(2)) * self.B_mix)
                term5 = log(term5_1/term5_2)

                fv = term1 - term2 - (term3 * term4 * term5)
                fv = exp(fv)
                fv = fv * yi * self.pressure * 10**(-5)
                self.fv[compound] = fv
            except ValueError:
                raise ValueError(f'Cannot calculate partial vapor fugacities of compounds at {self.temperature-273.15} '
                                 f'C and {self.pressure * 10**-5} Bar because mixture cannot exist as vapor at '
                                 f'given conditions')

    @staticmethod
    def __pull_z_value__(z_dict, phase):
        if phase in z_dict.keys():
            z = z_dict[phase]
        elif 'single state' in z_dict.keys():
            z = z_dict['single state']
        else:
            raise TypeError("Phase not found in z dict")
        return z

    def __calculate_Ai_Bi__(self, ai, bi):
        Ai = (ai * self.pressure) / ((self.R * self.temperature) ** 2)
        Bi = (bi * self.pressure) / (self.R * self.temperature)
        return Ai, Bi

    def __calculate_sum_yj_Aij__(self, compound_i):
        sum_yj_Aij = 0
        yi = self.compounds[compound_i]
        for compound_j, yj in self.compounds.items():
            k = self.__fetch_kij__(compound_i, compound_j)
            ai = self.a_data[compound_i]
            aj = self.a_data[compound_j]
            aij = (yi*yj) * (1 - k) * sqrt(ai*aj)
            Aij = aij*self.pressure/((self.R*self.temperature)**2)
            sum_yj_Aij += yj*Aij
        return sum_yj_Aij


# Gather compounds, make sure they are correct...
# compounds = {'methane': 0.65, 'ethane': 0.20, 'propane': 0.15}
# compounds = {'oxygen': 1.0}
compounds = {'ethane': 0.5, 'n-butane': 0.5}

calculator = PendRob(compounds)
calculator.calculate_fv(pressure=1, temperature=100)
print(calculator.z_mix)
print(calculator.fv)
