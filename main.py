import pandas as pd
import numpy as np
from math import log, sqrt, exp


class PendRob:

    def __init__(self, compounds):

        """
        The goal of this class is to be able to solve for Z, V, or partial fugacities of vapor for any mixture found
        in BOTH table 6.6-1 and 9.4-1, please see examples.py.

        All equations and tables referenced are from "Chemical, Biochemical, and Engineering Thermodynamics"
        by Stanly I. Sandler

        The csv files containing data are just copied from pg.424 table 9.4-1 and pg.241 table 6.6-1

        The general PREOS can be found in pg.208, equation 6.4-2

        This the init function, the main purpose of this function is to gather data in csv files and define
        variables to be used in other parts of the script.

        As a little note, a function that is formatted __something__(), is meant to only be used by the class,
        and not the user. Anything that is formatted something() is okay to use.
        Also PREOS stands for Peng-Robinson Equation of State

        :param compounds: A dict of compounds, where keys are compound names found in preso_Tc_Pc_w_data.csv and values
        are the vapor mole fraction
        """

        self.compounds = compounds

        self.R = 8.314  # In Pa * m^3 * K^-1 * mol^-1
        self.raw_Tc_Pc_w_data = pd.read_csv('preso_Tc_Pc_w_data.csv')
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
        self.V = {}

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
        Goal is to calculate the a_mix and b_mix. The equation for this function is in pg.423 equation 9.4-8 and
        9.4-9. The binary interaction parameter is fetched in __fetch_kij__()

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

        """
        This function will fetch kji that is used in equation 9.4-9. This is just retrieved from the
        binary_interaction.csv file.

        :param compound: compound 1
        :param sub_compound: compound 2
        :return: kij = kji
        """

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
        equation used is from pg.209 equation 6.4-4

        where:
        Z^3 + uZ^2 + wZ + v = 0

        u = -1 + B
        w = A - 3B^2 - 2B
        v = -AB + B^2 + B^3

        A = aP / (RT)^2
        B = bP / RT

        Z = PV / RT

        where a, b, A, B, and Z can all be mixture or single component values

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

        """
        Simple little helper function to convert pressure from bar to Pa, and temperture from C to K

        :param pressure: pressure in Bar
        :param temperature: temperature in C
        :return: pressure in Pa, temperature in K
        """

        pressure = pressure * 10**5  # Convert bar to Pa
        temperature = temperature + 273.15  # Convert C to K
        return pressure, temperature

    def __calculate_A_B__(self, a, b):

        """
        This uses equation pg.210, table 6.4-3. It is important to note that equation works for both
        individual components (a, b) and mixtures (a_mix, b_mix)

        :param a: a from PREOS
        :param b: b from PREOS
        :return: A, B for PREOS
        """

        A = (a * self.pressure) / ((self.R * self.temperature) ** 2)
        B = (b * self.pressure) / (self.R * self.temperature)
        return A, B

    @staticmethod
    def __calculate_u_w_v__(A, B):
        """
        This uses equation pg.210, table 6.4-3. It is important to note that equation works for both
        individual components (A, B) and mixtures (A_mix, B_mix)

        :param A: A from PREOS
        :param B: B from PREOS
        :return: u, w, v for PREOS
        """

        u = -1 + B
        w = A - 3 * (B ** 2) - 2 * B
        v = -(A * B) + B ** 2 + B ** 3
        return u, w, v

    @staticmethod
    def __calculate_z_roots__(u, w, v):

        """
        This function is a bit more complicated than the other equations in terms of reasoning for operations.
        The equation used was pg.209, equation 6.4-4

        This equation will solve for Z using numpy to get the roots of the equation, thus the different values
        of Z.

        If Z has only one real root, then it is in only a single state.

        If Z has three real roots, the highest root is vapor and lowest is liquid. The middle root is meta-stable and
        is ignored. If one of the Z roots is negative, that root is also ignored.

        There should not be only two real roots ever. If there is, then there is a bug in the code and an error is
        raised.

        :param u: u from PREOS
        :param w: w from PREOS
        :param v: v from PREOS
        :return: Z root(s) as a dict, defining which state each Z is for
        """

        roots = np.roots([1, u, v, w])
        roots = roots[np.isreal(roots)]
        roots = np.real(roots)
        if len(roots) == 3:
            roots = roots[roots > 0]
            roots = {'vapor': max(roots), 'liquid': min(roots)}
        elif len(roots) == 1:
            root = roots[0]
            if root < 0:
                raise ValueError("No roots for Z found at given conditions")
            roots = {'single state': root}
        elif len(roots) == 0:
            raise ValueError("No roots for Z found at given conditions")
        else:
            raise Exception("Something broke, there should not be only two real roots")
        return roots

    def calculate_V(self, pressure, temperature):

        """
        Calculate molar volume of mixture in mol/m^3 using PREOS

        :param pressure: absolute pressure of mixture in bar
        :param temperature: temperture of mixture in C
        :return: molar volume in m^3/mol for mixture at conditions given
        """

        self.calculate_Z(pressure, temperature)
        for state, z in self.z_mix.items():
            self.V[state] = z * self.R * self.temperature / self.pressure

    def calculate_fv(self, pressure, temperature):

        """
        Goal of this function is to calculate partial vapor fugacites of the different components in the mixture.
        The equation used is in pg.423, equation 9.4-10.

        Most of the values used in the equations are generated when calculating z_mix, but there are a few terms
        that need to be calculated separately.

        SumOf(yjAij) is not the same as SumOf(yjaij), and the SumOf(yjAij) is calculated for each compound
        in __calculate_sum_yj_Aij__.

        And Bi is a single competent value, so that is calculated as well.

        Also, the Z value used should from the vapor phase, and not the liquid phase

        :param pressure: pressure in bar
        :param temperature: temperature in C
        :return: fugactiy in bar
        """

        self.calculate_Z(pressure, temperature)
        z_mix = self.__pull_z_value__(self.z_mix, 'vapor')

        for compound, yi in self.compounds.items():
            ai = self.a_data[compound]
            bi = self.b_data[compound]
            Ai, Bi = self.__calculate_A_B__(ai, bi)
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

        """

        :param z_dict: dict of values containing phases as keys and Z's as values
        :param phase: phase of to get Z for
        :return: Z value in phase
        """

        if phase in z_dict.keys():
            z = z_dict[phase]
        elif 'single state' in z_dict.keys():
            z = z_dict['single state']
        else:
            raise Exception(f"Phase {phase} not found in z dict")
        return z

    def __calculate_sum_yj_Aij__(self, compound_i):

        """
        The equation used does not have a direct equation in the book. This variable is needed when
        calculating partial fugactiy of a compound.

        We know from pg.9.4-9 what aij is, and we can use aij in pg.210, equation 6.4-3 to get Aij.

        :param compound_i:
        :return:
        """

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

    def generate_table_for_fv(self, pressures, temperature, file=None):
        data = []
        for pressure in pressures:
            self.calculate_fv(pressure, temperature)
            row = {'Pressure (bar)': pressure, 'Z': self.__pull_z_value__(self.z_mix, 'vapor')}
            fv_dict_edited = {}
            for compound, value in self.fv.items():
                fv_dict_edited[f'Vapor Fugacity of {compound}'] = value
            row.update(fv_dict_edited)
            data.append(row)
        df = pd.DataFrame(data)
        if file is not None:
            df.to_csv(file, index=False)
        return df


compounds = {'methane': 0.65, 'ethane': 0.20, 'propane': 0.15}
calculator = PendRob(compounds)
calculator.generate_table_for_fv(pressures=[1, 5, 15], temperature=100, file='output.csv')

