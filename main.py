import pandas as pd
import numpy as np


class PendRob:

    def __init__(self):

        """
        Goal is to solve for Z and V in PREOS

        where:
        Z^3 + uZ^2 + wZ + v = 0

        u = -1 + B
        w = A - 3B^2 - 2B
        v = -AB + B^2 + B^3

        A = aP / (RT)^2
        B = bP / RT

        Z = PV / RT
        """

        # Gather compounds, make sure they are correct...
        # self.compounds = {'oxygen': 20, 'nitrogen': 60, 'water': 20}
        self.compounds = {'methane': 1}
        self.a_data = dict(zip(self.compounds.keys(), [np.nan] * len(self.compounds)))
        self.b_data = dict(zip(self.compounds.keys(), [np.nan] * len(self.compounds)))

        self.R = 8.314  # In Pa * m^3 * K^-1 * mol^-1

        self.raw_Tc_Pc_w_data = pd.read_csv('preso_ab.csv')
        self.raw_Tc_Pc_w_data = self.raw_Tc_Pc_w_data.rename(columns={'Tc (K)': 'Tc', 'Pc (MPa)': 'Pc'})

        self.a_mix = np.nan
        self.b_mix = np.nan

    def __calculate_a_b__(self, temperture, compound):
        compound_data = self.raw_Tc_Pc_w_data[self.raw_Tc_Pc_w_data['compound'] == compound]
        compound_data = compound_data.to_dict('records')[0]
        Tc = compound_data['Tc']
        Pc = compound_data['Pc'] * 10**6  # Convert MPa to Pa
        w = compound_data['w']

        k = 0.37464 + 1.54226*w - 0.26992*(w**2)
        alpha = (1 + k*(1 - (temperture/Tc)**(1/2)))**2

        a = 0.45724 * alpha * (self.R*Tc)**2/Pc
        b = 0.07780 * self.R*Tc/Pc

        return a, b

    def __calculate_a_b_mix__(self):

        k = 0
        a_mix = 0
        b_mix = 0

        for compound, y1 in self.compounds.items():
            b_mix += y1 * self.b_data[compound]
            for sub_compound, y2 in self.compounds.items():
                a1 = self.a_data[compound]
                a2 = self.a_data[sub_compound]
                a_mix += (a1*a2)**(1/2)

        return a_mix, b_mix

    def calculate_Z(self, pressure, temperature):

        """
        :param pressure: absolute pressure of gas in bar
        :param temperature: temperture of gas in C
        :return: z constant for gas at conditions given
        """

        pressure, temperature = self.__convert_variables__(pressure, temperature)
        if pressure < 0:
            raise ValueError("Absolute pressure can not less than 0")
        if temperature < 0:
            raise ValueError("Temperature cannot be less than 0K")

        for compound in self.compounds:
            self.a_data[compound], self.b_data[compound] = self.__calculate_a_b__(temperature, compound)

        self.a_mix, self.b_mix = self.__calculate_a_b_mix__()

        A, B = self.__calculate_A_B__(pressure, temperature)
        u, w, v = self.__calculate_u_w_v__(A, B)
        z = self.__calculate_z_roots__(u, w, v)
        return z

    @staticmethod
    def __convert_variables__(pressure, temperature):
        pressure = pressure * 10**5  # Convert bar to Pa
        temperature = temperature + 273.15  # Convert C to K
        return pressure, temperature

    def __calculate_A_B__(self, P, T):
        A = (self.a_mix * P) / ((self.R * T) ** 2)
        B = (self.b_mix * P) / (self.R * T)
        return A, B

    @staticmethod
    def __calculate_u_w_v__(A, B):
        u = -1 + B
        w = A - 3*(B**2) - 2*B
        v = -(A * B) + B**2 + B**3
        return u, w, v

    @staticmethod
    def __calculate_z_roots__(u, w, v):
        roots = np.roots([1, u, v, w])
        roots = roots[np.isreal(roots)]
        roots = np.real(roots)
        if len(roots) == 3:
            roots = {'gas': max(roots), 'liquid': min(roots)}
        elif len(roots) == 1:
            roots = {'single state': roots[0]}
        else:
            raise Exception("Something broke, there should not be only two real roots")
        return roots

    def calculate_V(self, pressure, temperature):

        """
        :param pressure: absolute pressure of gas in bar
        :param temperature: temperture of gas in C
        :return: molar volume in m^3/mol for gas at conditions given
        """

        all_z = self.calculate_Z(pressure, temperature)
        pressure, temperature = self.__convert_variables__(pressure, temperature)
        V = {}
        for state, z in all_z.items():
            V[state] = z * self.R * temperature / pressure
        return V


calculator = PendRob()
value = calculator.calculate_V(pressure=0.1, temperature=300)
print(value)
