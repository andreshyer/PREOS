import pandas as pd
import numpy as np


class PendRob:

    def __init__(self, compound):

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

        :param compound: compound to calculate V and Z for
        """

        self.compound = compound

        self.ab_data = pd.read_csv('preso_ab.csv')
        self.ab_data = self.ab_data.rename(columns={'a (Pa*m^6*mol^-2)': 'a',
                                                    'b [(m^3/mol)*10^5]': 'b'})

        if self.compound not in list(self.ab_data['compound']):
            compounds = ', '.join(list(self.ab_data['compound']))
            raise ValueError(
                f"Compound {self.compound} was not found in spreadsheet.\n"
                f"Available compounds are: {compounds}\n"
                "Please check spelling or add compound to spreadsheet."
                )

        self.R = 8.314  # In Pa * m^3 * K^-1 * mol^-1
        self.a, self.b = self.__fetch_a_b__()

    def __fetch_a_b__(self):
        ab_data = self.ab_data[self.ab_data['compound'] == self.compound]
        ab_data = ab_data.to_dict('records')[0]
        a = ab_data['a']
        b = ab_data['b']
        return a, b

    def calculate_z(self, pressure, temperature):

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

        A, B = self.__calculate_A_B__(pressure, temperature)
        u, w, v = self.__calculate_u_w_v__(A, B)
        z = self.__calculate_z__(u, w, v)
        return z

    def calculate_V(self, pressure, temperature):
        z = self.calculate_z(pressure, temperature)
        V = z * pressure * temperature / self.R
        return V

    @staticmethod
    def __convert_variables__(pressure, temperature):
        pressure = pressure * 10**5  # Convert bar to Pa
        temperature = temperature + 273.15  # Convert C to K
        return pressure, temperature

    def __calculate_A_B__(self, P, T):
        A = (self.a * P) / (self.R * T)**2
        B = (self.b * P) / (self.R * T)
        return A, B

    @staticmethod
    def __calculate_u_w_v__(A, B):
        u = -1 + B
        w = A - 3*(B**2) - 2*B
        v = -(A * B) + B**2 + B**3
        return u, w, v

    @staticmethod
    def __calculate_z__(u, w, v):
        roots = np.roots([1, u, v, w])
        roots = roots[np.isreal(roots)]
        roots = np.real(roots)
        return max(roots)


oxygen = PendRob(compound='helium')
value = oxygen.calculate_V(pressure=0.05, temperature=300)
print(value)
