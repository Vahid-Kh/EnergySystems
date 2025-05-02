
import math as mt
from CoolProp.CoolProp import PropsSI as PSI

""" TDN class is the thermodynamic state class which at list two non-zero(Except for vapor quality which can be 
0 <= Q <= 1) parameter to find the all other thermodynamic properties.
 For all unknown values enter 0 and for quality"""

"""Reference condition for exergy calculation : """
t0 = 293.15
p0 = 101325


class TDN:

    def __init__(self, p, h, t, s, d, r):

        if p != 0 and h != 0:
            self.p = p                            # pressure in pa
            self.h = h                            # enthalpy in j/kg
            self.t = PSI('T', 'H', h, 'P', p, r)  # temperature in k
            self.s = PSI('S', 'H', h, 'P', p, r)  # entropy in j/kg/k
            self.d = PSI('D', 'H', h, 'P', p, r)  # density in kg/m3
            self.r = r                            # Refrigerant

        elif p != 0 and t != 0:
            self.p = p                            # pressure in pa
            self.h = PSI('H', 'T', t, 'P', p, r)  # enthalpy in j/kg
            self.t = t                            # temperature in k
            self.s = PSI('S', 'T', t, 'P', p, r)  # entropy in j/kg/k
            self.d = PSI('D', 'T', t, 'P', p, r)  # density in kg/m3
            self.r = r                            # Refrigerant

        elif p != 0 and s != 0:
            self.p = p                            # pressure in pa
            self.h = PSI('H', 'S', s, 'P', p, r)  # enthalpy in j/kg
            self.t = PSI('T', 'S', s, 'P', p, r)  # temperature in k
            self.s = s                            # entropy in j/kg/k
            self.d = PSI('D', 'S', s, 'P', p, r)  # density in kg/m3
            self.r = r                            # Refrigerant

        elif p != 0 and d != 0:
            self.p = p                            # pressure in pa
            self.h = PSI('H', 'D', d, 'P', p, r)  # enthalpy in j/kg
            self.t = PSI('T', 'D', d, 'P', p, r)  # temperature in k
            self.s = PSI('S', 'D', d, 'P', p, r)  # entropy in j/kg/k
            self.d = d                            # density in kg/m3
            self.r = r                            # Refrigerant

        elif h != 0 and t != 0:
            self.p = PSI('P', 'H', h, 'T', t, r)  # pressure in pa
            self.h = h                            # enthalpy in j/kg
            self.t = t                            # temperature in k
            self.s = PSI('S', 'H', h, 'T', t, r)  # entropy in j/kg/k
            self.d = PSI('D', 'H', h, 'T', t, r)  # density in kg/m3
            self.r = r                            # Refrigerant

        elif h != 0 and s != 0:
            self.p = PSI('P', 'H', h, 'S', s, r)  # pressure in pa
            self.h = h                            # enthalpy in j/kg
            self.t = PSI('T', 'H', h, 'S', s, r)  # temperature in k
            self.s = s                            # entropy in j/kg/k
            self.d = PSI('D', 'H', h, 'S', s, r)  # density in kg/m3
            self.r = r                            # Refrigerant

        elif h != 0 and d != 0:
            self.p = PSI('P', 'H', h, 'D', d, r)  # pressure in pa
            self.h = h                            # enthalpy in j/kg
            self.t = PSI('T', 'H', h, 'D', d, r)  # temperature in k
            self.s = PSI('S', 'H', h, 'D', d, r)  # entropy in j/kg/k
            self.d = d                            # density in kg/m3
            self.r = r                            # Refrigerant

        elif t != 0 and s != 0:
            self.p = PSI('P', 'T', t, 'S', s, r)  # pressure in pa
            self.h = PSI('H', 'T', t, 'S', s, r)  # enthalpy in j/kg
            self.t = t                            # temperature in k
            self.s = s                            # entropy in j/kg/k
            self.d = PSI('D', 'T', t, 'S', s, r)  # density in kg/m3
            self.r = r                            # Refrigerant

        elif t != 0 and d != 0:
            self.p = PSI('P', 'T', t, 'D', d, r)  # pressure in pa
            self.h = PSI('H', 'T', t, 'D', d, r)  # enthalpy in j/kg
            self.t = t                            # temperatute in k
            self.s = PSI('S', 'T', t, 'D', d, r)  # entropy in j/kg/k
            self.d = d                            # density in kg/m3
            self.r = r                            # Refrigerant

        elif s != 0 and d != 0:
            self.p = PSI('P', 'S', s, 'D', d, r)  # pressure in pa
            self.h = PSI('H', 'S', s, 'D', d, r)  # enthalpy in j/kg
            self.t = PSI('T', 'S', s, 'D', d, r)  # temperature in k
            self.s = s                            # entropy in j/kg/k
            self.d = d                            # density in kg/m3
            self.r = r                            # Refrigerant

        else:
            self.p = p                            # pressure in pa
            self.h = h                            # enthalpy in j/kg
            self.t = PSI('T', 'H', h, 'P', p, r)  # temperature in k
            self.s = PSI('S', 'H', h, 'P', p, r)  # entropy in j/kg/k
            self.d = PSI('D', 'H', h, 'P', p, r)  # density in kg/m3
            self.r = r                            # Refrigerant

            """Not enough argument provided or more argument than required entered" + 
            "You should only provide 2 parameter and name of refrigerant as the last parameter"""

    def __copy__(self):
        return TDN(self.p, self.h, self.t, self.s, self.d, self.r)

    def dp(self, dp):                              # Pressure drop function
        return self.p-dp

    def dt(self, dt):                              # Temperature  variations function
        return self.t + dt

    def comp(self, eta, h_isen):                   # Compressor function
        return self.h + (h_isen-self.h) / eta

    def cp(self):                                  # Heat capacity at constant pressure
        return PSI('CPMASS', 'H', self.h, 'P', self.p, self.r)

    def cv(self):                                  # Heat capacity at constant Volume
        return PSI('CVMASS', 'H', self.h, 'P', self.p, self.r)

    def ex(self):                                  # Exergy
        return self.h - PSI('H', 'T', t0, 'P', p0, self.r) - t0 * (self.s - PSI('S', 'T', t0, 'P', p0, self.r))

    def z(self):                                   # Compressibility facto
        return PSI('Z', 'T', self.t, 'P', self.p, self.r)

    def m(self):                                   # Molecular weight
        return PSI('M', 'T', self.t, 'P', self.p, self.r)

    def u(self):                                   # Internal energy
        return PSI('U', 'T', self.t, 'P', self.p, self.r)

    def prop(self):                                # Prints all properties in list form
        return [round(self.p, 5), round(self.h/1e3, 5), round(self.t, 5), round(self.s/1e3, 5), round(1/self.d, 5), round(self.q(), 5), self.r, round(self.ex()/1e3, 5)]

    def propl(self):                                # Prints all properties in list form
        return [round(self.p, 5), round(self.h/1e3, 5), round(self.t, 5), round(self.s/1e3, 5), round(1/self.d, 5),  self.r]

    def propex(self):                                # Prints all properties in list form
        return [round(self.p, 5), round(self.h/1e3, 5), round(self.t, 5), round(self.s/1e3, 5), round(self.ex()/1e3, 5), round(1/self.d, 5),  self.r]

    def q(self):
        return PSI('Q', 'H', self.h, 'P', self.p, self.r)

    def propf(self):
        return [self.p, self.h/1e3, self.t, self.s/1e3, 1/self.d, self.r
        , PSI('Q', 'H', self.h, 'P', self.p, self.r)                  # Quality of vapor
        , PSI('C', 'H', self.h, 'P', self.p, self.r)                  # Heat capacity in (j/kg/k)
        , PSI('U', 'H', self.h, 'P', self.p, self.r)                  # Internal energy in (j/kg)
        , PSI('V', 'H', self.h, 'P', self.p, self.r)                  # Viscosity in (Pa.sec)
        , PSI('L', 'H', self.h, 'P', self.p, self.r)]                 # Thermal conductivity in (W/m/K)]


"""'All TDN states are defined by pressure as p, and enthalpy as h and the selected refrigerant'
      'Tref=293.15   Kpref=101325pa  href=0Jkg−1   sref=0Jkg−1K−1')
"""
""" TDNF class is the thermodynamic state class full parameters  C is heat capacity in (j/kg/k), U is internal energy 
in (j/kg),  V is viscosity in (Pa.sec) and L is thermal conductivity in (W/m/K) """


class TDNQ(TDN):

    def __init__(self, p, h, t, s, d, r, q):
        super().__init__(p, h, t, s, d, r)

        if q == '-':
            self.q = PSI('Q', 'H', h, 'P', p, r)  # Quality of vapor
            self.c = PSI('C', 'H', h, 'P', p, r)  # Heat capacity in (j/kg/k)
            self.u = PSI('U', 'H', h, 'P', p, r)  # Internal energy in (j/kg)
            self.v = PSI('V', 'H', h, 'P', p, r)  # Viscosity in (Pa.sec)
            self.l = PSI('L', 'H', h, 'P', p, r)  # Thermal conductivity in (W/m/K)

        else:
            self.q = q                            # Quality of vapor
            self.c = PSI('C', 'H', h, 'P', p, r)  # Heat capacity in (j/kg/k)
            self.u = PSI('U', 'H', h, 'P', p, r)  # Internal energy in (j/kg)
            self.v = PSI('V', 'H', h, 'P', p, r)  # Viscosity in (Pa.sec)
            self.l = PSI('L', 'H', h, 'P', p, r)  # Thermal conductivity in (W/m/K)


def mlog(a, b):

    return (a - b) / (mt.log(a / b))


def sort(sub):
    # reverse = None (Sorts in Ascending order)
    # key is set to sort using second element of
    # sublist lambda has been used
    sub.sort(key=lambda x: x[1], reverse=True)
    return sub


tdn700 = TDN(70000000, 0, 300,0, 0, "Hydrogen")

# print(tdn700.propl())
tdn20 = TDN(2000000, tdn700.h, 0,0, 0, "Hydrogen")
# print(tdn20.cp())
# print(tdn20.propl())
tdn3= TDN(300000,tdn20.h,0,0,0,"Hydrogen")
# print(tdn3.propl())
# print(tdn3.cp())

# print("Cooled down at 20 bar")
tdn20cool = tdn20 = TDN(2000000, 0, 300,0, 0, "Hydrogen")
tdn3= TDN(300000,tdn20cool.h,0,0,0,"Hydrogen")
# print(tdn3.propl())
# print(tdn3.cp())

print("-------Compression power calculated here ------")
p_in = 15e5 # bar
p_out = 100e5  # bar
eta_comp = 0.3
# eta_comp = 1
print("p_in  "  ,p_in,"   |  p_out   ",p_out)
# print(p_out)
tdn_in = TDN(p_in, 0, 300, 0, 0, "Hydrogen")
tdn_out_isen = TDN(p_out, 0, 0, tdn_in.s, 0, "Hydrogen")
print("Density difference " , tdn_in.d, tdn_out_isen.d)
print(tdn_in.h, (TDN(p_in, 0, 0, tdn_in.s, 0, "Hydrogen")).h)
# print((tdn_in.h - (TDN(p_in , 0, 0, tdn_in.s, 0, "Hydrogen")).h) / 0.8)
# print("This is a test: ", tdn_in.h, tdn_out_isen)
tdn_out_real = TDN(p_out, (tdn_in.h + (-tdn_in.h + (TDN(p_out, 0, 0, tdn_in.s, 0, "Hydrogen")).h) / eta_comp), 0, 0, 0, "Hydrogen")
# print("p_in  "  ,p_in,"   |  p_out   ",p_out)
# print("heat avilable from out comp t to 300k   ",(tdn_out_real.t - 300) * tdn_out_real.cp())
# print("Electricity supplied to compressor",(tdn_out_real.t - 300) * tdn_out_real.cp()*(1/(1-eta_comp)) )
# print("Temp out   ", tdn_out_real.t)
print("p_in="  ,p_in/1e5, "|  p_out=",p_out/1e5, " | heat",round(((tdn_out_real.t - 300) * tdn_out_real.cp())/1e6)," | Elec=",round(((tdn_out_real.t - 300) * tdn_out_real.cp()*(1/(1-eta_comp)) )/1e6), "| Eta_Isen = ", eta_comp)

tdn8 = TDN(800000,0,300,0,0,"Hydrogen")
# print((TDN(35000000,0,0,tdn20cool.s,0, "Hydrogen")).h -tdn8.h, TDN(35000000,0,0,tdn20cool.s,0, "Hydrogen").t)
# print(TDN(100000,0,273.15,0,0,"Hydrogen").d * 500)









