# encoding: utf-8

from TDN import *
import CoolProp
from CoolProp.CoolProp import Props as PS
from CoolProp.CoolProp import PropsSI as PSI
from ht import *   #  Dont use PIP but use add library and find "ht"
import math as mt
import pandas as pd
from tabulate import tabulate


CtoK = 273.15
"""Testing TDN based on CoolProp"""
# tdn40 = TDN(4000000,0,300,0,0,"Hydrogen")
# tdn700 = TDN(70000000, 0, 300,0, 0, "Hydrogen")
# print(TDNex(10e5, 3e5, 0, 0, 0, 'co2', 300).exergy())
# print(CoolProp.CoolProp.Props("V", "P", 1e7, "T", 350, "CO2"))



""" Critical Heat Flux for CO2 boiling and water at 30 degC 
#  Parameters
#  sigma [float] Surface tension of liquid [N/m]
#  Hvap [float] Heat of vaporization of the fluid at P, [J/kg]
#  rhol [float] Density of the liquid [kg/m^3]
#  rhog [float] Density of the produced gas [kg/m^3]
#  K [float] Constant []
#  Returns
#  q: float Critical heat flux [W/m^2]
"""

t_boil = 24

row = 1
for media in ["Water", "CO2", "R12", "R717"]:
    h_v         = PSI('H','T',t_boil+CtoK,'Q',1,media)   # Saturated vapor enthalpy of Water at 1 atm in J/kg
    h_l         = PSI('H','T',t_boil+CtoK,'Q',0,media)   # Saturated liquid enthalpy of Water at 1 atm in J/kg
    rho_v       = PSI('D','T',t_boil+CtoK,'Q',1,media)   # [float] Density of the gas [kg/m^3]
    rho_l       = PSI('D','T',t_boil+CtoK,'Q',0,media)   # [float] Density of the liquid [kg/m^3]
    h_vap       =  h_v - h_l # Latent heat of vaporization of Water at 1 atm in J/kg
    k_l         = PSI('L','T',t_boil+CtoK,'Q',0,media)   #  kl [float] Thermal conductivity of liquid [W/m/K]
    surf_ten    = PSI('I','T',t_boil+CtoK,'Q',0,media)  # 	Surface Tension [N/m]
    q_crit      = boiling_nucleic.Zuber(surf_ten, h_vap, rho_l, rho_v, K=0.18)    # Critical heat flux boiling
    lambda_hx   = 0.01   #  lambda is the wavelength of the corrugations
    # d_hyd       = area/perimeter   #  dh = hydraulic diameter (m
    cp_v        = PSI('C','T',t_boil+CtoK,'Q',1,media)   #
    cp_l        = PSI('C','T',t_boil+CtoK,'Q',0,media)   #
    mu          = PSI("V", "T", t_boil+CtoK, "Q", 0, media)  # [Pa-s]
    h_boiling   = h_boiling_Han_Lee_Kim(m=3E-5, x=.4, Dh=0.002, rhol=rho_l , rhog=rho_v , kl=0.086, mul=156E-6, Hvap=h_vap,
                                      Cpl=cp_l, q=1E5, A_channel_flow=0.0003, wavelength=3.7E-3,
                                      chevron_angle=45)  # h [float] Boiling heat transfer coefficient [W/m^2/K]
    #           =
    #           =
    #           =
    #           =
    #Nu           =h*l/k
    #Pr           =mu*c/k
    #           =
    #           =

    data = [{
        'Media':                media,
        'rho_v kg/m3':          rho_v,
        'rho_l kg/m3':          rho_l,
        'h_vap J/kg':           h_vap,
        'mu PaSec':             mu,
        # 'k_l W/m/K':          k_l,
        # 'surf_ten N/m':       surf_ten,
        'q_crit W/cm2':         q_crit*1e-4,
        'h_boil W/cm2/K':       h_boiling*1e-4,
        'h_vap_vol J/m3':       h_vap*rho_v, #
        # ' ': ,   #
        # ' ': ,   #
        # ' ': ,   #
        }]
    if row == 1:
        df_med = pd.DataFrame(columns=data[0].keys())
    values_list = list(data[0].values())
    new_row = pd.Series(values_list, index=data[0].keys())
    df_med_row = pd.DataFrame(data)
    df_med = pd.concat([df_med, new_row.to_frame().T], ignore_index=True)
    row += 1
    # print(" -- For ", media, " the Critical Heat Flux --> ", round(q_crit, 2))

row = 1
for material in ['stainless steel', 'extruded aluminium 6061', 'copper alloy', 'carbon']:
    near_mat    = nearest_material(material)
    k_mat       = k_material(near_mat)
    rho_mat     = rho_material(near_mat)
    cp_mat      = Cp_material(near_mat)
    #           =
    #           =
    #           =
    #Nu           =h*l/k
    #Pr           =mu*c/k
    #           =
    #           =

    data = [{
        'Material':             near_mat,
        'k_mat' :               k_mat,
        'rho_mat' :             rho_mat,
        'cp_mat' :              cp_mat,
        # ' ': ,   #
        # ' ': ,   #
        # ' ': ,   #
        }]
    if row == 1:
        df_mat = pd.DataFrame(columns=data[0].keys())
    values_list = list(data[0].values())

    new_row = pd.Series(values_list, index=data[0].keys())
    df_mat_row = pd.DataFrame(data)
    df_mat = pd.concat([df_mat, new_row.to_frame().T], ignore_index=True)
    row += 1
    # print(" -- For ", media, " the Critical Heat Flux --> ", round(q_crit, 2))

print(tabulate(df_med, df_med.columns, numalign="center", tablefmt="grid"))
print(tabulate(df_mat, df_mat.columns, numalign="center", tablefmt="grid"))





"""Testing HT(heat transfer)  library """

"""https://pypi.org/project/thermo/#what-is-thermo  """
from thermo.chemical import Chemical
tol = Chemical('toluene')
print(tol.rhog, tol.Cpg, tol.kg, tol.mug)

from thermo import Mixture
vodka = Mixture(['water', 'ethanol'], Vfls=[.6, .4], T=300, P=1E5)
print(vodka.Prl,vodka.Prg)
air = Mixture('air', T=400, P=1e5)
print(air.Cp)


"""https://pypi.org/project/thermochem/"""
import pytest
from thermochem.combustion import SimpleCombustor, Combustor
from thermochem.burcat import Elementdb

"""http://pyromat.org/features.html"""
# import pyromat as pm
# air = pm.get('ig.air')
# air.gam(T=[300, 400, 500, 600])
# print(air.gam(T=[300,400]))

from EESConnect import EESConnector
from tkinter import filedialog
import tkinter as tk

# Select the EES file path
# root = tk.Tk()
# root.withdraw()
# ees_file_path = filedialog.askopenfilename()

# with EESConnector() as ees:
#     ees.ees_file_path = ees_file_path
#     result = ees.calculate(["air_ha", 110, 1013.25])
#     print(result[1])
