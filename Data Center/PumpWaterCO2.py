# encoding: utf-8

from scipy import *
from CoolProp.CoolProp import PropsSI as PSI
from CoolProp.CoolProp import Props as PS
from ht import *   #  Dont use PIP but use add library and find "ht"
import CoolProp
import math as mt
from CoolProp.CoolProp import PropsSI as PSI
from TDN import *


tk = 273.15
g = 9.81  # Gravity constant
l = 5000  # meter
r = 0.5  # meter

"""Load assumption"""
q_heat = 100e6  # Heat from data center 100MW

"""Pressure definition"""
# p_co2, r_co2_in, r_co2_out = 0.69e7, 0.5, 0.5  # Pressure of CO2 [Pa]
p_co2, r_co2_in, r_co2_out = 1e7, 0.5, 0.5  # Pressure of CO2 [Pa]
p_water = 5e5  # Pressure of water [Pa]

r_water = 0.5  # meter

t_co2_in = 28  # temp degC
t_co2_out = 43  # temp degC

t_water_in = 28  # temp degC
t_water_out = 43  # temp degC

q_vol = 5700  # m3/h

co2_in = TDN(p_co2, 0, tk + t_co2_in, 0, 0, "co2")
co2_out = TDN(p_co2, 0, tk + t_co2_out, 0, 0, "co2")

water_in = TDN(p_water, 0, tk + t_water_in, 0, 0, "water")
water_out = TDN(p_water, 0, tk + t_water_out, 0, 0, "water")

m_dot_co2 = q_heat /(co2_out.h-co2_in.h)
m_dot_water = q_heat /(water_out.h-water_in.h)

""" kwH/lit for 2Ø CO2 6 1Ø Water"""
co2_spes_heat = ((co2_out.h-co2_in.h)/co2_out.d)*1000
water_spes_heat = ((water_out.h-water_in.h)/water_out.d)*1000

q_vol_co2_in = m_dot_co2 /  co2_in.d *3600     # m3/h
q_vol_co2_out = m_dot_co2 / co2_out.d *3600    # m3/h
q_vol_water  = m_dot_water / water_out.d *3600  # m3/h

mu_co2 = PSI("V", "P", p_co2, "T", tk+t_co2_in, "CO2")
mu_water = PSI("V", "P", p_water, "T", tk+t_water_in, "H2O")  # [Pa-s]

rho_co2 = co2_in.d
rho_water = water_in.d

"""https://chemenggcalc.com/hagen-poiseuille-equation-calculator/"""
"""The pressure drop in a pipe can be expressed by Hagen–Poiseuille equation:
Δ𝑃=(8𝜇𝐿𝑄)/(𝜋𝑅^(4))
Q: Volumetric flow rate (in cubic meters per second, m3/s)r: Pipe radius (in meters, m)ΔP: Pressure difference across the length of the pipe (in pascals, Pa)η: Dynamic viscosity of the fluid (in pascals-seconds, Pa·s)L: Pipe length (in meters, m)
"""
dp_co2_out = (8 * mu_co2 * l * q_vol_co2_out) / (mt.pi * r_co2_out ** (4))
dp_co2_in = (8 * mu_co2 * l * q_vol_co2_in) / (mt.pi * r_co2_in ** (4))
dp_water = (8 * mu_water * l * q_vol_water) / (mt.pi * r_water ** (4))

h_co2 = dp_co2_in/1e4
h_water = dp_water/1e4

h = 50  # meter
efficiency = 0.6  # Pump to shaft efficiency

"""https://www.engineeringtoolbox.com/pumps-power-d_505.html"""
p_pump_co2 = (q_vol_co2_out * h_co2 * rho_co2 * g) / (3.6 * 1e6)
p_pump_water = (q_vol_water * h_water * rho_water * g) / (3.6 * 1e6)

print("mass flow co2 & water   ",m_dot_co2, m_dot_water)
print("Spesific heat in [kWh/lit]   ", co2_spes_heat/ 3.6e6,water_spes_heat/ 3.6e6)
print("pressure drop ", dp_co2_out,dp_co2_in,dp_water )
print("Volumetric flow rate", q_vol_co2_in, q_vol_co2_out, q_vol_water  )
print("H CO2 vs Water : ", h_co2, h_water)
print("Pump power for CO2  ", p_pump_co2, "  [kW]  | Calculating shaft power    ",  p_pump_co2/ efficiency)
print("Pump power for Water  ", p_pump_water, "  [kW]  | Calculating shaft power   ", p_pump_water/ efficiency)
print("# Viscosity in (Pa.sec)", "Water -->  ", mu_water,  "      CO2 -->", mu_co2 )
