"""Importing required libraries"""
from TDN import TDN,PSI
from CoolProp.HumidAirProp import HAPropsSI
import numpy as np
from functions_SkyTemp import T_Sky_Average, dew_point_temperature, T_Sky_Building_Dymoala_BlBd
import math
# import OMPython
# help(OMPython)



"""
Calculations related to the assessment of Radiative Cool (SkyCool) concept

Background:
    - It is based on radiation of the heat to space while at the same time reflecting most of the irradiation received
    on the surface.
    - Ideal material would be a "Black body"(Emissivity of "1") with surface irradiance absorptivity of "0".

Assumptions:
    1-Heat interactions with environment include only:
        a- Sun irradiation (Over year average value)
        b- Heat dissipation to space due to radiation
        c- Convection of heat to ambient(wind)

    2- Calculations done for 1m2 panel
    3- Lumped model heat transfer is assumed ("0D")
    4- Evaporation model uses a simplified model that has dependency only on
        a- Wind velocoty
        b- Plate area
        c- HumRatMax & humRat (kg water/ kg air)

API CoolProp
http://www.coolprop.org/coolprop/HighLevelAPI.html
"""

"""Conversion factors, constants """
T_zero = 273.15

""" Constants """
AntioineA = 8.07131                     # Valid for range 0 t0 100
AntioineB = 1730.63                     # Valid for range 0 t0 100
AntioineC = 233.426                     # Valid for range 0 t0 100
BOLTZMANN_CONSTANT = 1.38064852e-23     # in J/K
const_Stefan_Boltzman = 5.6697e-8       # W/m2/K4
PLANCK_CONSTANT = 6.62607004e-34        # in J.s
SPEED_OF_LIGHT = 299792458              # in m/s


"""Ambient condition"""
T_amb_C = 30             # Ambient temperature [K]
P_amb = 101325.0            # Ambient pressure in [mBar]
V_wind = 1               # Wind velocity [m/s]
relHumP = 50             # Relative humidity (0 to 100) [-]
P_amb_mBar = 995.0          # Ambient pressure in [mBar]
nOpaque = 1              # Cloud coverage (0 to 1) [-] - The lower the more heat transfer
HorIR = 50                 # Horizontal Infrared Radiation [W/m2]
Time_24 = 12.0              # Time in 24-hour format [h]
"""Initialization values """
# # # Temperature of the surface is first assumed a value
# # # !!!!!!!!!!!!!!!! to solve for later!!!!!!!!!!!!!!
C_evap_resistance = 0.8
T_surf_C = 40 # Surface temperature [K]
Q_Rad_ = 0
Q_TRE_ = 0

# T_amb_C, V_wind, relHumP, nOpaque, HorIR, Time_24, T_surf_C, Q_Rad_, Q_TRE_ =  35, 1, 12, 0.05, 0.0, 12, 42, 210, 300
# T_amb_C, V_wind, relHumP, nOpaque, HorIR, Time_24, T_surf_C, Q_Rad_, Q_TRE_ = 36, 3, 27, 0.05, 0.0, 14, 47, 170, 300
# T_amb_C, V_wind, relHumP, nOpaque, HorIR, Time_24, T_surf_C, Q_Rad_, Q_TRE_ = 31, 3, 60, 0.05, 0.0, 12, 50, 150, 300

"""CALCULATED  Ambient condition"""
T_amb = T_amb_C + T_zero
T_surf = T_surf_C + T_zero
relHum = relHumP/100
T_dew = dew_point_temperature(T_amb-T_zero, relHumP)
# Calculate the saturation vapor pressure
PSat = 10 ** (AntioineA - (AntioineB / (AntioineC + T_amb - T_zero)))
Pw = ((relHum / 100) * PSat)
# T_Sky = T_Sky_Building_Dymoala_BlBd(nOpaque, HorIR, T_amb-T_zero, T_dew)


""" Location Input parameters """
elevation_angle = 30  # in degrees
atmospheric_pressure = 101325  # in Pa
wavelength = 9.6e-6  # in meters


"""Surface properties"""
# a_plt = 0.2 * 0.2  # Assuming unit surface area
a_plt = 1  # Assuming unit surface area


"""Coating material selection"""
absrb = 0.03  # Absorptivity of Manganese Oxide Ref: solarmirror.com
# emsvt = 0.1  # Emissivity of Manganese Oxide  Ref: solarmirror.com
emsvt = 0.8  # Emissivity random value
k = 42  # Thermal conductivity of Manganese Oxide Ref: crystran.co.uk

""" Heat absorption from the sun irradiation average has been extracted from SolarAtlas """
irr_dk = 255  # [W/m2] Total irradiance average yearly for Denmark Ref: SolarAtlas
# # # World average energy balance is 185.9 according to NASA earths energy budget
q_sun = irr_dk * absrb * a_plt

"""T_SKY FROM FUNCTIONS"""
T_Sky = T_Sky_Average(Time_24, P_amb_mBar, nOpaque, HorIR, T_amb-T_zero, T_dew, Pw, irr_dk) + T_zero  # Surrounding temperature from ambient [K] - based on yearly average Appl. Sci. 2020, 10, 8057; doi:10.3390/app10228057


""" Heat dissipation to space can e modeled using the radiative heat transfer based on the emissivity """
boltz = 5.67e-8 # Boltzman constant [W/m2/K4]
q_rad = boltz * emsvt * ((T_surf) ** 4 - (T_Sky) ** 4) * a_plt  # Radiation heat to space [W/m2]


""" Heat absorption from environment through forced/natural convection (wind speed 4 m/s) """
g = 9.8  # Gravitational acceleration [m/s2]
beta = 1 / T_amb  # Coefficiat of volume expansion [1/K] , (Ideal gas assumption)

t_inf = T_Sky  # For Grashof number calc
d_air = TDN(P_amb, 0, T_amb, 0, 0, 'Air').d
mu_air = PSI('V', 'T', T_amb, 'P', P_amb, 'Air')
c_air = PSI('C', 'T', T_amb, 'P', P_amb, 'Air')
k_air = PSI('L', 'T', T_amb, 'P', P_amb, 'Air')
l_re = 1  # Unit characteristic length [m]
v_k = k_air/d_air   # Kinematic viscosoty = mu/rho  []

reynolds = (d_air * V_wind * l_re) / mu_air
prantl = mu_air*c_air/k
nusselt = 0.664*reynolds**(0.5) *prantl**(1/3)  #  For plate with Reynolds<5e5
grashof = g * beta * (T_surf - T_amb) * l_re ** 3 / v_k ** 2
richarson =  grashof/reynolds**2
h_plt_air = k_air*nusselt/l_re
q_h_plt_air = h_plt_air * (T_amb - T_surf) * a_plt

""" Heat dissipation through evaporation """
humRat = HAPropsSI('W','P', P_amb, 'T', T_amb, 'R', relHum)  # kg water/kg dry air
humRatMax = HAPropsSI('W','P', P_amb, 'T', T_amb, 'R', 1)  # kg water/kg dry air
gramEvapW = (25 + 19 * V_wind) * a_plt * (humRatMax - humRat) / 3600 * (1-C_evap_resistance)
cEvapW = PSI('H', 'P', P_amb, 'Q', 1, "Water") - PSI('H', 'P', P_amb, 'Q', 0, "Water")
q_evap = gramEvapW * cEvapW


print("_______________________________________________________________________________________________________________\n")


print(
    " Ambient temperature [K]               :  ", T_amb,
    "\n Ambient pressure in [mBar]            :  ", P_amb_mBar,
    "\n Cloud coverage (0 to 1) [-]           :  ", nOpaque,
    "\n Relative humidity (0 to 100) [-]      :  ", relHumP,
    "\n Horizontal Infrared Radiation [W/m2]  :  ", HorIR,
    "\n Time in 24-hour format [h]            :  ", Time_24,
)
print("_______________________________________________________________________________________________________________\n")
print("Variables & Parameters: ")
print("       Relative humidity  ", relHum*100 , " %")
print("       Ambient Temp       ", T_amb - T_zero, " degC")
print("       Dew Temp           ", T_dew, " degC")
print("       Wind velocity      ", V_wind, " m/s")


print("_______________________________________________________________________________________________________________\n")
print("Heat transfer dimensionless numbers: ")
print("     - Reynolds number - for convection heat transfer =           ", round(reynolds,2))
print("     - Grashof number - for free OR forced convection check =     ", round(grashof,3))
print("     - Richardson number - for free OR forced convection check =  ", round(richarson,10))
print("     - Convection heat transfer coeficient                        ", round(h_plt_air,8))
print("     Q: Is it forced connections?                                 ", richarson<1e-1)

print("_______________________________________________________________________________________________________________\n")
print("Heat transfer under lumped assumption: ", "Units are in [W] for surface area of", round(a_plt,6), " [m^2] \n")

print("     - Q_rad_space - Dissipated radiation heat to space =                   ", round(-q_rad,4), " vs. ", -Q_Rad_ )
print("     - Q_conv_plt_air - Absorbed heat through convection from ambient air = ", round(q_h_plt_air*2,4))
print("     - Q_rad_sun - Absorbed radiation heat from sun(yearly average) =       ", round(q_sun,4))
print("     - Q_Evap_W - Heat dissipated due to evaporation =                      ", round(-q_evap,4), "vs. ", -(Q_TRE_-Q_Rad_))


print("     - Total heat balance of  =                                             ", round(q_h_plt_air*2+q_sun-q_rad-q_evap, 4)," vs. ", -Q_TRE_)



					