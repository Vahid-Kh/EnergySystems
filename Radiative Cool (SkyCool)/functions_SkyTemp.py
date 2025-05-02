import math
import numpy as np
import pandas as pd
import pvlib
from pvlib.tools import sind
from pvlib._deprecation import warn_deprecated
from pvlib.tools import _get_sample_intervals
import scipy
import scipy.constants
import warnings



# Input parameters
elevation_angle = 30  # in degrees
atmospheric_pressure = 101325  # in Pa
wavelength = 9.6e-6  # in meters



# Constants
BOLTZMANN_CONSTANT = 1.38064852e-23  # in J/K
const_Stefan_Boltzman = 5.6697e-8  # W/m2/K4
PLANCK_CONSTANT = 6.62607004e-34  # in J.s
SPEED_OF_LIGHT = 299792458  # in m/s
EMISSIVITY = 0.7
T_zero = 273.15

# Convert temperature to Kelvin


def dew_point_temperature(t, rh):
    """
    Calculates the dew point temperature in Celsius using the Arden Buck equation.
    Args:
    t: float, temperature in Celsius
    rh: float, relative humidity in percent
    Returns:
    float, dew point temperature in Celsius
    """
    # Calculate the saturation vapor pressure
    AntioineA = 8.07131  # Valid for range 0 t0 100
    AntioineB = 1730.63  # Valid for range 0 t0 100
    AntioineC = 233.426  # Valid for range 0 t0 100
    PSat = 10 ** (AntioineA - (AntioineB / (AntioineC + t)))   # in mmHg  (760 mmHg = 101.325 kPa = 1.000 atm = normal pressure)
    Pw = ((rh / 100) * PSat)
    a = 17.27
    b = 237.7
    alpha = ((a * t) / (b + t)) + np.log(rh/100.0)
    dew_point = (b * alpha) / (a - alpha)
    # print(rh, PSat/760,Pw/760)
    return round(dew_point,2)





def T_sky_Allen_Clark(T_amb , T_dew , nOpaque):
    emmi_sky = 0.787 + 0.0028 * T_dew
    emmi_sky_cloudy = emmi_sky * (1 + 0.0224 * nOpaque - 0.0035*nOpaque**2 + 0.00028*nOpaque**3)
    t_sky = emmi_sky_cloudy**0.25 *T_amb
    return t_sky


def T_sky_Marting_Berdahl(time, P_amb_mbar,T_amb,T_dew, nOpaque, emmi_thin = 0.4 , emmi_thick= 1.0,h_thin= 2, h_thick=8, h_0=8.2):
    """time in 24 hour scale [0-24]
    P_amb in milibar
    Cloud coverage ratio as a number between 0 & 1 for both thin/thick type of clouds
    Height of clods h in Km, default values used 2 & 8 , h_0 =8.2
    """
    emmi_m = 0.711 + 0.56 * (T_dew/100) + 0.73*(T_dew/100)**2
    emmi_h = 0.013 * math.cos(2*math.pi*time/24)
    emmi_e = 0.00012 * (P_amb_mbar-1000)
    emmi_clear = emmi_m + emmi_h + emmi_e
    C = nOpaque * emmi_thin * math.exp(-h_thin/h_0)  + nOpaque* emmi_thick * math.exp(-h_thick/h_0)
    # C = n_cloud_coverage_thin * emmi_thin * math.exp(-h_thin / h_0) + n_cloud_coverage_thick * emmi_thick * math.exp( -h_thick / h_0)

    emmi_sky = emmi_clear + (1-emmi_clear) * C
    t_sky = emmi_sky**0.25 *T_amb
    return t_sky


def T_sky_dew_cloud_effect(T_amb, T_dew, nOpaque):
    """ Input tempretures must be in Kelvin"""
    emsvt_sky = (0.787 + 0.764 * math.log((T_dew+273)/273)) * (1 + 0.0224*nOpaque -0.0035*nOpaque**2 + 0.00028*nOpaque**3)
    ir_horizontal = const_Stefan_Boltzman * emsvt_sky * (T_amb+273)**4
    T_sky = emsvt_sky**0.25 * T_amb

    return T_sky

def T_sky_grag_1982(T_amb):
    T_sky = T_amb - 20
    return T_sky

def T_sky_Swinbank_1963(T_amb):
    T_sky = 0.0552*T_amb**(1.5)
    return T_sky


def T_sky_Fuentus_1987(T_amb):
    T_sky = 0.037536*T_amb**1.5 + 0.32*T_amb
    return T_sky


# def T_sky_Daguenet_1985(T_amb=T_amb, Pw=Pw * 0.000145038, n_cloud_coverage= nOpaque):
#     A = 10.1 * math.log(Pw)-12.3
#     B = 1.7 * T_amb + 107
#     C = -0.22 * math.log(Pw) + 1.25
#     L_0 = 3.6 * T_amb + 231
#     N = 8 - n_cloud_coverage*8
#
#     L = -L_0*(1+1.01*A) + B * C * ((8-N)/8)
#     print(L)
#     T_sky = (L/const_Stefan_Boltzman)**0.25
#     print(Pw, T_sky, A  ,B  ,C  ,L_0 ,N  ,L)
#     return T_sky


def T_sky_Daguenet_Sweden(T_amb, Pw , nOpaque):

    emmi_sky_cloudy = (0.43 + 0.082 * Pw**0.5) * (1-0.1 * nOpaque) + 0.1 * nOpaque
    t_sky = emmi_sky_cloudy**0.25 *T_amb
    return t_sky


def T_Sky_Building_Dymoala_BlBd(nOpaque, HorIR, T_amb, T_dew):
    """ Computation of black-body sky temperature degC Dry bulb temperature at ground level
        Dew point temperature   degC
        Opaque sky cover [0, 1]
        Black-body sky temperature
        Horizontal infrared irradiation W/m2
    """
    T_amb+=T_zero
    T_dew+=T_zero
    nOpa10 =  10 * nOpaque  #"Input nOpa is scaled to [0,1] instead of [0,10]"
    epsSky =  (0.787 + 0.764 * math.log(T_dew/T_zero))*(1 + 0.0224*nOpa10 - 0.0035*(nOpa10**2) + 0.00028*(nOpa10**3))
    TBlaSky =  T_amb*(epsSky**0.25)

    return round(TBlaSky-T_zero,2)



def T_Sky_faiman_rad(poa_global, temp_air, wind_speed=1.0, ir_down=0.9,
               u0=25.0, u1=6.84, sky_view=1.0, emissivity=0.88):


    abs_zero = -273.15
    sigma = scipy.constants.Stefan_Boltzmann

    if ir_down is None:
        qrad_sky = 0.0
    else:
        ir_up = sigma * ((temp_air - abs_zero)**4)
        qrad_sky = emissivity * sky_view * (ir_up - ir_down)

    heat_input = poa_global - qrad_sky
    total_loss_factor = u0 + u1 * wind_speed
    temp_difference = heat_input / total_loss_factor
    return temp_air + temp_difference




# skytemp= pvlib.temperature.faiman_rad(0, 50, wind_speed=5.0, ir_down=100, u0=25.0, u1=6.84, sky_view=0.90, emissivity=0.88)
# print(skytemp)
# print(T_Sky_faiman_rad(0, 30,wind_speed=3, ir_down=100,u0=25,u1=6.84, sky_view=0.9, emissivity=EMISSIVITY))

def T_Sky_Average(time_24,P_amb_mBar,nOpaque , HorIR , T_amb , T_dew, Pw, irr_dk):
    T_Sky_Mean = [
        round(T_sky_Allen_Clark(T_amb,T_dew, nOpaque),2),
        round(T_sky_Marting_Berdahl(time_24, P_amb_mBar,T_amb,T_dew, nOpaque,),2),
        round(T_sky_dew_cloud_effect(T_amb, T_dew, nOpaque),2),
        round(T_sky_grag_1982(T_amb),2),
        round(T_sky_Swinbank_1963(T_amb),2),
        round(T_sky_Fuentus_1987(T_amb),2),
        round(T_sky_Daguenet_Sweden(T_amb, Pw, nOpaque),2),
        round(T_Sky_Building_Dymoala_BlBd(nOpaque, HorIR, T_amb, T_dew),2),
        round(T_Sky_faiman_rad(irr_dk, T_amb,wind_speed=3, ir_down=HorIR,u0=25,u1=6.84, sky_view=0.9, emissivity=EMISSIVITY),2)
    ]


    print(
        " T_sky_Allen_Clark                      :    ", round(T_sky_Allen_Clark(T_amb, T_dew, nOpaque), 2),
        "\n T_sky_Marting_Berdahl                 :    ",
        round(T_sky_Marting_Berdahl(time_24, P_amb_mBar, T_amb, T_dew, nOpaque), 2),
        "\n T_sky_dew_cloud_effect                :    ", round(T_sky_dew_cloud_effect(T_amb, T_dew, nOpaque), 2),
        # "\n T_sky_grag_1982                       :    ", round(T_sky_grag_1982(T_amb),2),
        "\n T_sky_Swinbank_1963                   :    ", round(T_sky_Swinbank_1963(T_amb), 2),
        "\n T_sky_Fuentus_1987                    :    ", round(T_sky_Fuentus_1987(T_amb), 2),
        "\n T_sky_Daguenet_Sweden                 :    ", round(T_sky_Daguenet_Sweden(T_amb, Pw, nOpaque), 2),
        "\n T_Sky_Building_Dymoala_BlBd           :    ",
        round(T_Sky_Building_Dymoala_BlBd(nOpaque, HorIR, T_amb, T_dew), 2),
        "\n T_Sky_PVLib_Faiman           :    ",
        round(T_Sky_faiman_rad(HorIR, T_amb, wind_speed=3, ir_down=100, u0=25, u1=6.84, sky_view=0.9, emissivity=EMISSIVITY),2)
    )
    print("\n T_Sky_Average                         :    ", round(sum(T_Sky_Mean)/len(T_Sky_Mean), 2))
    return round(sum(T_Sky_Mean)/len(T_Sky_Mean),2)


print("______________________________________FROM FUNCTIONS_____________________________________")





