# GCBooster.py
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np


def COP(T_amb_C,TtoPlot):
    import CoolProp.CoolProp as CP
    import matplotlib.pyplot as plt
    import numpy as np
    # Inputs
    t_amb_C = T_amb_C

    # Inputs
    P_gc_2 = 131.4 * (10 ** 5)
    eta_booster = 0.65

    T_amb = t_amb_C + 273.15
    if T_amb <= 283.15:
        Q_dot_MT = 100
    else:
        Q_dot_MT = 4 * (T_amb - 273.15) + 60
    Q_dot_LT = 30

    T_MT = -8 + 273.15
    T_LT = -32 + 273.15

    P_rec = 3625000
    Delta_SH = 10
    Delta_ext_SH = 10
    Delta_appr_gc = 2

    # Receiver
    T_rec = CP.PropsSI('T', 'P', P_rec, 'Q', 0, 'CO2')
    h_out_rec_liq = CP.PropsSI('H', 'P', P_rec, 'Q', 0, 'CO2')
    h_out_rec_vap = CP.PropsSI('H', 'P', P_rec, 'Q', 1, 'CO2')

    # MT evaporators
    h_in_evap_MT = h_out_rec_liq
    P_evap_MT = CP.PropsSI('P', 'Q', 0, 'T', T_MT, 'CO2')
    h_out_evap_MT = CP.PropsSI('H', 'P', P_evap_MT, 'T', T_MT + Delta_SH, 'CO2')
    m_dot_MT = Q_dot_MT * 1000 / (h_out_evap_MT - h_in_evap_MT)

    # LT evaporators
    h_in_evap_LT = h_out_rec_liq
    P_evap_LT = CP.PropsSI('P', 'Q', 0, 'T', T_LT, 'CO2')
    h_out_evap_LT = CP.PropsSI('H', 'P', P_evap_LT, 'T', T_LT + Delta_SH, 'CO2')
    m_dot_LT = Q_dot_LT * 1000 / (h_out_evap_LT - h_in_evap_LT)

    # LT compressors
    T_suction_LT = T_LT + Delta_SH + Delta_ext_SH
    h_suction_LT = CP.PropsSI('H', 'P', P_evap_LT, 'T', T_suction_LT, 'CO2')
    s_suction_LT = CP.PropsSI('S', 'P', P_evap_LT, 'T', T_suction_LT, 'CO2')
    eta_LT = -0.007 * ((P_evap_MT / P_evap_LT) ** 2) + 0.0153 * (P_evap_MT / P_evap_LT) + 0.6273
    s_discharge_LT_is = s_suction_LT
    h_discharge_LT_is = CP.PropsSI('H', 'P', P_evap_MT, 'S', s_discharge_LT_is, 'CO2')
    h_discharge_LT = (h_discharge_LT_is - h_suction_LT) / eta_LT + h_suction_LT
    W_dot_LT = m_dot_LT * (h_discharge_LT - h_suction_LT) / 1000

    # Condenser/gas cooler #1
    T_gc_out = Delta_appr_gc + T_amb

    # Function to calculate P_gc
    def FloatingGCref(Tamb, dTgc):
        A0 = 4.81864e7
        A1 = -4.21537e5
        A2 = 9.44047e2
        alpha = 0.1
        Pref = 67e5
        T100 = 39.0 + 273.15
        dTsub = 3
        HRMODE = 0
        if HRMODE < 2:
            T = min(50, max(7, Tamb)) + 273.15
            Tref = ((A1 * A1 - 4.0 * (A0 - Pref) * A2) ** 0.5 - A1) / (2.0 * A2)
            tmp1 = T100 + dTsub
            tmp2 = T100 - Tref
            beta = (1e7 - A0 - (A1 + A2 * tmp1) * tmp1) / (tmp2 + (tmp2 * tmp2 + alpha * alpha) ** 0.5)
            b0 = A0 + (A1 + A2 * dTsub) * dTsub - beta * Tref
            b1 = A1 + 2 * A2 * dTsub + beta
            b2 = A2
            delta = T - Tref
            root = (delta * delta + alpha * alpha) ** 0.5
            P = b0 + (b1 + b2 * T) * T + beta * root
            Sgc_ref = T + dTgc - 273.15
            Pgc_ref = P * 1e-5
            return Sgc_ref, Pgc_ref
        else:
            return None, None

    _, P_gc_bar = FloatingGCref(t_amb_C, Delta_appr_gc)
    P_gc_1 = P_gc_bar * 1e5

    h_gc_out_1 = CP.PropsSI('H', 'P', P_gc_1, 'T', T_gc_out, 'CO2')
    s_gc_out_1 = CP.PropsSI('S', 'P', P_gc_1, 'T', T_gc_out, 'CO2')

    # Booster
    h_suction_booster = h_gc_out_1
    s_discharge_booster_is = s_gc_out_1
    h_discharge_booster_is = CP.PropsSI('H', 'P', P_gc_2, 'S', s_discharge_booster_is, 'CO2')
    h_discharge_booster = (h_discharge_booster_is - h_suction_booster) / eta_booster + h_suction_booster

    # Condenser/gas cooler #2
    h_gc_out_2 = CP.PropsSI('H', 'P', P_gc_2, 'T', T_gc_out, 'CO2')
    h_in_rec = h_gc_out_2
    x_in_rec = CP.PropsSI('Q', 'P', P_rec, 'H', h_in_rec, 'CO2')

    # Vapour by-pass valve
    h_out_vbv = h_out_rec_vap
    m_dot_vbv = (m_dot_MT + m_dot_LT) * x_in_rec / (1 - x_in_rec)

    # Mixing points
    h_mix_1 = (m_dot_MT * h_out_evap_MT + m_dot_LT * h_discharge_LT) / (m_dot_MT + m_dot_LT)
    T_mix_1 = CP.PropsSI('T', 'P', P_evap_MT, 'H', h_mix_1, 'CO2')
    h_mix_2 = ((m_dot_MT + m_dot_LT) * h_mix_1 + m_dot_vbv * h_out_vbv) / (m_dot_MT + m_dot_LT + m_dot_vbv)
    T_mix_2 = CP.PropsSI('T', 'P', P_evap_MT, 'H', h_mix_2, 'CO2')

    # MT compressors
    T_suction_MT = T_mix_2 + Delta_ext_SH
    h_suction_MT = CP.PropsSI('H', 'P', P_evap_MT, 'T', T_suction_MT, 'CO2')
    s_suction_MT = CP.PropsSI('S', 'P', P_evap_MT, 'T', T_suction_MT, 'CO2')
    eta_MT = -0.0144 * ((P_gc_1 / P_evap_MT) ** 2) + 0.0719 * (P_gc_1 / P_evap_MT) + 0.5826
    s_discharge_MT_is = s_suction_MT
    h_discharge_MT_is = CP.PropsSI('H', 'P', P_gc_1, 'S', s_discharge_MT_is, 'CO2')
    h_discharge_MT = (h_discharge_MT_is - h_suction_MT) / eta_MT + h_suction_MT
    W_dot_MT = (m_dot_MT + m_dot_LT + m_dot_vbv) * (h_discharge_MT - h_suction_MT) / 1000

    W_dot_booster = (m_dot_MT + m_dot_LT + m_dot_vbv) * (h_discharge_booster - h_suction_booster) / 1000

    Q_dot_GC_1 = (m_dot_MT + m_dot_LT + m_dot_vbv) * (h_discharge_MT - h_gc_out_1) / 0.97 / 1000
    W_dot_fan_1 = Q_dot_GC_1 * 0.03

    Q_dot_GC_2 = (m_dot_MT + m_dot_LT + m_dot_vbv) * (h_discharge_booster - h_gc_out_2) / 0.97 / 1000
    W_dot_fan_2 = Q_dot_GC_2 * 0.03

    COP = (Q_dot_LT + Q_dot_MT) / (W_dot_MT + W_dot_LT + W_dot_fan_1 + W_dot_booster + W_dot_fan_2)


    if T_amb_C==TtoPlot:
        # Function to calculate saturation properties
        def saturation_curve(fluid, t_min, t_max, num_points=100):
            t_range = np.linspace(t_min, t_max, num_points)
            p_sat = [CP.PropsSI('P', 'T', t, 'Q', 0, fluid) for t in t_range]
            h_liq = [CP.PropsSI('H', 'T', t, 'Q', 0, fluid) for t in t_range]
            h_vap = [CP.PropsSI('H', 'T', t, 'Q', 1, fluid) for t in t_range]
            return h_liq, h_vap, p_sat

        # Calculate saturation curve
        t_min = CP.PropsSI('Tmin', 'CO2')
        t_crit = CP.PropsSI('Tcrit', 'CO2')
        h_liq, h_vap, p_sat = saturation_curve('CO2', t_min, t_crit)

        # Create the plot
        plt.figure(figsize=(12, 8))
        plt.plot(h_liq, p_sat, 'b', h_vap, p_sat, 'r')
        plt.yscale('log')
        plt.xlabel('Enthalpy (J/kg)')
        plt.ylabel('Pressure (Pa)')
        plt.title('Pressure-Enthalpy Diagram for CO2 Refrigeration System')

        # Function to plot a point
        def plot_point(h, p, label):
            plt.plot(h, p, 'ko', markersize=6)
            plt.text(h, p, label, fontsize=8, ha='left', va='bottom')


        # Print Booster cooling contribiution
        plot_point(h_suction_booster, P_gc_1, 'Booster_suction')
        plot_point(h_discharge_booster, P_gc_2, 'Booster_discharge')
        plot_point(h_gc_out_2, P_gc_2, 'GC2_out')


        # Plot the points
        plot_point(h_out_rec_liq, P_rec, 'Rec_liq')
        plot_point(h_in_evap_MT, P_evap_MT, 'MT_evap_in')
        plot_point(h_out_rec_vap, P_rec, 'Rec_vap')
        plot_point(h_out_evap_MT, P_evap_MT, 'MT_evap_out')
        plot_point(h_in_evap_LT, P_evap_LT, 'LT_evap_in')
        plot_point(h_out_evap_LT, P_evap_LT, 'LT_evap_out')
        plot_point(h_suction_LT, P_evap_LT, 'LT_suction')
        plot_point(h_discharge_LT, P_evap_MT, 'LT_discharge')
        plot_point(h_gc_out_1, P_gc_1, 'GC1_out')
        plot_point(h_suction_booster, P_gc_1, 'Booster_suction')
        plot_point(h_discharge_booster, P_gc_2, 'Booster_discharge')
        plot_point(h_gc_out_2, P_gc_2, 'GC2_out')
        plot_point(h_in_rec, P_rec, 'Rec_in')
        plot_point(h_out_vbv, P_evap_MT, 'VBV_out')
        plot_point(h_mix_1, P_evap_MT, 'Mix1')
        plot_point(h_mix_2, P_evap_MT, 'Mix2')
        plot_point(h_suction_MT, P_evap_MT, 'MT_suction')
        plot_point(h_discharge_MT, P_gc_1, 'MT_discharge')

        plt.grid(True)

    Booster_effect = (h_suction_booster-h_gc_out_2)/(h_discharge_booster-h_suction_booster)
    return [COP, Booster_effect]

