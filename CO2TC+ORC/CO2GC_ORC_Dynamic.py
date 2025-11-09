import numpy as np
import pandas as pd
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# --- 0) Floating gas‐cooler reference ---
def FloatingGCref(Tamb, dTgc):
    A0, A1, A2 = 4.81864e7, -4.21537e5, 9.44047e2
    alpha, Pref = 0.1, 67e5
    T100, dTsub, HRMODE = 39+273.15, 3, 0
    if HRMODE < 2:
        T = min(50, max(7, Tamb)) + 273.15
        Tref = ((A1*A1 - 4*(A0-Pref)*A2)**0.5 - A1)/(2*A2)
        tmp1, tmp2 = T100 + dTsub, T100 - Tref
        beta = (1e7 - A0 - (A1+A2*tmp1)*tmp1) / (tmp2 + np.sqrt(tmp2*tmp2 + alpha*alpha))
        b0 = A0 + (A1 + A2*dTsub)*dTsub - beta*Tref
        b1 = A1 + 2*A2*dTsub + beta
        b2 = A2
        delta = T - Tref
        root  = np.sqrt(delta*delta + alpha*alpha)
        P = b0 + (b1 + b2*T)*T + beta*root
        Sgc_ref = T + dTgc - 273.15
        return Sgc_ref, P*1e-5  # returns (gas cooler out T °C, P in bar)
    else:
        return None, None

# --- 1) Time series ---
time_index = pd.date_range('2025-05-01', periods=24, freq='H')
df = pd.DataFrame({
    'ambient_temp'   : 5 + 40*np.sin(np.linspace(0,2*np.pi,24)),
    'co2_inlet_temp' : 40 + 80*np.sin(np.linspace(0,2*np.pi,24)),
    'co2_outlet_temp': 30 + 5*np.sin(np.linspace(0,2*np.pi,24))+5
}, index=time_index)

# --- 2) System parameters ---
params = {
    'co2_mass_flow': 1.0,    # kg/s (for ORC heat input)
    'orc_fluid'    : 'R245fa',
    'turbine_eff'  : 0.8,
    'pump_eff'     : 0.7,
    'min_pinch'    : 4      # °C
}
Delta_appr_gc = 2          # approach temperature for gas cooler (°C)
Delta_ext_SH = 10          # external superheat for compressors (°C)

# --- 3) Saturation‐curve helper ---
def sat_curve(fluid):
    Tmin = CP.PropsSI('TMIN', fluid)
    Tcr  = CP.PropsSI('TCRIT',fluid)
    Ts   = np.linspace(Tmin, Tcr*0.99, 200)
    h_l, h_v, P_s = [], [], []
    for T in Ts:
        try:
            P_s.append(CP.PropsSI('P','T',T,'Q',0,fluid)/1e5)
            h_l.append(CP.PropsSI('H','T',T,'Q',0,fluid)/1e3)
            h_v.append(CP.PropsSI('H','T',T,'Q',1,fluid)/1e3)
        except ValueError:
            pass
    return np.array(h_l), np.array(h_v), np.array(P_s)

hliq_orc, hvap_orc, Psat_orc = sat_curve(params['orc_fluid'])
hliq_co2, hvap_co2, Psat_co2 = sat_curve('CO2')

# --- 4) Precompute all states & curves for each timestep ---
data = []
for ts, row in df.iterrows():
    Tamb = row['ambient_temp']
    Tin = row['co2_inlet_temp']
    Tout= row['co2_outlet_temp']

    # --- ORC Cycle P–h ---
    P_cond_orc = CP.PropsSI('P','T',Tamb+273.15,'Q',0, params['orc_fluid'])
    T_evap_orc= Tin - params['min_pinch']
    P_evap_orc= CP.PropsSI('P','T',T_evap_orc+273.15,'Q',0, params['orc_fluid'])

    # State 1: saturated liquid
    h1 = CP.PropsSI('H','P',P_cond_orc,'Q',0, params['orc_fluid'])/1e3
    s1 = CP.PropsSI('S','P',P_cond_orc,'Q',0, params['orc_fluid'])/1e3
    # State 2: after pump
    h2s= CP.PropsSI('H','P',P_evap_orc,'S',s1*1e3, params['orc_fluid'])/1e3
    h2 = h1 + (h2s - h1)/params['pump_eff']
    # State 3: saturated vapor
    h3 = CP.PropsSI('H','P',P_evap_orc,'Q',1, params['orc_fluid'])/1e3
    # State 4: after turbine
    s3 = CP.PropsSI('S','P',P_evap_orc,'H',h3*1e3, params['orc_fluid'])/1e3
    h4s= CP.PropsSI('H','P',P_cond_orc,'S',s3*1e3, params['orc_fluid'])/1e3
    h4 = h3 - params['turbine_eff']*(h3 - h4s)

    orc_h = [h1, h2, h3, h4, h1]
    orc_P = [P_cond_orc/1e5, P_evap_orc/1e5, P_evap_orc/1e5, P_cond_orc/1e5, P_cond_orc/1e5]

    # --- Heat‐Exchanger T–h profile between ORC and CO₂ ---
    # CO₂ side: from Tout → Tin (enthalpy space)
    P_gc_temp = Tamb + 273.15
    # we'll overwrite P_gc in full-cycle below, but for hex, use same FloatingGCref
    _, Pgc_bar_hex = FloatingGCref(Tamb, Delta_appr_gc)
    Pgc_hex = Pgc_bar_hex*1e5
    h_co2_curve = np.linspace(
        CP.PropsSI('H','T',Tout+273.15,'P',Pgc_hex,'CO2')/1e3,
        CP.PropsSI('H','T',Tin+273.15,'P',Pgc_hex,'CO2')/1e3,
        50
    )
    T_co2_curve = [CP.PropsSI('T','H',h*1e3,'P',Pgc_hex,'CO2')-273.15
                   for h in h_co2_curve]
    h_orc_curve = np.linspace(h2, h3, 50)
    T_orc_curve = [CP.PropsSI('T','H',h*1e3,'P',P_evap_orc,params['orc_fluid'])-273.15
                   for h in h_orc_curve]

    # --- CO₂ Full Booster Refrigeration Cycle P–h ---
    # 1) Receiver
    P_rec    = 3_625_000
    h_rec_li = CP.PropsSI('H','P',P_rec,'Q',0,'CO2')/1e3
    h_rec_vp = CP.PropsSI('H','P',P_rec,'Q',1,'CO2')/1e3

    # 2) MT Evaporator
    T_MT      = -8 + 273.15
    Q_MT_kW   = 100 if Tamb<=10 else (4*Tamb + 60)
    P_ev_MT   = CP.PropsSI('P','T',T_MT,'Q',0,'CO2')
    h_ev_out_MT = CP.PropsSI('H','P',P_ev_MT,'T',T_MT+10,'CO2')/1e3
    m_dot_MT    = Q_MT_kW / (h_ev_out_MT - h_rec_li)

    # 3) LT Evaporator
    T_LT      = -32 + 273.15
    Q_LT_kW   = 30
    P_ev_LT   = CP.PropsSI('P','T',T_LT,'Q',0,'CO2')
    h_ev_out_LT = CP.PropsSI('H','P',P_ev_LT,'T',T_LT+10,'CO2')/1e3
    m_dot_LT    = Q_LT_kW / (h_ev_out_LT - h_rec_li)

    # 4) LT Compressor
    T_suc_LT = T_LT + 10 + Delta_ext_SH
    h_suc_LT = CP.PropsSI('H','P',P_ev_LT,'T',T_suc_LT,'CO2')/1e3
    s_suc_LT = CP.PropsSI('S','P',P_ev_LT,'T',T_suc_LT,'CO2')/1e3
    eta_LT   = -0.007*(P_ev_MT/P_ev_LT)**2 + 0.0153*(P_ev_MT/P_ev_LT) + 0.6273
    h_dis_LT_is = CP.PropsSI('H','P',P_ev_MT,'S',s_suc_LT*1e3,'CO2')/1e3
    h_dis_LT    = h_suc_LT + (h_dis_LT_is - h_suc_LT)/eta_LT

    # 5) Gas‐Cooler via FloatingGCref
    Sgc_out, Pgc_bar = FloatingGCref(Tamb, Delta_appr_gc)
    P_gc = Pgc_bar * 1e5
    T_gc = Sgc_out + 273.15
    h_gc_out = CP.PropsSI('H','P',P_gc,'T',T_gc,'CO2')/1e3

    # 6) Bypass valve & mixing
    x_in_rec = CP.PropsSI('Q','P',P_rec,'H',h_gc_out*1e3,'CO2')
    m_dot_vbv = (m_dot_MT + m_dot_LT) * x_in_rec / (1 - x_in_rec)
    h_mix1 = (m_dot_MT*h_ev_out_MT + m_dot_LT*h_dis_LT) / (m_dot_MT + m_dot_LT)
    P_mix1 = P_ev_MT
    h_mix2 = ((m_dot_MT+m_dot_LT)*h_mix1 + m_dot_vbv*h_rec_vp) / (m_dot_MT+m_dot_LT+m_dot_vbv)
    P_mix2 = P_ev_MT

    # 7) MT Compressor
    T_suc_MT = CP.PropsSI('T','P',P_mix2,'H',h_mix2*1e3,'CO2') + Delta_ext_SH
    h_suc_MT = CP.PropsSI('H','P',P_mix2,'T',T_suc_MT,'CO2')/1e3
    s_suc_MT = CP.PropsSI('S','P',P_mix2,'T',T_suc_MT,'CO2')/1e3
    eta_MT   = -0.0144*(P_gc/P_ev_MT)**2 + 0.0719*(P_gc/P_ev_MT) + 0.5826
    h_dis_MT_is = CP.PropsSI('H','P',P_gc,'S',s_suc_MT*1e3,'CO2')/1e3
    h_dis_MT    = h_suc_MT + (h_dis_MT_is - h_suc_MT)/eta_MT

    co2_h = [
        h_rec_li, h_rec_vp,
        h_ev_out_MT, h_mix1, h_mix2,
        h_suc_LT, h_dis_LT,
        h_gc_out,
        h_suc_MT, h_dis_MT,
        h_rec_li
    ]
    co2_P = [
        P_rec/1e5, P_rec/1e5,
        P_ev_MT/1e5, P_mix1/1e5, P_mix2/1e5,
        P_ev_LT/1e5, P_ev_MT/1e5,
        P_gc/1e5,
        P_mix2/1e5, P_gc/1e5,
        P_rec/1e5
    ]

    data.append({
        'time': ts,
        'orc_h': orc_h,       'orc_P': orc_P,
        'h_co2_curve': h_co2_curve, 'T_co2_curve': T_co2_curve,
        'h_orc_curve': h_orc_curve, 'T_orc_curve': T_orc_curve,
        'co2_h': co2_h,       'co2_P': co2_P
    })

# --- 5) Animate ORC P–h diagram ---
fig1, ax1 = plt.subplots()
ax1.semilogy(hliq_orc, Psat_orc, 'b-', label='liq sat')
ax1.semilogy(hvap_orc, Psat_orc, 'r-', label='vap sat')
line1, = ax1.plot([],[], 'g-'); pts1, = ax1.plot([],[], 'ko')
ax1.set_xlabel('h (kJ/kg)'); ax1.set_ylabel('P (bar)'); ax1.legend()
def update_orc(i):
    d = data[i]
    line1.set_data(d['orc_h'], d['orc_P'])
    pts1.set_data(d['orc_h'][:-1], d['orc_P'][:-1])
    ax1.relim(); ax1.autoscale_view()
    ax1.set_title(f"ORC P–h @ {d['time']}")
    return line1, pts1
ani1 = FuncAnimation(fig1, update_orc, frames=len(data), interval=500, blit=True)

# --- 6) Animate HEX T–h profile ---
fig2, ax2 = plt.subplots()
l1, = ax2.plot([],[], 'r-', label='CO₂'); l2, = ax2.plot([],[], 'b-', label=params['orc_fluid'])
pt, = ax2.plot([],[], 'ko')
ax2.set_xlabel('h (kJ/kg)'); ax2.set_ylabel('T (°C)'); ax2.legend()
def update_hex(i):
    d = data[i]
    l1.set_data(d['h_co2_curve'], d['T_co2_curve'])
    l2.set_data(d['h_orc_curve'], d['T_orc_curve'])
    pt.set_data([d['h_orc_curve'][-1]], [d['T_orc_curve'][-1]])
    ax2.relim(); ax2.autoscale_view()
    ax2.set_title(f"HEX T–h @ {d['time']}")
    return l1, l2, pt
ani2 = FuncAnimation(fig2, update_hex, frames=len(data), interval=500, blit=True)

# --- 7) Animate CO₂ full booster cycle P–h (markers only) ---
fig3, ax3 = plt.subplots()
ax3.semilogy(hliq_co2, Psat_co2, 'b-', label='liq sat')
ax3.semilogy(hvap_co2, Psat_co2, 'r-', label='vap sat')
# Plot markers only—no line
pts3, = ax3.plot([], [], marker='o', linestyle='', color='k', label='cycle points')
ax3.set_xlabel('h (kJ/kg)')
ax3.set_ylabel('P (bar)')
ax3.legend()

def update_co2(i):
    d = data[i]
    # Only update the marker positions
    pts3.set_data(d['co2_h'], d['co2_P'])
    ax3.relim()
    ax3.autoscale_view()
    ax3.set_title(f"CO₂ Full Cycle P–h @ {d['time']}")
    return (pts3,)

ani3 = FuncAnimation(fig3, update_co2, frames=len(data), interval=500, blit=True)

plt.show()
