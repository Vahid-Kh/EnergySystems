
import time
import pandas as pd                 # Dataframe library
import numpy as np                  # Scientific computing with nD object support
from datetime import datetime
from TDN import TDN, PSI, mlog, mt, sort
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from functions import mov_ave, plot_2,plot_3



"""----------------- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%------------- """
"""----------------- %%%%%%%%%%%--- CALCULATIONS ---%%%%%%%%%------------- """
"""----------------- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%------------- """

""" Refrigerant """
r = 'CO2'

""" Time """
dt = 1/60   # Time steps (seconds)

# t1, t2 = 1, 4
# t1, t2 = 2, 9
# t1, t2 = 3, 13
t1, t2 = 4, 17
# t1, t2 = 5, 21
# t1, t2 = 9, 42

tt = t1+t2    # Reliability time (seconds)

"""Cooling load on MT cabinets"""
load_pwm = 60e3
load_pc = 80e3

"""_________ Reference Condition (Ambient used in 2nd law) __________ : """
tamb = 293.15
pamb = 101325
rgas = 8.31446261815324  # Gas constant J⋅K−1⋅mol−1

"""Pressure levels::::"""
p_lt = 16.83e5
pmt = 26.5e5
pit = 40e5
pgc = 90e5

""" Opt GS Temperature: """
tgc = 273.15 + 35

""" Volumes  [m^3] """
vd = 500e-3
v60 = 30.9 / 3600   # m3/s
v50 = 25.6 / 3600   # m3/s
v30 = v60-(v60-v50)*3

""" T evap and P gc for Bitzer polynom"""
to1 = PSI('T', 'Q', 1, 'P', pit, r)      # T evap at Mt level
p_HP1 = 90                                  # Pgc in bar

to2 = PSI('T', 'Q', 1, 'P', pmt, r)                                # T evap at Mt level
p_HP2 = 90                                  # Pgc in bar

"""___Temperature [K]___ : """
dtssh = 10  # DT of superheat in (k)
"""   P drops               """
# dp_suc = 2e4  # Pressure drop at the suction side of compressor (pa)
# dp_dis = 2e4  # Pressure drop at the discharge side of compressor (pa)
# dp_evap = 7e4  # Pressure drop at the Evaporator in (pa)
# dp_cond = 2e4  # Pressure drop at the condenser in (pa)
# dp_liq = 1e4  # Pressure drop at the liquid line (pa)

"""------------------------------------"""

"""___Number of Data Points___ : """
num_data_1 = 6  # Number of data points for Refrigerant
num_data_2 = 8  # Number of data points for Refrigerant

""" Efficiencies :::: """
eta_vol_mt = 0.9

"""________________________________________________________________________________"""
"""      MT COMPRESSOR GROUP

Input Values:										

Compressor model			4CTC-30K							
Mode			Refrigeration and Air conditioning							
Refrigerant			R744							
Reference temperature			Dew point temp.							
Gas cooling outlet			25,0 °C							
Suct. gas superheat			10,00 K							
Operating mode			Transcritical							
Power supply			400V-3-50Hz							
Capacity control			100%							
Useful superheat			100%							

------------------------------------------------------------										
Polynomial:										
y = c1 + c2*to + c3*p_HP + c4*to^2 + c5*to*p_HP + c6*p_HP^2 + c7*to^3 + c8*p_HP*to^2 + c9*to*p_HP^2 + c10*p_HP^3										

Coefficients:										
	c1	c2	c3	c4	c5	c6	c7	c8	c9	c10
Q [W]	0,00000000000000000000	0,00000000000000000000	0,00000000000000000000	0,00000000000000000000	0,00000000000000000000	0,00000000000000000000	0,00000000000000000000	0,00000000000000000000	0,00000000000000000000	0,00000000000000000000
P [W]	-1,089,201,779,535,750,000,000,000	-72,294,659,113,107,600,000,000	60,333,996,681,101,100,000,000	-612,308,664,344,127,000,000	1,018,534,579,823,220,000,000	-0,94555636305316300000	0,04733890604629470000	0,01295718415594570000	-0,00988455081063097000	-0,00354475538872107000
m [kg/h]	251,756,974,213,016,000,000,000	7,490,666,951,329,370,000,000	-890,364,355,501,862,000,000	0,91933920026369500000	-0,13593927823451300000	0,02042525932763880000	0,00623565082052914000	-0,00055656997425671500	-0,00002308991356673690	-0,00004912171980273300
I [A]	-930,259,250,099,061,000,000	-102,977,881,948,575,000,000	0,85644027637668600000	-0,00990311670803154000	0,01336367508179710000	-0,00091712826978452300	0,00006728186844906470	0,00002477743497250250	0,00000064926809632384	-0,00000614632125452918

------------------------------------------------------------										
Validity range of Polynomials										
Evaporating SST:			-20°C	...	0,1°C					
High pressure:			73,8bar(a)	...	106bar(a)					
Attention: Consider also application range of compressor!										
"""

"""  CONSTANTS OF POLYNOMIALS :  [W] [kh/h]"""
m_4ctc = [2517.56974213016000000000, 74.90666951329370000000, -8.90364355501862000000, 0.91933920026369500000, -0.13593927823451300000, 0.02042525932763880000, 0.00623565082052914000, -0.00055656997425671500, -0.00002308991356673690, -0.00004912171980273300]
p_4ctc = [-10892.01779535750000000000, -722.94659113107600000000, 603.33996681101100000000, -6.12308664344127000000, 10.18534579823220000000, -0.94555636305316300000, 0.04733890604629470000, 0.01295718415594570000, -0.00988455081063097000, -0.00354475538872107000]
""" Used for IT in Parallel compression """
m_4mtc = [657.01624068665400000000, 26.10576311752670000000, -4.23401414159409000000, 0.59904521603669500000, -0.09928305379032660000, 0.02405203148063410000, 0.00843630458267610000, -0.00131558951278288000, 0.00023077486161974100, -0.00007192875085841990]
p_4mtc = [-4382.67634199997000000000, -132.05744640000000000000, 185.89802219999900000000, -1.01561159700000000000, -0.26110268909999600000, -0.54099135089999100000, 0.00344413347300000000, -0.02375192358000000000, 0.01719785826000000000, -0.00000000000000002818]
"""______________________________VCC REFRIGERANT TDN STATES__________________________________________________"""

"""________________________________________________________________________________"""

""" _____________   PWM MODE 1 ::::::::::    _____________   """

PWM1 = [TDN for ii in range(num_data_1)]
PWM1_C = [TDN for ii in range(num_data_1)]
"""
SAmple for TDN states :::::  

PWM1[] = TDN(, , , , , r)  # 
PWM1_C[] = PWM1[].__copy__()

"""
PWM1[0] = TDN(pit, 0, PSI('T', 'Q', 1, 'P', pit, r) + dtssh, 0, 0, r)  # TDN state in suction side of compressor
PWM1[1] = TDN(pit, 0, PSI('T', 'Q', 1, 'P', pit, r) + dtssh, 0, 0, r)  # TDN state in suction side of compressor
PWM1[2] = TDN(pgc, 0, 0, PWM1[1].s, 0, r)  # TDN state in suction side of compressor
PWM1[3] = TDN(pgc, 0, tgc, 0, 0, r)
PWM1[4] = TDN(pit, PWM1[3].h, 0, 0, 0, r)
q1rec =   PSI('Q', 'H', PWM1[4].h, 'P', PWM1[4].p, r)
PWM1[5] = TDN(pit, PSI('H', 'Q', 0, 'P', pit, r), 0, 0, 0, r)

"""________________________________________________________________________________"""
"""  _____________  PWM MODE 2 ::::::::::  _______________    """

PWM2 = [TDN for ii in range(num_data_2)]
PWM2_C = [TDN for ii in range(num_data_2)]

PWM2[0] = TDN(pmt, 0, PSI('T', 'Q', 1, 'P', pmt, r) + dtssh, 0, 0, r)  # Pressure drop in suction side entrance of compressor
PWM2[1] = TDN(pmt, 0, PSI('T', 'Q', 1, 'P', pmt, r) + dtssh, 0, 0, r)  # Pressure drop in suction side entrance of compressor
PWM2[2] = TDN(pgc, 0, 0, PWM2[1].s, 0, r)  # TDN state in suction side of compressor
PWM2[3] = TDN(pgc, 0, tgc, 0, 0, r)
PWM2[4] = TDN(pit, PWM2[3].h, 0, 0, 0, r)
q2rec =   PSI('Q', 'H', PWM2[4].h, 'P', PWM2[4].p, r)
PWM2[5] = TDN(pit, PSI('H', 'Q', 0, 'P', pit, r), 0, 0, 0, r)
PWM2[6] = TDN(pit, PSI('H', 'Q', 1, 'P', pit, r), 0, 0, 0, r)
PWM2[7] = TDN(pmt, PWM2[5].h, 0, 0, 0, r)


"""________________________________________________________________________________"""


mmte = load_pwm / (PWM2[1].h - PWM2[7].h)
mmte_pc = load_pc/(PWM2[1].h-PWM2[7].h)
mfg = mmte * q2rec
mfg_pc = mmte_pc * q2rec

print('Accumulated mass [Kg/s],  Evap and Rec : ' ,mmte, mfg)
print('Total Accumulated mass,   Evap and Rec : ' ,mmte*t1, mfg*t2)

print('Total Accumulated Volume, Evap and Rec : ' ,mmte*t1/PWM2[1].d*1e3, mfg*t2/PWM1[1].d*1e3)

"""________________________________________________________________________________"""

""" ASSUMPTION :::   COMP RUN AT 60 HZ WHEN SUCTION IS FROM 10 Liter DEAD VOLUME (WORST CASE SCENARIO)"""
"""  PWM  """
pomt_l = []
mmt_l = []

"""  Parallel Compression  """
pomt_l_pc = []
mmt_l_pc  = []
poit_l_pc = []
mit_l_pc  = []

""" Booster System  """
pomt_l_b  = []
mmt_l_b   = []

pnew = pit

for i in range(int(tt/dt)):
    pomt_l_pc.append(pmt)
    poit_l_pc.append(pit)
    pomt_l_b.append(pmt)

    if i <= t1/dt:
        pomt_l.append(pit)
        # print('Mode 1')

    elif int(t1/dt) < i < int((t1+t2)/dt) and pnew>pmt:
        dnew = (TDN(pomt_l[i-1], 0, PSI('T', 'Q', 1, 'P', pomt_l[i-1], r) + dtssh, 0, 0, r).d * vd -
                TDN(pomt_l[i-1], 0, PSI('T', 'Q', 1, 'P', pomt_l[i-1], r) + dtssh, 0, 0, r).d*v60/60)\
               /vd

        pnew = TDN(0,0,PSI('T', 'Q', 1, 'D', dnew, r) + dtssh,0,dnew,r).p
        pomt_l.append(pnew)
        # print("Mode 1 to mode 2 transition ")

    elif int(t1 / dt) < i < int((t1 + t2) / dt) and pnew < pmt:
        pomt_l.append(pmt)
        # print('Mode 2')

    else:
        print("Pressure condition not defined")


"""Bitzer polynomial"""
vmt_l_b  = []
vmt_l_pc = []
vit_l_pc = []
vmt_l    = []

for i in range(tt):
    mmt_l_b.append (mmte+mfg)
    mmt_l_pc.append(mmte_pc)
    mit_l_pc.append(mfg_pc)


    vmt_l_b.append((mmte + mfg)/TDN(pmt, 0, PSI('T', 'Q', 1, 'P', pmt, r) + dtssh, 0, 0, r).d)
    vmt_l_pc.append((mmte_pc)/TDN(pmt, 0, PSI('T', 'Q', 1, 'P', pmt, r) + dtssh, 0, 0, r).d)
    vit_l_pc.append((mfg_pc)/TDN(pit, 0, PSI('T', 'Q', 1, 'P', pit, r) + dtssh, 0, 0, r).d)



    if i <= t1:
        mmt_l.append(mfg * tt / t1)
        vmt_l.append((mfg * tt / t1)/TDN(pit, 0, PSI('T', 'Q', 1, 'P', pit, r) + dtssh, 0, 0, r).d)

    elif t1 < i < (t1 + t2) :
        mmt_l.append(mmte * tt / t2)
        vmt_l.append((mmte * tt / t2)/TDN(pmt, 0, PSI('T', 'Q', 1, 'P', pmt, r) + dtssh, 0, 0, r).d)

    else:
        print("Mass balance condition not defined")



    """________________________________________________________________________________"""

""" PWM SYSTEM   :::  """
pwrmt = []
time = []
p_HP = pgc/1e5
for i in range(tt):
    to = PSI('T', 'Q', 1, 'P', pomt_l[int(i / dt)], r) - 273.15
    time.append(i)
    c = m_4ctc
    mdot_4ctc = (c[0] + c[1] * to + c[2] * p_HP + c[3] * to ** 2 + c[4] * to * p_HP + c[5] * p_HP ** 2 + c[
        6] * to ** 3 + c[7] * p_HP * to ** 2 + c[8] * to * p_HP ** 2 + c[9] * p_HP ** 3) / 3600

    rc = (mmt_l[i] / mdot_4ctc * 60 - 30) / 30 * 100

    # print(pomt_l[int(i/dt)])
    # print(mdot_4ctc)
    if i == 0 :
        print('Running capacity  PWM IT mode' ,rc)
    if  i == t1+1:
        print('Running capacity  PWM MT mode' ,rc)

    c = p_4ctc
    pwr_4ctc = (c[0] + c[1] * to + c[2] * p_HP + c[3] * to ** 2 + c[4] * to * p_HP + c[5] * p_HP ** 2 + c[6] * to ** 3 +
                c[7] * p_HP * to ** 2 + c[8] * to * p_HP ** 2 + c[9] * p_HP ** 3)

    pmt_bitzer = (rc / 100 * 30 + 30) / 60 * pwr_4ctc
    pwrmt.append(pmt_bitzer)

""" PARALLEL COMPRESSION SYSTEM  ::: """
pwr_mt_pc = []
pwr_it_pc = []

for i in range(tt):
    """ MT comp Parallel compression """
    to = PSI('T', 'Q', 1, 'P', pomt_l_pc[int(i/dt)], r) - 273.15

    c = m_4ctc
    mdot_4ctc = (c[0] + c[1] * to + c[2] * p_HP + c[3] * to ** 2 + c[4] * to * p_HP + c[5] * p_HP ** 2 + c[
        6] * to ** 3 + c[7] * p_HP * to ** 2 + c[8] * to * p_HP ** 2 + c[9] * p_HP ** 3) / 3600

    rc = (mmt_l_pc[i] / mdot_4ctc * 60-30)/30*100
    if i == 0:
        print('Running capacity  PC MT comp ' ,rc)

    c = p_4ctc
    pwr_4ctc = (c[0] + c[1] * to + c[2] * p_HP + c[3] * to ** 2 + c[4] * to * p_HP + c[5] * p_HP ** 2 + c[6] * to ** 3 +
                c[7] * p_HP * to ** 2 + c[8] * to * p_HP ** 2 + c[9] * p_HP ** 3)

    pmt_bitzer = (rc / 100 * 30 + 30) / 60 * pwr_4ctc
    pwr_mt_pc.append(pmt_bitzer)

    """ IT comp Parallel compression """
    to = PSI('T', 'Q', 1, 'P', poit_l_pc[int(i / dt)], r) - 273.15

    c = m_4mtc
    mdot_4ctc = (c[0] + c[1] * to + c[2] * p_HP + c[3] * to ** 2 + c[4] * to * p_HP + c[5] * p_HP ** 2 + c[
        6] * to ** 3 + c[7] * p_HP * to ** 2 + c[8] * to * p_HP ** 2 + c[9] * p_HP ** 3) / 3600

    rc = (mit_l_pc[i] / mdot_4ctc * 60 - 30) / 30 * 100

    # print(pomt_l[int(i/dt)])
    # print(mdot_4ctc)
    if i == 0:
        print('Running capacity  PC IT comp ' ,rc)
    c = p_4mtc
    pwr_4ctc = (c[0] + c[1] * to + c[2] * p_HP + c[3] * to ** 2 + c[4] * to * p_HP + c[5] * p_HP ** 2 + c[6] * to ** 3 +
                c[7] * p_HP * to ** 2 + c[8] * to * p_HP ** 2 + c[9] * p_HP ** 3)

    pmt_bitzer = (rc / 100 * 30 + 30) / 60 * pwr_4ctc
    pwr_it_pc.append(pmt_bitzer)


""" BOOSTER SYSTEM  ::: """
pwr_b = []
for i in range(tt):
    to = PSI('T', 'Q', 1, 'P', pomt_l_b[int(i / dt)], r) - 273.15

    c = m_4ctc
    mdot_4ctc = (c[0] + c[1] * to + c[2] * p_HP + c[3] * to ** 2 + c[4] * to * p_HP + c[5] * p_HP ** 2 + c[
        6] * to ** 3 + c[7] * p_HP * to ** 2 + c[8] * to * p_HP ** 2 + c[9] * p_HP ** 3) / 3600

    rc = (mmt_l_b[i] / mdot_4ctc * 60 - 30) / 30 * 100

    # print(pomt_l[int(i/dt)])
    # print(mdot_4ctc)
    if i == 0:
        print('Running capacity  Booster    ' ,rc)
    c = p_4ctc
    pwr_4ctc = (c[0] + c[1] * to + c[2] * p_HP + c[3] * to ** 2 + c[4] * to * p_HP + c[5] * p_HP ** 2 + c[6] * to ** 3 +
                c[7] * p_HP * to ** 2 + c[8] * to * p_HP ** 2 + c[9] * p_HP ** 3)

    pmt_bitzer = (rc / 100 * 30 + 30) / 60 * pwr_4ctc
    pwr_b.append(pmt_bitzer)


cop_pwm = list(load_pwm/np.array(pwrmt))
cop_pwm_b = list(load_pwm/np.array(pwr_b))
cop_pwm_pc = list(load_pc/(np.array(pwr_it_pc)+np.array(pwr_mt_pc)))

"""_________________________________________________________________________________________________________________
________________________________________________PRINT_______________________________________________________________
_________________________________________________________________________________________________________________"""

print("INSTANTANEOUS COP PWM in a period :                  ", np.mean(cop_pwm))
print("INSTANTANEOUS COP Booster in a period :              ", np.mean(cop_pwm_b))
print("INSTANTANEOUS COP Parallel Compression in a period : ", np.mean(cop_pwm_pc))

print(['Pressure(pa)', 'Enthalpy(kJ/kg)', 'Temperature(K)', 'Entropy(kJ/kg/K)', 'Volume(M3/kg)', 'Refrigerant'])
print('-------------------Single stage VCC------------')
i = 0
for til in PWM1:
    print(i, TDN.prop(til))
    i += 1
i = 0

print(['Pressure(pa)', 'Enthalpy(kJ/kg)', 'Temperature(K)', 'Entropy(kJ/kg/K)', 'Volume(M3/kg)', 'Refrigerant'])
print('-------------------Single stage VCC------------')
i = 0
for til in PWM2:
    print(i, TDN.prop(til))
    i += 1
i = 0

"""_________________________________________________________________________________________________________________
____________________________________________________PLOTS___________________________________________________________
_________________________________________________________________________________________________________________"""

"""  INSTANTANEOUS COP comparison  """
plt.figure('COP comparison ')
plt.plot(cop_pwm_b)
plt.plot(cop_pwm_pc)
plt.plot(cop_pwm)


plt.xlabel('Time [Sec]')
plt.ylabel("COP [-]")
plt.grid(c='grey', linestyle='-', linewidth=0.3, alpha=0.4)
plt.legend(['Booster system', 'Parallel compression', 'PWM system'])

"""  Compessor suction pressure comparison  """
plt.figure('Compessor suction pressure comparison ')

plt.plot(pomt_l_b)
plt.plot(pomt_l_pc)
plt.plot(poit_l_pc)
plt.plot(pomt_l)

plt.xlabel('Time [Sec]')
plt.ylabel("Compessor suction pressure [Pa]")
plt.grid(c='grey', linestyle='-', linewidth=0.3, alpha=0.4)
plt.legend(['Booster system', 'Parallel compression MT', 'Parallel compression IT', 'PWM system'])


"""  Compessor Mass flow rate comparison [Kg/sec]   """
plt.figure('Compessor Mass flow rate comparison [Kg/sec]')

plt.plot(mmt_l_b)
plt.plot(mmt_l_pc)
plt.plot(mit_l_pc)
plt.plot(mmt_l)

plt.xlabel('Time [Sec]')
plt.ylabel("Compessor Mass flow rate [Kg/sec]")
plt.grid(c='grey', linestyle='-', linewidth=0.3, alpha=0.4)
plt.legend(['Booster system', 'Parallel compression MT', 'Parallel compression IT', 'PWM system'])

"""  Compessor Volumetric flow rate [m3/sec]   """
plt.figure('Compessor Volumetric flow rate [m3/sec]')

plt.plot(vmt_l_b)
plt.plot(vmt_l_pc)
plt.plot(vit_l_pc)
plt.plot(vmt_l)

plt.xlabel('Time [Sec]')
plt.ylabel("Compessor Volumetric flow rate [m3/sec]")
plt.grid(c='grey', linestyle='-', linewidth=0.3, alpha=0.4)
plt.legend(['Booster system', 'Parallel compression MT', 'Parallel compression IT', 'PWM system'])

"""___________________________Thermodynamic states plot____________________________________________________ : """
plt.figure('Thermodynamic states plot ')
h_plot = []
h_plot_ro = []  # First plots the lines and then 'ro' spesifies red dots as the TDN states
p_plot = []
p_plot_ro = []

for i in range(len(PWM1)):
    h_plot.append(PWM1[i].h)
    h_plot_ro.append(PWM1[i].h)
    p_plot.append(PWM1[i].p)
    p_plot_ro.append(PWM1[i].p)
h_p_ro_1 = plt.semilogy(h_plot_ro, p_plot_ro, 'g*', markersize=15)

for i in range(len(PWM2)):
    h_plot.append(PWM2[i].h)
    h_plot_ro.append(PWM2[i].h)
    p_plot.append(PWM2[i].p)
    p_plot_ro.append(PWM2[i].p)


# h_p = plt.semilogy(h_plot, p_plot)
h_p_ro_2 = plt.semilogy(h_plot_ro, p_plot_ro, 'ro')
plt.title('Vapor Compressionn cycle TDN states of refrigerant' + r)
plt.xlabel('Entalpy [j/kg]')
plt.ylabel('Log Pressure, Log(P) [Pa]')

"""--------------------------------------------------------------------------------------------------------------"""

"""    Liquid & Vapor saturation plot"""
sat_liq_h = []
sat_vap_h = []
sat__p = []

"""______________Two-Phase plot from T_min to P_Critical for Refrigerant__________________"""
for p in range(int(round(PSI('P', 'Q', 0, 'T', PSI(r, 'Tmin'), r), -4) + 1e4), int(round(PSI(r, 'pcrit'), -4) - 1e4),
               1000):
    sat_liq_h.append(PSI('H', 'Q', 0, 'P', p, r))  # At vapor quality of 0 the sat_liquid line is obtained
    sat_vap_h.append(PSI('H', 'Q', 1, 'P', p, r))  # At vapor quality of 1 the sat_vapor line is obtained
    sat__p.append(p)

plot_liq_sat = plt.semilogy(sat_liq_h, sat__p)
plot_vap_sat = plt.semilogy(sat_vap_h, sat__p)
plt.grid()

plt.show()

