

"""
Electric Pressure Switch Reliability
"""

import numpy as np


def rlblt_or (r1, r2, r3=1, r4=1, r5=1):
    " Can take min 2 & max 4 arguments as reliability "
    if max(r1, r2, r3, r4, r5) > 1 or min(r1, r2, r3, r4, r5) < 0:
        print( " !!! WARNING !!! : Reliability values must be between 0 & 1 ")
    return r1 * r2 * r3 * r4 * r5


def rlblt_and (r1, r2, r3=0, r4=0, r5=0):
    " Can take min 2 & max 4 arguments as reliability "
    if max(r1, r2, r3, r4, r5) > 1 or min(r1, r2, r3, r4, r5) < 0:
        print( " !!! WARNING !!! : Reliability values must be between 0 & 1 ")
    return 1 - (1 - r1) * (1 - r2) * (1 - r3) * (1 - r4) * (1 - r5)



r_all = ((1e6 - 1) / 1e6)

"""  Pressure Safety Switch "A"  """

""" HW --> Customized IC | GreenPAK 
    MTTF = 1.88e10  ; Failure rate = 0.05 FIT
    Source: SLG46855V_Reliability Report.pdf
    Source: https://miro.com/app/board/o9J_lJ__Sfg=/?moveToWidget=3458764520499430525&cot=14
"""
mttf_hw = 1.88e10
r_hw_a = (mttf_hw-1)/mttf_hw
r_hw_b = (mttf_hw-1)/mttf_hw

"""--------------------------------------------------------------------------------------------------------"""

""" RY --> Relay
    MTTFd = (B10d/0.1) * number of operations component take in a year
    @ 3 A 230 V AC: B10D = 1400000 & @ At 3 A 24 V DC: B10D = 400000
    Assumption of 20 occurances of needing relay to act per year 
    Source: https://www.gt-engineering.it/en/en-iso-standards/en-ISO-13849-1-performance-level-estimation/ciAO 
    Source: https://www.elesta-gmbh.com/en/helpers/download/?file=reliability_b10d_sisf_3.pdf """
mttf_ry_a = (400000/0.1) * 20
r_ry_a = (mttf_ry_a-1)/mttf_ry_a
# r_ps_a = r_rb_a * r_as_a # Pressure Sensor  - A
"""--------------------------------------------------------------------------------------------------------"""

# """ PS --> Pressure sensor
#     Source: https://www.sintef.no/projectweb/pds-main-page/pds-handbooks/pds-data-handbook/
#
#     """
# r_ps_a = (1e6 - 1.98) / 1e6  # Pressure Sensor  - A
# r_ps_b = (1e6 - 1.98) / 1e6  # Pressure Sensor  - B


""" PS --> Pressure sensor 
    Source: Sensing solutions  | John Hansen  | Product Manager |  Application & Product Support 
    The MTTFd value for AKS 2050 is > 100 years at 125 C
    What is assumed here is 85 degC as ambient temperature which gives MTTFd of 1000 years
    """


r_ps_a = ((200*365*24 - 1) )/ (200*365*24)  # Pressure Sensor  - A
r_ps_b = ((200*365*24 - 1)) / (200*365*24)  # Pressure Sensor  - B
# r_ps_b = r_rb_b * r_as_b # Pressure Sensor  - B
print(r_ps_a)
r_rb_a = np.sqrt(r_ps_a)
r_as_a = np.sqrt(r_ps_a)
r_rb_b = np.sqrt(r_ps_b)
r_as_b = np.sqrt(r_ps_b)
"""--------------------------------------------------------------------------------------------------------"""

""" !!! AVERAGE VALUE TAKEN FOR START !!! """
r_pt_a = 1  # Potentiometer setting - A
r_pt_b = 1  # Potentiometer setting - B
r_ajs = r_all   # Adjustment screw - One used in all system for giving set point
r_ry_b = r_all  #

"""_____________________________________________________________________________________"""
"""  SCENARIO 1 : ADJUSTABLE PRESSURE SAFETY - Full redundancy
https://miro.com/app/board/o9J_lJ__Sfg=/?moveToWidget=3458764518301431404&cot=14
 """


r_ss_a = rlblt_or(r_rb_a
                  , r_as_a
                  , r_pt_a
                  , r_hw_a
                  , r_ry_a
                  )
r_ss_b = rlblt_or(r_rb_b
                  , r_as_b
                  , r_pt_b
                  , r_hw_b
                  , r_ry_b
                  )

r_pss = rlblt_or(rlblt_and(r_ss_a, r_ss_b), r_ajs)
print("  -----------------------------------------------------  ")
print("SCENARIO 1 : ADJUSTABLE PRESSURE SAFETY - Full REDUNDANCY")
print("Reliability value: ", r_pss)
print("Failure rate: ", round((1-r_pss)*1e7,3)," Failiure in 10 Million")



"""_____________________________________________________________________________________"""

"""  SCENARIO 2 : NON -ADJUSTABLE PRESSURE SAFETY - FULL REDUNDANCY
https://miro.com/app/board/o9J_lJ__Sfg=/?moveToWidget=3458764519132390194&cot=14
"""


r_ss_a = rlblt_or(  r_rb_a
                  , r_as_a
                  , r_hw_a
                  , r_ry_a
                  )
r_ss_b = rlblt_or(  r_rb_b
                  , r_as_b
                  , r_hw_b
                  , r_ry_b
                  )

r_pss = rlblt_and(r_ss_a, r_ss_b)

print("  -----------------------------------------------------  ")
print("SCENARIO 2 : NON-ADJUSTABLE PRESSURE SAFETY - FULL REDUNDANCY")
print("Reliability value: ",r_pss)
print("Failure rate: ", round((1-r_pss)*1e7,4)," Failure in 10 Million")

"""_____________________________________________________________________________________"""

"""  SCENARIO 3 : NON-ADJUSTABLE PRESSURE SAFETY - PARTIAL REDUNDANCY
https://miro.com/app/board/o9J_lJ__Sfg=/?moveToWidget=3458764519135725598&cot=14
"""

r_ps_a = rlblt_or(  r_rb_a
                  , r_as_a
                  )
r_ps_b = rlblt_or(  r_rb_b
                  , r_as_b
                  )
r_pss = rlblt_or(rlblt_and(r_ps_a,r_ps_b), r_hw_a, r_ry_a)

print("  -----------------------------------------------------  ")
print("SCENARIO 3 : NON-ADJUSTABLE PRESSURE SAFETY - PARTIAL REDUNDANCY")
print("Reliability value: ",r_pss)
print("Failure rate: ", round((1-r_pss)*1e7,3)," Failure in 10 Million")

"""_____________________________________________________________________________________"""


"""_____________________________________________________________________________________"""

"""  SCENARIO 4 : NON-ADJUSTABLE PRESSURE SAFETY - NO REDUNDANCY
https://miro.com/app/board/o9J_lJ__Sfg=/?moveToWidget=3458764519155100294&cot=14
"""

r_pss = rlblt_or(r_rb_a, r_as_a, r_hw_a, r_ry_a)
print("  -----------------------------------------------------  ")
print("SCENARIO 4 : NON-ADJUSTABLE PRESSURE SAFETY - NO REDUNDANCY")
print("Reliability value: ",r_pss)
print("Failure rate: ", round((1-r_pss)*1e7,3)," Failure in 10 Million")
print("  -----------------------------------------------------  ")
"""_____________________________________________________________________________________"""


"""_____________________________________________________________________________________"""

"""  SCENARIO 5 : NON-ADJUSTABLE PRESSURE SAFETY - REDUNDANCY
https://miro.com/app/board/o9J_lJ__Sfg=/?moveToWidget=3458764519155100294&cot=14
"""
r_ps_tot = rlblt_and(r_ps_a,r_ps_b)

r_ai = ((1e7 - 1) / 1e7)  # Reliability of adjustment interface is assumed to be in 10 M failure

r_sic = rlblt_or(r_ai,r_ps_tot,r_hw_a) # Reliability switch IC
r_pss =  rlblt_or(rlblt_and(r_ry_b,r_ry_b),r_hw_a, r_sic)
print("  -----------------------------------------------------  ")
print("SCENARIO 5 : NON-ADJUSTABLE PRESSURE SAFETY - NO REDUNDANCY")
print("Reliability value: ",r_pss)
print("Failure rate: ", round((1-r_pss)*1e7,3)," Failure in 10 Million")
print("  -----------------------------------------------------  ")
"""_____________________________________________________________________________________"""


