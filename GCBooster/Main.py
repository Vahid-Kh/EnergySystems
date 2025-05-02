# main.py
import GCBaseline
import GCBooster
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import PlotLiveBaseline
import PlotLiveBooster

T_amb_C_max = 55
TtoPlot = 40
i=0

COP_Base_l = []
COP_Boost_l = []
COP_Boost_gain_l = []
for i in range(T_amb_C_max):

    COP_Base = GCBaseline.COP(i,TtoPlot)
    COP_Base_l.append(COP_Base)

    COP_Booster = GCBooster.COP(i,TtoPlot)
    COP_Boost_l.append(COP_Booster[0])
    COP_Boost_gain_l.append(COP_Booster[1])

    i+=1

    print(COP_Booster,COP_Base)


# Create a DataFrame from the lists
df = pd.DataFrame({
    'COP_Base_l': COP_Base_l,
    'COP_Boost_l': COP_Boost_l,
    'COP_Boost_gain_l': COP_Boost_gain_l,

})

# Plot the data
plt.figure(figsize=(12, 8))
plt.plot(df['COP_Base_l'], label='COP Base')
plt.plot(df['COP_Boost_l'], label='COP Booster')
plt.plot(df['COP_Boost_gain_l'], label='Booster cooling/comp. power')

plt.xlabel('Index')
plt.ylabel('COP Value')
plt.title('COP Booster vs COP Base')
plt.legend()

# import subprocess
# subprocess.run(["python", "PlotLiveBaseline.py"])
# import subprocess
# subprocess.run(["python", "PlotLiveBooster.py"])

# import os
# os.system("PlotLiveBaseline.py")

# with open("PlotLiveBaseline.py") as file:
#     exec(file.read())

plt.show()


