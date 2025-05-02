import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import Main


# main.py
import LiveGCBaseline
import LiveGCBooster
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

try:
    T_amb_C_max = Main.T_amb_C_max
except:
    T_amb_C_max = 50
i=0

S_Base_l = []
S_Boost_l = []


for i in range(T_amb_C_max):

    COP_Base = LiveGCBaseline.COP(T_amb_C_max-i)
    S_Base_l.append(COP_Base[1])
    # print(S_Base_l)

    # COP_Booster = LiveGCBooster.COP(i,TtoPlot)
    # S_Boost_l.append(COP_Booster[0])

    i+=1


# Create a DataFrame from the lists

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import CoolProp.CoolProp as CP


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

# Sample data: List of lists, each containing [H, P, Label] triples
# Replace this with your actual data
data = [
    [[2e5, 1e6, "1"], [3e5, 2e6, "2"], [4e5, 3e6, "3"], [5e5, 4e6, "4"]],
    [[2.5e5, 1.5e6, "A"], [3.5e5, 2.5e6, "B"], [4.5e5, 3.5e6, "C"], [5.5e5, 4.5e6, "D"]],
    [[3e5, 2e6, "X"], [4e5, 3e6, "Y"], [5e5, 4e6, "Z"], [6e5, 5e6, "W"]],
]
data = S_Base_l
# Set up the figure and axis
fig, ax = plt.subplots(figsize=(12, 8))

# Plot the saturation curve
ax.plot(h_liq, p_sat, 'b', h_vap, p_sat, 'r')
ax.set_yscale('log')
ax.set_xlabel('Enthalpy (J/kg)')
ax.set_ylabel('Pressure (Pa)')
ax.set_title('Pressure-Enthalpy Diagram for CO2 Refrigeration System')

# Initialize empty lists to store the lines and text annotations
lines = []
texts = []

scatter_points = []
texts = []


# Animation update function
def update(frame):
    # Clear previous scatter points and text annotations
    for scatter in scatter_points:
        scatter.remove()
    for text in texts:
        text.remove()
    scatter_points.clear()
    texts.clear()

    # Get the current data set
    current_data = data[frame]

    # Separate H, P coordinates, and labels
    h_coords, p_coords, labels = zip(*current_data)

    # Plot new data as scatter points
    new_scatter = ax.scatter(h_coords, p_coords, c='g', s=50)  # Green circles
    scatter_points.append(new_scatter)

    # Add labels to each point
    for h, p, label in current_data:
        text = ax.text(h, p, label, fontsize=10, ha='right', va='bottom')
        texts.append(text)

    return scatter_points + texts


# Create the animation
anim = animation.FuncAnimation(fig, update, frames=len(data), interval=200, blit=True, repeat=False)

plt.tight_layout()
plt.show()
