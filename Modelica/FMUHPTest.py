from fmpy import simulate_fmu
import os
import matplotlib.pyplot as plt

# Path to your FMU file
fmu_filename = 'CircuitInertiaFMU.fmu'

# Define simulation options
start_time = 0.0
stop_time = 1.0  # Adjust based on your model
input_value = 5.0  # Example input

import numpy as np

# Create structured input array
input_data = np.array([
    (start_time, input_value),
    (stop_time, input_value)
], dtype=[('time', np.float64), ('valueInput', np.float64)])


# Run the simulation
result = simulate_fmu(
    filename=fmu_filename,
    start_time=start_time,
    stop_time=stop_time,
    input=input_data,
    output=['valueOutput']
)

# Access the output
for row in result:
    time = row[0]
    output_value = row[1]
    print(f"Time: {time:.2f}, valueOutput: {output_value}")
