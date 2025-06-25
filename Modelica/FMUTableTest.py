import numpy as np
from fmpy import simulate_fmu

# Path to your FMU file
fmu_filename = 'CircuitInertiaFMU.fmu'

# Example input table: time and valueInput
# Replace these with your actual data points
time_input_table = np.array([
    (0.0, 1.0),
    (0.5, 3.5),
    (1.0, 5.0),
    (1.5, 2.5),
    (2.0, 4.0),
], dtype=[('time', np.float64), ('valueInput', np.float64)])

# Define simulation stop time (based on last timestamp in input)
stop_time = time_input_table['time'][-1]

# Simulate FMU
result = simulate_fmu(
    filename=fmu_filename,
    start_time=0.0,
    stop_time=stop_time,
    input=time_input_table,
    output=['valueOutput']
)

# Print the results
for row in result:
    time, output_value = row[0], row[1]
    print(f"Time: {time:.2f}, valueOutput: {output_value}")
