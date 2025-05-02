import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

# Create figure and axis
fig, ax = plt.subplots()
line, = ax.plot([], [])

# Set up plot parameters
ax.set_xlim(0, 100)
ax.set_ylim(0, 1)
ax.set_title("Real-time Data Animation")
ax.set_xlabel("Time")
ax.set_ylabel("Value")

# Initialize data
data = []
max_points = 100


def animate(frame):
    # Generate new data point
    new_point = np.random.random()
    data.append(new_point)

    # Remove old data if we exceed max_points
    if len(data) > max_points:
        data.pop(0)

    # Update the line data
    line.set_data(range(len(data)), data)

    return line,


# Create animation
ani = animation.FuncAnimation(fig, animate, frames=200, interval=50, blit=True)

plt.show()
