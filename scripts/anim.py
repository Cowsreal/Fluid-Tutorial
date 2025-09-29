import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob

# Path to CSV files
files = sorted(glob.glob("data/density_step*.csv"))

# Read the first file to get grid size
first_file = files[0]
data = np.loadtxt(first_file, comments='#', delimiter=',')
ny, nx = data.shape

# Set up the figure
fig, ax = plt.subplots()
im = ax.imshow(data, origin='lower', cmap='viridis', interpolation='nearest')
cbar = fig.colorbar(im)
cbar.set_label('Density')

def update(frame):
    filename = files[frame]
    data = np.loadtxt(filename, comments='#', delimiter=',')
    im.set_data(data)
    ax.set_title(f"Step {frame}, File: {filename}")
    return [im]

ani = animation.FuncAnimation(fig, update, frames=len(files), interval=200, blit=True)

plt.show()

ani.save('density_animation.mp4', fps=5, dpi=150)
