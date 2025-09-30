import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import os

# Get file lists for each variable
files_rho = sorted(glob.glob("../data/rho/step*.csv"))
files_u   = sorted(glob.glob("../data/u/step*.csv"))
files_v   = sorted(glob.glob("../data/v/step*.csv"))
files_p   = sorted(glob.glob("../data/pressure/step*.csv"))

nframes = min(len(files_rho), len(files_u), len(files_v), len(files_p))

# Load first file to get grid size
first_file = np.loadtxt(files_rho[0], comments='#', delimiter=',')
ny, nx = first_file.shape

# Set up figure with 2x2 subplots
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
titles = ["Density (ρ)", "Velocity u", "Velocity v", "Pressure p"]

# Initialize images
rho = np.loadtxt(files_rho[0], comments='#', delimiter=',')
u   = np.loadtxt(files_u[0], comments='#', delimiter=',')
v   = np.loadtxt(files_v[0], comments='#', delimiter=',')
p   = np.loadtxt(files_p[0], comments='#', delimiter=',')

fields = [rho, u, v, p]
ims = []
for ax, title, field in zip(axes.flat, titles, fields):
    im = ax.imshow(field, origin="lower", cmap="viridis", interpolation="nearest")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, shrink=0.7)
    ims.append(im)

def update(frame):
    rho = np.loadtxt(files_rho[frame], comments='#', delimiter=',')
    u   = np.loadtxt(files_u[frame], comments='#', delimiter=',')
    v   = np.loadtxt(files_v[frame], comments='#', delimiter=',')
    p   = np.loadtxt(files_p[frame], comments='#', delimiter=',')

    new_fields = [rho, u, v, p]
    for im, field in zip(ims, new_fields):
        im.set_data(field)

    fig.suptitle(f"Step {frame}")
    return ims

ani = animation.FuncAnimation(fig, update, frames=nframes, interval=200, blit=False)

plt.tight_layout()
plt.show()

ani.save("flow_fields.mp4", fps=5, dpi=150)
