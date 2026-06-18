# MASTER LOG FILE FOR ALL EMAILS: /home/hmalik342/Working/Email/MahediMasum24_Workspace/Mahedi_Email/Text_Logs/sent_emails.txt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.spatial.distance import cdist

# --- Parameters ---
R = 100.0             # Cell radius
N_enzymes = 15        # Number of enzymes
N_substrates = 300    # Initial number of substrates
dt = 0.5              # Time step
D_sub = 8.0           # Diffusion coeff for substrate/product
D_enz = 1.5           # Diffusion coeff for enzyme
reaction_radius = 5.0 # Distance for reaction to occur
total_frames = 300    # Duration of video (300 frames @ 20fps = 15 seconds)

# --- Initialization ---
np.random.seed(42) # For reproducibility

# Random positions within circle
def random_in_circle(n, radius):
    r = radius * np.sqrt(np.random.rand(n))
    theta = 2 * np.pi * np.random.rand(n)
    return np.column_stack((r * np.cos(theta), r * np.sin(theta)))

pos_e = random_in_circle(N_enzymes, R)
pos_s = random_in_circle(N_substrates, R)
pos_p = np.empty((0, 2))

# Arrays to keep track of counts for plotting
time_arr = []
sub_count = []
prod_count = []

# --- Figure setup ---
fig, (ax_cell, ax_graph) = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={'width_ratios': [1, 1]})
plt.subplots_adjust(wspace=0.3)

# Setup Cell Axis
ax_cell.set_xlim(-R-15, R+15)
ax_cell.set_ylim(-R-15, R+15)
ax_cell.set_aspect('equal')
ax_cell.set_title("2D Virtual Cell: Enzymatic Reaction", fontsize=16, fontweight='bold', pad=15)
ax_cell.axis('off')
membrane = plt.Circle((0, 0), R, color='#2c3e50', fill=False, linewidth=4)
ax_cell.add_patch(membrane)
ax_cell.text(0, -R-20, 'Cell Membrane', horizontalalignment='center', fontsize=12, color='#2c3e50')

# Plot Particles
scatter_e = ax_cell.scatter([], [], c='#3498db', s=120, label='Enzyme', edgecolors='white', zorder=3, alpha=0.9)
scatter_s = ax_cell.scatter([], [], c='#e74c3c', s=30, label='Substrate', alpha=0.7, zorder=2)
scatter_p = ax_cell.scatter([], [], c='#2ecc71', s=30, label='Product', alpha=0.7, zorder=2)
ax_cell.legend(loc='upper right', bbox_to_anchor=(1.25, 1.05), frameon=True)

# Setup Graph Axis
ax_graph.set_xlim(0, total_frames)
ax_graph.set_ylim(0, N_substrates + 20)
ax_graph.set_title("Real-time Reaction Kinetics", fontsize=16, fontweight='bold', pad=15)
ax_graph.set_xlabel("Time (frames)", fontsize=12)
ax_graph.set_ylabel("Number of Molecules", fontsize=12)
ax_graph.grid(True, linestyle='--', alpha=0.5)

line_s, = ax_graph.plot([], [], c='#e74c3c', lw=3, label='Substrate')
line_p, = ax_graph.plot([], [], c='#2ecc71', lw=3, label='Product')
ax_graph.legend(loc='center right')

def reflect_boundary(pos, radius):
    dist = np.linalg.norm(pos, axis=1)
    out_of_bounds = dist > radius
    if np.any(out_of_bounds):
        # Push particle slightly inside the boundary and reverse its direction
        # Here we just project it back proportionally inside the circle
        pos[out_of_bounds] = pos[out_of_bounds] * (2 * radius / dist[out_of_bounds, None] - 1)
    return pos

def update(frame):
    global pos_e, pos_s, pos_p
    
    # 1. Brownian motion (Random Walk)
    pos_e += np.sqrt(2 * D_enz * dt) * np.random.randn(N_enzymes, 2)
    
    if len(pos_s) > 0:
        pos_s += np.sqrt(2 * D_sub * dt) * np.random.randn(len(pos_s), 2)
        
    if len(pos_p) > 0:
        pos_p += np.sqrt(2 * D_sub * dt) * np.random.randn(len(pos_p), 2)
        
    # 2. Reflective boundary
    pos_e = reflect_boundary(pos_e, R)
    if len(pos_s) > 0:
        pos_s = reflect_boundary(pos_s, R)
    if len(pos_p) > 0:
        pos_p = reflect_boundary(pos_p, R)
        
    # 3. Reaction Kinetics (Substrate -> Product if near Enzyme)
    if len(pos_s) > 0 and len(pos_e) > 0:
        dists = cdist(pos_s, pos_e) # NxM distance matrix between all substrates and enzymes
        min_dists = np.min(dists, axis=1) # Distance to closest enzyme for each substrate
        reacted_idx = np.where(min_dists < reaction_radius)[0]
        
        if len(reacted_idx) > 0:
            new_products = pos_s[reacted_idx]
            pos_p = np.vstack([pos_p, new_products]) if len(pos_p) > 0 else new_products
            pos_s = np.delete(pos_s, reacted_idx, axis=0)
            
    # 4. Update data for plotting
    time_arr.append(frame)
    sub_count.append(len(pos_s))
    prod_count.append(len(pos_p))
    
    # 5. Update graphical objects
    scatter_e.set_offsets(pos_e)
    if len(pos_s) > 0:
        scatter_s.set_offsets(pos_s)
    else:
        scatter_s.set_offsets(np.empty((0,2)))
        
    if len(pos_p) > 0:
        scatter_p.set_offsets(pos_p)
    else:
        scatter_p.set_offsets(np.empty((0,2)))
        
    line_s.set_data(time_arr, sub_count)
    line_p.set_data(time_arr, prod_count)
    
    return scatter_e, scatter_s, scatter_p, line_s, line_p

# Create the animation
print("Generating Virtual Cell Animation... This might take a minute.")
ani = animation.FuncAnimation(fig, update, frames=total_frames, interval=50, blit=True)

# Save as MP4 using FFmpeg
Writer = animation.writers['ffmpeg']
writer = Writer(fps=20, metadata=dict(artist='Virtual Cell Lab'), bitrate=1800)
ani.save('virtual_cell_simulation.mp4', writer=writer)
print("Simulation complete. Video saved successfully!")
