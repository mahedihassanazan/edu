# MASTER LOG FILE FOR ALL EMAILS: /home/hmalik342/Working/Email/MahediMasum24_Workspace/Mahedi_Email/Text_Logs/sent_emails.txt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.distance import cdist

# --- Parameters ---
R = 50.0              # Cell radius
R_ext = R + 20.0      # External boundary
N_ext_glc = 150       # Initial external glucose
N_glut = 12           # Number of Glucose Transporters (Membrane)
N_hk = 8              # Number of Hexokinases (Inside)
N_pk = 8              # Number of Pyruvate Kinases (Inside)
dt = 0.5
D_fast = 6.0          # Diffusion for small molecules (Glc, G6P, Pyr)
D_slow = 1.0          # Diffusion for enzymes
rxn_radius_glut = 10.0 # Reaction distance for GLUT
rxn_radius_enz = 10.0  # Reaction distance for enzymes
total_frames = 250    # Duration of video

# --- Initialization ---
np.random.seed(42)

def random_surface(n, radius):
    phi = np.random.uniform(0, 2*np.pi, n)
    costheta = np.random.uniform(-1, 1, n)
    theta = np.arccos(costheta)
    x = radius * np.sin(theta) * np.cos(phi)
    y = radius * np.sin(theta) * np.sin(phi)
    z = radius * np.cos(theta)
    return np.column_stack((x, y, z))

def random_shell(n, r_min, r_max):
    pos = np.random.randn(n, 3)
    pos /= np.linalg.norm(pos, axis=1)[:, np.newaxis]
    r = np.random.uniform(r_min**3, r_max**3, n)**(1/3)
    return pos * r[:, np.newaxis]

def random_inside(n, radius):
    return random_shell(n, 0, radius)

pos_glut = random_surface(N_glut, R)
pos_hk = random_inside(N_hk, R)
pos_pk = random_inside(N_pk, R)

pos_ext_glc = random_shell(N_ext_glc, R+2, R_ext)
pos_int_glc = np.empty((0, 3))
pos_g6p = np.empty((0, 3))
pos_pyr = np.empty((0, 3))

time_arr = []
count_ext = []
count_int = []
count_pyr = []

# --- Figure Setup ---
fig = plt.figure(figsize=(16, 7))
ax_cell = fig.add_subplot(121, projection='3d')
ax_graph = fig.add_subplot(122)

# Wireframe sphere for the membrane
u = np.linspace(0, 2 * np.pi, 20)
v = np.linspace(0, np.pi, 10)
x_sphere = R * np.outer(np.cos(u), np.sin(v))
y_sphere = R * np.outer(np.sin(u), np.sin(v))
z_sphere = R * np.outer(np.ones(np.size(u)), np.cos(v))
ax_cell.plot_wireframe(x_sphere, y_sphere, z_sphere, color='#34495e', alpha=0.15)

ax_cell.set_xlim([-R_ext, R_ext])
ax_cell.set_ylim([-R_ext, R_ext])
ax_cell.set_zlim([-R_ext, R_ext])
ax_cell.set_title("3D Cell: Glucose Transport & Glycolysis", fontsize=16, fontweight='bold', pad=20)
ax_cell.axis('off')

# Plot Particles
# Using dummy points [0] just to initialize the scatter objects, they will be updated in frame 0
scat_glut = ax_cell.scatter(*pos_glut.T, c='yellow', s=120, label='GLUT (Transporter)', marker='s', edgecolors='black', depthshade=False)
scat_hk = ax_cell.scatter([0], [0], [0], c='blue', s=80, label='Hexokinase', alpha=0.8)
scat_pk = ax_cell.scatter([0], [0], [0], c='purple', s=80, label='Pyruvate Kinase', alpha=0.8)

scat_ext = ax_cell.scatter([0], [0], [0], c='white', edgecolors='gray', s=30, label='Extracellular Glucose')
scat_int = ax_cell.scatter([0], [0], [0], c='cyan', s=30, label='Intracellular Glucose')
scat_g6p = ax_cell.scatter([0], [0], [0], c='orange', s=30, label='G6P / Intermediate')
scat_pyr = ax_cell.scatter([0], [0], [0], c='red', s=30, label='Pyruvate')

ax_cell.legend(loc='upper right', bbox_to_anchor=(1.15, 1.05), fontsize=10, frameon=True)

# Setup Graph
ax_graph.set_xlim(0, total_frames)
ax_graph.set_ylim(0, N_ext_glc + 20)
ax_graph.set_title("Metabolic Flux over Time", fontsize=16, fontweight='bold', pad=15)
ax_graph.set_xlabel("Time (frames)", fontsize=12)
ax_graph.set_ylabel("Number of Molecules", fontsize=12)
ax_graph.grid(True, linestyle='--', alpha=0.5)

line_ext, = ax_graph.plot([], [], c='gray', lw=3, label='Ext Glucose')
line_int, = ax_graph.plot([], [], c='orange', lw=3, label='Int Glucose + Intermediates')
line_pyr, = ax_graph.plot([], [], c='red', lw=3, label='Pyruvate (x2)')
ax_graph.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=3)
plt.subplots_adjust(bottom=0.2)

def reflect_inside(pos, radius):
    if len(pos) == 0: return pos
    dist = np.linalg.norm(pos, axis=1)
    out = dist > radius
    if np.any(out):
        pos[out] = pos[out] * (2*radius/dist[out, None] - 1)
    return pos

def reflect_shell(pos, r_min, r_max):
    if len(pos) == 0: return pos
    dist = np.linalg.norm(pos, axis=1)
    out_outer = dist > r_max
    out_inner = dist < r_min
    if np.any(out_outer):
        pos[out_outer] = pos[out_outer] * (2*r_max/dist[out_outer, None] - 1)
    if np.any(out_inner):
        pos[out_inner] = pos[out_inner] * (2*r_min/dist[out_inner, None] - 1)
    return pos

def update_scatter3d(scatter, pos):
    if len(pos) > 0:
        scatter._offsets3d = (pos[:,0], pos[:,1], pos[:,2])
    else:
        scatter._offsets3d = ([], [], [])

def update(frame):
    global pos_hk, pos_pk, pos_ext_glc, pos_int_glc, pos_g6p, pos_pyr
    
    # Smooth camera rotation
    ax_cell.view_init(elev=20., azim=frame * (360/total_frames))
    
    # 1. Diffusion
    if len(pos_hk) > 0: pos_hk += np.sqrt(2*D_slow*dt)*np.random.randn(N_hk, 3)
    if len(pos_pk) > 0: pos_pk += np.sqrt(2*D_slow*dt)*np.random.randn(N_pk, 3)
    if len(pos_ext_glc) > 0: pos_ext_glc += np.sqrt(2*D_fast*dt)*np.random.randn(len(pos_ext_glc), 3)
    if len(pos_int_glc) > 0: pos_int_glc += np.sqrt(2*D_fast*dt)*np.random.randn(len(pos_int_glc), 3)
    if len(pos_g6p) > 0: pos_g6p += np.sqrt(2*D_fast*dt)*np.random.randn(len(pos_g6p), 3)
    if len(pos_pyr) > 0: pos_pyr += np.sqrt(2*D_fast*dt)*np.random.randn(len(pos_pyr), 3)
        
    # 2. Reactions & Boundaries
    # a) Ext Glucose hitting GLUT
    if len(pos_ext_glc) > 0:
        dist_ext = np.linalg.norm(pos_ext_glc, axis=1)
        crossing = dist_ext < R + 5 # close to membrane
        if np.any(crossing):
            cross_pos = pos_ext_glc[crossing]
            dists = cdist(cross_pos, pos_glut)
            min_d = np.min(dists, axis=1)
            transported_local_idx = np.where(min_d < rxn_radius_glut)[0]
            
            if len(transported_local_idx) > 0:
                global_crossing_idx = np.where(crossing)[0]
                transported_global_idx = global_crossing_idx[transported_local_idx]
                
                new_int = pos_ext_glc[transported_global_idx]
                # Push slightly inside
                new_int = new_int * ((R - 4) / np.linalg.norm(new_int, axis=1)[:, None])
                pos_int_glc = np.vstack([pos_int_glc, new_int]) if len(pos_int_glc)>0 else new_int
                pos_ext_glc = np.delete(pos_ext_glc, transported_global_idx, axis=0)
        
    # Reflect remaining ext_glc
    pos_ext_glc = reflect_shell(pos_ext_glc, R, R_ext)
    
    # b) Int Glucose + HK -> G6P
    if len(pos_int_glc) > 0 and len(pos_hk) > 0:
        dists = cdist(pos_int_glc, pos_hk)
        min_d = np.min(dists, axis=1)
        reacted_idx = np.where(min_d < rxn_radius_enz)[0]
        if len(reacted_idx) > 0:
            new_g6p = pos_int_glc[reacted_idx]
            pos_g6p = np.vstack([pos_g6p, new_g6p]) if len(pos_g6p)>0 else new_g6p
            pos_int_glc = np.delete(pos_int_glc, reacted_idx, axis=0)
            
    # c) G6P + PK -> 2 Pyruvate
    if len(pos_g6p) > 0 and len(pos_pk) > 0:
        dists = cdist(pos_g6p, pos_pk)
        min_d = np.min(dists, axis=1)
        reacted_idx = np.where(min_d < rxn_radius_enz)[0]
        if len(reacted_idx) > 0:
            new_pyr_1 = pos_g6p[reacted_idx]
            new_pyr_2 = new_pyr_1 + np.random.randn(len(reacted_idx), 3) * 2.0
            new_pyr = np.vstack([new_pyr_1, new_pyr_2])
            pos_pyr = np.vstack([pos_pyr, new_pyr]) if len(pos_pyr)>0 else new_pyr
            pos_g6p = np.delete(pos_g6p, reacted_idx, axis=0)

    # Internal reflections
    pos_hk = reflect_inside(pos_hk, R)
    pos_pk = reflect_inside(pos_pk, R)
    pos_int_glc = reflect_inside(pos_int_glc, R)
    pos_g6p = reflect_inside(pos_g6p, R)
    pos_pyr = reflect_inside(pos_pyr, R)

    # 3. Graph updates
    time_arr.append(frame)
    count_ext.append(len(pos_ext_glc))
    count_int.append(len(pos_int_glc) + len(pos_g6p))
    count_pyr.append(len(pos_pyr))
    
    line_ext.set_data(time_arr, count_ext)
    line_int.set_data(time_arr, count_int)
    line_pyr.set_data(time_arr, count_pyr)
    
    # Dynamically adjust y-axis if Pyruvate exceeds initial limit
    if len(pos_pyr) > ax_graph.get_ylim()[1]:
        ax_graph.set_ylim(0, len(pos_pyr) + 20)
    
    # 4. Scatter updates
    update_scatter3d(scat_hk, pos_hk)
    update_scatter3d(scat_pk, pos_pk)
    update_scatter3d(scat_ext, pos_ext_glc)
    update_scatter3d(scat_int, pos_int_glc)
    update_scatter3d(scat_g6p, pos_g6p)
    update_scatter3d(scat_pyr, pos_pyr)
    
    return scat_hk, scat_pk, scat_ext, scat_int, scat_g6p, scat_pyr, line_ext, line_int, line_pyr

print("Generating 3D Virtual Cell Animation... This will take a moment.")
ani = animation.FuncAnimation(fig, update, frames=total_frames, interval=50, blit=False)

Writer = animation.writers['ffmpeg']
writer = Writer(fps=20, metadata=dict(artist='Virtual Cell Lab'), bitrate=1800)
ani.save('glycolysis_3d_simulation.mp4', writer=writer)
print("Simulation complete. Video saved successfully!")
