# MASTER LOG FILE FOR ALL EMAILS: /home/hmalik342/Working/Email/MahediMasum24_Workspace/Mahedi_Email/Text_Logs/sent_emails.txt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.distance import cdist

# --- Parameters ---
R = 50.0
N_cox_init = 25
N_aa_init = 200
N_asp_dose = 60
dt = 0.5
D_fast = 6.0
D_slow = 1.0
rxn_radius = 12.0
total_frames = 250
intervention_frame = 100

np.random.seed(42)

def random_inside(n, radius):
    pos = np.random.randn(n, 3)
    pos /= np.linalg.norm(pos, axis=1)[:, np.newaxis]
    r = np.random.uniform(0, radius**3, n)**(1/3)
    return pos * r[:, np.newaxis]

def reflect_inside(pos, radius):
    if len(pos) == 0: return pos
    dist = np.linalg.norm(pos, axis=1)
    out = dist > radius
    if np.any(out):
        pos[out] = pos[out] * (2*radius/dist[out, None] - 1)
    return pos

pos_cox = random_inside(N_cox_init, R)
pos_cox_inact = np.empty((0, 3))
pos_aa = random_inside(N_aa_init, R)
pos_pg = np.empty((0, 3))
pos_asp = np.empty((0, 3))

time_arr = []
count_pg = []
count_cox_active = []
count_cox_inact = []

fig = plt.figure(figsize=(16, 7))
ax_cell = fig.add_subplot(121, projection='3d')
ax_graph = fig.add_subplot(122)

# Sphere
u = np.linspace(0, 2 * np.pi, 20)
v = np.linspace(0, np.pi, 10)
x_sphere = R * np.outer(np.cos(u), np.sin(v))
y_sphere = R * np.outer(np.sin(u), np.sin(v))
z_sphere = R * np.outer(np.ones(np.size(u)), np.cos(v))
ax_cell.plot_wireframe(x_sphere, y_sphere, z_sphere, color='#34495e', alpha=0.1)
ax_cell.set_xlim([-R, R])
ax_cell.set_ylim([-R, R])
ax_cell.set_zlim([-R, R])
ax_cell.set_title("3D Cell: Aspirin Inhibition of COX", fontsize=16, fontweight='bold', pad=20)
ax_cell.axis('off')

# Scatter
scat_cox = ax_cell.scatter([0], [0], [0], c='blue', s=150, label='Active COX', alpha=0.9, edgecolors='white')
scat_cox_inact = ax_cell.scatter([0], [0], [0], c='gray', s=150, label='Inactivated COX', alpha=0.9, edgecolors='black')
scat_aa = ax_cell.scatter([0], [0], [0], c='white', edgecolors='gray', s=30, label='Arachidonic Acid')
scat_pg = ax_cell.scatter([0], [0], [0], c='red', s=40, label='Prostaglandin (Inflammation)')
scat_asp = ax_cell.scatter([0], [0], [0], c='yellow', s=60, label='Aspirin (Drug)', marker='D', edgecolors='orange')

ax_cell.legend(loc='upper right', bbox_to_anchor=(1.25, 1.05), fontsize=10, frameon=True)

# Graph
ax_graph.set_xlim(0, total_frames)
ax_graph.set_ylim(0, N_aa_init + 10)
ax_graph.set_title("Inflammation Levels vs Time", fontsize=16, fontweight='bold', pad=15)
ax_graph.set_xlabel("Time (frames)", fontsize=12)
ax_graph.set_ylabel("Molecule Count", fontsize=12)
ax_graph.grid(True, linestyle='--', alpha=0.5)

line_pg, = ax_graph.plot([], [], c='red', lw=4, label='Prostaglandin (Inflammation)')
line_cox_a, = ax_graph.plot([], [], c='blue', lw=2, linestyle='--', label='Active COX')
line_cox_i, = ax_graph.plot([], [], c='gray', lw=2, linestyle=':', label='Inactivated COX')
ax_graph.axvline(x=intervention_frame, color='orange', linestyle='-', lw=3, label='Aspirin Administered', alpha=0.7)
ax_graph.legend(loc='upper left', fontsize=10)

def update_scatter3d(scatter, pos):
    if len(pos) > 0:
        scatter._offsets3d = (pos[:,0], pos[:,1], pos[:,2])
    else:
        scatter._offsets3d = ([], [], [])

def update(frame):
    global pos_cox, pos_cox_inact, pos_aa, pos_pg, pos_asp
    
    # Introduce drug
    if frame == intervention_frame:
        pos_asp = random_inside(N_asp_dose, R)
        
    ax_cell.view_init(elev=20., azim=frame * (360/total_frames))
    
    # Diffusion
    if len(pos_cox) > 0: pos_cox += np.sqrt(2*D_slow*dt)*np.random.randn(len(pos_cox), 3)
    if len(pos_cox_inact) > 0: pos_cox_inact += np.sqrt(2*D_slow*dt)*np.random.randn(len(pos_cox_inact), 3)
    if len(pos_aa) > 0: pos_aa += np.sqrt(2*D_fast*dt)*np.random.randn(len(pos_aa), 3)
    if len(pos_pg) > 0: pos_pg += np.sqrt(2*D_fast*dt)*np.random.randn(len(pos_pg), 3)
    if len(pos_asp) > 0: pos_asp += np.sqrt(2*D_fast*dt)*np.random.randn(len(pos_asp), 3)
        
    # Reactions
    # 1. Aspirin + Active COX -> Inactive COX (Irreversible)
    if len(pos_asp) > 0 and len(pos_cox) > 0:
        dists = cdist(pos_asp, pos_cox)
        min_d = np.min(dists, axis=1)
        reacted_asp = np.where(min_d < rxn_radius)[0]
        if len(reacted_asp) > 0:
            hit_cox = np.argmin(dists[reacted_asp], axis=1)
            hit_cox = np.unique(hit_cox)
            
            new_inact = pos_cox[hit_cox]
            pos_cox_inact = np.vstack([pos_cox_inact, new_inact]) if len(pos_cox_inact)>0 else new_inact
            pos_cox = np.delete(pos_cox, hit_cox, axis=0)
            
            # Aspirin is consumed when it binds to COX
            pos_asp = np.delete(pos_asp, reacted_asp, axis=0)

    # 2. AA + Active COX -> PG
    if len(pos_aa) > 0 and len(pos_cox) > 0:
        dists = cdist(pos_aa, pos_cox)
        min_d = np.min(dists, axis=1)
        reacted_aa = np.where(min_d < rxn_radius)[0]
        if len(reacted_aa) > 0:
            new_pg = pos_aa[reacted_aa]
            pos_pg = np.vstack([pos_pg, new_pg]) if len(pos_pg)>0 else new_pg
            pos_aa = np.delete(pos_aa, reacted_aa, axis=0)
            
    # Reflections
    pos_cox = reflect_inside(pos_cox, R)
    pos_cox_inact = reflect_inside(pos_cox_inact, R)
    pos_aa = reflect_inside(pos_aa, R)
    pos_pg = reflect_inside(pos_pg, R)
    pos_asp = reflect_inside(pos_asp, R)
    
    # Graph updates
    time_arr.append(frame)
    count_pg.append(len(pos_pg))
    count_cox_active.append(len(pos_cox))
    count_cox_inact.append(len(pos_cox_inact))
    
    line_pg.set_data(time_arr, count_pg)
    line_cox_a.set_data(time_arr, count_cox_active)
    line_cox_i.set_data(time_arr, count_cox_inact)
    
    # Scatters
    update_scatter3d(scat_cox, pos_cox)
    update_scatter3d(scat_cox_inact, pos_cox_inact)
    update_scatter3d(scat_aa, pos_aa)
    update_scatter3d(scat_pg, pos_pg)
    update_scatter3d(scat_asp, pos_asp)
    
    return scat_cox, scat_cox_inact, scat_aa, scat_pg, scat_asp, line_pg, line_cox_a, line_cox_i

print("Generating Aspirin Inhibition Animation... This will take a moment.")
ani = animation.FuncAnimation(fig, update, frames=total_frames, interval=50, blit=False)
Writer = animation.writers['ffmpeg']
writer = Writer(fps=20, metadata=dict(artist='Virtual Cell Lab'), bitrate=1800)
ani.save('aspirin_inhibition_3d.mp4', writer=writer)
print("Simulation complete. Video saved as aspirin_inhibition_3d.mp4")
