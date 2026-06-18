# MASTER LOG FILE FOR ALL EMAILS: /home/hmalik342/Working/Email/MahediMasum24_Workspace/Mahedi_Email/Text_Logs/sent_emails.txt
import matplotlib.pyplot as plt
import numpy as np

def read_xvg(filename):
    x, y = [], []
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith(('@', '#')):
                parts = line.split()
                if len(parts) >= 2:
                    x.append(float(parts[0]))
                    y.append(float(parts[1]))
    return np.array(x), np.array(y)

# 1. Protein RMSD
x, y = read_xvg('rmsd_protein.xvg')
plt.figure(figsize=(8, 5))
plt.plot(x, y, color='blue', linewidth=1.5)
plt.title('Protein Backbone RMSD (50 ns)')
plt.xlabel('Time (ns)')
plt.ylabel('RMSD (nm)')
plt.grid(True, linestyle='--', alpha=0.7)
plt.savefig('protein_rmsd.png', dpi=300, bbox_inches='tight')
plt.close()

# 2. Ligand RMSD
try:
    x, y = read_xvg('rmsd_ligand.xvg')
    plt.figure(figsize=(8, 5))
    plt.plot(x, y, color='red', linewidth=1.5)
    plt.title('Ligand RMSD (50 ns)')
    plt.xlabel('Time (ns)')
    plt.ylabel('RMSD (nm)')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig('ligand_rmsd.png', dpi=300, bbox_inches='tight')
    plt.close()
except FileNotFoundError:
    pass

# 3. Protein RMSF
x, y = read_xvg('rmsf_protein.xvg')
plt.figure(figsize=(8, 5))
plt.plot(x, y, color='green', linewidth=1.5)
plt.title('Protein C-alpha RMSF (50 ns)')
plt.xlabel('Residue Index')
plt.ylabel('RMSF (nm)')
plt.grid(True, linestyle='--', alpha=0.7)
plt.savefig('protein_rmsf.png', dpi=300, bbox_inches='tight')
plt.close()

# 4. Radius of Gyration
x, y = read_xvg('gyrate_protein.xvg')
plt.figure(figsize=(8, 5))
plt.plot(x, y, color='purple', linewidth=1.5)
plt.title('Radius of Gyration (50 ns)')
plt.xlabel('Time (ps)') # gyrate outputs time in ps by default
plt.ylabel('Rg (nm)')
plt.grid(True, linestyle='--', alpha=0.7)
plt.savefig('protein_gyrate.png', dpi=300, bbox_inches='tight')
plt.close()

print("Plots successfully generated!")
