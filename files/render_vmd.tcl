display resize 800 600
display projection Orthographic
color Display Background white
axes location Off

mol new md_0_1.tpr waitfor all
mol addfile md_cluster.xtc waitfor all

set num_frames [molinfo top get numframes]
puts "Total frames: $num_frames"

mol delrep 0 top

# Protein
mol representation Cartoon
mol color ColorID 7
mol selection "protein"
mol material Opaque
mol addrep top

# Ligand
mol representation Licorice 0.3 12 12
mol color Name
mol selection "resname UNL"
mol material Opaque
mol addrep top

display resetview
scale by 1.5

file mkdir vmd_frames

set frame_idx 0
for {set i 0} {$i < $num_frames} {incr i 100} {
    animate goto $i
    display update
    set filename [format "vmd_frames/frame_%04d.dat" $frame_idx]
    render Tachyon $filename
    incr frame_idx
}

quit
