from mdanalysis import read_residue_pairs
from mdanalysis import compute_distances
from mdanalysis import mod_csv
from mdanalysis import calc_zcoords
from mdanalysis import compute_inward_open_contact
from mdanalysis import compute_outward_open_contact
from mdanalysis import calc_contact_mean
from mdanalysis import calc_contact_mean
from mdanalysis import calc_contact_delta
from mdanalysis import filt_contact
from mdanalysis import convert_to_csv

# Extract distance data from MD trajectories
topologies = [
    "../MSM_FESred/round0/7XK3/system.pdb", 
    ]

dcds = [
    "../MSM_FESred/round0/7XK3/prod1.dcd", 
    ]

residue_pairs = read_residue_pairs("input_features.csv")
outs = [
    "../MSM_FESred/round0/7XK3/_pair_distances.csv", 
    ]
outs2 = [
    "../MSM_FESred/round0/7XK3/pair_distances.csv", 
    ]

reference_file = '../MSM_FESred/round0/7XK3/system.pdb'

for topology, dcd, out, out2 in zip(topologies, dcds, outs, outs2):
    compute_distances(topology, dcd, residue_pairs, out)
    mod_csv(out, out2)
    calc_zcoords(reference_file, topology, dcd)


# Invistigate inward-open- or outward-open- specific contacts in NqrD/E
topology_file = "../MSM_FESred/round0/7XK3/system.pdb"
trajectory_file = "../MSM_FESred/round0/7XK3/prod1.dcd"

out1 = "../MSM_FESred/round0/7XK3/contact_inward_open.npy"
frame1 = "../MSM_FESred/round0/7XK3/frames_inward_open.txt"
compute_inward_open_contact(topology_file, trajectory_file, out1, frame1, 50, 100)

out2 = "../MSM_FESred/round0/7XK3/contact_outward_open.npy"
frame2 = "../MSM_FESred/round0/7XK3/frames_outward_open.txt"
compute_outward_open_contact(topology_file, trajectory_file, out2, frame2, 50, 10)

inward_opens = [
    "../MSM_FESred/round0/7XK3/contact_inward_open.npy", 
]
out_io = "../MSM_FESred/msm_data/ave_contact_inward_open.npy"
calc_contact_mean(inward_opens, out_io)

outward_opens = [
    "../MSM_FESred/round0/7XK3/contact_outward_open.npy", 
]
out_oo = "../MSM_FESred/msm_data/ave_contact_outward_open.npy"
calc_contact_mean(outward_opens, out_oo)

out_delta = "../MSM_FESred/msm_data/delta_contact.npy"
calc_contact_delta(out_oo, out_io, out_delta)

thresh = 0.4
out_thresh = "../MSM_FESred/msm_data/specific_contacts.dat"
filt_contact(out_delta, thresh, out_thresh)

out_csv = "../script/input_features.csv"
convert_to_csv(out_thresh, out_csv)

