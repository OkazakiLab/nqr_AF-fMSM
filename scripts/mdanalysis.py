import numpy as np
import pandas as pd
import glob
import re
import os
import csv
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.analysis.align import AlignTraj
import mdtraj as md
import numpy as np

# Make residue_pairs from a csv file (:input_features.csv)
def read_residue_pairs(filename):
    residue_pairs = []
    with open(filename, newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        for row in reader:
            seg1, res1, seg2, res2 = row
            residue_pairs.append((seg1, int(res1), seg2, int(res2)))
    return residue_pairs


def compute_distances(topology_file, trajectory_files, residue_pairs, output_csv):
    u = mda.Universe(topology_file, trajectory_files)
    data = []
    
    # Calculate distances in each frame
    for ts in u.trajectory:
        row = {'Frame': ts.frame}
        for segid1, resid1, segid2, resid2 in residue_pairs:

            group1 = u.select_atoms(f'segid {segid1} and resid {resid1} and not name H*')
            group2 = u.select_atoms(f'segid {segid2} and resid {resid2} and not name H*')
            
            if len(group1) > 0 and len(group2) > 0:
                distances = distance_array(group1.positions, group2.positions)
                min_distance = np.min(distances)  # Closest heavy atoms
            else:
                min_distance = np.nan
            
            row[f'{segid1}{resid1}-{segid2}{resid2}'] = min_distance
        
        data.append(row)
    
    df = pd.DataFrame(data).round(2)
    df.to_csv(output_csv, index=False)
    

def mod_csv(input, output):
    df = pd.read_csv(input)
    df = df.drop(columns=["Frame"], errors="ignore")
    df.to_csv(output, index=False)


# Calculate sodium z-coodinate after alignment to a reference structure(PDB ID: 7XK3)
def calc_zcoords(reference_file, topology_file, trajectory_file):
    reference_u = mda.Universe(reference_file)
    output_dir = os.path.dirname(os.path.abspath(topology_file))

    u = mda.Universe(topology_file, trajectory_file)
    AlignTraj(u, reference_u, select='(segid NQRD or segid NQRE) and name CA', in_memory=False).run()
    frame_data = []

    for ts in u.trajectory:
        sel1 = u.select_atoms("name FE1 and segid N1")
        sel2 = u.select_atoms("element Na and around 20 (segid N1)")

        if len(sel1) == 0 or len(sel2) == 0:
            continue

        pos1 = sel1.positions
        res1 = sel1.resids
        pos2 = sel2.positions
        res2 = sel2.resids
        frame = ts.frame

        frame_data.extend(
            {'Frame': frame, 'Z_coord_n1': pos11[2], 'na_res': resid2, 'Z_coord_na': pos22[2]}
            for pos11, resid2, pos22 in zip(pos1, res2, pos2)
        )

    df = pd.DataFrame(frame_data)
    output_path = os.path.join(output_dir, f'zcoords.csv')    
    df.to_csv(output_path, index=False)


#--------------------------------#
"""

The following scripts: for investigating input features
; determing inward-open- or outward-open- specific contacts in NqrD/E

"""
# Calculate contact-formation ratio in the inward-open state
def compute_inward_open_contact(topology_file, trajectory_file, out, frame_out, chunk_size=100, stride=10):
    topology = md.load_topology(topology_file)

    protein_residues = [res for res in topology.residues if res.is_protein]
    nres = len(protein_residues)
    print(f"Detected protein residues: {nres}")

    inward1 = topology.select("resid 1143 and not element H")
    inward2 = topology.select("resid 1442 and not element H")
    outward1 = topology.select("resid 1224 and not element H")
    outward2 = topology.select("resid 1350 and not element H")

    pair1 = np.array([(i, j) for i in inward1 for j in inward2])
    pair2 = np.array([(i, j) for i in outward1 for j in outward2])

    avg = np.zeros((nres, nres), dtype=np.float32)
    ntot = 0
    total_selected_frames = 0

    for traj in md.iterload(trajectory_file, top=topology, chunk=chunk_size, stride=stride):
        all_pairs = np.vstack((pair1, pair2))
        alldist = md.compute_distances(traj, all_pairs)

        n1, n2 = len(pair1), len(pair2)
        min_dist1 = alldist[:, :n1].min(axis=1)
        min_dist2 = alldist[:, n1:].min(axis=1)

        selected_frames = np.where((min_dist1 > 0.45) & (min_dist2 < 0.45))[0]   # inward>4.5Å, outward<4.5Å 
        total_selected_frames += len(selected_frames)
        print(f"Selected frames in chunk: {len(selected_frames)}")

        if len(selected_frames) > 0:
            filtered_traj = traj.slice(selected_frames)

            dist, pair = md.compute_contacts(
                filtered_traj,
                contacts="all",
                scheme="closest-heavy",
                ignore_nonprotein=True,
                periodic=False
            )

            contact = (dist < 0.45).astype(np.float32)

            frame_sums = contact.sum(axis=0)
            avg[pair[:, 0], pair[:, 1]] += frame_sums
            avg[pair[:, 1], pair[:, 0]] += frame_sums

            ntot += contact.shape[0]

    if ntot > 0:
        avg /= float(ntot)

    np.save(out, avg)

    with open(frame_out, "w") as f:
        f.write(f"Total selected frames: {total_selected_frames}\n")


# Calculate contact-formation ratio in the outward-open state
def compute_outward_open_contact(topology_file, trajectory_file, out, frame_out, chunk_size=100, stride=10):
    topology = md.load_topology(topology_file)

    protein_residues = [res for res in topology.residues if res.is_protein]
    nres = len(protein_residues)
    print(f"Detected protein residues: {nres}")

    inward1 = topology.select("resid 1143 and not element H")
    inward2 = topology.select("resid 1442 and not element H")
    outward1 = topology.select("resid 1224 and not element H")
    outward2 = topology.select("resid 1350 and not element H")

    pair1 = np.array([(i, j) for i in inward1 for j in inward2])
    pair2 = np.array([(i, j) for i in outward1 for j in outward2])

    avg = np.zeros((nres, nres), dtype=np.float32)
    ntot = 0
    total_selected_frames = 0

    for traj in md.iterload(trajectory_file, top=topology, chunk=chunk_size, stride=stride):
        all_pairs = np.vstack((pair1, pair2))
        alldist = md.compute_distances(traj, all_pairs)

        n1, n2 = len(pair1), len(pair2)
        min_dist1 = alldist[:, :n1].min(axis=1)
        min_dist2 = alldist[:, n1:].min(axis=1)

        selected_frames = np.where((min_dist1 < 0.45) & (min_dist2 > 0.45))[0]   # inward<4.5Å, outward>4.5Å 
        total_selected_frames += len(selected_frames)
        print(f"Selected frames in chunk: {len(selected_frames)}")

        if len(selected_frames) > 0:
            filtered_traj = traj.slice(selected_frames)

            dist, pair = md.compute_contacts(
                filtered_traj,
                contacts="all",
                scheme="closest-heavy",
                ignore_nonprotein=True,
                periodic=False
            )

            contact = (dist < 0.45).astype(np.float32)

            frame_sums = contact.sum(axis=0)
            avg[pair[:, 0], pair[:, 1]] += frame_sums
            avg[pair[:, 1], pair[:, 0]] += frame_sums

            ntot += contact.shape[0]

    if ntot > 0:
        avg /= float(ntot)

    np.save(out, avg)

    with open(frame_out, "w") as f:
        f.write(f"Total selected frames: {total_selected_frames}\n")


def calc_contact_mean(contact_files, output):
    conc = np.stack(
        [np.load(f) for f in contact_files],
        axis=0
    )
    avg = np.mean(conc, axis=0)
    np.save(output, avg)


def calc_contact_delta(contact_outward_open, contact_inward_open, out_delta):
    outward_open = np.load(contact_outward_open)
    inward_open = np.load(contact_inward_open)
    delta = outward_open - inward_open
    np.save(out_delta, delta)



def filt_contact(delta_contact, thresh, out_thresh):
    delta = np.load(delta_contact, allow_pickle=True)
    if delta.ndim != 2:
        raise ValueError(f"delta must be 2D, got shape {delta.shape}")

    nres = delta.shape[0]

    istart = 862
    iend = min(1654, nres)

    with open(out_thresh, "w") as f:
        for i in range(istart, iend):
            for j in range(i + 1, iend):
                val = delta[i, j]

                if val > thresh:
                    f.write(f":{i+1}@CA :{j+1}@CA blue\n")
                elif val < -thresh:
                    f.write(f":{i+1}@CA :{j+1}@CA red\n")


# Make input_features.csv
def convert_to_csv(input, out_csv):
    with open(input, 'r') as file:
        data = file.read()

    def adjust_number(number):
        if 862 <= number <= 1118:
            return number - 861, 'NQRC'
        elif 1119 <= number <= 1328:
            return number - 1118, 'NQRD'
        elif 1329 <= number <= 1526:
            return number - 1328, 'NQRE'
        elif 1527 <= number <= 1934:
            return number - 1526, 'NQRF'
        else:
            return number, ''

    lines = data.strip().split('\n')
    extracted_data = []

    for line in lines:
        numbers = re.findall(r':(\d+)@', line)
        adjusted_numbers = [adjust_number(int(num)) for num in numbers]

        seg1, res1 = adjusted_numbers[0][1], adjusted_numbers[0][0]
        seg2, res2 = adjusted_numbers[1][1], adjusted_numbers[1][0]
    
        extracted_data.append([seg1, res1, seg2, res2])

    df = pd.DataFrame(extracted_data, columns=['seg1', 'res1', 'seg2', 'res2'])
    df.to_csv(out_csv, index=False)    

