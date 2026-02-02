import pyemma
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import math
from pyemma.plots.markovtests import _add_ck_subplot
from contextlib import redirect_stdout

np.NAN = np.nan
np.alltrue = np.all

"""
Index
#01. tICA from 43 distance data (input features)
#02. Convert experimental data to tICA data
#03. Make dataset including input features, tics, and Na+ z-coords
#04. Plot MD trajectories along with tIC1-2
#05. Load dataset from csv files
#06. Calculate VAMP2 score
#07. Save VAMP2 score
#08. Plot VAMP2 score
#09. Calculate implied timescale (its)
#10. Plot stationary distribution & free energy
#11. Plot free energy along with gate sizes in NqrD/E
#12. Plot free energy along with Na+ z-coords and inward gate size
#13. Plot free energy using Na+-bound data
#14. Plot free energy using Na+-unbound data
#15. Plot coarse-grained states (PCCA+ result)
#16. Chapman-Kolmogorov test (ck_test)
#17-1,2,3. Transition path analysis

"""


#01. tICA from 43 distance data (input features)
def perform_tica(inputs, dim, outputs):
    data = [pd.read_csv(file).values for file in inputs]
    tica = pyemma.coordinates.tica(data, lag=10, dim=dim)
    tica_outputs = tica.get_output()

    for tica_output, output in zip(tica_outputs, outputs):
        df = pd.DataFrame(tica_output, columns=[f'tIC{i+1}' for i in range(dim)])
        df = df.round(2)
        df.to_csv(output, index=False)

    eigenvectors = tica.eigenvectors
    for i in range(min(dim, eigenvectors.shape[0])): 
        np.savetxt(f'../MSM_FESred/msm_data/tic{i+1}_eigenvector.csv', eigenvectors[:, i], delimiter=',')

    mean_data = tica.mean
    np.savetxt('../MSM_FESred/msm_data/tica_mean.csv', mean_data, delimiter=',')

    eigenvalues = tica.eigenvalues 
    df_contributions = pd.DataFrame({
        'tics': [f'tic{i+1}' for i in range(dim)],
        'eigenvalues': eigenvalues[:dim]
    })
    df_contributions.to_csv('../MSM_FESred/msm_data/tica_eigenvalues.csv', index=False)


#02. Convert experimental data to tICA data
def apply_tica(input_csv, tica_dim, output_csv):
    X_new = np.loadtxt(input_csv, delimiter=',', skiprows=1)
    # Load tICA mean
    mean_data = np.loadtxt('../MSM_FESred/msm_data/tica_mean.csv').reshape(1, -1)

    tic_results = {}
    dim = tica_dim + 1
    for i in range(1, dim):
        vector_path = f'../MSM_FESred/msm_data/tic{i}_eigenvector.csv'
        tic_vector = np.loadtxt(vector_path).reshape(-1, 1)
        tic_projection = np.dot(X_new - mean_data, tic_vector)
        tic_results[f'tIC{i}'] = tic_projection.flatten()

    df = pd.DataFrame(tic_results)
    df = df.round(2)
    df.to_csv(output_csv, index=False)


#03. Make dataset including input features, tics, and Na+ z-coords
def make_dataset(z_coords, distances, tics, output):
    df3 = pd.read_csv(distances, header=0)
    n_frames = len(df3)

    try:
        df1 = pd.read_csv(z_coords, header=0)
    except pd.errors.EmptyDataError:
        df1 = pd.DataFrame()

    if df1.empty:
        df2 = pd.DataFrame({
            "Frame": range(0, n_frames),
            "Z_coord_na": ["none"] * n_frames,
            "na_res": ["none"] * n_frames,
            "Z_coord_n1": ["none"] * n_frames,
        })
    else:
        df1["Frame"] = pd.to_numeric(df1["Frame"], errors="coerce")
        df1["Z_coord_na"] = pd.to_numeric(df1["Z_coord_na"], errors="coerce")
        df1["Z_coord_n1"] = pd.to_numeric(df1["Z_coord_n1"], errors="coerce")

        # keep only the row with the smallest |Z_coord_na| for each frame
        idx_min = df1.groupby("Frame")["Z_coord_na"].apply(lambda x: x.abs().idxmin())
        df1 = df1.loc[idx_min]

        full_Frame = pd.DataFrame({"Frame": range(0, n_frames)})
        df2 = full_Frame.merge(df1, on="Frame", how="left")

        # convert "NaN" into "none"
        df2["Z_coord_na"] = df2["Z_coord_na"].fillna("none")
        df2["na_res"] = df2["na_res"].fillna("none")
        df2["Z_coord_n1"] = df2["Z_coord_n1"].fillna("none")

    z_na = pd.to_numeric(df2["Z_coord_na"], errors="coerce")
    z_n1 = pd.to_numeric(df2["Z_coord_n1"], errors="coerce")

    df2["na_z"] = (z_na - z_n1).round(2)
    df2["na_z"] = df2["na_z"].fillna("none")

    df4 = pd.read_csv(tics, header=0)

    if len(df4) != n_frames:
        raise ValueError(f"Row count mismatch: distances={n_frames}, tics={len(df4)}")

    df5 = df2.copy()

    df5["NQRD26-NQRE115"] = df3["NQRD26-NQRE115"].values
    df5["NQRD107-NQRE23"] = df3["NQRD107-NQRE23"].values

    for i in range(1, 44):
        col = f"tIC{i}"
        if col not in df4.columns:
            raise KeyError(f"Missing column in tics: {col}")
        df5[col] = df4[col].values

    df5.to_csv(output, index=False)


#04. Plot MD trajectories along with tIC1-2
def plot_tics(exp_tics, tica_outputs, out_fig):
    df = pd.read_csv(exp_tics)
    plt.figure(figsize=(5, 4))

    for file in tica_outputs:
        df1 = pd.read_csv(file)
        label = os.path.basename(os.path.dirname(file))
    
        plt.plot(df1['tIC1'], df1['tIC2'], linewidth=0.4, alpha=0.3, label=label)

    plt.scatter(df['tIC1'], df['tIC2'], color='red')

    plt.xlabel('tIC1', fontsize=18)
    plt.ylabel('tIC2', fontsize=18)
    plt.xticks([-1,0,1,2], fontsize=16)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.savefig(out_fig, bbox_inches='tight')    


#05. Load dataset from csv files
def load_columns_from_csv(csv_files, columns):
    data = []
    for f in csv_files:
        df = pd.read_csv(f, header=0)
        dff = df[[c for c in columns if c in df.columns]].copy()
        # 数値化（"None" や文字列は NaN へ）
        for col in dff.columns:
            dff[col] = pd.to_numeric(dff[col], errors="coerce")

        data.append(dff.values.astype(float))
    return data


#06. Calculate VAMP2 score
def cluster_centers(n_clustercenters, tica_output):
    scores = np.zeros((len(n_clustercenters), 5))
    for n, k in enumerate(n_clustercenters):
        for m in range(5):
            cl = pyemma.coordinates.cluster_kmeans(
                tica_output, k=k, max_iter=50, stride=50
            )
            msm = pyemma.msm.estimate_markov_model(cl.dtrajs, lag=10)
            scores[n, m] = msm.score_cv(
                cl.dtrajs, n=1, score_method='VAMP2', score_k=min(10, k)
            )

    try:
        lower, upper = pyemma.util.statistics.confidence_interval(scores.T, conf=0.9)
    except AttributeError:
        mean = np.mean(scores, axis=1)
        std = np.std(scores, axis=1)
        lower = mean - std
        upper = mean + std

    return scores, lower, upper


#07. Save VAMP2 score
def save_plot_data_csv(n_clustercenters, results_dict, output_csv):
    rows = []
    for label, (scores, lower, upper) in results_dict.items():
        mean = np.mean(scores, axis=1)
        error = upper - mean
        for k, m, e in zip(n_clustercenters, mean, error):
            rows.append({'n_cluster': k, 'mean_score': m, 'error': e, 'label': label})
    
    df = pd.DataFrame(rows)
    df.to_csv(output_csv, index=False)


#08. Plot VAMP2 score
def plot_from_csv(csv_file, fig_out):
    df = pd.read_csv(csv_file)

    fig, ax = plt.subplots(figsize=(5,4))

    for label in df['label'].unique():
        sub = df[df['label'] == label]
        ax.errorbar(
            sub['n_cluster'], sub['mean_score'], yerr=sub['error'],
            fmt='-o', capsize=4, label=f'{label}'
        )

#    ax.set_xscale('log')
    ax.set_xlabel('Number of cluster centers', fontsize=18)
    ax.set_ylabel('VAMP-2 score', fontsize=18)
    ax.tick_params(labelsize=16)
    fig.tight_layout()
    plt.savefig(fig_out)


#09. Calculate implied timescale (its)
def calc_its(dataset, tics, k, seed, lags, fig_out):
    columns = [f"tIC{i+1}" for i in range(tics + 1)]
    tica_output = load_columns_from_csv(dataset, columns)
    cluster = pyemma.coordinates.cluster_kmeans(tica_output, k=k, max_iter=50, stride=10, fixed_seed=seed)
    its = pyemma.msm.its(cluster.dtrajs, lags=lags, nits=10)

    fig, ax = plt.subplots(figsize=(5,4))
    pyemma.plots.plot_implied_timescales(its, units='ns', dt=0.1)
    ax.tick_params(labelsize=10)
    fig.tight_layout()
    plt.savefig(fig_out)


#10. Plot stationary distribution & free energy
def plot_distrib_energy(dataset, tics, k, seed, lag, fig_out):
    columns = [f"tIC{i+1}" for i in range(tics + 1)]
    tica_output = load_columns_from_csv(dataset, columns)
    cluster = pyemma.coordinates.cluster_kmeans(tica_output, k=k, max_iter=50, stride=10, fixed_seed=seed)
    msm = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=lag)

    tica_concatenated = np.concatenate(tica_output)
    dtrajs_concatenated = np.concatenate(cluster.dtrajs)

    # MSM-active cluster indexes
    active_set = msm.active_set
    mapping = np.full(np.max(dtrajs_concatenated)+1, -1)
    mapping[active_set] = np.arange(len(active_set))

    # extract only existing-state in MSM
    mapped_dtrajs = mapping[dtrajs_concatenated]
    mask = mapped_dtrajs != -1
    # masked tica coordinates
    tica_masked = tica_concatenated[mask, :2].T  # shape: (2, N)
    z_values = msm.pi[mapped_dtrajs[mask]]       # shape: (N,)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)
    pyemma.plots.plot_contour(
        *tica_masked,
        z_values,
        ax=axes[0],
        mask=True,
        cbar_label='stationary distribution')
    pyemma.plots.plot_free_energy(
        *tica_concatenated[:, :2].T,
        weights=np.concatenate(msm.trajectory_weights()),
        ax=axes[1],
        legacy=False)
    for ax in axes.flat:
        ax.set_xlabel('tIC1', fontsize=18)
        ax.set_xticks([-1,0,1,2])
        ax.set_xticklabels([-1,0,1,2], fontsize=16)
        ax.set_yticks([-2,-1,0,1,2,3])
        ax.set_yticklabels([-2,-1,0,1,2,3], fontsize=16)
    axes[0].set_ylabel('tIC2', fontsize=18)
    axes[0].set_title('Stationary distribution', fontweight='bold', fontsize=16)
    axes[1].set_title('Reweighted free energy surface', fontweight='bold', fontsize=16)
    fig.tight_layout()

    plt.savefig(fig_out)


#11. Plot free energy along with gate sizes in NqrD/E
def plot_free_energy_gate_all(
    dataset, tics, k, seed, lag, fig_out, kT=2.479, nbins=100,
    display_max=12.0, show_ticks=(0, 2, 4, 6, 8, 10, 12)
):
    tic_columns = [f"tIC{i+1}" for i in range(tics + 1)]
    tica_output = load_columns_from_csv(dataset, tic_columns)
    gate_columns = ["NQRD26-NQRE115", "NQRD107-NQRE23"]  # inward gate = "NQRD26-NQRE115", outward gate = "NQRD107-NQRE23"
    gate_data = load_columns_from_csv(dataset, gate_columns)
    cluster = pyemma.coordinates.cluster_kmeans(tica_output, k=k, max_iter=50, stride=10, fixed_seed=seed)
    msm = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=lag)

    gate_data_concat = np.concatenate(gate_data)
    weights = np.concatenate(msm.trajectory_weights())    
    
    # Remove row data "NaN"
    valid_mask = np.isfinite(gate_data_concat[:, 0]) & np.isfinite(gate_data_concat[:, 1]) & np.isfinite(weights)
    x = gate_data_concat[valid_mask, 0]
    y = gate_data_concat[valid_mask, 1]
    w = weights[valid_mask]

    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(111)

    pyemma.plots.plot_free_energy(
        x, y,
        weights=w, kT=kT, cmap='nipy_spectral',
        nbins=nbins, ax=ax, cbar=False
    )

    if ax.collections:
        mappable = ax.collections[-1]
    else:
        mappable = ax.images[-1]

    vmin, vmax = mappable.get_clim()
    data_min, data_max = vmin, vmax

    cbar = fig.colorbar(mappable, ax=ax)
    cbar.set_label('Free energy / kT', fontsize=12)

    tick_positions = [data_min + (data_max - data_min) * (t / display_max) for t in show_ticks]
    cbar.set_ticks(tick_positions)
    cbar.set_ticklabels([str(t) for t in show_ticks])

    ax.set_xlim(2, 14)
    ax.set_ylim(2, 12)
    ax.set_xticks(np.linspace(2, 14, 7))
    ax.set_yticks(np.linspace(2, 12, 6))
    ax.tick_params(labelsize=12)
    ax.set_xlabel('Cytoplasmic gate / Å', fontsize=14)
    ax.set_ylabel('Periplasmic gate / Å', fontsize=14)

    plt.tight_layout()
    plt.savefig(fig_out, dpi=200)


#12. Plot free energy along with Na+ z-coords and inward gate size
def plot_free_energy_zna_inward_all(
    dataset, fig_out,
    kT=2.479, nbins=150, display_max=12.0,
    show_ticks=(0, 2, 4, 6, 8, 10, 12),
    xlim=(2, 14), ylim=(-25, 25),
    interpolation="bilinear",
):
    
    gate_columns = ["NQRD26-NQRE115", "na_z"]
    gate_data = load_columns_from_csv(dataset, gate_columns)
    gate_concatenated = np.concatenate(gate_data)

    x = gate_concatenated[:, 0]
    y = gate_concatenated[:, 1]   # na_z

    # remove na_z = "None"
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    w = np.ones_like(x, dtype=float)

    x_edges = np.linspace(xlim[0], xlim[1], nbins + 1)
    y_edges = np.linspace(ylim[0], ylim[1], nbins + 1)
    H, xe, ye = np.histogram2d(x, y, bins=[x_edges, y_edges], weights=w)

    # convert probability into free energy
    H_sum = H.sum()
    P = H / H_sum if H_sum > 0 else H
    with np.errstate(divide="ignore", invalid="ignore"):
        F = -kT * np.log(P)

    if np.isfinite(F).any():
        F = F - np.nanmin(F)

    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(111)

    im = ax.imshow(
        F.T,
        origin="lower",
        extent=[xe[0], xe[-1], ye[0], ye[-1]],
        cmap="nipy_spectral",
        norm=Normalize(vmin=0.0, vmax=display_max),
        interpolation=interpolation,
        aspect="auto",
    )

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Free energy / kT", fontsize=12)
    cbar.set_ticks(list(show_ticks))
    cbar.set_ticklabels([str(t) for t in show_ticks])

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.invert_yaxis()
    ax.set_xticks(np.linspace(xlim[0], xlim[1], 7))
    ax.set_yticks([-20, -10, 0, 10, 20])
    ax.tick_params(labelsize=12)
    ax.set_xlabel("Cytoplasmic gate / Å", fontsize=14)
    ax.set_ylabel("Relative sodium Z-coordinates / Å", fontsize=14)

    plt.tight_layout()
    plt.savefig(fig_out, dpi=200)


#13. Plot free energy using Na+-bound data
def plot_free_energy_na_bound(
    dataset, tics, k, seed, lag, fig_out,
    kT=2.479, nbins=150, display_max=12.0,
    show_ticks=(0, 2, 4, 6, 8, 10, 12),
    xlim=(2, 14), ylim=(2, 12),
    interpolation='bilinear', 
    gaussian_sigma=None):

    tic_columns = [f"tIC{i+1}" for i in range(tics + 1)]
    tica_output = load_columns_from_csv(dataset, tic_columns)
    gate_columns = ["NQRD26-NQRE115", "NQRD107-NQRE23"]
    gate_data = load_columns_from_csv(dataset, gate_columns)
    zcoords_columns = ["na_z"]
    zcoords_data = load_columns_from_csv(dataset, zcoords_columns)

    cluster = pyemma.coordinates.cluster_kmeans(tica_output, k=k, max_iter=50, stride=10, fixed_seed=seed)
    msm = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=lag)

    gate_concatenated = np.concatenate(gate_data)
    zcoords_concatenated = np.concatenate(zcoords_data).flatten()
    weights = np.concatenate(msm.trajectory_weights())
    

    # define Na+-bound state
    mask = (zcoords_concatenated >= -8) & (zcoords_concatenated <= 8)
    x = gate_concatenated[mask, 0]
    y = gate_concatenated[mask, 1]
    w = weights[mask]

    x_edges = np.linspace(xlim[0], xlim[1], nbins + 1)
    y_edges = np.linspace(ylim[0], ylim[1], nbins + 1)
    H, xe, ye = np.histogram2d(x, y, bins=[x_edges, y_edges], weights=w)

    # convert probability into free energy
    H_sum = H.sum()
    P = H / H_sum if H_sum > 0 else H
    with np.errstate(divide='ignore', invalid='ignore'):
        F = -kT * np.log(P)
    # minimum energy == 0
    if np.isfinite(F).any():
        F = F - np.nanmin(F)

    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(111)

    im = ax.imshow(
        F.T, 
        origin='lower',
        extent=[xe[0], xe[-1], ye[0], ye[-1]],
        cmap='nipy_spectral',
        norm=Normalize(vmin=0.0, vmax=display_max),
        interpolation=interpolation 
    )

    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Free energy / kT', fontsize=12)
    cbar.set_ticks(list(show_ticks))
    cbar.set_ticklabels([str(t) for t in show_ticks])

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xticks(np.linspace(*xlim, 7))
    ax.set_yticks(np.linspace(*ylim, 6))
    ax.tick_params(labelsize=12)
    ax.set_xlabel('Cytoplasmic gate / Å', fontsize=14)
    ax.set_ylabel('Periplasmic gate / Å', fontsize=14)

    plt.tight_layout()
    plt.savefig(fig_out, dpi=200)


#14. Plot free energy using Na+-unbound data
def plot_free_energy_na_unbound(
    dataset, tics, k, seed, lag, fig_out,
    kT=2.479, nbins=150, display_max=12.0, 
    show_ticks=(0, 2, 4, 6, 8, 10, 12),
    xlim=(2, 14), ylim=(2, 12),
    interpolation='bilinear', 
    gaussian_sigma=None
):
    tic_columns = [f"tIC{i+1}" for i in range(tics + 1)]
    tica_output = load_columns_from_csv(dataset, tic_columns)
    gate_columns = ["NQRD26-NQRE115", "NQRD107-NQRE23"]
    gate_data = load_columns_from_csv(dataset, gate_columns)
    zcoords_columns = ["na_z"]
    zcoords_data = load_columns_from_csv(dataset, zcoords_columns)

    cluster = pyemma.coordinates.cluster_kmeans(tica_output, k=k, max_iter=50, stride=10, fixed_seed=seed)
    msm = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=lag)

    gate_concatenated = np.concatenate(gate_data)
    zcoords_concatenated = np.concatenate(zcoords_data).flatten()
    weights = np.concatenate(msm.trajectory_weights())
    

    # define Na+-unbound state (not Na+-bound state)
    mask = (zcoords_concatenated >= -8) & (zcoords_concatenated <= 8)
    x = gate_concatenated[~mask, 0]
    y = gate_concatenated[~mask, 1]
    w = weights[~mask]

    x_edges = np.linspace(xlim[0], xlim[1], nbins + 1)
    y_edges = np.linspace(ylim[0], ylim[1], nbins + 1)
    H, xe, ye = np.histogram2d(x, y, bins=[x_edges, y_edges], weights=w)

    # convert probability into free energy
    H_sum = H.sum()
    P = H / H_sum if H_sum > 0 else H
    with np.errstate(divide='ignore', invalid='ignore'):
        F = -kT * np.log(P)
    # minimum energy == 0
    if np.isfinite(F).any():
        F = F - np.nanmin(F)

    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(111)

    im = ax.imshow(
        F.T, 
        origin='lower',
        extent=[xe[0], xe[-1], ye[0], ye[-1]],
        cmap='nipy_spectral',
        norm=Normalize(vmin=0.0, vmax=display_max),
        interpolation=interpolation
    )

    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Free energy / kT', fontsize=12)
    cbar.set_ticks(list(show_ticks))
    cbar.set_ticklabels([str(t) for t in show_ticks])

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xticks(np.linspace(*xlim, 7))
    ax.set_yticks(np.linspace(*ylim, 6))
    ax.tick_params(labelsize=12)
    ax.set_xlabel('Cytoplasmic gate / Å', fontsize=14)
    ax.set_ylabel('Periplasmic gate / Å', fontsize=14)

    plt.tight_layout()
    plt.savefig(fig_out, dpi=200)


#15. Plot coarse-grained states (PCCA+ result)
def plot_free_energy_with_pcca_states(
        dataset, tics, k, seed, lag, nstates, fig_out,
        kT=2.479, nbins=100, metastates=None):
    
    tic_columns = [f"tIC{i+1}" for i in range(tics + 1)]
    tica_output = load_columns_from_csv(dataset, tic_columns)
    gate_columns = ["NQRD26-NQRE115", "NQRD107-NQRE23"]
    gate_data = load_columns_from_csv(dataset, gate_columns)

    cluster = pyemma.coordinates.cluster_kmeans(tica_output, k=k, max_iter=50, stride=10, fixed_seed=seed)
    msm = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=lag)
    gate_data_concat = np.concatenate(gate_data)
    weights = np.concatenate(msm.trajectory_weights())

    msm.pcca(nstates)
    meta_assignments = msm.metastable_assignments
    dtrajs_concat = np.concatenate(cluster.dtrajs)
    meta_traj = meta_assignments[dtrajs_concat]    

    # remove "NaN"
    valid_mask = np.isfinite(gate_data_concat[:, 0]) & np.isfinite(gate_data_concat[:, 1]) & np.isfinite(weights)
    x_all = gate_data_concat[valid_mask, 0]
    y_all = gate_data_concat[valid_mask, 1]
    w_all = weights[valid_mask]
    meta_traj = meta_traj[valid_mask]

    if metastates is None:
        metastates = np.unique(meta_traj)

    n_states = len(metastates)
    fig, axes = plt.subplots(1, n_states, figsize=(4*n_states, 4),
                             sharex=True, sharey=True)

    if n_states == 1:
        axes = [axes]

    for ax, m_id in zip(axes, metastates):
        # background: free energy
        pyemma.plots.plot_free_energy(
            x_all, y_all, weights=w_all, kT=kT,
            cmap='nipy_spectral', nbins=nbins,
            ax=ax, cbar=(m_id == metastates[-1]), cbar_label='Free energy / kT'
        )

        mask = (meta_traj == m_id)
        ax.scatter(x_all[mask], y_all[mask],
                   c='red', s=5, alpha=0.6)

        ax.set_xlim(2, 14)
        ax.set_ylim(2, 12)
        ax.set_xticks(np.linspace(2, 14, 7))
        ax.set_yticks(np.linspace(2, 12, 6))
        ax.tick_params(labelsize=12)
        ax.set_xlabel('Cytoplasmic gate / Å', fontsize=14)
        if ax == axes[0]:
            ax.set_ylabel('Periplasmic gate / Å', fontsize=14)
        ax.legend(loc='best', frameon=False)

    plt.tight_layout()
    plt.savefig(fig_out)


#16. Chapman-Kolmogorov test (ck_test)
def plot_ck_test(
    dataset, tics, k, seed, lag, nstates, outfile,
    *, diag=False, y01=True, layout=None,
    figsize=None, units="ns", dt=0.1,
    padding_between=0.1, padding_top=0.075,
    dpi=300, bbox_inches="tight", close=True,
    tick_labelsize=10,
    **plot_kwargs
):
    # ---- MSM + CK test ----
    cols = [f"tIC{i+1}" for i in range(tics + 1)]
    X = load_columns_from_csv(dataset, cols)
    cl = pyemma.coordinates.cluster_kmeans(X, k=k, max_iter=50, stride=10, fixed_seed=seed)
    msm = pyemma.msm.estimate_markov_model(cl.dtrajs, lag=lag, dt_traj="1 ns")
    ck = msm.cktest(nstates)
    n = ck.nsets

    # ---- layout / figsize ----
    if diag:
        if layout in (None, "wide"):
            ncol = min(4, n)
            layout = (math.ceil(n / ncol), ncol)
        elif layout == "tall":
            nrow = min(4, n)
            layout = (nrow, math.ceil(n / nrow))
        # (rows, cols) が渡されてきたらそのまま
    else:
        layout = (n, n)

    if figsize is None:
        s = min(3.0, 10.0 / max(layout))
        figsize = (s * layout[1], s * layout[0])

    fig, axes = plt.subplots(*layout, sharex=True, sharey=y01, figsize=figsize)
    axs = np.atleast_1d(axes).ravel()

    # ---- plot ----
    last_lest = last_lpred = None
    for kk, ax in enumerate(axs):
        if diag and kk >= n:
            ax.axis("off")
            continue
        i, j = (kk, kk) if diag else divmod(kk, n)
        extra = {"ipos": kk // layout[1], "jpos": kk % layout[1]} if diag else {}
        last_lest, last_lpred = _add_ck_subplot(
            ck, 0, ax, i, j, y01=y01, units=units, dt=dt, **extra, **plot_kwargs
        )

    # ---- legend / styling / save ----
    if last_lest and last_lpred:
        label = "predict" + (f"     conf. {100.0 * ck.conf:3.1f}%" if ck.has_errors else "")
        fig.legend([last_lest[0], last_lpred[0]], [label, "estimate"],
                   loc="upper center", ncol=2, frameon=False)

    for ax in fig.axes:
        ax.tick_params(labelsize=tick_labelsize)

    fig.subplots_adjust(top=1.0 - padding_top, wspace=padding_between, hspace=padding_between)
    fig.savefig(outfile, dpi=dpi, bbox_inches=bbox_inches)
    if close:
        plt.close(fig)
    return fig, axes


#17-1. Transition path analysis
def assign_by_membership(M):
    return np.argmax(M, axis=1)


#17-2. Transition path analysis
def clusters_by_membership(M):
    a = np.argmax(M, axis=1)
    m = np.shape(M)[1]
    res = []
    for i in range(m):
        res.append(np.where(a == i)[0])
    return res

#17-3. Transition path analysis
def transition_path_analysis(dataset, tics, k, seed, lag, nstates, start, final, out_path=None):
    """
    out_path: str or None
        None: print
    """
    def _run():
        columns = [f"tIC{i+1}" for i in range(tics + 1)]
        tica_output = load_columns_from_csv(dataset, columns)
        cluster = pyemma.coordinates.cluster_kmeans(
            tica_output, k=k, max_iter=50, stride=10, fixed_seed=seed
        )
        M = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=lag)
        M.pcca(nstates)
        M2 = M.metastable_memberships
        pcca_clusters = clusters_by_membership(M2)

        print("PCCA clusters", pcca_clusters)

        tpt2 = pyemma.msm.tpt(M, start, final)

        (tpt_sets, tpt2_coarse) = tpt2.coarse_grain(pcca_clusters)
        print(tpt_sets)
        print(tpt2_coarse.A)
        print(tpt2_coarse.B)
        print("pi coarse: ")
        print(tpt2_coarse.stationary_distribution)
        print("sum = ", np.sum(tpt2_coarse.stationary_distribution))
        print("committors : ")
        print(tpt2_coarse.committor)
        print(tpt2_coarse.backward_committor)
        print("F coarse : ")
        print(tpt2_coarse.gross_flux)
        print("F net coarse : ")
        print(tpt2_coarse.net_flux)

        (paths, pathfluxes) = tpt2_coarse.pathways()
        cumflux = 0
        total_flux = float(np.sum(pathfluxes))  # sum(pathfluxes) でもOKだが明示
        print("Path flux\t\t%path\t%of total\tpath")
        for i in range(len(paths)):
            cumflux += pathfluxes[i]
            print(
                pathfluxes[i], "\t",
                "%3.1f" % (100.0 * pathfluxes[i] / total_flux), "%\t",
                "%3.1f" % (100.0 * cumflux / total_flux), "%\t\t",
                paths[i]
            )

    if out_path is None:
        _run()
    else:
        with open(out_path, "w", encoding="utf-8") as f:
            with redirect_stdout(f):
                _run()

