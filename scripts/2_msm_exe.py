from msm import load_columns_from_csv
from msm import cluster_centers
from msm import save_plot_data_csv
from msm import plot_from_csv
from msm import calc_its
from msm import plot_distrib_energy
from msm import plot_free_energy_gate_all
from msm import plot_free_energy_zna_inward_all
from msm import plot_free_energy_na_bound
from msm import plot_free_energy_na_unbound
from msm import plot_free_energy_with_pcca_states
from msm import plot_ck_test
from msm import transition_path_analysis

dataset = [
    "../MSM_FESred/round0/7XK3/dataset.csv",
    "../MSM_FESred/round0/8A1U/dataset.csv",
    "../MSM_FESred/round0/9UD2/dataset.csv",
    "../MSM_FESred/round0/9UDF/dataset.csv",
    "../MSM_FESred/round0/AF2/dataset.csv",
    "../MSM_FESred/round1/cluster09_1/dataset.csv", 
    "../MSM_FESred/round1/cluster09_2/dataset.csv", 
    "../MSM_FESred/round1/cluster09_3/dataset.csv", 
    "../MSM_FESred/round1/cluster09_4/dataset.csv", 
    "../MSM_FESred/round1/cluster09_5/dataset.csv", 
    "../MSM_FESred/round1/cluster12_1/dataset.csv", 
    "../MSM_FESred/round1/cluster12_2/dataset.csv", 
    "../MSM_FESred/round1/cluster12_3/dataset.csv", 
    "../MSM_FESred/round1/cluster12_4/dataset.csv", 
    "../MSM_FESred/round1/cluster12_5/dataset.csv", 
    "../MSM_FESred/round1/cluster17_1/dataset.csv", 
    "../MSM_FESred/round1/cluster17_2/dataset.csv", 
    "../MSM_FESred/round1/cluster17_3/dataset.csv", 
    "../MSM_FESred/round1/cluster17_4/dataset.csv", 
    "../MSM_FESred/round1/cluster17_5/dataset.csv", 
    "../MSM_FESred/round1/cluster21_1/dataset.csv", 
    "../MSM_FESred/round1/cluster21_2/dataset.csv", 
    "../MSM_FESred/round1/cluster21_3/dataset.csv", 
    "../MSM_FESred/round1/cluster21_4/dataset.csv", 
    "../MSM_FESred/round1/cluster21_5/dataset.csv", 
    "../MSM_FESred/round1/cluster08_1/dataset.csv", 
    "../MSM_FESred/round1/cluster08_2/dataset.csv", 
    "../MSM_FESred/round1/cluster08_3/dataset.csv", 
    "../MSM_FESred/round1/cluster08_4/dataset.csv", 
    "../MSM_FESred/round1/cluster08_5/dataset.csv", 
    "../MSM_FESred/round2/cluster12_1/dataset.csv", 
    "../MSM_FESred/round2/cluster12_2/dataset.csv", 
    "../MSM_FESred/round2/cluster12_3/dataset.csv", 
    "../MSM_FESred/round2/cluster12_4/dataset.csv", 
    "../MSM_FESred/round2/cluster12_5/dataset.csv", 
    "../MSM_FESred/round2/cluster22_1/dataset.csv", 
    "../MSM_FESred/round2/cluster22_2/dataset.csv", 
    "../MSM_FESred/round2/cluster22_3/dataset.csv", 
    "../MSM_FESred/round2/cluster22_4/dataset.csv", 
    "../MSM_FESred/round2/cluster22_5/dataset.csv", 
    "../MSM_FESred/round2/cluster23_1/dataset.csv", 
    "../MSM_FESred/round2/cluster23_2/dataset.csv", 
    "../MSM_FESred/round2/cluster23_3/dataset.csv", 
    "../MSM_FESred/round2/cluster23_4/dataset.csv", 
    "../MSM_FESred/round2/cluster23_5/dataset.csv", 
    "../MSM_FESred/round2/cluster26_1/dataset.csv", 
    "../MSM_FESred/round2/cluster26_2/dataset.csv", 
    "../MSM_FESred/round2/cluster26_3/dataset.csv", 
    "../MSM_FESred/round2/cluster26_4/dataset.csv", 
    "../MSM_FESred/round2/cluster26_5/dataset.csv", 
    "../MSM_FESred/round2/cluster28_1/dataset.csv", 
    "../MSM_FESred/round2/cluster28_2/dataset.csv", 
    "../MSM_FESred/round2/cluster28_3/dataset.csv", 
    "../MSM_FESred/round2/cluster28_4/dataset.csv", 
    "../MSM_FESred/round2/cluster28_5/dataset.csv", 
    "../MSM_FESred/round2/cluster52_1/dataset.csv", 
    "../MSM_FESred/round2/cluster52_2/dataset.csv", 
    "../MSM_FESred/round2/cluster52_3/dataset.csv", 
    "../MSM_FESred/round2/cluster52_4/dataset.csv", 
    "../MSM_FESred/round2/cluster52_5/dataset.csv", 
    "../MSM_FESred/round2/cluster53_1/dataset.csv", 
    "../MSM_FESred/round2/cluster53_2/dataset.csv", 
    "../MSM_FESred/round2/cluster53_3/dataset.csv", 
    "../MSM_FESred/round2/cluster53_4/dataset.csv", 
    "../MSM_FESred/round2/cluster53_5/dataset.csv", 
    "../MSM_FESred/round2/cluster56_1/dataset.csv", 
    "../MSM_FESred/round2/cluster56_2/dataset.csv", 
    "../MSM_FESred/round2/cluster56_3/dataset.csv", 
    "../MSM_FESred/round2/cluster56_4/dataset.csv", 
    "../MSM_FESred/round2/cluster56_5/dataset.csv", 
    "../MSM_FESred/round2/cluster57_1/dataset.csv", 
    "../MSM_FESred/round2/cluster57_2/dataset.csv", 
    "../MSM_FESred/round2/cluster57_3/dataset.csv", 
    "../MSM_FESred/round2/cluster57_4/dataset.csv", 
    "../MSM_FESred/round2/cluster57_5/dataset.csv", 
    "../MSM_FESred/round2/cluster58_1/dataset.csv", 
    "../MSM_FESred/round2/cluster58_2/dataset.csv", 
    "../MSM_FESred/round2/cluster58_3/dataset.csv", 
    "../MSM_FESred/round2/cluster58_4/dataset.csv", 
    "../MSM_FESred/round2/cluster58_5/dataset.csv", 
    "../MSM_FESred/AF3/AF3_806_v1/dataset.csv", 
    "../MSM_FESred/AF3/AF3_806_v2/dataset.csv", 
    "../MSM_FESred/AF3/AF3_806_v3/dataset.csv", 
    "../MSM_FESred/AF3/AF3_806_v4/dataset.csv", 
    "../MSM_FESred/AF3/AF3_806_v5/dataset.csv", 
    "../MSM_FESred/AF3/AF3_832_v1/dataset.csv", 
    "../MSM_FESred/AF3/AF3_832_v2/dataset.csv", 
    "../MSM_FESred/AF3/AF3_832_v3/dataset.csv", 
    "../MSM_FESred/AF3/AF3_832_v4/dataset.csv", 
    "../MSM_FESred/AF3/AF3_832_v5/dataset.csv", 
    "../MSM_FESred/AF3/AF3_885_v1/dataset.csv", 
    "../MSM_FESred/AF3/AF3_885_v2/dataset.csv", 
    "../MSM_FESred/AF3/AF3_885_v3/dataset.csv", 
    "../MSM_FESred/AF3/AF3_885_v4/dataset.csv", 
    "../MSM_FESred/AF3/AF3_885_v5/dataset.csv", 
    "../MSM_FESred/AF3/AF3_898_v1/dataset.csv", 
    "../MSM_FESred/AF3/AF3_898_v2/dataset.csv", 
    "../MSM_FESred/AF3/AF3_898_v3/dataset.csv", 
    "../MSM_FESred/AF3/AF3_898_v4/dataset.csv", 
    "../MSM_FESred/AF3/AF3_898_v5/dataset.csv", 
    "../MSM_FESred/round3/cluster005_1/dataset.csv", 
    "../MSM_FESred/round3/cluster005_2/dataset.csv", 
    "../MSM_FESred/round3/cluster005_3/dataset.csv", 
    "../MSM_FESred/round3/cluster005_4/dataset.csv", 
    "../MSM_FESred/round3/cluster005_5/dataset.csv", 
    "../MSM_FESred/round3/cluster064_1/dataset.csv", 
    "../MSM_FESred/round3/cluster064_2/dataset.csv", 
    "../MSM_FESred/round3/cluster064_3/dataset.csv", 
    "../MSM_FESred/round3/cluster064_4/dataset.csv", 
    "../MSM_FESred/round3/cluster064_5/dataset.csv", 
    "../MSM_FESred/round3/cluster096_1/dataset.csv", 
    "../MSM_FESred/round3/cluster096_2/dataset.csv", 
    "../MSM_FESred/round3/cluster096_3/dataset.csv", 
    "../MSM_FESred/round3/cluster096_4/dataset.csv", 
    "../MSM_FESred/round3/cluster096_5/dataset.csv", 
    "../MSM_FESred/round3/cluster099_1/dataset.csv", 
    "../MSM_FESred/round3/cluster099_2/dataset.csv", 
    "../MSM_FESred/round3/cluster099_3/dataset.csv", 
    "../MSM_FESred/round3/cluster099_4/dataset.csv", 
    "../MSM_FESred/round3/cluster099_5/dataset.csv", 
    "../MSM_FESred/round3/cluster166_1/dataset.csv", 
    "../MSM_FESred/round3/cluster166_2/dataset.csv", 
    "../MSM_FESred/round3/cluster166_3/dataset.csv", 
    "../MSM_FESred/round3/cluster166_4/dataset.csv", 
    "../MSM_FESred/round3/cluster166_5/dataset.csv", 
]


#6. Calculate VAMP2 score
n_clustercenters = [10,50,100,200,400,600]
tic_ranges = [2,3,4,5,6,7,8,9,10,20,30,43] 
results_dict = {}
for n_tic in tic_ranges:
    columns = [f"tIC{i+1}" for i in range(n_tic)]
    tica_output = load_columns_from_csv(dataset, columns)
    scores, lower, upper = cluster_centers(n_clustercenters, tica_output)
    results_dict[f"tIC1~{n_tic}"] = (scores, lower, upper)


#7. Save VAMP2 score
csv_vamp2 = "../MSM_FESred/msm_data/vamp2_cluster_num.csv"
save_plot_data_csv(n_clustercenters, results_dict, csv_vamp2)


#8. Plot VAMP2 score
fig_vamp2 = "../MSM_FESred/msm_data/vamp2_cluster_num.png"
plot_from_csv(csv_vamp2, fig_vamp2)


#9. Calculate implied timescale (its)
tics = 20
k = 100
seed = 1
lags = [1,2,4,10,20,30,40,50,60,70,80,90,100]  # steps
fig_its = "../MSM_FESred/msm_data/its_k100.png"
calc_its(dataset, tics, k, seed, lags, fig_its)


#10. Plot stationary distribution & free energy
lag = 40
fig_distrib_energy = "../MSM_FESred/msm_data/distrib_energy_k100_lag4ns.png"
plot_distrib_energy(dataset, tics, k, seed, lag, fig_distrib_energy)


#11. Plot free energy along with gate sizes in NqrD/E
kT = 1
fig_free_energy_gate_all = "../MSM_FESred/msm_data/free_energy_k100_lag4ns_gate_all.png"
plot_free_energy_gate_all(
    dataset, tics, k, seed, lag, fig_free_energy_gate_all, kT = kT, 
    nbins=100, display_max=12.0, show_ticks=(0, 2, 4, 6, 8, 10, 12))


#12. Plot free energy along with Na+ z-coords and inward gate size
fig_free_energy_zna_inward_all = "../MSM_FESred/msm_data/free_energy_k100_lag4ns_gate_all.png"
plot_free_energy_zna_inward_all(, fig_free_energy_zna_inward_all, kT=kT)


#13. Plot free energy using Na+-bound data
fig_free_energy_gate_na_bound = "../MSM_FESred/msm_data/free_energy_k100_lag4ns_gate_na_bound.png"
plot_free_energy_na_bound(dataset, tics, k, seed, lag, fig_free_energy_gate_na_bound, kT=kT)


#14. Plot free energy using Na+-unbound data
fig_free_energy_gate_na_unbound = "../MSM_FESred/msm_data/free_energy_k100_lag4ns_gate_na_unbound.png"
plot_free_energy_na_unbound(dataset, tics, k, seed, lag, fig_free_energy_gate_na_unbound, kT=kT)


#15. Plot coarse-grained states (PCCA+ result)
nstates_ck = 3
fig_CG = "../MSM_FESred/msm_data/3CG_states.png"
plot_free_energy_with_pcca_states(dataset, tics, k, seed, lag, nstates_ck, fig_CG, kT=kT)


#16. Chapman-Kolmogorov test (ck_test)
fig_ck_test = "../MSM_FESred/msm_data/ck_test.png"
plot_ck_test(
    dataset, tics, k, seed, lag, nstates_ck, fig_ck_test, 
    diag=False, y01=True)


#17-1,2,3. Transition path analysis
nstates_tpt = 30
start = [59,67]
final = [30,36]
out_tpt = "../MSM_FESred/msm_data/tpt.dat"
transition_path_analysis(dataset, tics, k, seed, lag, nstates_tpt, start, final, out_tpt)
