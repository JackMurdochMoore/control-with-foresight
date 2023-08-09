# alite
Code implementing local linearisation-based control methods for nonlinear systems. Developed in MATLAB R2020a.

__Associated with the paper__\
"Foresight and relaxation enable efficient control of nonlinear complex systems"\
by\
Xin Zhang, Xiaozhu Zhang, Gang Yan,∗ and Jack Murdoch Moore.

__About__\
This code should allow reproduction of Fig. 2-3 of the paper via the following steps.
1. To generate the data, run the script _gen_data_comparison_ with "comparisonType = 1" and "comparisonType = 2" in lines 11-12 of the script.
    * This will produce the .MAT files of data "comparison-1e-06,0.001,0.1,1_bist_cont-3,5_strat-43,43,43,43_energy-1,1,1,1.mat" and "comparison-0.0001,0.0001,0.0001_bist_cont-3,4_strat-45,46,43_energy-100000,1,1.mat".
2. Run _plot_relaxation_comparison_ and _plot_strategy_comparison_ to use the .MAT files to produce the image files _comparison-1e-06,0.001,0.1,1_bist_cont-3,5_strat-43,43,43,43.png_ (Fig. 3 of paper) and _comparison-0.0001,0.0001,0.0001_bist_cont-3,4_strat-45,46,43.png_ (Fig. 2 of paper).

__Functions and scripts:__\
_gen_data_comparison.m:_ Script which defines dynamical system, finds its fixed points and runs function _control_func.m_ to impl;ement control.\
_control_func.m:_ Function which implements different control strategies.\
_plot_strategy_comparison.m:_ Script which uses data _comparison-0.0001,0.0001,0.0001_bist_cont-3,4_strat-45,46,43_energy-100000,1,1.mat_ to produce image _comparison-1e-06,0.001,0.1,1_bist_cont-3,5_strat-43,43,43,43.png_ (Fig. 3 of paper).\
_plot_relaxation_comparison.m:_ Script which uses data _comparison-0.0001,0.0001,0.0001_bist_cont-3,4_strat-45,46,43_energy-100000,1,1.mat_ to produce image _comparison-0.0001,0.0001,0.0001_bist_cont-3,4_strat-45,46,43_energy-100000,1,1.mat_ (Fig. 2 of paper).

__Data and images:__\
_comparison-0.0001,0.0001,0.0001_bist_cont-3,4_strat-45,46,43_energy-100000,1,1.mat:_ Data, produced by script _gen_data_comparison.m_ via function _control_func.m_, which allows comparison of different control strategies.\
_comparison-1e-06,0.001,0.1,1_bist_cont-3,5_strat-43,43,43,43_energy-1,1,1,1.mat:_ Data, produced by script _gen_data_comparison.m_ via function _control_func.m_, which allows comparison of different values of the relaxation parameter r.\
_est_corr_dim_4.m:_ Estimate correlation dimension and scaling interval of a networks using different methods and model C(s) = s^D.\
_example.m:_ An example to illustrate generation of a synthetic network and estimation of its correlation dimension.\
_find_local_minima.m:_ Find local minima in a vector.\
_load_network.m:_ Load an empirical network from data in folder "networks".\
_log_like_3.m:_ Calculate log-likelihood per observation for model c(s) = s^(D - 1).\
_log_like_4.m:_ Calculate log-likelihood per observation for model C(s) = s^D.\
_plot_corr_dim_estimates.m:_ Plot results of benchmarking, previously saved in folders "results-est-dim-3" and "results-est-dim-4".\
_show_fit_func.m:_ Fit power-law to an interval and illustrate fit and objective function (negative log-likelihood per observation).\
_show_fit_script.m:_ Illustrate fit and fitting process for empirical or synthetic networks by calling function show_fit_func.m.\
_small_world_manhattan.m:_ Generate lattice* or small world network*.\
_small_world_manhattan_lcc.m:_ Generate lattice* or small world network* and retain only its largest connected component.

\* Lattices/small world networks are/are derived from regular $D$-dimensional toroidal lattices defined using a periodic version of the city block (or Manhattan or taxi cab) metric mentioned but not explored in "Epidemic dynamics on higher-dimensional small world networks", _Applied Mathematics and Computation_
421, 126911, by H. Wang, J. M. Moore, M. Small, J. Wang, H. Yang and C. Gu (2022) (associated code at https://github.com/JackMurdochMoore/small-world).
