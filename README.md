# control-with-foresight
Code implementing local linearisation-based control methods for nonlinear systems. Developed in MATLAB R2020a.

__Associated with the paper__\
"Foresight and relaxation enable efficient control of nonlinear complex systems"\
by\
Xin Zhang, Xiaozhu Zhang, Gang Yan,∗ and Jack Murdoch Moore.

__This code__ allows reproduction of Fig. 2-3 of the paper via the following steps (if you like, go straight to Step 2).
1. Generate data.\
   Run the script _gen_data_comparison_.
    * This will produce the .MAT files of data _comparison-1e-06,0.001,0.1,1_bist_cont-3,5_strat-43,43,43,43_energy-1,1,1,1.mat_ and _comparison-0.0001,0.0001,0.0001_bist_cont-3,4_strat-45,46,43_energy-100000,1,1.mat_.
2. Make plots from data.\
   Run _plot_comparison_.
    * This will use the .MAT files to produce the image files _comparison-1e-06,0.001,0.1,1_bist_cont-3,5_strat-43,43,43,43.png_ (Fig. 3 of paper) and _comparison-0.0001,0.0001,0.0001_bist_cont-3,4_strat-45,46,43.png_ (Fig. 2 of paper).

__Functions and scripts:__
* _gen_data_comparison.m:_ Script which defines dynamical system, finds its fixed points and runs function _control_func.m_ to implement control.
* _control_func.m:_ Function which implements different control strategies.
* _plot_comparison.m:_ Script which uses .MAT data files to produce .PNG image files (Fig. 2-3 of paper).

__Data and images:__
* _comparison-0.0001,0.0001,0.0001_bist_cont-3,4_strat-45,46,43_energy-100000,1,1.mat:_ Data, produced by script _gen_data_comparison.m_ via function _control_func.m_, which allows comparison of different control strategies.
* _comparison-1e-06,0.001,0.1,1_bist_cont-3,5_strat-43,43,43,43_energy-1,1,1,1.mat:_ Data, produced by script _gen_data_comparison.m_ via function _control_func.m_, which allows comparison of different values of the relaxation parameter r.
* _comparison-0.0001,0.0001,0.0001_bist_cont-3,4_strat-45,46,43.png:_ Image, produced by script _plot_strategy_comparison.m_ using data _comparison-0.0001,0.0001,0.0001_bist_cont-3,4_strat-45,46,43_energy-100000,1,1.mat_, which juxtaposes outcomes of different control strategies.
* _comparison-1e-06,0.001,0.1,1_bist_cont-3,5_strat-43,43,43,43.png:_ Image, produced by script _plot_relaxation_comparison.m_ using data _comparison-1e-06,0.001,0.1,1_bist_cont-3,5_strat-43,43,43,43_energy-1,1,1,1.mat_, which juxtaposes outcomes of different values of the relaxation parameter r.
