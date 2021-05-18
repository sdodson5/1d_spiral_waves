# 1d_spiral_waves
Continue 1D spiral waves and heteroclinic bifurcation point as described in Dodson &amp; Lewis, 2021.


**Main Driver Files:**
1. **continue_1d_spiral.m:** Continues an initial condition/guess for the 1D spiral in any system parameter.
2. **continue_1d_spiral_unstable_eigenvalue.m:** Continues an initial condition/guess for the 1D spiral in any system parameter and also computes the corresponding usntable eigenvalue-eigenvector pair
3. **continue_heteroclinic_bifur_pt.m:** Continues an initial condition/guess for the heteroclinic bifucation point at which the 1D spiral appears. 



**Initial Conditions:**
1. **MD_oneDspiral_soln.mat:** Initial condition for 1D spiral continuation
2. **ML_1dspiral_efcn_soln.mat:** Initial condition for 1D spiral + unstable eigenvalue continuation
3. **ML_SsF_soln.mat:** Initial condition for heteroclinic bifurcation point continuation


Folder **Util** contains auxiliary functions called by the driver files. 








