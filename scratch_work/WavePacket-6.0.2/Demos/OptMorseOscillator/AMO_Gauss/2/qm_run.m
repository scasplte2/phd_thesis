%% Solve the TDSE in PDE setting
qm_setup; 
state=wave;
qm_init;
qm_bound(state);    
qm_propa(state);
qm_cleanup;


