%% Solve TISE
qm_setup;  
state=wave;
state.sav_export = true; % Save wavefunctions etc
qm_init; 
qm_bound(state);            
qm_cleanup;

%% Matrix representations
qm_setup; 
state=wave; 
qm_init; 
qm_matrix (state); 
qm_cleanup;

%% Switch to ABNCD control language
qm_setup;  
qm_init; 
qm_abncd ('tdse');      
qm_cleanup;

%% Optimal control
qm_setup;  
qm_init; 
qm_optimal ('tdse'); 
qm_cleanup;