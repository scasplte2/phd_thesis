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

%% LvNE: Get ABNCD matrices
qm_setup(); 
qm_init; 
qm_abncd('lvne');
qm_cleanup();

%% LvNE: H2 optimal model reduction
dim=170;
qm_setup();
qm_init();
qm_H2model('lvne',dim); 
qm_cleanup();

%% LvNE: Optimal control theory in reduced dimensionality
qm_setup();
qm_init();
qm_optimal('lvne','h',dim);
qm_cleanup();