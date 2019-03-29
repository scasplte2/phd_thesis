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

%% LvNE: Balancing transformation
qm_setup();
qm_init();
qm_balance('lvne'); 
qm_cleanup();

%% LvNE: Perform truncation
dim = 170;
qm_setup();
qm_init;
qm_truncate('lvne','t',dim);
qm_cleanup;

%% LvNE: Optimal control theory in reduced dimensionality
qm_setup();
qm_init();
qm_optimal('lvne','t',dim);
qm_cleanup();