%% Solve TISE in position/momentum representation 
qm_setup(); 
state=wave; 
state.sav_export = true; % Save wavefunctions etc
qm_init(state); 
qm_bound(state); 
qm_cleanup();  

%% Get matrix elements
qm_setup(); 
state=wave; 
qm_init(state); 
qm_matrix(state); 
qm_cleanup();

%% Transform to energy representation
qm_setup(); 
qm_init(state); 
qm_abncd('tdse'); 
qm_cleanup(); 

%% Solve TDSE in energy representation 
qm_setup(); 
qm_init(state); 
qm_control('tdse'); 
qm_cleanup();  
saveas(7,'uxy.jpg')