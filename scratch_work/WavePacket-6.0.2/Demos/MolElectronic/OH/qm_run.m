%% Bare state representation
cd 0-1/BareState;   
qm_setup(); 
state=wave;
qm_init(state); 
qm_propa(state); 
qm_cleanup();   
cd ../..;

%% Floquet dressed state representation
cd 0-1/Dressed02;   
qm_setup(); 
state=wave;
qm_init(state); 
qm_propa(state); 
qm_cleanup();   
cd ../..;
