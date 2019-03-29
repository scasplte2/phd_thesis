%% Short-time simulation
cd 1;   
qm_setup(); 
state=wave;
qm_init(state); 
qm_propa(state); 
qm_cleanup();   
cd ..;

%% Long-time simulation
cd 2;   
qm_setup(); 
state=wave;
qm_init(state); 
qm_propa(state); 
qm_cleanup();   
cd ..;
