% Script to run all demos in all subdirectories. Usually, you do not need this.
% The idea is that, together with the prerun and postrun scripts, we can
% automatically run all the demos, produce all the output files, in particular
% movies, images, etc., and end up with a directory tree we can directly use for
% the our wiki webpages maintained at sourceforge.net. Use this to upload
%
% scp -r * username,wavepacket@web.sourceforge.net:htdocs/Demos/ 
%
% Should you want, for some reason, to run a subset of the demos, you can find
% similar scripts in each subdirectory.

cd ChemReaction;        qm_run;    cd ..;
cd ConicalInter;        qm_run;    cd ..;
cd CrossingBerlin;      qm_run;    cd ..;
cd CrossingTully;       qm_run;    cd ..;
cd DoubleWell;          qm_run;    cd ..;
cd FemtoChem;           qm_run;    cd ..;
cd FreeParticle;        qm_run;    cd ..;
cd HarmOscillator;      qm_run;    cd ..;
cd HenonHeiles;         qm_run;    cd ..;
cd LinearRamp;          qm_run;    cd ..;
cd MolElectronic;       qm_run;    cd ..;
cd MolRotation;         qm_run;    cd ..;
cd MolTorsion;          qm_run;    cd ..;
cd MolVibration;        qm_run;    cd ..;
cd MorseOscillator;     qm_run;    cd ..;
cd PlanePendulum;       qm_run;    cd ..;
cd SpinBoson;           qm_run;    cd ..;
