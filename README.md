Author: Simon Tarboush                                      
Last Modified: March, 2024

If you use this code or any (modified) part of it in any publication, please cite the paper: 
Simon Tarboush, Anum Ali, Tareq Y. Al-Naffouri, 
"Cross-Field Channel Estimation for Ultra Massive-MIMO THz Systems", IEEE Transactions on Wireless Communications.
(https://ieeexplore.ieee.org/document/10410228)

You may also refer to the conference version of this work, which specifically delves into 
Hybrid Spherical Planar Wave Model (HSPWM) channel estimation and introduces the associated
reduced dictionary technique.
Simon Tarboush, Anum Ali, Tareq Y. Al-Naffouri, 
"Compressive Estimation of Near Field Channels for Ultra Massive-MIMO Wideband THz Systems", 
ICASSP 2023 - 2023 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP).
(Codes are also available on IEEExplore https://ieeexplore.ieee.org/document/10096832 
and GitHub https://github.com/SimonTarboush/Compressive-Estimation-of-Near-Field-Channels-for-Ultra-Massive-Mimo-Wideband-THz-Systems)

If you use the channel simulator code "TeraMIMO" or any (modified) part of it in any publication, please cite 
the paper: Simon Tarboush, Hadi Sarieddeen, Hui Chen, Mohamed Habib Loukil, Hakim Jemaa, Mohamed-Slim Alouini, Tareq Y. Al-Naffouri
"TeraMIMO: A Channel Simulator for Wideband Ultra-Massive MIMO Terahertz Communications",
IEEE Transactions on Vehicular Technology.

Contact person email: simon.w.tarboush@gmail.com
You can also watch a Youtube video where we explained the work: https://www.youtube.com/watch?v=sxT4gYmqbaI&t=1s
# Brief Description:

Fig7_OfflineTraining.m : the script to obtain results of Fig. 7
Fig9_Fig8_Proposed_CrossField_ReducedDictionary.m : the script to obtain results of Fig. 9 or 8

All other functions are called by the these two scripts.

The main file is "run.sh" which calls only one script (because of simulation time)
Please, try first to run Fig9_Fig8_Proposed_CrossField_ReducedDictionary alone with a small
number of iterations and rotations, since this function requires heavy computation resources.
# Notes: 

This set of code was tested on MATLAB_R2021_a
