%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Simon Tarboush                                      
% Last Modified: March, 2024
%
% If you use this code or any (modified) part of it in any publication, please cite the paper: 
% Simon Tarboush, Anum Ali, Tareq Y. Al-Naffouri, 
% "Cross-Field Channel Estimation for Ultra Massive-MIMO THz Systems", IEEE Transactions on Wireless Communications.
% (https://ieeexplore.ieee.org/document/10410228)
%
% You may also refer to the conference version of this work, which specifically delves into 
% Hybrid Spherical Planar Wave Model (HSPWM) channel estimation and introduces the associated
% reduced dictionary technique.
% Simon Tarboush, Anum Ali, Tareq Y. Al-Naffouri, 
% "Compressive Estimation of Near Field Channels for Ultra Massive-MIMO Wideband THz Systems", 
% ICASSP 2023 - 2023 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP).
% (Codes are also available on IEEExplore https://ieeexplore.ieee.org/document/10096832 
% and GitHub https://github.com/SimonTarboush/Compressive-Estimation-of-Near-Field-Channels-for-Ultra-Massive-Mimo-Wideband-THz-Systems)
%
% If you use the channel simulator code "TeraMIMO" or any (modified) part of it in any publication, please cite 
% the paper: Simon Tarboush, Hadi Sarieddeen, Hui Chen, Mohamed Habib Loukil, Hakim Jemaa, Mohamed-Slim Alouini, Tareq Y. Al-Naffouri
% "TeraMIMO: A Channel Simulator for Wideband Ultra-Massive MIMO Terahertz Communications",
% IEEE Transactions on Vehicular Technology.
%
% Contact person email: simon.w.tarboush@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function implements Algorithm 2 in the manuscript to construct the proposed Reduced Dictionary (RD) method by selecting the most correlated
% columns from an oversampled dictionary. More details in Sec. III-C1 SWM-based channel estimation with RD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% Pbar_OVS: The oversampled Polar-domain dictionary (defined in Step 1 of Algorithm 2)
% Labels_Pbar_OVS: The labels (the angular (linear) and distance (non-linear) samplings points) of the oversampled Polar-domain dictionary
% Pbar_Est: The selected columns of the Polar-domain dictionary based on the first "reference" T-R SA channel estimation
% G_RD: The size of the reduced dictionary
% Output Arguments:
% Pbar_RD: The RD Polar-domain dictionary
% Labels_Pbar_RD: The labels of the RD Polar-domain dictionary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pbar_RD, Labels_Pbar_RD] = get_Proposed_RD_NearField(Pbar_OVS, Labels_Pbar_OVS, Pbar_Est, G_RD)
U_corr = Pbar_OVS'*Pbar_Est;
U_corr_norm = sqrt(diag(U_corr*U_corr'));
[~ , sel_ind] = sort(U_corr_norm,'descend');
Pbar_RD = Pbar_OVS(:,sel_ind(1:G_RD));
Labels_Pbar_RD = Labels_Pbar_OVS(:,sel_ind(1:G_RD));