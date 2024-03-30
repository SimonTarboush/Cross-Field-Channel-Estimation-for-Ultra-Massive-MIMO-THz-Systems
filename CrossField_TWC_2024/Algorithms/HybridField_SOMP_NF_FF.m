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
% This function implements the Hybrid-field SOMP algorithm by first estimating the near-field (NF) components and substract their effect 
% after that start estimating the far-field (FF) components
% Since we  do not assume any prior knowledge of the number of far-field and near-field components,
% we assume the mixing ratio of the fields is half i.e: Lbar_Far = 1/2*Lbar and Lbar_Near = 1/2*Lbar
% This implementation is inspried from [1] where the authors define the OMP CS-estimator
% [1] X. Wei and L. Dai, “Channel estimation for extremely large-scale massive MIMO: Far-field, near-field, or hybrid-field?” IEEE Commun. Lett., vol. 26, no. 1, pp. 177-181, 2021.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% y: Measurements vector
% Upsilon_Far: The sensing matrix based on far-field dictionaries
% Upsilon_Near: The sensing matrix based on near-field dictionaries
% Pilot_pwr: The transmission power in the training phase
% Lbar: Estimation of the sparsity level
% K: Number of subcarriers
% Output Arguments:
% H_Est_Far: The far-field components estimation of the channel in beamspace
% Support_Far: The index selected by the SOMP algorithm (known as the support) for the far-field components
% H_Est_Near: The near-field components estimation of the channel in beamspace
% Support_Near: The index selected by the SOMP algorithm (known as the support) for the near-field components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H_est_Far, Support_Far, H_est_Near, Support_Near] = HybridField_SOMP_NF_FF(y, Upsilon_Far, Upsilon_Near, Pilot_pwr, Lbar, K)
R = y;
G_Far = size(Upsilon_Far,2);
G_Near = size(Upsilon_Near,2);
Index_Far = [];
Index_Near = [];
norm_Upsilon_Far = sqrt(sum(abs(Upsilon_Far).^2, 1)).';
norm_Upsilon_Near = sqrt(sum(abs(Upsilon_Near).^2, 1)).';
for ell =1: Lbar/2 % threshold = sparsity level
    % Find AoD/AoA pair
    [~, indx_c_Near] = max(sum(abs(Upsilon_Near'*R./norm_Upsilon_Near).^2,2)); 
    % Update AoD/AoA pair index
    Index_Near = [Index_Near indx_c_Near];
    H_est_Near = zeros(G_Near,K);
    % Estimate channel gains by solving the Least Square problem
    H_est_Near(Index_Near,:) = 1/sqrt(Pilot_pwr)*pinv(Upsilon_Near(:,Index_Near)'*Upsilon_Near(:,Index_Near))*Upsilon_Near(:,Index_Near)'*y;
    % Update residual
    R = y - sqrt(Pilot_pwr)*Upsilon_Near(:,Index_Near)*H_est_Near(Index_Near,:);
end
Support_Near = Index_Near;
for ell =1: Lbar/2 % threshold = sparsity level
    % Find AoD/AoA pair
    [~, indx_c_Far] = max(sum(abs(Upsilon_Far'*R./norm_Upsilon_Far).^2,2));
    % Update AoD/AoA pair index
    Index_Far = [Index_Far indx_c_Far]; 
    H_est_Far = zeros(G_Far,K);
    % Estimate channel gains by solving the Least Square problem
    H_est_Far(Index_Far,:) = 1/sqrt(Pilot_pwr)*pinv(Upsilon_Far(:,Index_Far)'*Upsilon_Far(:,Index_Far))*Upsilon_Far(:,Index_Far)'*y;
    % Update residual
    if ~isempty(Support_Near)
        R = y - sqrt(Pilot_pwr)*Upsilon_Near(:,Index_Near)*H_est_Near(Index_Near,:)- sqrt(Pilot_pwr)*Upsilon_Far(:,Index_Far)*H_est_Far(Index_Far,:);
    else
        R = y - sqrt(Pilot_pwr)*Upsilon_Far(:,Index_Far)*H_est_Far(Index_Far,:);
    end
end
Support_Far = Index_Far;
end
