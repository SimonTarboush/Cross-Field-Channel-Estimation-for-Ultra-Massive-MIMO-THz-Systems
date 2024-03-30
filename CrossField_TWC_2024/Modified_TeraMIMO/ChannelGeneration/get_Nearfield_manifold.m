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
% This function generates a near-field array response vector, mainly implementing Eq. (9) in the manuscript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% p: Channel struct that contains the channel parameters.
% Qbar: number of AEs inside a SA (ULA-based)
% d: the distance of the ell-th scatterer from the qbar-th Tx/Rx AE
% d0: communication distance of the  ell-th path (Rx for ell = 1 and scatterer for ell > 1) from the Tx center
% lambda: wavelength
% Output Arguments:
% b: near-field array response vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  b= get_Nearfield_manifold(Qbar, d, d0, lambda)
    b = exp(-1j*2*pi*(d-d0)/lambda)/sqrt(Qbar);
end