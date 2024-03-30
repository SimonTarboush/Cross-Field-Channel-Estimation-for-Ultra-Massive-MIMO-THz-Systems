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
% This function generates the Tx/Rx dictionaries \bar{A}_T/\bar{A}_R which
% will be used for far field (PWM) and intermediate field (HSPWM) channel estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% Quantized_Phase_Shifts: The actual quantized phase shift values
% Qbar:    The number of the ULA antennas
% G:          Grid size to quantize the angular region
% useQPH: Choose to use the quantized phase shift values or not; 'True' or 'False'
% Output Arguments:
% A_bar: The dictionary \bar{A}_T/\bar{A}_R used in compressed sensing (Eq. 24 in the manuscript) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A_bar = get_DictionaryULA(Quantized_Phase_Shifts, Qbar, G, useQPH)
OVS = ceil(G/Qbar);
OVSQbar = OVS*Qbar;
row = (-(Qbar-1)/2:(Qbar-1)/2)' ;
col = -1+2/OVSQbar:2/OVSQbar:1;
A_ULA  =  exp(1j*pi*row*col) / sqrt(Qbar);
% Alternative method to generate DFT dictionary without oversampling
% A_tmp = dftmtx(Qbar)/sqrt(Qbar); 
if strcmp(useQPH,'True')
    % Obtain the closest phase in the quantized phase set due to finite resolution phase-shifters following Eq. 24 in the manuscript
    Ang_A_tmp = wrapTo2Pi(angle(A_ULA));
    Modified_Quan_PhaseShifts = [Quantized_Phase_Shifts Quantized_Phase_Shifts(end)+(Quantized_Phase_Shifts(end)-Quantized_Phase_Shifts(end-1))];
    % Quantizing based on phase shifters' resolution
    QuantAng_A_tmp = interp1(Modified_Quan_PhaseShifts, Modified_Quan_PhaseShifts, Ang_A_tmp, 'nearest');
    % Define the Dictionary matrix
    A_bar = 1/sqrt(Qbar)*exp(1j*QuantAng_A_tmp);
else % useQPH = 'False'
    A_bar = A_ULA;
end
end

