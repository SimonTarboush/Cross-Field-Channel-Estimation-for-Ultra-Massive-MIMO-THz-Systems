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
% This function construct the quantized phase Dictionary based on an unquantized phase Dictionary 
% by obtaining the closest phase in the quantized phase set due to finite resolution phase-shifters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% Dict: a near or far field dictionary
% Quantized_Phase_Shifts: The actual quantized phase shift values
% Output Arguments:
% Dict_quant_phase: a near or far field dictionary following quantized phase constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Dict_quant_phase = get_Quantized_Dict(Dict, Quantized_Phase_Shifts)
Dict_quant_phase = complex(zeros(size(Dict)));
Q = size(Dict_quant_phase,1);
Dict_quant_col = size(Dict_quant_phase,2);
for indx_c = 1:Dict_quant_col
    Dict_col = Dict(:,indx_c);
    % Obtain the closest phase in the quantized phase set due to finite resolution phase-shifters
    Ang_Dict_tmp = wrapTo2Pi(angle(Dict_col));
    Modified_Quan_PhaseShifts = [Quantized_Phase_Shifts Quantized_Phase_Shifts(end)+(Quantized_Phase_Shifts(end)-Quantized_Phase_Shifts(end-1))];
    % Quantizing based on phase shifters' resolution
    QunatAng_Dict_tmp = interp1(Modified_Quan_PhaseShifts, Modified_Quan_PhaseShifts, Ang_Dict_tmp, 'nearest');
    % Define the Dictionary matrix
    Dict_quant_phase(:,indx_c) = 1/sqrt(Q)*exp(1j*QunatAng_Dict_tmp);
end