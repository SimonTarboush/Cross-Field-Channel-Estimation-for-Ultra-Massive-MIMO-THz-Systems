%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Simon Tarboush                                      
% Last Modified: March, 2024
%
% The original code is from [1] and the authors' public available MATLAB codes:
% [1] M. Cui and L. Dai, “Channel estimation for extremely large-scale MIMO: Far-field or near-field?” IEEE Trans. Commun., vol. 70, no. 4, pp. 2663–2677, 2022.
% However, I made some necessary modifications to make the manifold compatible with our system model and assumptions
% I added the quantization option also for achievable rate calculations
%
% Contact person email: simon.w.tarboush@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates the near-field array response vector based on Eq. 9 in the manuscript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% Qbar: The number of the ULA antennas
% d_ae:  Antenna elements spacing
% lambda: Wavelength
% dist: Sampled distance for Polar-domain dictionary
% theta: Sampled angle for Polar-domain dictionary
% Quantized_Phase_Shifts: The actual quantized phase shift values
% useQPH: Choose to use the quantized phase shift values or not; 'True' or 'False'
% Output Arguments:
% b: The near-field array response vector (Eq. 9 in the manuscript) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  b = get_Polar_Domain_Manifold(Qbar, d_ae, lambda, dist, theta, Quantized_Phase_Shifts, useQPH)
qbar_indx = -(Qbar-1)/2:1:(Qbar-1)/2;
d_exact = dist-qbar_indx*d_ae*cos(theta)+(qbar_indx.^2*d_ae^2*(1-cos(theta)^2)/(2*dist));
b_ULA = exp(-1j*2*pi*(d_exact-dist)/lambda)/sqrt(Qbar);
if strcmp(useQPH,'True')
    Ang_P_tmp = wrapTo2Pi(angle(b_ULA));
    Modified_Quan_PhaseShifts = [Quantized_Phase_Shifts Quantized_Phase_Shifts(end)+(Quantized_Phase_Shifts(end)-Quantized_Phase_Shifts(end-1))];
    QuantAng_P_tmp = interp1(Modified_Quan_PhaseShifts, Modified_Quan_PhaseShifts, Ang_P_tmp, 'nearest');
    b = 1/sqrt(Qbar)*exp(1j*QuantAng_P_tmp);
else
    b = b_ULA;
end
end