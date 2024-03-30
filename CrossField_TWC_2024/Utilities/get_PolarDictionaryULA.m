%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Simon Tarboush                                      
% Last Modified: March, 2024
%
% The original code is from [1] and the authors' public available MATLAB codes:
% [1] M. Cui and L. Dai, “Channel estimation for extremely large-scale MIMO: Far-field or near-field?” IEEE Trans. Commun., vol. 70, no. 4, pp. 2663–2677, 2022.
% However, I made some necessary modifications to make the manifold compatible with our system model and assumptions
%
% Contact person email: simon.w.tarboush@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates the Polar-domain dictionary proposed in [1] which is suitable for near-field channel estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% Qbar: The number of the ULA antennas
% OVS: Angular domain oversampling factor
% AE_spac:  Antenna elements spacing
% lambda: Wavelength
% dist: Sampled distance for Polar-domain dictionary
% theta: Sampled angle for Polar-domain dictionary
% beta, rho_min, rho_max: Parameters to define the Polar-domain dictionary [1]
% Quantized_Phase_Shifts: The actual quantized phase shift values
% useQPH: Choose to use the quantized phase shift values or not; 'True' or 'False'
% Output Arguments:
% PolarDictionary: The Polar-domain dictionary \bar{P}_T/\bar{P}_R used in near-field compressed sensing based on [1]
% Label: the angular (linear) and distance (non-linear) samplings for the 2D grid based on Polar-domain dictionary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PolarDictionary, Label] = get_PolarDictionaryULA(Qbar, OVS, AE_spac, lambda, beta, rho_min, rho_max, Quantized_Phase_Shifts, useQPH)
QbarOVS = OVS*Qbar;
spatial_ang = -1 + 2/QbarOVS : 2/QbarOVS : 1;
PolarDictionary = cell( QbarOVS, 1 );
Label = cell( QbarOVS, 1);
Zmax = (Qbar * AE_spac)^2 / 2 / lambda / beta^2;
kmax = floor(Zmax/rho_min);
for indx = 1:QbarOVS
    Z = (Qbar * AE_spac)^2 * ( 1 - spatial_ang(indx)^2) / 2 / lambda / beta^2;
    kmax = floor(Z/rho_min);
    kmin = floor(Z/rho_max) + 1;
    dist = zeros(1, kmax - kmin + 2);
    dist(:,1) = (Qbar * AE_spac)^2 * 2 / lambda;
    dist(:,2:end) = Z./(kmin:kmax);
    PolarDictionary{indx} = zeros(Qbar, kmax + 1);
    Label{indx} = zeros(2, kmax + 1);
    for indx_dist = 1 : kmax - kmin + 2
        PolarDictionary{indx}(:, indx_dist) = get_Polar_Domain_Manifold(Qbar, AE_spac, lambda, dist(indx_dist), acos(spatial_ang(indx)), Quantized_Phase_Shifts, useQPH);
        Label{indx}(:, indx_dist) = [spatial_ang(indx), dist(indx_dist)]';
    end
end
PolarDictionary = merge(PolarDictionary, QbarOVS, Qbar); % converted to matrix
Label = merge(Label, QbarOVS, 2); % converted to matrix
end

function B = merge(A, D, Q)
S = zeros(1, D);
for idx = 1 : D
    S(idx) = size(A{idx}, 2);
end
B = zeros(Q, sum(S));
for idx = 1 : D
    B(:, sum(S(1:idx)) - S(idx) + 1: sum(S(1:idx))) = A{idx};
end
end