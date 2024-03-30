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
% This function updates the THz channel parameters based on the input q struct (q is the channel information struct)
% This is a modified version of the original code:
% https://github.com/Simon-Tarboush/TeraMIMO || https://github.com/hasarieddeen/TeraMIMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% q: Channel struct that contains the initial THz channel paramters
% Output Arguments:
% p: The updated channel struct that contains all of the final simulation parameters to generate the THz channel using TeraMIMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p  = update_channel_param_TIV(q)
%% Scenario Parameters
% Define channel type
if isfield(q,'channelType')
    p.channelType = q.channelType;
else
    p.channelType = 'LoS';
end
% Method for calucating the molecular absorption coefficient
if isfield(q,'absorptionType')
    p.absorptionType = q.absorptionType;
else
    p.absorptionType = 'Hitran';
    % Options are: /'Hitran' /'Approx1' /'Approx2'
    % Hitran: Exact absorption coefficient, valid in the frequency band: [0.1-10] THz
    % Approx1: First approximation of absorption coefficient, valid in the frequency band: [275-400] GHz
    % Approx2: Second approximation of absorption coefficient, valid in the frequency band: [100-450] GHz
end
% Spherical and planar wave model, HSPM
if isfield(q,'WaveModelSA') 
    p.WaveModelSA = q.WaveModelSA; 
else
    p.WaveModelSA = 'Sphere'; %'Sphere'/'Plane'
end
if isfield(q,'WaveModelAE') 
    p.WaveModelAE = q.WaveModelAE;
else
    p.WaveModelAE = 'Sphere';  %'Sphere'/'Plane'
end
% Antenna Model
if isfield(q,'AntennaModel')
    p.AntennaModel = q.AntennaModel;
else
    p.AntennaModel = 'Isotropic';
end
%% Constants
p.c = physconst('LightSpeed');                         % Speed of light in vacuum (m/sec)
p.T0 = 296;                                                        % Reference temperature (Kelvin)
p.p0 = 1;                                                            % Standard pressure (atm)
p.Kb = physconst('Boltzmann');                       % Boltzmann constant (Joule/Kelvin)
% System Constants
if isfield(q,'T')
    p.T = q.T;              % System temperature (Kelvin) [25°C Celsius], 0 Kelvin - 273.15 = -273.1°C
else 
    p.T = p.T0;
end
if isfield(q,'p')
    p.p = q.p;              % System pressure (atm)
else
    p.p = p.p0;
end
if isfield(q,'PLE')
    p.PLE = q.PLE;
else
    p.PLE = 2;              % Path loss exponent
end
%% Transmission Parameters (Fc & BW)
if isfield(q,'Fc')
    p.Fc = q.Fc;                        % Center frequency of the transmission bandwidth (Hz)
else
    p.Fc = 0.3e12;                   % Center frequency of the transmission bandwidth (Hz)
end
if isfield(q,'BW')
    p.BW = q.BW;                    % Total channel bandwidth (Hz)
else
    p.BW = 0.01e12;                % Total channel bandwidth (Hz)
end
if isfield(q,'Nsub_c')
    p.Nsub_c = q.Nsub_c;        % Number of subcarriers to divide the total Bandwidth (K-subcarriers)  
else
    p.Nsub_c = 16;                      % Number of subcarriers to divide the total Bandwidth (K-subcarriers)
end
if isfield(q,'Nsub_b')
    p.Nsub_b = q.Nsub_b;        % Number of sub-bands in each subcarrier
else
    p.Nsub_b = 1;                      % Number of sub-bands in each subcarrier
end
p = get_FrequencyRng_param(p);
%% Absorption coefficient calculation
if strcmp(p.absorptionType,'Hitran')
    if (p.Fc - p.BW/2) < 0.1e12 || (p.Fc + p.BW/2) > 10e12
        error('Error: HITRAN is valid in: [0.1-10] THz');
    end
    p = get_HITRAN_param(p);
elseif strcmp(p.absorptionType,'Approx1')
    if (p.Fc - p.BW/2) < 0.275e12 || (p.Fc + p.BW/2) > 0.4e12
        error('Error: approx1, approximation of abs_coef is valid in: [275-400] GHz');
    end
    if isfield(q,'rel_humidity')
        p.rel_humidity = q.rel_humidity;              % Relative humidity
    else
        p.rel_humidity = 50;                                  % Relative humidity
    end
    
elseif strcmp(p.absorptionType,'Approx2') 
    if (p.Fc - p.BW/2) < 0.1e12 || (p.Fc + p.BW/2) > 0.45e12
        error('Error: approx2, approximation of abs_coef is valid in: [100-450] GHz');
    end
    if isfield(q,'rel_humidity')
        p.rel_humidity = q.rel_humidity;              % Relative humidity
    else
        p.rel_humidity = 50;                                  % Relative humidity
    end
else
    error('This method for computing absorption coefficient isn''t implemented, options are: Hitran, Approx1, Approx2');
end
%% UM-MIMO transceiver design (AoSA structure)
% This part includes the following definitions for each component 'i' in {BS, UE}:
% 1- Number of subarrays (SAs): total number is [Q_i = M_i*N_i where M_i (# of rows) and N_i (# of columns)]
% 2- Number of antenna elements (AEs) inside each SA: total number is [Qbar_i = Mbar_i*Nbar_i where Mbar_i (# of rows) and Nbar_i (# of columns)]
% 3- SA Spacing: Deltaxx is defined from center of SA(q) to center of SA(q+1) (default spacD*p.lambda_c/2, where spacD = sM x Mbar_i or sN x Nbar_i), we only define spacD
% 4- AE Spacing: deltaxx is defined from center of AE(qb) to center of AE(qb+1) (default spacd*p.lambda_c/2, spacd = 1), we fix this to half-wavelength

% Tx AoSA parameters
if isfield(q,'Mt')
    p.Mt = q.Mt; % Number of transmitter SAs (row)
else
    p.Mt = 1;
end
if isfield(q,'Nt')
    p.Nt = q.Nt; % Number of transmitter SAs (column)
else
    p.Nt = 1;
end
if isfield(q,'Mat')
    p.Mat = q.Mat; % Number of transmitter AEs (row) inside each SA
else
    p.Mat = 1;
end
if isfield(q,'Nat')
    p.Nat = q.Nat; % Number of transmitter AEs (column) inside each SA
else
    p.Nat = 1;
end
if isfield(q,'Sep_factTx')
    p.Sep_factTx = q.Sep_factTx;
else
    p.Sep_factTx = 1;
end

p.DeltaMt = p.Sep_factTx*p.Mat;
p.DeltaNt = p.Sep_factTx*p.Nat;
p.Tx_AoSA = get_AoSA_param(p, p.Mt, p.Nt, p.Mat, p.Nat, p.DeltaMt, p.DeltaNt, 1, 1);

% Rx AoSA parameters
if isfield(q,'Mr')
    p.Mr = q.Mr; % Number of Receiver SAs (row)
else
    p.Mr = 1;
end
if isfield(q,'Nr')
    p.Nr = q.Nr; % Number of Receiver SAs (column)
else
    p.Nr = 1; 
end
if isfield(q,'Mar')
    p.Mar = q.Mar; % Number of receiver AEs (row) inside each SA
else
    p.Mar = 1;
end
if isfield(q,'Nar')
    p.Nar = q.Nar; % Number of receiver AEs (column) inside each SA
else
    p.Nar = 1;
end
if isfield(q,'Sep_factRx')
    p.Sep_factRx = q.Sep_factRx;
else
    p.Sep_factRx = 1;
end

p.DeltaMr = p.Sep_factRx*p.Mar;
p.DeltaNr = p.Sep_factRx*p.Nar;
p.Rx_AoSA = get_AoSA_param(p, p.Mr, p.Nr, p.Mar, p.Nar, p.DeltaMr, p.DeltaNr, 1, 1);
%% Design geometry
% Define local/global position and Euler angles

% Tx center 3D positions (global coordinates)
% Tx Euler rotation angles, following ZYX intrinsic rotation
% [Px Py Pz];[alphadot betadot gammadot] (degrees)

if isfield(q,'positionTx') 
    p.positionTx = q.positionTx;     % Tx center 3D positions (global coordinates)
else
    p.positionTx = [0; 0; 0];
end
if isfield(q,'eulerTx')
    p.eulerTx = q.eulerTx;           % Tx Euler rotation angles, following ZYX intrinsic rotation
else
    p.eulerTx = [0; 0; 0];
end
p.Tx = get_Geometry_param(p.positionTx(1), p.positionTx(2), p.positionTx(3), p.eulerTx(1), p.eulerTx(2), p.eulerTx(3));

% Rx center 3D positions (global coordinates)
% Rx Euler rotation angles, following ZYX intrinsic rotation
if isfield(q,'positionRx')
    p.positionRx = q.positionRx;     % Rx center 3D positions (global coordinates)
else
    p.positionRx = [1; 0; 0];
end
if isfield(q,'eulerRx')
    p.eulerRx = q.eulerRx;           % Rx Euler rotation angles, following ZYX intrinsic rotation
else
    p.eulerRx = [rad2deg(pi); 0; 0];
end
p.Rx = get_Geometry_param(p.positionRx(1), p.positionRx(2), p.positionRx(3), p.eulerRx(1), p.eulerRx(2), p.eulerRx(3)); 

% Transmission distance (m)(from center to center)
p.d_tx_rx = norm(p.Tx.Pos-p.Rx.Pos);
p.d_rx_tx = norm(p.Rx.Pos-p.Tx.Pos);

% to compute the far-field && near-field regions (based on some definitions of the boundaries)
% Maximum array size, here we ignore the physical dimension of each AE
p.D_max = max([p.Tx_AoSA.Dmax p.Rx_AoSA.Dmax]);

p.Dist_Fresnel_MISO_SIMO = 0.62*sqrt(p.D_max^3/p.lambda_c0);
p.Dist_Fraunhofer_MISO_SIMO = 2*(p.D_max)^2/p.lambda_c0;
p.Dist_Fraunhofer_MIMO = 2*(p.Tx_AoSA.Dmax+p.Rx_AoSA.Dmax)^2/p.lambda_c0;
p.Dist_Fraunhofer_ARD_MIMO = 4*p.Tx_AoSA.Dmax*p.Rx_AoSA.Dmax/p.lambda_c0; % see [24] in the manuscript for the definition
% % Maximum array size on the SA-level
p.TxSA_D_max = max([p.Mat-1 p.Nat-1]*p.lambda_c0/2);p.RxSA_D_max = max([p.Mar-1 p.Nar-1]*p.lambda_c0/2); 
p.SA_D_max = max([p.TxSA_D_max p.RxSA_D_max]);

p.SA_Dist_Fresnel_MISO_SIMO = 0.62*sqrt(p.SA_D_max^3/p.lambda_c0);
p.SA_Dist_Fraunhofer_MISO_SIMO = 2*(p.SA_D_max)^2/p.lambda_c0;
p.SA_Dist_Fraunhofer_MIMO = 2*(p.TxSA_D_max+p.RxSA_D_max)^2/p.lambda_c0;
p.SA_Dist_Fraunhofer_ARD_MIMO = 4*p.TxSA_D_max*p.RxSA_D_max/p.lambda_c0;
%% Tx, Rx antenna gain (frequency-independent)

switch p.AntennaModel
    case 'Isotropic'
        p.Tx.GaindBi = 0;    % Gain (dBi)
        p.Tx.Gain = 10.^(p.Tx.GaindBi/10); % Gain (scalar value)
        p.Rx.GaindBi = 0;    % Gain (dBi)
        p.Rx.Gain = 10.^(p.Rx.GaindBi/10); % Gain (scalar value)   
    case 'Directional'
        % Ref [1]: Constantine A. Balanis, "Antenna Theory: Analysis and Design" 4th edition Eq.(2-26)
        % assuming Directivity = Gain i.e. e_cd = e_c * e_d = 1
        % No losses, i.e. perfect conduction efficiency & perfect dielectric efficiency
        % Also, perfect reflection (mismatch) efficiency e_r (i.e. e_r=1)
        % Reflection coefficient equal to zero (matching between source & antenna)
        
        % half-power beamwidth in azimuth-plane, examples: rad2deg(2*sqrt(pi)) 27.7 60 120
        % half-power beamwidth in elevation-plane, examples: rad2deg(2*sqrt(pi)) 27.7 30 60
        % input value in (degree), output in (Rad)
        
        % Tx
        p.Tx.psi_azi = deg2rad(rad2deg(60));
        p.Tx.psi_elev = deg2rad(rad2deg(30));
        p.Tx.GaindBi = pow2db(4*pi/(p.Tx.psi_azi*p.Tx.psi_elev));    % Gain (dBi)
        % Eq. (51) in TeraMIMO IEEE_TVT
        p.Tx.Gain = 10.^(p.Tx.GaindBi/10);                           % Gain (scalar value)
        % Rx
        p.Rx.psi_azi = deg2rad(rad2deg(60));
        p.Rx.psi_elev = deg2rad(rad2deg(30));
        p.Rx.GaindBi = pow2db(4*pi/(p.Rx.psi_azi*p.Rx.psi_elev));    % Gain (dBi)
        % Eq. (51) in TeraMIMO IEEE_TVT
        p.Rx.Gain = 10.^(p.Rx.GaindBi/10);                           % Gain (scalar value)        
    otherwise
        error('This Antenna Model isn''t implemented !!');
end

%% MPC parameters based on Multi-ray channel modeling
% This code is inspired by the following papers:
% [1] C. Han, A. O. Bicen, and I. F. Akyildiz, “Multi-ray channel modeling and wideband characterization for wireless communications in the terahertz band,” 
%       IEEE Trans. Wireless Commun., vol. 14, no. 5, pp. 2402–2412, 2014.
% [2] K. Dovelos et al., “Channel estimation and hybrid combining for wideband terahertz massive MIMO systems,” 
%       IEEE J. Sel. Areas Commun., vol. 39, no. 6, pp. 1604–1620, 2021.
if isfield(q,'L')
    p.L = q.L;
else
    p.L = 3;   % number of paths (the first component will be the LoS and the remaining (L-1) are NLoS)
end
% Required parameters to define Eqs. (17) and (18) in the manuscript
% for details check: Sec.II-C THz Specific Channel Characteristics and TABLE II: Simulation parameters in the manuscript
if isfield(q,'kappa_T')
    p.kappa_T = q.kappa_T;
else
    p.kappa_T = 2.24 - 0.025i;      % refractive index
end
if isfield(q,'sigma_rough')
    p.sigma_rough= q.sigma_rough;
else
    p.sigma_rough = 0.088e-3;   % roughness factor
end
p.varphi_in = (pi/2)*rand(p.L,1);                                   % Angle of incidence (for ell > 1). Random Variable ~ U(0,pi/2)
p.varphi_ref = asin((1/p.kappa_T)*sin(p.varphi_in));   % Angle of refracted wave (for ell > 1)
p.d_min = p.d_tx_rx+p.d_tx_rx/1000;              % Scatter distance range: [d_min,d_max]
p.d_max = p.d_tx_rx+6*p.d_tx_rx;                    % Scatter distance range: [d_min,d_max]
p.d_ell = unifrnd(p.d_min,p.d_max,[p.L 1]);    % Distance from ell-scatterer to the center of the AoSA. Random Variable ~ U(d_min,d_max)
p.d_ell(1) = p.d_tx_rx;        % The first component is the LoS always
p.loc_scat = rand(p.L,1);    % The location of any scatterer is also a random variable
p.d_ell_Tx = p.loc_scat.*p.d_ell;
p.d_ell_Rx = (1-p.loc_scat).*p.d_ell;
p.tauLoS = p.d_tx_rx/p.c;            % Delay of LoS path
p.tau_LoS_NLoS = p.d_ell/p.c;   % Delay of NLoS paths
p.tau_LoS_NLoS(1) = p.tauLoS;
% Generate random variable for AoA/AoD
p.theta_TxNLoS = pi*rand(p.L,1);    % Tx AoD (ULA-based implementation). Random Variable ~ U(0,pi)
p.phi_RxNLoS = pi*rand(p.L,1);       % Rx AoA (ULA-based implementation). Random Variable ~ U(0,pi)
end