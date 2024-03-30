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
% This function generates the channel estimation for a UM-MIMO THz channel based on the proposed algorithm "Cross-field Reduced-Dictionary (RD)" 
% and the baseline algorithms to ensure a fair comparison.
% The proposed algorithms 2 and 3 in the manuscript are implemented in this function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% Dist_val: Communication distance (m)
% BetaDot_val: Rotation angle (degree), for current implementation of the ULA, the rotation is around the y-axis (Sec II-B)
% Output Arguments:
% .mat file that contains the NMSE and AR for all methods including the proposed one.
% After running this function over multiple communication distances and array rotations,
% we can use the function "Fig9_Fig8_Proposed_CrossField_ReducedDictionary" to get Figs. 9 or 8 of the manuscript 
% (and Fig. 9 or 8 - after adjusting the simulation parameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function get_Proposed_CrossField_RD_Parallel(Dist_val, BetaDot_val)
%% Initialize Channel Parameters
p_ch = generate_channel_param_TIV();
p_ch.L = 3; % Number of channel paths: L=1 LoS only, you can use L= 3, 5, ...
p_ch.Sep_factTx = 2; % 1 for scenario 1 and 2 for scenario 2
p_ch.Sep_factRx = 10; % 1 for scenario 1 and 10 for scenario 2
% To change the aperture as per the paper and get Fig. 9 use: p_ch.Sep_factTx = 2; p_ch.Sep_factRx = 10
% To change the aperture as per the paper and get Fig. 9 use: p_ch.Sep_factTx = 2; p_ch.Sep_factRx = 10
p_ch = update_channel_param_TIV(p_ch);
%% Calculation of Absorption Coefficient 
K_abs = get_Abs_Coef(p_ch);
%% Main Simulation Parameters
RequiredSNR = 6;
E = 10; % 45 is the used value in the manuscript % Number of offline trials/runs/iterations
Lbar = 10; % Estimation of the sparsity level value (for SOMP based on for-loop)
gamma_S_H = 0.4;   % for Fig. 9 use 0.4, for Fig. 8 0.25
gamma_H_P = 0.12; % for Fig. 9 use 0.12, for Fig. 8 0.125
SaveSimulationResults = 1; % Decide to save or skip the simulation file
%% Parameters related to the transmitter and receiver uniform linear array
% Antennas Parameters
d_ae = 0.5; % antenna elements spacing, for DFT matrix d_ae normalized to lambda/2
% Number of SAs (For this work, since we use a ULA lies on Z-axis, Q_T_v > 1 and Q_T_h= 1)
Q_T_v = p_ch.Tx_AoSA.Qdim(1); % Number of transmit SAs on the z-axis of planar array
Q_T_h = p_ch.Tx_AoSA.Qdim(2); % Number of transmit SAs on the y-axis of planar array
Q_R_v = p_ch.Rx_AoSA.Qdim(1); % Number of receiver SAs on the z-axis of planar array
Q_R_h = p_ch.Rx_AoSA.Qdim(2); % Number of receiver SAs on the y-axis of planar array
% Number of AEs
Qbar_T_v = p_ch.Tx_AoSA.Qbardim(1); % Number of transmit antenna elements (AEs) per SA on the z-axis of planar array
Qbar_T_h = p_ch.Tx_AoSA.Qbardim(2); % Number of transmit antenna elements (AEs) per SA on the y-axis of planar array
Qbar_R_v = p_ch.Rx_AoSA.Qbardim(1); % Number of receiver antenna elements (AEs) per SA on the z-axis of planar array
Qbar_R_h = p_ch.Rx_AoSA.Qbardim(2); % Number of receiver antenna elements (AEs) per SA on the y-axis of planar array
% Total Number of SAs and AEs
Q_T = p_ch.Tx_AoSA.Q;  % Equal to: Q_T_h*Q_T_v; % Number of total SA Tx antennas
Q_R = p_ch.Rx_AoSA.Q;  % Equal to: Q_R_h*Q_R_v; % Number of total SA Rx antennas
Qbar_T = p_ch.Tx_AoSA.Qbar;% Equal to: Qbar_T_h*Qbar_T_v; % Number of total Tx antenna elements (AEs) per SA
Qbar_R = p_ch.Rx_AoSA.Qbar;% Equal to: Qbar_R_h*Qbar_R_v; % Number of total Rx antenna elements (AEs) per SA
%% Measurements/Beams during the training transmission phase
Ns_train = 1; % Number of spatial streams during training, always equal to 1
M_T = Qbar_T/2;         % Number of exhaustive-search training Tx (Full Training)/(Partial Training)
M_R = Qbar_R;            % Number of exhaustive-search training Rx (Full Training)/(Partial Training)
Compression_RatioTx = 1; Compression_RatioRx = 1;
Sel_SAs_IndTx = 1:Q_T;Sel_SAs_IndRx = 1:Q_R;
Num_of_UsedSATx = length(Sel_SAs_IndTx);
Num_of_UsedSARx = length(Sel_SAs_IndRx);
M_T_measMAT = zeros(Q_T,1);  % The number of measurements at Tx
M_R_measMAT = zeros(Q_R,1); % The number of measurements at Rx
M_T_measMAT(Sel_SAs_IndTx,1) = Compression_RatioTx*M_T; % The number of measurements at Tx
M_R_measMAT(Sel_SAs_IndRx,1) = Compression_RatioRx*M_R; % The number of measurements at Rx
%% Parameters related to the Quantization of Phase Shifters for the Codebooks
Q_T_quant = 2; % Number of Tx PSs quantization bits
Q_R_quant = 2; % Number of Rx PSs quantization bits
% Quantization of Phase for the Codebooks and RF Beamformers/Combiners
QPhaseShifts_tx_sector_start = 0;QPhaseShifts_tx_sector_end = 2*pi;
QPhaseShifts_rx_sector_start = 0;QPhaseShifts_rx_sector_end = 2*pi;
% We add hardware constraints on the Analog Phase Shifters
QPhaseShifts_tx_sector_rng = QPhaseShifts_tx_sector_end-QPhaseShifts_tx_sector_start;
QPhaseShifts_rx_sector_rng = QPhaseShifts_rx_sector_end-QPhaseShifts_rx_sector_start;
% Possible values of Tx and Rx Quantized phase-shifts following finite resolution phase-shifters (PSs)
Quant_PS_T = QPhaseShifts_tx_sector_start:QPhaseShifts_tx_sector_rng/2^Q_T_quant:QPhaseShifts_tx_sector_end-QPhaseShifts_tx_sector_rng/2^Q_T_quant;
Quant_PS_R = QPhaseShifts_rx_sector_start:QPhaseShifts_rx_sector_rng/2^Q_R_quant:QPhaseShifts_rx_sector_end-QPhaseShifts_rx_sector_rng/2^Q_R_quant;
%% Define the Tx power based on the required SNR
% Subcarriers and bandwidth
K = p_ch.Nsub_c;     % Number of subcarriers
B_sys = p_ch.BW;     % System Bandwidth
% Noise Power
N0_dBm = -173.8+10*log10(B_sys); % Noise power (dBm) % 10*log10(T*Kb) + 30 = -173.8
N0_sc = 10.^((N0_dBm-30)/10); % in Watts
% PL
lambda  = p_ch.lambda_c0; % wavelength
PLE_val = p_ch.PLE; % path loss exponent
d_TxRx3D = Dist_val; % 3D distance between transmitter and receiver
PL = ((4*pi*d_TxRx3D)./lambda).^PLE_val.*exp(K_abs(ceil(K/2),1)*d_TxRx3D); % Path loss check TeraMIMO IEEE_TVT for the definition
PL_dB = 10*log10(PL);
% Tx Power Scalar
Tx_power_tot_dBm = 10*log10(K) + RequiredSNR + N0_dBm + PL_dB;
Tx_power_tot = 10.^((Tx_power_tot_dBm-30)/10); % Total Tx power in Watt
%% Construct the Tx and Rx Dictionaries
% Parameters related to the grid/ovresampling of used dictionaries for SOMP
G_OVS = 2; % Dictionary Oversampling factor for proposed RD method
% Define Parameters for the Tx/Rx intermediate-/far-field dictionaries
num_multG = 1; % grid oversampling (=1 : no oversampling in the angular domain)
% Note that for efficient on-grid algorithm: G >= max(Qbar_T,Qbar_R)
% number of grids G_T and G_R
G_T_v = num_multG*Qbar_T_v; % Number of Quantization levels for Tx vertical dimension dictionary of UPA per SA
G_T_h = num_multG*Qbar_T_h; % Number of Quantization levels for Tx horizontal dimension dictionary of UPA per SA
G_R_v = num_multG*Qbar_R_v; % Number of Quantization levels for Rx vertical dimension dictionary of UPA per SA
G_R_h = num_multG*Qbar_R_h; % Number of Quantization levels for Rx horizontal dimension dictionary of UPA per SA
G_T = G_T_h*G_T_v;          % Total Number of Quantization levels for Tx dictionary
G_R = G_R_h*G_R_v;         % Total Number of Quantization levels for Rx dictionary
% Define the Tx/Rx intermediate-/far-field dictionaries
% for channel estimation:
Abar_T = get_DictionaryULA(Quant_PS_T, Qbar_T, G_T,'False');
Abar_R = get_DictionaryULA(Quant_PS_R, Qbar_R, G_R,'False');
% for AR Analysis
Abar_T_quant = get_DictionaryULA(Quant_PS_T, Qbar_T, G_T,'True');
Abar_R_quant = get_DictionaryULA(Quant_PS_R, Qbar_R, G_R,'True');
% Get the oversampled dictionaries used in the proposed RD method 
Abar_T_OVS = get_DictionaryULA(Quant_PS_T, Qbar_T, G_OVS*G_T,'False');
Abar_R_OVS = Abar_R; % Since the number of AEs of a SA at the Rx is less than 64 % get_DictionaryULA(Quant_PS_R, Qbar_R, G_OVS*G_R,'False'); 
% for AR Analysis
Abar_T_OVS_quant = get_DictionaryULA(Quant_PS_T, Qbar_T, G_OVS*G_T,'True');
Abar_R_OVS_quant = Abar_R_quant; % get_DictionaryULA(Quant_PS_R, Qbar_R, G_OVS*G_R,'True');
% Define the Tx/Rx near-field dictionaries (Polar-domain)
% Polar Dictionary Parameters: for more details, please read [1]
% [1] M. Cui and L. Dai, “Channel estimation for extremely large-scale MIMO: Far-field or near-field?” IEEE Trans. Commun., vol. 70, no. 4, pp. 2663–2677, 2022.
Beta_PolDict = 1.2;
Dist_min_PolDict = 0.5;
Dist_max_PolDict = 64;
% Define the Tx/Rx Polar-domain near-field dictionary 
% for channel estimation:
[Pbar_T, Labels_Pbar_T] = get_PolarDictionaryULA(Qbar_T, num_multG, d_ae*lambda, lambda, Beta_PolDict, Dist_min_PolDict, Dist_max_PolDict, Quant_PS_T, 'False');
[Pbar_R, Labels_Pbar_R] = get_PolarDictionaryULA(Qbar_R, num_multG, d_ae*lambda, lambda, Beta_PolDict, Dist_min_PolDict, Dist_max_PolDict, Quant_PS_R, 'False');
% for AR Analysis
[Pbar_T_quant, Labels_Pbar_T_quant] = get_PolarDictionaryULA(Qbar_T, num_multG, d_ae*lambda, lambda, Beta_PolDict, Dist_min_PolDict, Dist_max_PolDict, Quant_PS_T, 'True');
[Pbar_R_quant, Labels_Pbar_R_quant] = get_PolarDictionaryULA(Qbar_R, num_multG, d_ae*lambda, lambda, Beta_PolDict, Dist_min_PolDict, Dist_max_PolDict, Quant_PS_R, 'True');
% Get the oversampled dictionaries used in the proposed RD method 
[Pbar_T_OVS, Labels_Pbar_T_OVS] = get_PolarDictionaryULA(Qbar_T, G_OVS, d_ae*lambda, lambda, Beta_PolDict, Dist_min_PolDict, Dist_max_PolDict, Quant_PS_T, 'False');
Pbar_R_OVS = Pbar_R;
Labels_Pbar_R_OVS = Labels_Pbar_R;
% get_PolarDictionaryULA(Qbar_R, G_OVS, d_ae*lambda, lambda, Beta_PolDict, Dist_min_PolDict, Dist_max_PolDict, Quant_PS_R, 'False');
%% Define the Dictionary Matrix
Thetabar = kron(Abar_T,Abar_R); % is the matrix \bar{\Theta} (Eqs. 34 and 35 in the manuscript, can be used for PWM- and HSPWM-based estimation)
Thetabar_SWM = kron(Pbar_T,Pbar_R); % is the matrix \bar{\Theta} (Eq. 33 in the manuscript, can be used for the SWM-based estimation)
% TeraMIMO channel implementation uses transpose for the Tx array response vectors (ARVs), if the channel model uses Hermitian instead of transpose, we should use conj(Abar_T)
%% Reduced Dictionary Initial Simulation Parameters
% Initial values will be changed inside the proposed RD method because these values are linked to the dictionary sizes
G_T_RD = Qbar_T/2;
G_R_RD = Qbar_R;
G_T_Pol_RD = floor(size(Pbar_T,2)/4);
G_R_Pol_RD = size(Pbar_R,2);
% step size to find the required phase shifts for PWM-RD estimation 
StepSize_Phase = 1000;
%% Initialize Arrays for Results: NMSE, AR, Proposed Metric, Complexity
% NMSE
NMSE_SOMP_NF = zeros(E,K);
NMSE_SOMP_FF = zeros(E,K);
NMSE_SOMP_Hybrid_NF_FF = zeros(E,K);
NMSE_SOMP_Hybrid_FF_NF = zeros(E,K);
NMSE_SOMP_Cross = zeros(E,K);
% AR
AR_PCSI = zeros(E,K);
AR_SOMP_NF = zeros(E,K);
AR_SOMP_FF = zeros(E,K);
AR_SOMP_Hybrid_NF_FF = zeros(E,K);
AR_SOMP_Hybrid_FF_NF = zeros(E,K);
AR_SOMP_Cross = zeros(E,K);
% Proposed metric
if Q_R>1
    Num_Poss_EtaComb = nchoosek(Q_R,2);
else
    Num_Poss_EtaComb =1;
end
Eta_Online = zeros(E, Num_Poss_EtaComb+1);
% Selected Channel Model
SelectedModel = strings(E,1);
% Complexity Analysis
Complexity_NF = zeros(E,1);
Complexity_FF = zeros(E,1);
Complexity_Cross = zeros(E,1);
%% Main Loop
tic;
% Update simulation and channel parameters
Tx_pwr = Tx_power_tot/Ns_train;
% Tx_pwr is "P_T" defined in Eq. 29 in the manuscript, i.e., is the total Tx power used per transmission during the training phase.
% during the training there is only one SA to be used, i.e., number of streams during training is Ns_train = 1
DataTx_pwr = Tx_pwr;
noise_pwr = N0_sc;
p_ch.positionRx = [Dist_val; 0; 0];
p_ch.eulerRx = [rad2deg(pi); BetaDot_val; 0];
p_ch = update_channel_param_TIV(p_ch);
parfor indx_iter = 1:E
    % AoSAs THz channel generation
    SWM = channel_TIV(p_ch, K_abs);
    % The spherical wave model (SWM) channel model is used in all simulations as it is the most accurate model for all communication distances
    H_AoSA = cell2mat(SWM.H);
    % Arrays/Cells to save UM-MIMO channel estimation, support index, beamforming/combining vectors
    % Reference and Estimated Channels
    H_AoSA_Cell = SWM.H;
    H_Est_SOMP_UM_NF = cell(Q_R,Q_T,K);H_Est_SOMP_UM_FF = cell(Q_R,Q_T,K);
    H_Est_SOMP_UM_Hybrid_NF_FF = cell(Q_R,Q_T,K);H_Est_SOMP_UM_Hybrid_FF_NF = cell(Q_R,Q_T,K);
    H_Est_SOMP_UM_Cross = cell(Q_R,Q_T,K);
    % Beamforming/Combining Vectors
    F_Hat_SOMP_NF = cell(Q_R,Q_T);W_Hat_SOMP_NF = cell(Q_R,Q_T);
    F_Hat_SOMP_FF = cell(Q_R,Q_T);W_Hat_SOMP_FF = cell(Q_R,Q_T);
    F_Hat_SOMP_Hybrid_NF_FF = cell(Q_R,Q_T);W_Hat_SOMP_Hybrid_NF_FF = cell(Q_R,Q_T);
    F_Hat_SOMP_Hybrid_FF_NF = cell(Q_R,Q_T);W_Hat_SOMP_Hybrid_FF_NF = cell(Q_R,Q_T);
    F_Hat_SOMP_Cross = cell(Q_R,Q_T);W_Hat_SOMP_Cross = cell(Q_R,Q_T);
    % Save grid size for complexity calculations
    Num_Grids_Exact_NF = zeros(Q_T,Q_R);
    Num_Grids_Exact_FF = zeros(Q_T,Q_R);
    % Support (AoD and AoA indices)
    AoDInd_SOMP_Cross_cand = cell(Q_R,Q_T);AoAInd_SOMP_Cross_cand = cell(Q_R,Q_T);
    % Best N Beams
    L_dom_Beams_SOMP_UM_Cross = zeros(Q_R,Q_T);
    % Proposed Metric Initial Matrix
    Chi_Online_Tmp = zeros(M_R*M_T*K,Q_R);
    % Save random beamformers, combiners, and noise to ensure a fair comparison
    Z_AoSA = complex(zeros(Qbar_T, Qbar_T, Num_of_UsedSATx));
    C_AoSA = complex(zeros(Qbar_R, Qbar_R, Num_of_UsedSATx));
    N_AoSA = cell(Num_of_UsedSATx, Num_of_UsedSARx);
    N_AoSA(:,:,:) = {complex(zeros(Qbar_R, Qbar_T, K))};
    for indx_UsedSATx = 1:Num_of_UsedSATx
        % RF Beamforming Matrix
        Z_AoSA(:,:,indx_UsedSATx) = get_RandomCodebook(Quant_PS_T, Q_T_quant, Qbar_T, Qbar_T); % Random beam training
        %  RF Combining Matrix 
        C_AoSA(:,:,indx_UsedSATx) = get_RandomCodebook(Quant_PS_R, Q_R_quant, Qbar_R, Qbar_R); % Random beam training
        for indx_UsedSARx = 1:Num_of_UsedSARx
            %%%%%%%%%%%%%%%%%% Noise Generation %%%%%%%%%%%%%%%%%%
            N_AoSA{indx_UsedSATx,indx_UsedSARx} = sqrt(K)*sqrt(noise_pwr)/sqrt(2)*(randn(Qbar_R, Qbar_T, K)+1j*randn(Qbar_R, Qbar_T, K));
        end
    end
    %%%%%%%%%%%%%%%% Conventional Baseline Methods %%%%%%%%%%%%%%%%
    %
    %%%%%%%%%%%%%%%%%% Begin loops over Tx/Rx SAs %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% loop over the Tx SAs %%%%%%%%%%%%%%%%%%
    for indx_UsedSATx = 1:Num_of_UsedSATx
        %%%%%%%%%%%%%%%%%% Select the Tx SA %%%%%%%%%%%%%%%%%%
        [qth,qtv] = ind2sub([Q_T_h Q_T_v],Sel_SAs_IndTx(indx_UsedSATx));
        %%%%%%%%%%%%%%%%%% RF Beamforming Matrix %%%%%%%%%%%%%%%%%%
        Z = squeeze(Z_AoSA(:,:,indx_UsedSATx));
        %%%%%%%%%%%%%%%%%%  RF Combining Matrix  %%%%%%%%%%%%%%%%%%
        % following our proposal in Sec. III-A, the combining weights are applied to all Rx SAs
        C = squeeze(C_AoSA(:,:,indx_UsedSATx));
        %%%%%%%%%%%%%%%%%% loop over the Rx SAs %%%%%%%%%%%%%%%%%%
        for indx_UsedSARx = 1:Num_of_UsedSARx
            %%%%%%%%%%%%%%%%%% Select the Rx SA %%%%%%%%%%%%%%%%%%
            [qrh,qrv] = ind2sub([Q_R_h Q_R_v],Sel_SAs_IndRx(indx_UsedSARx));
            %%%%%%%%%%%%%%%%%% Define the number of measurements %%%%%%%%%%%%%%%%%
            M_T_meas = M_T_measMAT((qtv-1)*Q_T_h+qth,1);
            M_R_meas = M_R_measMAT((qrv-1)*Q_R_h+qrh, 1);
            %%%%%%%%%%%%%%%%%% THz channel between q_T(th) Tx SA and q_R(th) Rx SA %%%%%%%%%%%%%%%%%%
            H_qRqT = H_AoSA(((qrv-1)*Q_R_h+qrh-1)*Qbar_R_v*Qbar_R_h+1:((qrv-1)*Q_R_h+qrh)*Qbar_R_v*Qbar_R_h,((qtv-1)*Q_T_h+qth-1)*Qbar_T_v*Qbar_T_h+1:((qtv-1)*Q_T_h+qth)*Qbar_T_v*Qbar_T_h,:);
            % This is H_{q_R,q_T} defined in Eq. (5) in the manuscript (for SWM, Eqs. (7) and (8))
            %%%%%%%%%%%%%%%%%% Noise Generation %%%%%%%%%%%%%%%%%%
            N = N_AoSA{indx_UsedSATx,indx_UsedSARx};
            %%%%%%%%%%%%%%%%%% The vectorized version of the received signal of Eq. (29) for all subcarriers based on Tx and Rx signals %%%%%%%%%%%%%%%%%%
            y = get_Channel_Output(H_qRqT,N(:,1:M_T_meas,:),Z(:,1:M_T_meas),C(:,1:M_R_meas),K,Tx_pwr);
            %%%%%%%%%%%%%%%%%% Compute Measurement Matrix %%%%%%%%%%%%%%%%%%
            Psi = kron(transpose(Z(:,1:M_T_meas)),(C(:,1:M_R_meas))');
            % following Eq. (32) in the manuscript Psi = kron(Z^Transpose, C^Herm) is the measurement matrix
            %%%%%%%%%%%%%%%%%% Compute Sensing Matrix %%%%%%%%%%%%%%%%%%
            % The Sensing matrix will be denoted as Upsilon = Psi*Thetabar (Measurement Matrix * Quantized Dictionary Matrix), 
            % check for example Eq. (33), (34), and (35)
            % for FF, Hybrid
            Upsilon_SOMP_FF = Psi*Thetabar;
            % for NF, Hybrid
            Upsilon_SOMP_NF = Psi*Thetabar_SWM;
            %%%%%%%%%%%%%%%%%% Save Sensing Matrix Size for Complexity Analysis %%%%%%%%%%%%%%%%%%
            Num_Grids_Exact_FF(indx_UsedSATx, indx_UsedSARx) = numel(Upsilon_SOMP_FF);
            Num_Grids_Exact_NF(indx_UsedSATx, indx_UsedSARx) = numel(Upsilon_SOMP_NF);
            %%%%%%%%%%%%%%%%%% SOMP CS-Estimator %%%%%%%%%%%%%%%%%%
            % bs here stands for BeamSpace domain 
            % SOMP NF  
            [hest_bs_SOMP_NF, support_SOMP_NF] = SOMP(y,Upsilon_SOMP_NF,Tx_pwr,Lbar,K);
            H_bs_SOMP_NF = reshape(hest_bs_SOMP_NF,size(Pbar_R,2),size(Pbar_T,2),K);
            % SOMP FF
            [hest_bs_SOMP_FF, support_SOMP_FF] = SOMP(y,Upsilon_SOMP_FF,Tx_pwr,Lbar,K);
            H_bs_SOMP_FF = reshape(hest_bs_SOMP_FF,size(Abar_R,2),size(Abar_T,2),K);            
            % SOMP Hybrid NF then FF 
            [hest_bs_SOMP_HybNFFF_FF, support_SOMP_HybNFFF_FF, hest_bs_SOMP_HybNFFF_NF, support_SOMP_HybNFFF_NF] = ...
                HybridField_SOMP_NF_FF(y,Upsilon_SOMP_FF,Upsilon_SOMP_NF,Tx_pwr,Lbar,K);
            H_bs_SOMP_HybNFFF_NF = reshape(hest_bs_SOMP_HybNFFF_NF,size(Pbar_R,2),size(Pbar_T,2),K);
            H_bs_SOMP_HybNFFF_FF = reshape(hest_bs_SOMP_HybNFFF_FF,size(Abar_R,2),size(Abar_T,2),K);
            % SOMP Hybrid FF then NF 
            [hest_bs_SOMP_HybFFNF_FF, support_SOMP_HybFFNF_FF, hest_bs_SOMP_HybFFNF_NF, support_SOMP_HybFFNF_NF] = ...
                HybridField_SOMP_FF_NF(y,Upsilon_SOMP_FF,Upsilon_SOMP_NF,Tx_pwr,Lbar,K);
            H_bs_SOMP_HybFFNF_NF = reshape(hest_bs_SOMP_HybFFNF_NF,size(Pbar_R,2),size(Pbar_T,2),K);
            H_bs_SOMP_HybFFNF_FF = reshape(hest_bs_SOMP_HybFFNF_FF,size(Abar_R,2),size(Abar_T,2),K); 
            %%%%%%%%%%%%%%%%%% Extract the support and select angles for beamforming/combining %%%%%%%%%%%%%%%%%%
            % SOMP NF            
            [rowMaxBeamRx_SOMP_NF,colMaxBeamTx_SOMP_NF] = ind2sub(size(H_bs_SOMP_NF),support_SOMP_NF(1:Lbar));
            F_Hat_SOMP_NF{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Pbar_T_quant(:,colMaxBeamTx_SOMP_NF(1));
            W_Hat_SOMP_NF{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Pbar_R_quant(:,rowMaxBeamRx_SOMP_NF(1));
            % SOMP FF
            [rowMaxBeamRx_SOMP_FF,colMaxBeamTx_SOMP_FF] = ind2sub(size(H_bs_SOMP_FF),support_SOMP_FF(1:Lbar));
            F_Hat_SOMP_FF{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Abar_T_quant(:,colMaxBeamTx_SOMP_FF(1));
            W_Hat_SOMP_FF{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Abar_R_quant(:,rowMaxBeamRx_SOMP_FF(1));
            % SOMP Hybrid NF then FF
            [rowMaxBeamRx_SOMP_HybridNFFF,colMaxBeamTx_SOMP_HybridNFFF] = ind2sub(size(H_bs_SOMP_HybNFFF_NF),support_SOMP_HybNFFF_NF(1:Lbar/2));            
            F_Hat_SOMP_Hybrid_NF_FF{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Pbar_T_quant(:,colMaxBeamTx_SOMP_HybridNFFF(1));
            W_Hat_SOMP_Hybrid_NF_FF{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Pbar_R_quant(:,rowMaxBeamRx_SOMP_HybridNFFF(1));
            % SOMP Hybrid FF then NF
            [rowMaxBeamRx_SOMP_HybridFFNF,colMaxBeamTx_SOMP_HybridFFNF] = ind2sub(size(H_bs_SOMP_HybFFNF_FF),support_SOMP_HybFFNF_FF(1:Lbar/2));
            F_Hat_SOMP_Hybrid_FF_NF{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Abar_T_quant(:,colMaxBeamTx_SOMP_HybridFFNF(1));
            W_Hat_SOMP_Hybrid_FF_NF{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Abar_R_quant(:,rowMaxBeamRx_SOMP_HybridFFNF(1));
            %%%%%%%%%%%%%%%%%% Fill the SA (MIMO channel) estimation into the UM-MIMO channels  %%%%%%%%%%%%%%%%%%
            H_Est_SOMP_NF = zeros(Qbar_R,Qbar_T,K);
            H_Est_SOMP_FF = zeros(Qbar_R,Qbar_T,K);
            H_Est_SOMP_Hybrid_NF_FF = zeros(Qbar_R,Qbar_T,K);
            H_Est_SOMP_Hybrid_FF_NF = zeros(Qbar_R,Qbar_T,K);
            for indx_subc = 1:K
                % Create Angle-Domain Channels
                % Relation with beamspace domain representation is presented in Eqs. (21), (25), and (27) in the manuscript
                % SOMP NF
                H_Est_SOMP_NF(:,:,indx_subc) = Pbar_R*H_bs_SOMP_NF(:,:,indx_subc)*Pbar_T.';
                % SOMP FF
                H_Est_SOMP_FF(:,:,indx_subc) = Abar_R*H_bs_SOMP_FF(:,:,indx_subc)*Abar_T.';
                % SOMP Hybrid NF then FF
                if ~isempty(support_SOMP_HybNFFF_NF)
                    if ~isempty(support_SOMP_HybNFFF_FF)
                        H_Est_SOMP_Hybrid_NF_FF(:,:,indx_subc) = Pbar_R*H_bs_SOMP_HybNFFF_NF(:,:,indx_subc)*Pbar_T.' + ...
                            Abar_R*H_bs_SOMP_HybNFFF_FF(:,:,indx_subc)*Abar_T.';
                    else
                        H_Est_SOMP_Hybrid_NF_FF(:,:,indx_subc) = Pbar_R*H_bs_SOMP_HybNFFF_NF(:,:,indx_subc)*Pbar_T.';
                    end
                else
                    H_Est_SOMP_Hybrid_NF_FF(:,:,indx_subc) = Abar_R*H_bs_SOMP_HybNFFF_FF(:,:,indx_subc)*Abar_T.';
                end
                % SOMP Hybrid FF then NF
                if  ~isempty(support_SOMP_HybFFNF_FF)
                    if ~isempty(support_SOMP_HybFFNF_NF)
                        H_Est_SOMP_Hybrid_FF_NF(:,:,indx_subc) = Pbar_R*H_bs_SOMP_HybFFNF_NF(:,:,indx_subc)*Pbar_T.' + ...
                            Abar_R*H_bs_SOMP_HybFFNF_FF(:,:,indx_subc)*Abar_T.';
                    else
                        H_Est_SOMP_Hybrid_FF_NF(:,:,indx_subc) = Abar_R*H_bs_SOMP_HybFFNF_FF(:,:,indx_subc)*Abar_T.';
                    end
                else
                    H_Est_SOMP_Hybrid_FF_NF(:,:,indx_subc) = Pbar_R*H_bs_SOMP_HybFFNF_NF(:,:,indx_subc)*Pbar_T.';
                end
                % Fill the UM-MIMO arrays by estimations
                H_Est_SOMP_UM_NF{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,indx_subc} = H_Est_SOMP_NF(:,:,indx_subc);
                H_Est_SOMP_UM_FF{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,indx_subc} = H_Est_SOMP_FF(:,:,indx_subc);
                H_Est_SOMP_UM_Hybrid_NF_FF{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,indx_subc} = H_Est_SOMP_Hybrid_NF_FF(:,:,indx_subc);
                H_Est_SOMP_UM_Hybrid_FF_NF{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,indx_subc} = H_Est_SOMP_Hybrid_FF_NF(:,:,indx_subc);
            end
            %%%%%%%%%%%%%%%%%% End of loop on Tx-Rx SA channel Estimation  %%%%%%%%%%%%%%%%%%
        end
    end
    %%%%%%%%%%%%%%%%%% End of loops over UM-MIMO channels  %%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%% Proposed Cross-field Reduced Dictionary Solution %%%%%%%%%%%%%%%%
    %
    % Decide the appropriate estimation strategy based on the proposed model-selection metric in the online operation
    %
    % Compute the proposed metric
    %
    %%%%%%%%%%%%%%%%%% Select the Tx SA %%%%%%%%%%%%%%%%%%
    % Start training only from the 1st "reference" SA. The received signals only from the reference Tx SA as it is sufficient for determining the region
    [qth,qtv] = ind2sub([Q_T_h Q_T_v],Sel_SAs_IndTx(1));
    %%%%%%%%%%%%%%%%%% RF Beamforming Matrix %%%%%%%%%%%%%%%%%%
    Z = squeeze(Z_AoSA(:,:,1)); % Random beam training 
    %%%%%%%%%%%%%%%%%%  RF Combining Matrix  %%%%%%%%%%%%%%%%%%
    % following our proposal in Sec. III-A, the combining weights are applied to all Rx SAs
    C = squeeze(C_AoSA(:,:,1)); % Random beam training
    for indx_UsedSARx = 1:Num_of_UsedSARx
        %%%%%%%%%%%%%%%%%% Select the Rx SA %%%%%%%%%%%%%%%%%%
        [qrh,qrv] = ind2sub([Q_R_h Q_R_v],Sel_SAs_IndRx(indx_UsedSARx));
        %%%%%%%%%%%%%%%%%% Define the number of measurements %%%%%%%%%%%%%%%%%
        M_T_meas = M_T_measMAT((qtv-1)*Q_T_h+qth,1);
        M_R_meas = M_R_measMAT((qrv-1)*Q_R_h+qrh,1);
        %%%%%%%%%%%%%%%%%% THz channel between q_T(th) Tx SA and q_R(th) Rx SA %%%%%%%%%%%%%%%%%%
        H_qRqT = H_AoSA(((qrv-1)*Q_R_h+qrh-1)*Qbar_R_v*Qbar_R_h+1:((qrv-1)*Q_R_h+qrh)*Qbar_R_v*Qbar_R_h,((qtv-1)*Q_T_h+qth-1)*Qbar_T_v*Qbar_T_h+1:((qtv-1)*Q_T_h+qth)*Qbar_T_v*Qbar_T_h,:);
        % This is H_{q_R,q_T} defined in Eq. (5) in the manuscript (for SWM, Eqs. (7) and (8))
        %%%%%%%%%%%%%%%%%% Noise Generation %%%%%%%%%%%%%%%%%%
        N = N_AoSA{1,indx_UsedSARx};
        %%%%%%%%%%%%%%%%%% The vectorized version of the received signal of Eq. (29) for all subcarriers based on Tx and Rx signals %%%%%%%%%%%%%%%%%%
        Chibar_qr_1 = get_Channel_Output(H_qRqT,N(:,1:M_T_meas,:),Z(:,1:M_T_meas),C(:,1:M_R_meas),K,Tx_pwr);
        %%%%%%%%%%%%%%%%%% Compute the normalized vector \boldsymbol{\chi} used later in Eq. (30) %%%%%%%%%%%%%%%%%%
        Chi_Online_Tmp(:,indx_UsedSARx) = abs(Chibar_qr_1(:))/sqrt(sum(abs(Chibar_qr_1(:)).^2));
    end
    %%%%%%%%%%%%%%%%%% Compute all possible Combinations for the proposed model selection metric using Eqs. (30) and (31) %%%%%%%%%%%%%%%%%%
    indx_cnt = 0;
    Eta_Online_r_c = zeros(Num_Poss_EtaComb,1);
    for indx_pwr_row = 1:Q_R-1
        for indx_pwr_col = indx_pwr_row+1:Q_R
            indx_cnt = indx_cnt + 1;
            Chi_diff = Chi_Online_Tmp(:,indx_pwr_row)-Chi_Online_Tmp(:,indx_pwr_col);
            % Eq. (30) in the manuscript
            Eta_Online_r_c(indx_cnt,1) = sum(abs(Chi_diff).^2);
        end
    end
    % Eq. (30) and (31) in the manuscript, the max value (last element) is the proposed metric, however, we keep all possible values for analysis purposes
    Eta_Online(indx_iter,:) =  [Eta_Online_r_c(1:Num_Poss_EtaComb); max(Eta_Online_r_c(1:Num_Poss_EtaComb))];
    % Apply the decision rule based on the Eta_Online_Max and pre-computed offline thresholds
    Eta_Online_Max = max(Eta_Online_r_c(1:Num_Poss_EtaComb));
    if Eta_Online_Max >= gamma_S_H
        Est_Ch_Model = "SWM";
    elseif Eta_Online_Max >=  gamma_H_P && Eta_Online_Max < gamma_S_H
        Est_Ch_Model = "HSPWM";
    else
        Est_Ch_Model = "PWM";
    end
    SelectedModel(indx_iter,1) = Est_Ch_Model;
    %%%%%%%%%%%%%%%%%% Begin loops over Tx/Rx SAs %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% loop over the Tx SAs %%%%%%%%%%%%%%%%%%
    switch Est_Ch_Model
        case "SWM"
            % Initialize parameters for the proposed method and complexity analysis
            Pbar_T_RD_quant = []; Pbar_R_RD_quant = [];       
            Num_Grids_Exact_Cross = zeros(Q_T, Q_R);
            %%%%%%%%%%%%%%%%%% Begin loops over Tx/Rx SAs %%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%% loop over the Tx SAs %%%%%%%%%%%%%%%%%%
            for indx_UsedSATx = 1:Num_of_UsedSATx
                %%%%%%%%%%%%%%%%%% Select the Tx SA %%%%%%%%%%%%%%%%%%
                [qth,qtv] = ind2sub([Q_T_h Q_T_v],Sel_SAs_IndTx(indx_UsedSATx));
                %%%%%%%%%%%%%%%%%% RF Beamforming Matrix %%%%%%%%%%%%%%%%%%
                Z = squeeze(Z_AoSA(:,:,indx_UsedSATx));
                %%%%%%%%%%%%%%%%%%  RF Combining Matrix  %%%%%%%%%%%%%%%%%%
                % following our proposal in Sec. III-A, the combining weights are applied to all Rx SAs
                C = squeeze(C_AoSA(:,:,indx_UsedSATx));
                %%%%%%%%%%%%%%%%%% loop over the Rx SAs %%%%%%%%%%%%%%%%%%
                for indx_UsedSARx = 1:Num_of_UsedSARx
                    %%%%%%%%%%%%%%%%%% Select the Rx SA %%%%%%%%%%%%%%%%%%
                    [qrh,qrv] = ind2sub([Q_R_h Q_R_v],Sel_SAs_IndRx(indx_UsedSARx));
                    %%%%%%%%%%%%%%%%%% Define the number of measurements %%%%%%%%%%%%%%%%%
                    M_T_meas = M_T_measMAT((qtv-1)*Q_T_h+qth,1);
                    M_R_meas = M_R_measMAT((qrv-1)*Q_R_h+qrh, 1);
                    %%%%%%%%%%%%%%%%%% Proposed Algorithm: Cross-Field RD (SWM case) %%%%%%%%%%%%%%%%%%
                    % based on "Algorithm 2 Reduced Dictionary Method" in the manuscript and
                    % Sec. III-C1 SWM-based channel estimation with RD.
                    if indx_UsedSARx > 1
                        % Detect the best previously estimated candidate support based on the first "reference" T-R SA channel estimation (we added -1 to the index of Rx-plane)
                        Supp_T_Est = AoDInd_SOMP_Cross_cand{(qrv-1)*Q_R_h+qrh -1,(qtv-1)*Q_T_h+qth};                 
                        % Construct the Reduced Dictionary based on the previous support by selecting the most correlated
                        % columns from an oversampled dictionary following Algorithm 2 in the manuscript
                        if indx_UsedSARx == 2
                            Pbar_T_Est = Pbar_T(:,Supp_T_Est);
                            Pbar_T_RD = get_Proposed_RD_NearField(Pbar_T_OVS, Labels_Pbar_T_OVS, Pbar_T_Est, G_T_Pol_RD);
                            Pbar_T_RD_quant = get_Quantized_Dict(Pbar_T_RD, Quant_PS_T);
                        else
                            Pbar_T_Est = Pbar_T_RD(:,Supp_T_Est);
                            Pbar_T_RD = get_Proposed_RD_NearField(Pbar_T_OVS, Labels_Pbar_T_OVS, Pbar_T_Est, G_T_Pol_RD);
                            Pbar_T_RD_quant = get_Quantized_Dict(Pbar_T_RD, Quant_PS_T);
                        end
                        % Apply the same method for the receiver side if the MIMO array is massive (per SA: AEs >= 64)
                        % For now keep the receiver without applying the RD 
                        Pbar_R_RD = Pbar_R;
                        Pbar_R_RD_quant = Pbar_R_quant;
                        % Compute RD Dictionary Matrix
                        Thetabar_SWM_RD = kron(Pbar_T_RD, Pbar_R_RD);
                    end
                    %%%%%%%%%%%%%%%%%% THz channel between q_T(th) Tx SA and q_R(th) Rx SA %%%%%%%%%%%%%%%%%%
                    H_qRqT = H_AoSA(((qrv-1)*Q_R_h+qrh-1)*Qbar_R_v*Qbar_R_h+1:((qrv-1)*Q_R_h+qrh)*Qbar_R_v*Qbar_R_h,((qtv-1)*Q_T_h+qth-1)*Qbar_T_v*Qbar_T_h+1:((qtv-1)*Q_T_h+qth)*Qbar_T_v*Qbar_T_h,:);
                    % This is H_{q_R,q_T} defined in Eq. (5) in the manuscript (for SWM, Eqs. (7) and (8))
                    %%%%%%%%%%%%%%%%%% Noise Generation %%%%%%%%%%%%%%%%%%
                    N = N_AoSA{indx_UsedSATx,indx_UsedSARx};
                    %%%%%%%%%%%%%%%%%% The vectorized version of the received signal of Eq. (29) for all subcarriers based on Tx and Rx signals %%%%%%%%%%%%%%%%%%
                    y = get_Channel_Output(H_qRqT,N(:,1:M_T_meas,:),Z(:,1:M_T_meas),C(:,1:M_R_meas),K,Tx_pwr);
                    %%%%%%%%%%%%%%%%%% Compute Measurement Matrix %%%%%%%%%%%%%%%%%%
                    Psi = kron(transpose(Z(:,1:M_T_meas)),(C(:,1:M_R_meas))');
                    % following Eq. (32) in the manuscript Psi = kron(Z^Transpose, C^Herm) is the measurement matrix
                    %%%%%%%%%%%%%%%%%% Compute Sensing Matrix %%%%%%%%%%%%%%%%%%
                    % The SWM RD sensing matrix 
                    if indx_UsedSARx > 1
                        Upsilon_SOMP_Cross = Psi*Thetabar_SWM_RD;
                    else
                        Upsilon_SOMP_Cross = Psi*Thetabar_SWM;
                    end
                    %%%%%%%%%%%%%%%%%% Save Sensing Matrix Size for Complexity Analysis %%%%%%%%%%%%%%%%%%
                    Num_Grids_Exact_Cross(indx_UsedSATx, indx_UsedSARx) = numel(Upsilon_SOMP_Cross);
                    %%%%%%%%%%%%%%%%%% Cross Field SWM-RD SOMP CS-Estimator %%%%%%%%%%%%%%%%%%
                    % bs here stands for BeamSpace domain
                    % Cross Field SWM-RD + SOMP
                    [hest_bs_SOMP_Cross, support_SOMP_Cross] = SOMP(y,Upsilon_SOMP_Cross,Tx_pwr,Lbar,K);
                    if indx_UsedSARx > 1
                        H_bs_SOMP_Cross = reshape(hest_bs_SOMP_Cross,size(Pbar_R_RD,2),size(Pbar_T_RD,2),K);
                    else
                        H_bs_SOMP_Cross = reshape(hest_bs_SOMP_Cross,size(Pbar_R,2),size(Pbar_T,2),K);
                    end
                    %%%%%%%%%%%%%%%%%% Extract the number of dominant/best beams %%%%%%%%%%%%%%%%%%
                    % Accumulate until you get 95% of the estimated channel power
                    L_dom_Beam_SOMP = get_NumBestBeams(hest_bs_SOMP_Cross, support_SOMP_Cross);
                    L_dom_Beams_SOMP_UM_Cross((qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth) = L_dom_Beam_SOMP;
                    %%%%%%%%%%%%%%%%%% Extract the support and select angles for beamforming/combining %%%%%%%%%%%%%%%%%%
                    % SOMP Cross RD SWM
                    [rowMaxBeamRx_SOMP_Cross,colMaxBeamTx_SOMP_Cross] = ind2sub(size(H_bs_SOMP_Cross),support_SOMP_Cross(1:Lbar));
                    AoDInd_SOMP_Cross_cand{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = colMaxBeamTx_SOMP_Cross;
                    AoAInd_SOMP_Cross_cand{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = rowMaxBeamRx_SOMP_Cross;
                    if indx_UsedSARx > 1
                        F_Hat_SOMP_Cross{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Pbar_T_RD_quant(:,colMaxBeamTx_SOMP_Cross(1));
                        W_Hat_SOMP_Cross{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Pbar_R_RD_quant(:,rowMaxBeamRx_SOMP_Cross(1));
                    else
                        F_Hat_SOMP_Cross{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Pbar_T_quant(:,colMaxBeamTx_SOMP_Cross(1));
                        W_Hat_SOMP_Cross{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Pbar_R_quant(:,rowMaxBeamRx_SOMP_Cross(1));
                    end
                    %%%%%%%%%%%%%%%%%% Fill the SA (MIMO channel) estimation into the UM-MIMO channels  %%%%%%%%%%%%%%%%%%
                    H_Est_SOMP_Cross = zeros(Qbar_R,Qbar_T,K);
                    for indx_subc = 1:K
                        % Create Angle-Domain Channels
                        % SOMP Cross RD SWM
                        if indx_UsedSARx > 1
                            H_Est_SOMP_Cross(:,:,indx_subc) = Pbar_R_RD*H_bs_SOMP_Cross(:,:,indx_subc)*Pbar_T_RD.';
                        else
                            H_Est_SOMP_Cross(:,:,indx_subc) = Pbar_R*H_bs_SOMP_Cross(:,:,indx_subc)*Pbar_T.';
                        end
                        % Fill the UM-MIMO arrays by estimations
                        H_Est_SOMP_UM_Cross{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,indx_subc} = H_Est_SOMP_Cross(:,:,indx_subc);
                    end
                    %%%%%%%%%%%%%%%%%% End of loop on Tx-Rx SA channel Estimation  %%%%%%%%%%%%%%%%%%
                end
            end
            %%%%%%%%%%%%%%%%%% End of loops over UM-MIMO channels  %%%%%%%%%%%%%%%%%%
        case "HSPWM"
            % Initialize parameters for the proposed method and complexity analysis
            Abar_T_RD_quant = []; Abar_R_RD_quant = [];       
            Num_Grids_Exact_Cross = zeros(Q_T, Q_R);
            %%%%%%%%%%%%%%%%%% Begin loops over Tx/Rx SAs %%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%% loop over the Tx SAs %%%%%%%%%%%%%%%%%%
            for indx_UsedSATx = 1:Num_of_UsedSATx
                %%%%%%%%%%%%%%%%%% Select the Tx SA %%%%%%%%%%%%%%%%%%
                [qth,qtv] = ind2sub([Q_T_h Q_T_v],Sel_SAs_IndTx(indx_UsedSATx));
                %%%%%%%%%%%%%%%%%% RF Beamforming Matrix %%%%%%%%%%%%%%%%%%
                Z = squeeze(Z_AoSA(:,:,indx_UsedSATx));
                %%%%%%%%%%%%%%%%%%  RF Combining Matrix  %%%%%%%%%%%%%%%%%%
                % following our proposal in Sec. III-A, the combining weights are applied to all Rx SAs
                C = squeeze(C_AoSA(:,:,indx_UsedSATx));
                %%%%%%%%%%%%%%%%%% loop over the Rx SAs %%%%%%%%%%%%%%%%%%
                for indx_UsedSARx = 1:Num_of_UsedSARx
                    %%%%%%%%%%%%%%%%%% Select the Rx SA %%%%%%%%%%%%%%%%%%
                    [qrh,qrv] = ind2sub([Q_R_h Q_R_v],Sel_SAs_IndRx(indx_UsedSARx));
                    %%%%%%%%%%%%%%%%%% Define the number of measurements %%%%%%%%%%%%%%%%%
                    M_T_meas = M_T_measMAT((qtv-1)*Q_T_h+qth,1);
                    M_R_meas = M_R_measMAT((qrv-1)*Q_R_h+qrh, 1);
                    %%%%%%%%%%%%%%%%%% Proposed Algorithm: Cross-Field RD (HSPWM case) %%%%%%%%%%%%%%%%%%
                    % based on "Algorithm 2 Reduced Dictionary Method" in the manuscript and
                    % Sec. III-C3 HSPWM-based channel estimation with RD.
                    if indx_UsedSARx > 1
                        % Detect the best previously estimated candidate support based on the first "reference" T-R SA channel estimation (we added -1 to the index of Rx-plane)
                        Supp_T_Est = AoDInd_SOMP_Cross_cand{(qrv-1)*Q_R_h+qrh -1,(qtv-1)*Q_T_h+qth};                 
                        % Construct the Reduced Dictionary based on the previous support by selecting the most correlated
                        % columns from an oversampled dictionary following Algorithm 2 in the manuscript
                        if indx_UsedSARx == 2
                            Abar_T_Est = Abar_T(:,Supp_T_Est);
                        else
                            Abar_T_Est = Abar_T_RD(:,Supp_T_Est);
                        end
                        Abar_T_RD = get_Proposed_RD_IntermediateField(Abar_T_OVS, Abar_T_Est, G_T_RD);
                        Abar_T_RD_quant = get_Quantized_Dict(Abar_T_RD, Quant_PS_T);
                        % Apply the same method for the receiver side if the MIMO array is massive (per SA: AEs >= 64)
                        % For now keep the receiver without applying the RD 
                        Abar_R_RD = Abar_R;
                        Abar_R_RD_quant = Abar_R_quant;
                        % Compute RD Dictionary Matrix
                        Thetabar_HSPWM_RD = kron(Abar_T_RD, Abar_R_RD);
                    end
                    %%%%%%%%%%%%%%%%%% THz channel between q_T(th) Tx SA and q_R(th) Rx SA %%%%%%%%%%%%%%%%%%
                    H_qRqT = H_AoSA(((qrv-1)*Q_R_h+qrh-1)*Qbar_R_v*Qbar_R_h+1:((qrv-1)*Q_R_h+qrh)*Qbar_R_v*Qbar_R_h,((qtv-1)*Q_T_h+qth-1)*Qbar_T_v*Qbar_T_h+1:((qtv-1)*Q_T_h+qth)*Qbar_T_v*Qbar_T_h,:);
                    % This is H_{q_R,q_T} defined in Eq. (5) in the manuscript (for SWM, Eqs. (7) and (8))
                    %%%%%%%%%%%%%%%%%% Noise Generation %%%%%%%%%%%%%%%%%%
                    N = N_AoSA{indx_UsedSATx,indx_UsedSARx};
                    %%%%%%%%%%%%%%%%%% The vectorized version of the received signal of Eq. (29) for all subcarriers based on Tx and Rx signals %%%%%%%%%%%%%%%%%%
                    y = get_Channel_Output(H_qRqT,N(:,1:M_T_meas,:),Z(:,1:M_T_meas),C(:,1:M_R_meas),K,Tx_pwr);
                    %%%%%%%%%%%%%%%%%% Compute Measurement Matrix %%%%%%%%%%%%%%%%%%
                    Psi = kron(transpose(Z(:,1:M_T_meas)),(C(:,1:M_R_meas))');
                    % following Eq. (32) in the manuscript Psi = kron(Z^Transpose, C^Herm) is the measurement matrix
                    %%%%%%%%%%%%%%%%%% Compute Sensing Matrix %%%%%%%%%%%%%%%%%%
                    % The HSPWM RD sensing matrix 
                    if indx_UsedSARx > 1
                        Upsilon_SOMP_Cross = Psi*Thetabar_HSPWM_RD;
                    else
                        Upsilon_SOMP_Cross = Psi*Thetabar;
                    end
                    %%%%%%%%%%%%%%%%%% Save Sensing Matrix Size for Complexity Analysis %%%%%%%%%%%%%%%%%%
                    Num_Grids_Exact_Cross(indx_UsedSATx, indx_UsedSARx) = numel(Upsilon_SOMP_Cross);
                    %%%%%%%%%%%%%%%%%% Cross Field HSPWM-RD SOMP CS-Estimator %%%%%%%%%%%%%%%%%%
                    % bs here stands for BeamSpace domain
                    % Cross Field HSPWM-RD + SOMP
                    [hest_bs_SOMP_Cross, support_SOMP_Cross] = SOMP(y,Upsilon_SOMP_Cross,Tx_pwr,Lbar,K);
                    if indx_UsedSARx > 1
                        H_bs_SOMP_Cross = reshape(hest_bs_SOMP_Cross,size(Abar_R_RD,2),size(Abar_T_RD,2),K);
                    else
                        H_bs_SOMP_Cross = reshape(hest_bs_SOMP_Cross,size(Abar_R,2),size(Abar_T,2),K);
                    end
                    %%%%%%%%%%%%%%%%%% Extract the number of dominant/best beams %%%%%%%%%%%%%%%%%%
                    % Accumulate until you get 95% of the estimated channel power
                    L_dom_Beam_SOMP = get_NumBestBeams(hest_bs_SOMP_Cross, support_SOMP_Cross);
                    L_dom_Beams_SOMP_UM_Cross((qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth) = L_dom_Beam_SOMP;
                    %%%%%%%%%%%%%%%%%% Extract the support and select angles for beamforming/combining %%%%%%%%%%%%%%%%%%
                    % SOMP Cross RD HSPWM
                    [rowMaxBeamRx_SOMP_Cross,colMaxBeamTx_SOMP_Cross] = ind2sub(size(H_bs_SOMP_Cross),support_SOMP_Cross(1:L_dom_Beam_SOMP));
                    AoDInd_SOMP_Cross_cand{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = colMaxBeamTx_SOMP_Cross;
                    AoAInd_SOMP_Cross_cand{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = rowMaxBeamRx_SOMP_Cross;
                    if indx_UsedSARx > 1
                        F_Hat_SOMP_Cross{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Abar_T_RD_quant(:,colMaxBeamTx_SOMP_Cross(1));
                        W_Hat_SOMP_Cross{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Abar_R_RD_quant(:,rowMaxBeamRx_SOMP_Cross(1));
                    else
                        F_Hat_SOMP_Cross{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Abar_T_quant(:,colMaxBeamTx_SOMP_Cross(1));
                        W_Hat_SOMP_Cross{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Abar_R_quant(:,rowMaxBeamRx_SOMP_Cross(1));
                    end
                    %%%%%%%%%%%%%%%%%% Fill the SA (MIMO channel) estimation into the UM-MIMO channels  %%%%%%%%%%%%%%%%%%
                    H_Est_SOMP_Cross = zeros(Qbar_R,Qbar_T,K);
                    for indx_subc = 1:K
                        % Create Angle-Domain Channels
                        % SOMP Cross RD HSPWM
                        if indx_UsedSARx > 1
                            H_Est_SOMP_Cross(:,:,indx_subc) = Abar_R_RD*H_bs_SOMP_Cross(:,:,indx_subc)*Abar_T_RD.';
                        else
                            H_Est_SOMP_Cross(:,:,indx_subc) = Abar_R*H_bs_SOMP_Cross(:,:,indx_subc)*Abar_T.';
                        end
                        % Fill the UM-MIMO arrays by estimations
                        H_Est_SOMP_UM_Cross{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,indx_subc} = H_Est_SOMP_Cross(:,:,indx_subc);
                    end
                    %%%%%%%%%%%%%%%%%% End of loop on Tx-Rx SA channel Estimation  %%%%%%%%%%%%%%%%%%
                end
            end
            %%%%%%%%%%%%%%%%%% End of loops over UM-MIMO channels  %%%%%%%%%%%%%%%%%%
        case "PWM"
            %%%%%%%%%%%%%%%%%% Proposed Algorithm: Cross-Field RD (PWM case) %%%%%%%%%%%%%%%%%%
            % based on Sec. III-C2 PWM-based channel estimation with RD.
            %
            % Drawing inspiration from the analysis of the far field PWM channel model and the corresponding beamspace representation, 
            % it becomes evident that the channels across different SAs are the same. Thus, we introduce a more efficient
            % approach by estimating only the reference Tx-Rx SA channel, which we refer to as PWM-RD
            %
            % Find min NMSE by concluding the phase-shift correction factor
            Phase_Search_Interval = -pi: 2*pi/StepSize_Phase:pi-pi/StepSize_Phase;
            Required_PhaseShift = complex(zeros(Q_R, Q_T, K));
            for indx_subc = 1:K
                H_11_ref = H_AoSA_Cell{1,1,indx_subc};
                for indx_UsedSATx = 1:Num_of_UsedSATx
                    Min_NMSE_AoSA= zeros(Q_R,StepSize_Phase);
                    for indx_UsedSARx = 1:Num_of_UsedSARx
                        if indx_UsedSARx ==1 && indx_UsedSATx == 1
                            Required_PhaseShift(1,1,indx_subc) = exp(1j*0);
                        else
                            H_xx_AoSA = H_AoSA_Cell{indx_UsedSARx,indx_UsedSATx,indx_subc};
                            for indx_tries = 1:StepSize_Phase
                                Min_NMSE_AoSA(indx_UsedSARx,indx_tries) = ...
                                    10*log10(norm(H_11_ref*exp(1j*Phase_Search_Interval(indx_tries))-H_xx_AoSA, 'fro')/norm(H_xx_AoSA, 'fro'));
                            end
                            [~, Ind_min] = min(Min_NMSE_AoSA(indx_UsedSARx,:));
                            Required_PhaseShift(indx_UsedSARx,indx_UsedSATx,indx_subc) = exp(1j*Phase_Search_Interval(Ind_min));
                        end
                    end
                end
            end
            % Estimation of H_11 only (H reference channel from 1st Tx SA to 1st Rx SA)
            %%%%%%%%%%%%%%%%%% Select the Tx SA %%%%%%%%%%%%%%%%%%
            [qth,qtv] = ind2sub([Q_T_h Q_T_v],Sel_SAs_IndTx(1));
            %%%%%%%%%%%%%%%%%% Select the Rx SA %%%%%%%%%%%%%%%%%%
            [qrh,qrv] = ind2sub([Q_R_h Q_R_v],Sel_SAs_IndRx(1));
            %%%%%%%%%%%%%%%%%% Define the number of measurements %%%%%%%%%%%%%%%%%
            M_T_meas = M_T_measMAT((qtv-1)*Q_T_h+qth, 1);
            M_R_meas = M_R_measMAT((qrv-1)*Q_R_h+qrh, 1);
            %%%%%%%%%%%%%%%%%% RF Beamforming Matrix %%%%%%%%%%%%%%%%%%
            Z = squeeze(Z_AoSA(:,:,1));
            %%%%%%%%%%%%%%%%%%  RF Combining Matrix  %%%%%%%%%%%%%%%%%%
            C = squeeze(C_AoSA(:,:,1));
            %%%%%%%%%%%%%%%%%% Noise Generation %%%%%%%%%%%%%%%%%%
            N = N_AoSA{1,1};
            %%%%%%%%%%%%%%%%%% THz channel between q_T(th) Tx SA and q_R(th) Rx SA %%%%%%%%%%%%%%%%%%
            H_qRqT = H_AoSA(((qrv-1)*Q_R_h+qrh-1)*Qbar_R_v*Qbar_R_h+1:((qrv-1)*Q_R_h+qrh)*Qbar_R_v*Qbar_R_h,((qtv-1)*Q_T_h+qth-1)*Qbar_T_v*Qbar_T_h+1:((qtv-1)*Q_T_h+qth)*Qbar_T_v*Qbar_T_h,:);
            % This is H_{q_R,q_T} defined in Eq. (5) in the manuscript (for SWM, Eqs. (7) and (8))
            %%%%%%%%%%%%%%%%%% The vectorized version of the received signal of Eq. (29) for all subcarriers based on Tx and Rx signals %%%%%%%%%%%%%%%%%%
            y = get_Channel_Output(H_qRqT,N(:,1:M_T_meas,:),Z(:,1:M_T_meas),C(:,1:M_R_meas),K,Tx_pwr);
            %%%%%%%%%%%%%%%%%% Compute Measurement Matrix %%%%%%%%%%%%%%%%%%
            Psi = kron(transpose(Z(:,1:M_T_meas)),(C(:,1:M_R_meas))');
            % following Eq. (32) in the manuscript Psi = kron(Z^Transpose, C^Herm) is the measurement matrix
            %%%%%%%%%%%%%%%%%% Compute Sensing Matrix %%%%%%%%%%%%%%%%%%
            Upsilon_SOMP_Cross = Psi*Thetabar;
            %%%%%%%%%%%%%%%%%% Save Sensing Matrix Size for Complexity Analysis %%%%%%%%%%%%%%%%%%
            Num_Grids_Exact_Cross = numel(Upsilon_SOMP_Cross);
            %%%%%%%%%%%%%%%%%% Cross Field PWM-RD SOMP CS-Estimator %%%%%%%%%%%%%%%%%%
            % bs here stands for BeamSpace domain
            % Cross Field PWM-RD + SOMP: we need only one estimation
            [hest_bs_SOMP_Cross, support_SOMP_Cross] = SOMP(y,Upsilon_SOMP_Cross,Tx_pwr,Lbar,K);
            H_bs_SOMP_Cross = reshape(hest_bs_SOMP_Cross,size(Abar_R,2),size(Abar_T,2),K);
            %%%%%%%%%%%%%%%%%% Extract the number of dominant/best beams %%%%%%%%%%%%%%%%%%
            % Accumulate until you get 95% of the estimated channel power
            L_dom_Beam_SOMP = get_NumBestBeams(hest_bs_SOMP_Cross, support_SOMP_Cross);
            L_dom_Beams_SOMP_UM_Cross((qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth) = L_dom_Beam_SOMP;
            %%%%%%%%%%%%%%%%%% Extract the support and select angles for beamforming/combining %%%%%%%%%%%%%%%%%%
            % SOMP Cross RD PWM
            [rowMaxBeamRx_SOMP_Cross,colMaxBeamTx_SOMP_Cross] = ind2sub(size(H_bs_SOMP_Cross),support_SOMP_Cross(1:Lbar));
            AoDInd_SOMP_Cross_cand{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = colMaxBeamTx_SOMP_Cross;
            AoAInd_SOMP_Cross_cand{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = rowMaxBeamRx_SOMP_Cross;
            F_Hat_SOMP_Cross{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Abar_T_quant(:,colMaxBeamTx_SOMP_Cross(1));
            W_Hat_SOMP_Cross{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = Abar_R_quant(:,rowMaxBeamRx_SOMP_Cross(1));
            %%%%%%%%%%%%%%%%%% Fill the SA (MIMO channel) estimation into the UM-MIMO channels  %%%%%%%%%%%%%%%%%%            
            % Complete the Estimation
            H_Est_SOMP_11_Cross = zeros(Qbar_R,Qbar_T,K);
            for indx_subc = 1:K
                % Create Angle-Domain Channels
                % SOMP Cross RD PWM
                H_Est_SOMP_11_Cross(:,:,indx_subc) = Abar_R*H_bs_SOMP_Cross(:,:,indx_subc)*Abar_T.';
            end
            %%%%%%%%%%%%%%%%%% Begin loops over Tx/Rx SAs %%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%% loop over the Tx SAs %%%%%%%%%%%%%%%%%%
            for indx_UsedSATx = 1:Num_of_UsedSATx
                for indx_UsedSARx = 1:Num_of_UsedSARx
                    [qth,qtv] = ind2sub([Q_T_h Q_T_v],Sel_SAs_IndTx(indx_UsedSATx));
                    [qrh,qrv] = ind2sub([Q_R_h Q_R_v],Sel_SAs_IndRx(indx_UsedSARx));
                    % Fill the UM-MIMO arrays by estimations
                    for indx_subc = 1:K    
                        H_Est_SOMP_UM_Cross{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth,indx_subc} = ...
                            H_Est_SOMP_11_Cross(:,:,indx_subc)*Required_PhaseShift(indx_UsedSARx,indx_UsedSATx,indx_subc);
                    end
                    %%%%%%%%%%%%%%%%%% Fill Precoding and Combining Matrices  %%%%%%%%%%%%%%%%%%
                    F_Hat_SOMP_Cross{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = F_Hat_SOMP_Cross{1,1}*Required_PhaseShift(indx_UsedSARx,indx_UsedSATx,floor(K/2));
                    W_Hat_SOMP_Cross{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth} = W_Hat_SOMP_Cross{1,1}*Required_PhaseShift(indx_UsedSARx,indx_UsedSATx,floor(K/2));
                    %%%%%%%%%%%%%%%%%% End of loop on Tx-Rx SA channel Estimation  %%%%%%%%%%%%%%%%%%
                end
            end
            %%%%%%%%%%%%%%%%%% End of loops over UM-MIMO channels  %%%%%%%%%%%%%%%%%%
        otherwise
            error('We only have three models: SWM, HSPWM, and PWM !');
    end
    %%%%%%%%%%%%%%%%% Performance Evaluation Metrics %%%%%%%%%%%%%%%%%%
    % 1- NMSE Computations for the UM-MIMO system
    for indx_subc = 1:K
        NMSE_SOMP_NF(indx_iter,indx_subc) = ...
            norm(cell2mat(H_Est_SOMP_UM_NF(:,:,indx_subc))-cell2mat(H_AoSA_Cell(:,:,indx_subc)),'fro')^2/norm(cell2mat(H_AoSA_Cell(:,:,indx_subc)),'fro')^2;
       NMSE_SOMP_FF(indx_iter,indx_subc) = ...
            norm(cell2mat(H_Est_SOMP_UM_FF(:,:,indx_subc))-cell2mat(H_AoSA_Cell(:,:,indx_subc)),'fro')^2/norm(cell2mat(H_AoSA_Cell(:,:,indx_subc)),'fro')^2;
        NMSE_SOMP_Hybrid_NF_FF(indx_iter,indx_subc) = ...
            norm(cell2mat(H_Est_SOMP_UM_Hybrid_NF_FF(:,:,indx_subc))-cell2mat(H_AoSA_Cell(:,:,indx_subc)),'fro')^2/norm(cell2mat(H_AoSA_Cell(:,:,indx_subc)),'fro')^2;
        NMSE_SOMP_Hybrid_FF_NF(indx_iter,indx_subc) = ...
            norm(cell2mat(H_Est_SOMP_UM_Hybrid_FF_NF(:,:,indx_subc))-cell2mat(H_AoSA_Cell(:,:,indx_subc)),'fro')^2/norm(cell2mat(H_AoSA_Cell(:,:,indx_subc)),'fro')^2;
        NMSE_SOMP_Cross(indx_iter,indx_subc) = ...
            norm(cell2mat(H_Est_SOMP_UM_Cross(:,:,indx_subc))-cell2mat(H_AoSA_Cell(:,:,indx_subc)),'fro')^2/norm(cell2mat(H_AoSA_Cell(:,:,indx_subc)),'fro')^2;
    end
    % 2- Achievable Rate Computations for the UM-MIMO system
    % Perfect channel state information (CSI)
    F_hat_bmf_PCSI = cell(Q_R, Q_T, K);
    W_hat_cmb_PCSI = cell(Q_R, Q_T, K);
    F_hat_bmf_PCSI(:,:) = {complex(zeros(Qbar_T,1))};
    W_hat_cmb_PCSI(:,:) = {complex(zeros(Qbar_R,1))};
    for indx_UsedSATx = 1:Num_of_UsedSATx
        %%%%%%%%%%%%%%%%%% Select the Tx SA %%%%%%%%%%%%%%%%%%
        [qth,qtv] = ind2sub([Q_T_h Q_T_v],Sel_SAs_IndTx(indx_UsedSATx));
        for indx_UsedSARx = 1:Num_of_UsedSARx
            %%%%%%%%%%%%%%%%%% Select the Rx SA %%%%%%%%%%%%%%%%%%
            [qrh,qrv] = ind2sub([Q_R_h Q_R_v],Sel_SAs_IndRx(indx_UsedSARx));
            for indx_subc = 1:K
                % Search for the best beam focusing vectors
                Dist_SWM_Ch = SWM.Dist_LoS{(qrv-1)*Q_R_h+qrh,(qtv-1)*Q_T_h+qth, indx_subc};
                F_focus = get_NearField_ARV(Qbar_T, Dist_SWM_Ch(floor(Qbar_R/2),:).',Dist_val,p_ch.lambda_kc(indx_subc));
                F_hat_bmf_PCSI{indx_UsedSARx,indx_UsedSATx,indx_subc} = conj(F_focus);
                W_focus = get_NearField_ARV(Qbar_R,Dist_SWM_Ch(:,floor(Qbar_T)/2),Dist_val,p_ch.lambda_kc(indx_subc));
                W_hat_cmb_PCSI{indx_UsedSARx,indx_UsedSATx,indx_subc} = W_focus;
            end
        end
    end
    AR_PCSI(indx_iter,:) = get_AchievableRate(p_ch, K, DataTx_pwr, noise_pwr, H_AoSA_Cell, H_AoSA_Cell, F_hat_bmf_PCSI, W_hat_cmb_PCSI);
    % SOMP NF
    F_hat_bmf_SOMP_EstCh_NF = repmat(cellfun(@conj, F_Hat_SOMP_NF, 'UniformOutput' ,false),[1 1 K]);
    W_hat_cmb_SOMP_EstCh_NF = repmat(W_Hat_SOMP_NF,[1 1 K]);
    AR_SOMP_NF(indx_iter,:) = get_AchievableRate(p_ch, K, DataTx_pwr, noise_pwr, H_Est_SOMP_UM_NF, H_AoSA_Cell, F_hat_bmf_SOMP_EstCh_NF, W_hat_cmb_SOMP_EstCh_NF);
    % SOMP FF
    F_hat_bmf_SOMP_EstCh_FF = repmat(cellfun(@conj, F_Hat_SOMP_FF, 'UniformOutput' ,false),[1 1 K]);
    W_hat_cmb_SOMP_EstCh_FF = repmat(W_Hat_SOMP_FF,[1 1 K]);
    AR_SOMP_FF(indx_iter,:) = get_AchievableRate(p_ch, K, DataTx_pwr, noise_pwr, H_Est_SOMP_UM_FF, H_AoSA_Cell, F_hat_bmf_SOMP_EstCh_FF, W_hat_cmb_SOMP_EstCh_FF);
    % SOMP Hybrid NF then FF
    F_hat_bmf_SOMP_EstCh_Hybrid_NF_FF = repmat(cellfun(@conj, F_Hat_SOMP_Hybrid_NF_FF, 'UniformOutput' ,false),[1 1 K]);
    W_hat_cmb_SOMP_EstCh_Hybrid_NF_FF = repmat(W_Hat_SOMP_Hybrid_NF_FF,[1 1 K]);
    AR_SOMP_Hybrid_NF_FF(indx_iter,:) = get_AchievableRate(p_ch, K, DataTx_pwr, noise_pwr, H_Est_SOMP_UM_Hybrid_NF_FF, H_AoSA_Cell, F_hat_bmf_SOMP_EstCh_Hybrid_NF_FF, W_hat_cmb_SOMP_EstCh_Hybrid_NF_FF);
    % SOMP Hybrid FF then NF
    F_hat_bmf_SOMP_EstCh_Hybrid_FF_NF = repmat(cellfun(@conj, F_Hat_SOMP_Hybrid_FF_NF, 'UniformOutput' ,false),[1 1 K]);
    W_hat_cmb_SOMP_EstCh_Hybrid_FF_NF = repmat(W_Hat_SOMP_Hybrid_FF_NF,[1 1 K]);
    AR_SOMP_Hybrid_FF_NF(indx_iter,:) = get_AchievableRate(p_ch, K, DataTx_pwr, noise_pwr, H_Est_SOMP_UM_Hybrid_FF_NF, H_AoSA_Cell, F_hat_bmf_SOMP_EstCh_Hybrid_FF_NF, W_hat_cmb_SOMP_EstCh_Hybrid_FF_NF);
    % SOMP Cross-field RD
    F_hat_bmf_SOMP_EstCh_Cross = repmat(cellfun(@conj, F_Hat_SOMP_Cross, 'UniformOutput' ,false),[1 1 K]);
    W_hat_cmb_SOMP_EstCh_Cross = repmat(W_Hat_SOMP_Cross,[1 1 K]);
    AR_SOMP_Cross(indx_iter,:) = get_AchievableRate(p_ch, K, DataTx_pwr, noise_pwr, H_Est_SOMP_UM_Cross, H_AoSA_Cell, F_hat_bmf_SOMP_EstCh_Cross, W_hat_cmb_SOMP_EstCh_Cross);
    % 3- Complexity Computations for the UM-MIMO
    Complexity_NF(indx_iter,1) = Lbar*K*sum(Num_Grids_Exact_NF,'all');
    Complexity_FF(indx_iter,1) = Lbar*K*sum(Num_Grids_Exact_FF,'all');
    Complexity_Cross(indx_iter,1) = Lbar*K*sum(Num_Grids_Exact_Cross,'all');
end
% Keep track of simulation time
Sim_Duration = toc;
hr = floor(Sim_Duration/3600);mint = floor((Sim_Duration - hr*3600)/60);sec = Sim_Duration - hr*3600 - mint*60;
fprintf('The simulation time is: %d hr %d min %f sec\n',hr,mint,sec);
%% Online Results
if SaveSimulationResults
    cd Res_Online
    simulationname = strcat('OnEst',num2str(Dist_val),'m',num2str(BetaDot_val),'RotAng',num2str(E),'Iter',...
        'L',num2str(p_ch.L),'Paths','SNR',num2str(RequiredSNR),'dB',num2str(Q_T_v),'by',num2str(Q_T_h),'x',num2str(Q_R_v),'by',num2str(Q_R_h),'SAs',...
        num2str(Qbar_T_v),'by',num2str(Qbar_T_h),'x',num2str(Qbar_R_v),'by',num2str(Qbar_R_h),'AEs',...
        num2str(M_T),'by',num2str(M_R),'T_R_meas',num2str(p_ch.Sep_factTx),'TxSep', num2str(p_ch.Sep_factRx),'RxSep');
    filename=strcat(simulationname,'.mat');
    save(sprintf('%s',filename), 'Dist_val', 'BetaDot_val', 'E', 'p_ch', 'Tx_pwr', 'noise_pwr', 'RequiredSNR', 'M_T', 'M_R', ...
        'NMSE_SOMP_NF', 'NMSE_SOMP_FF', 'NMSE_SOMP_Hybrid_NF_FF', 'NMSE_SOMP_Hybrid_FF_NF', 'NMSE_SOMP_Cross', ...
        'AR_PCSI', 'AR_SOMP_NF', 'AR_SOMP_FF', 'AR_SOMP_Hybrid_NF_FF', 'AR_SOMP_Hybrid_FF_NF', 'AR_SOMP_Cross', ...
        'Eta_Online', 'SelectedModel', 'Complexity_NF', 'Complexity_FF', 'Complexity_Cross', '-v7.3');
    cd ..
end
end