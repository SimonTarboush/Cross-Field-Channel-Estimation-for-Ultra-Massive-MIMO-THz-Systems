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
% This function generates the offline training curves/values for the proposed model selection metric
% (called Eta in the manuscript, Eq (31) in Sec. III-B.1) to solve the cross-field problem
% The metric "Eta" is a function of many parameters: communication distance, array rotation, Rx SNR, and number of measurements
% Inside the function, you can control the required Rx SNR and number of measurements, but as the implementation is performed in a parallel way,
% we can compute multiple values in the same run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% Dist_val: Communication distance (m)
% BetaDot_val: Rotation angle (degree), for current implementation of the ULA, the rotation is around the y-axis (Sec II-B)
% Output Arguments:
% .mat file that contains the training values of the proposed model selection metric "Eta".
% After running this function over multiple communication distances and array rotations,
% we can use the function "Fig7_OfflineTraining" to get Fig. 7 of the manuscript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function get_ProposedMetric_OfflineTraining_Parallel(Dist_val, BetaDot_val)
%% Initialize Channel Parameters
p_ch = generate_channel_param_TIV();
p_ch.L = 3; % Number of channel paths: L=1 LoS only, you can use L= 3, 5, ...
p_ch.Sep_factTx = 2; % 1 for scenario 1 and 2 for scenario 2
p_ch.Sep_factRx = 10; % 1 for scenario 1 and 10 for scenario 2
% To change the aperture as per the paper and get Figs. 7c and 7d use
% p_ch.Sep_factTx = 2; p_ch.Sep_factRx = 10
p_ch = update_channel_param_TIV(p_ch);
%% Calculation of Absorption Coefficient 
K_abs = get_Abs_Coef(p_ch);
%% Main Simulation Parameters
RequiredSNR = [0 6 8 14 16]; % This can be a vector to perform multiple simulation at every run
LenSNR = length(RequiredSNR);
E_offline = 10; % 45 is the used value in the manuscript % Number of offline trials/runs/iterations
SaveSimulationResults = 1; % Decide to save or skip the simulation file
%% Parameters related to the transmitter and receiver uniform linear array
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
M_T = [Qbar_T/2 Qbar_T/4 Qbar_T/8];  % Number of exhaustive-search training Tx (Full Training)/(Partial Training)
M_R = [Qbar_R Qbar_R Qbar_R];            % Number of exhaustive-search training Rx (Full Training)/(Partial Training)
Compression_RatioTx = 1; Compression_RatioRx = 1;
Sel_SAs_IndTx = 1:Q_T;Sel_SAs_IndRx = 1:Q_R;
Num_of_UsedSATx = length(Sel_SAs_IndTx);Num_of_UsedSARx = length(Sel_SAs_IndRx);
MeasTx_vec = length(M_T);
MeasRx_vec = length(M_R);
M_T_measMAT = zeros(Q_T,MeasTx_vec);  % The number of measurements at Tx
M_R_measMAT = zeros(Q_R,MeasRx_vec); % The number of measurements at Rx
M_T_measMAT(Sel_SAs_IndTx,:) = repmat(Compression_RatioTx*M_T,Num_of_UsedSATx,1); % The number of measurements at Tx
M_R_measMAT(Sel_SAs_IndRx,:) = repmat(Compression_RatioRx*M_R,Num_of_UsedSARx,1); % The number of measurements at Rx
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
%% Initialize Arrays for Results: Proposed metric
if Q_R>1
    Num_Poss_EtaComb = nchoosek(Q_R,2);
else
    Num_Poss_EtaComb =1;
end
Eta_Offline = zeros(LenSNR, MeasTx_vec, E_offline, Num_Poss_EtaComb+1);
% Based on Eq. (31) in the manuscript, we need to save the max value hence, we add one to the possible number of the proposed metric values "Num_Poss_EtaComb"
%% Main Loop
tic;
noise_pwr = N0_sc;
% Update channel parameters
p_ch.positionRx = [Dist_val; 0; 0];
p_ch.eulerRx = [rad2deg(pi); BetaDot_val; 0];
p_ch = update_channel_param_TIV(p_ch);
parfor indx_iter = 1:E_offline
    % AoSAs THz channel generation
    SWM = channel_TIV(p_ch, K_abs);
    % The spherical wave model (SWM) channel model is used in all simulations as it is the most accurate model for all communication distances
    H_AoSA =cell2mat(SWM.H);
    Chi_Offline_TmpCell = cell(LenSNR, MeasTx_vec);
    for indx_tr_pwr = 1:LenSNR
        for indx_meas = 1:MeasTx_vec
            Chi_Offline_TmpCell{indx_tr_pwr,indx_meas} = zeros(M_R(indx_meas)*M_T(indx_meas)*K,Q_R);
        end
    end
    % Start training only from the 1st "reference" SA. The received signals only from the reference Tx SA as it is sufficient for determining the region
    indx_UsedSATx = 1;
    %%%%%%%%%%%%%%%%%% Select the Tx SA %%%%%%%%%%%%%%%%%%
    [qth,qtv] = ind2sub([Q_T_h Q_T_v],Sel_SAs_IndTx(indx_UsedSATx));
    %%%%%%%%%%%%%%%%%% RF Beamforming Matrix %%%%%%%%%%%%%%%%%%
    Z = get_RandomCodebook(Quant_PS_T, Q_T_quant, Qbar_T, Qbar_T); % Random beam training
    %%%%%%%%%%%%%%%%%%  RF Combining Matrix  %%%%%%%%%%%%%%%%%%
    % following our proposal in Sec. III-A, the combining weights are applied to all Rx SAs
    C = get_RandomCodebook(Quant_PS_R, Q_R_quant, Qbar_R, Qbar_R); % Random beam training
    for indx_UsedSARx = 1:Num_of_UsedSARx
        %%%%%%%%%%%%%%%%%% Select the Rx SA %%%%%%%%%%%%%%%%%%
        [qrh,qrv] = ind2sub([Q_R_h Q_R_v],Sel_SAs_IndRx(indx_UsedSARx));
        %%%%%%%%%%%%%%%%%% THz channel between q_T(th) Tx SA and q_R(th) Rx SA %%%%%%%%%%%%%%%%%%
        H_qRqT = H_AoSA(((qrv-1)*Q_R_h+qrh-1)*Qbar_R_v*Qbar_R_h+1:((qrv-1)*Q_R_h+qrh)*Qbar_R_v*Qbar_R_h,((qtv-1)*Q_T_h+qth-1)*Qbar_T_v*Qbar_T_h+1:((qtv-1)*Q_T_h+qth)*Qbar_T_v*Qbar_T_h,:);
        % This is H_{q_R,q_T} defined in Eq. (5) in the manuscript (for SWM, Eqs. (7) and (8))
        %%%%%%%%%%%%%%%%%% Noise Generation %%%%%%%%%%%%%%%%%%
        N = sqrt(K)*sqrt(noise_pwr)/sqrt(2)*(randn(Qbar_R,Qbar_T,K)+1j*randn(Qbar_R,Qbar_T,K));
        %%%%%%%%%%%%%%%%%% Define the number of measurements %%%%%%%%%%%%%%%%%%
        for indx_tr_pwr = 1:LenSNR
            Tx_pwr = Tx_power_tot(indx_tr_pwr)/Ns_train;
            % Tx_pwr is "P_T" defined in Eq. 29 in the manuscript, i.e., is the total Tx power used per transmission during the training phase.
            % during the training there is only one SA to be used, i.e., number of streams during training is Ns_train = 1
            for indx_meas = 1:MeasTx_vec
                Tx_meas_indx = indx_meas;Rx_meas_indx = indx_meas;
                M_T_meas = M_T_measMAT((qtv-1)*Q_T_h+qth,Tx_meas_indx);
                M_R_meas = M_R_measMAT((qrv-1)*Q_R_h+qrh, Rx_meas_indx);
                %%%%%%%%%%%%%%%%%% The vectorized version of the received signal of Eq. (29) for all subcarriers based on Tx and Rx signals %%%%%%%%%%%%%%%%%%
                Chibar_qr_1 = get_Channel_Output(H_qRqT,N(:,1:M_T_meas,:),Z(:,1:M_T_meas),C(:,1:M_R_meas),K,Tx_pwr);
                %%%%%%%%%%%%%%%%%% Compute the normalized vector \boldsymbol{\chi} used later in Eq. (30) %%%%%%%%%%%%%%%%%%
                Chi_Offline_TmpCell{indx_tr_pwr,indx_meas}(:,indx_UsedSARx) = abs(Chibar_qr_1(:))/sqrt(sum(abs(Chibar_qr_1(:)).^2));
            end
        end
    end
    %%%%%%%%%%%%%%%%%% Compute all possible Combinations for the proposed model selection metric using Eqs. (30) and (31) %%%%%%%%%%%%%%%%%%
    for indx_tr_pwr = 1:LenSNR
        for indx_meas = 1:MeasTx_vec
            indx_cnt = 0;
            Chi_Offline_r_c=Chi_Offline_TmpCell{indx_tr_pwr,indx_meas};
            Eta_Offline_r_c = zeros(Num_Poss_EtaComb,1);
            for indx_pwr_row = 1:Q_R-1
                for indx_pwr_col = indx_pwr_row+1:Q_R
                    indx_cnt = indx_cnt + 1;
                    Chi_diff = Chi_Offline_r_c(:,indx_pwr_row)-Chi_Offline_r_c(:,indx_pwr_col);
                    % Eq. (30) in the manuscript
                    Eta_Offline_r_c(indx_cnt,1) = sum(abs(Chi_diff).^2);
                end
            end
            % Eq. (30) and (31) in the manuscript, the max value (last element) is the proposed metric, however, we keep all possible values for analysis purposes
            Eta_Offline(indx_tr_pwr,indx_meas,indx_iter,:) =  [Eta_Offline_r_c(1:Num_Poss_EtaComb); max(Eta_Offline_r_c(1:Num_Poss_EtaComb))];
        end
    end    
end
% Keep track of simulation time
Sim_Duration = toc;
hr = floor(Sim_Duration/3600);mint = floor((Sim_Duration - hr*3600)/60);sec = Sim_Duration - hr*3600 - mint*60;
fprintf('The simulation time is: %d hr %d min %f sec\n',hr,mint,sec);
%% Offline Results
% Proposed metric for determining the region based on Rx power difference across SAs
Eta_MAT_Max = Eta_Offline(:,:,:,end);
if SaveSimulationResults
    cd Res_Offline
    simulationname = strcat('OffTr',num2str(Dist_val),'m',num2str(BetaDot_val),'RotAng',num2str(E_offline),'Iter',...
        'L',num2str(p_ch.L),'Paths',num2str(Q_T_v),'by',num2str(Q_T_h),'x',num2str(Q_R_v),'by',num2str(Q_R_h),'SAs',...
        num2str(Qbar_T_v),'by',num2str(Qbar_T_h),'x',num2str(Qbar_R_v),'by',num2str(Qbar_R_h),'AEs',...
        num2str(p_ch.Sep_factTx),'TxSep', num2str(p_ch.Sep_factRx),'RxSep');
    filename=strcat(simulationname,'.mat');
    save(sprintf('%s',filename),'-v7.3');
    cd ..
end
end