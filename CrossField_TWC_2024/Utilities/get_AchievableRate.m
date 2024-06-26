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
% This function computes the achievable rate following Eq. (38) in the manuscript 
% and considering no pilot overhead (\varrho = 1) i.e. no loss in achievable rate due to the training
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% p: Channel struct that contains main simulation parameters 
% K: Number of subcarriers
% Tx_pwr_tot: The pilot transmission power
% noise_pwr: Noise power
% H_UM_Est: UM-MIMO AoSA channel estimation
% H_UM_Perf: UM-MIMO AoSA perfect channel
% F_hat_bmf: Beamforming vectors of UM-MIMO AoSA based on the compressed sensing channel estimation
% W_hat_cmb: Combining vectors of UM-MIMO AoSA based on the compressed sensing channel estimation
% Output Arguments:
% AR : Achievable rate vector (size of 1xK)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AR = get_AchievableRate(p, K, Tx_pwr_tot, noise_pwr, H_UM_Est, H_UM_Perf, F_hat_bmf, W_hat_cmb)
AR = zeros(1,K);
Q_T = p.Tx_AoSA.Q;
Q_R = p.Rx_AoSA.Q;  
Qbar_T = p.Tx_AoSA.Qbar;
Qbar_R = p.Rx_AoSA.Qbar;
% The number of spatial streams Ns
Ns_data = min(Q_T,Q_R);
% min(N_T_tot,N_R_tot);% (Ns = N_RF)  <= min(N_T_tot,N_R_tot);
N_RF_T = Ns_data; % will be equal to the number of Tx SAs during data transmission stage
N_RF_R = Ns_data; % will be equal to the number of Rx SAs during data transmission stage
% Number of SAs and AEs
Q_T_v = p.Tx_AoSA.Qdim(1);    % Number of transmit antennas on the z-axis of planar array
Q_T_h = p.Tx_AoSA.Qdim(2);    % Number of transmit antennas on the y-axis of planar array
Q_R_v = p.Rx_AoSA.Qdim(1);    % Number of receiver antennas on the z-axis of planar array
Q_R_h = p.Rx_AoSA.Qdim(2);    % Number of receiver antennas on the y-axis of planar array
% Used SAs
Sel_SAs_IndTx = 1:Q_T;
Sel_SAs_IndRx = 1:Q_R;
% Cells to save RF precoders and combiners 
F_RFtmp = cell(1,Ns_data,K);
W_RFtmp = cell(1,Ns_data,K);
F_RFtmp(:,:) = {complex(zeros(Qbar_T,1))};
W_RFtmp(:,:) = {complex(zeros(Qbar_R,1))};
% Loop over the subcarriers
for indx_subc = 1:K
    H_AoSAs = cell2mat(H_UM_Est(:,:,indx_subc));
    H_AoSAs_Perf = cell2mat(H_UM_Perf(:,:,indx_subc));
    for indx_strm = 1:Ns_data
        % Detect only the SAs facing each others
        indx_UsedSATx = indx_strm;
        indx_UsedSARx = indx_strm;
        [qth, qtv] = ind2sub([Q_T_h Q_T_v], Sel_SAs_IndTx(indx_UsedSATx));
        [qrh, qrv] = ind2sub([Q_R_h Q_R_v], Sel_SAs_IndRx(indx_UsedSARx));
        Tx_indx = (qtv-1)*Q_T_h+qth;
        Rx_indx = qrv*Q_R_h-qrh+1;
        F_RFtmp{1,indx_strm,indx_subc} = F_hat_bmf{Rx_indx,Tx_indx,indx_subc};
        W_RFtmp{1,indx_strm,indx_subc} = W_hat_cmb{Rx_indx,Tx_indx,indx_subc};
    end
    % RF precoding and combining for AoSA
    F_RF_tot = blkdiag(F_RFtmp{:,:,indx_subc});
    W_RF_tot = blkdiag(W_RFtmp{:,:,indx_subc});
    H_BB_to_BB = W_RF_tot'*H_AoSAs*F_RF_tot;
    % Extract BB digital precoding and combining based on SVD
    [U_svd_BB_to_BB, ~, V_svd_BB_to_BB] = svd(H_BB_to_BB);
    % Get the baseband precoder and combiner
    F_BB_opt = V_svd_BB_to_BB(:,1:N_RF_T);
    W_BB_opt = U_svd_BB_to_BB(:,1:N_RF_R);
    % Follow the power constraints
    W_BB_opt = W_BB_opt/norm(W_BB_opt,'fro');
    % Form the Hybrid precoding and combining matrices (RF analog + BB digital)
    F_est = F_RF_tot*F_BB_opt;
    W_est = W_RF_tot*W_BB_opt;
    % Compute Noise Power after combining (Matrix R in Eq. (38) in the manuscript)
    R_n = noise_pwr*(W_est'*W_est);
    % Compute Effective channel after hybrid precoding and combining (Matrix C in Eq. (38) in the manuscript)
    C = W_est'*H_AoSAs_Perf*F_est;
    % Compute the achievable rate, check Eq. (38) in the manuscript (but without averaging and using varrho)
    AR(1,indx_subc) = abs(log2(det(eye(N_RF_R) + (Tx_pwr_tot/(K*Ns_data)*(R_n^-1)*(C*C')))));
end
end