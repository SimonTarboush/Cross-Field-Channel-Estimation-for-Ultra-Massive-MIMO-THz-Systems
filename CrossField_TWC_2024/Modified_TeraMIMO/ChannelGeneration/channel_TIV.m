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
% This function generates a time-invariant THz channel
% This is a modified version of the original code:
% https://github.com/Simon-Tarboush/TeraMIMO || https://github.com/hasarieddeen/TeraMIMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% p: Channel struct that contains the channel parameters.
% K_abs: Molecular absorption coefficient, a matrix of size(p.Nsub_c, p.Nsub_b),
%              (number of subcarriers, number of subbands @ each subcarrier)
% Output Arguments:
% CH_Response: A struct contains the channel response, CH_Response.H: H(f) time-invariant frequency domain response, 
%                          the LoS only channel, the NLoS channel, and other important variables such as the Tx/Rx array response vectors and the magnitude of channel paths 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CH_Response = channel_TIV(p, K_abs)
% Initialize parameters and channnel matrix
num_subcarries = p.nFreq(1);
num_freqs_per_subcarr = p.nFreq(2);
H = cell(p.Rx_AoSA.Q,p.Tx_AoSA.Q,num_subcarries);
H_LoS = cell(p.Rx_AoSA.Q,p.Tx_AoSA.Q,num_subcarries);
H_NLoS = cell(p.Rx_AoSA.Q,p.Tx_AoSA.Q,num_subcarries);

B_T_tot = cell(p.Rx_AoSA.Q,p.Tx_AoSA.Q,num_subcarries);
B_T_tot(:,:,:) = {complex(zeros(p.Tx_AoSA.Qbar,p.L))};
B_R_tot = cell(p.Rx_AoSA.Q,p.Tx_AoSA.Q,num_subcarries);
B_R_tot(:,:,:) = {complex(zeros(p.Rx_AoSA.Qbar,p.L))};
AoD_LoS_loc = cell(p.Rx_AoSA.Q,p.Tx_AoSA.Q,num_subcarries);
AoA_LoS_loc = cell(p.Rx_AoSA.Q,p.Tx_AoSA.Q,num_subcarries);
Dist_LoS = cell(p.Rx_AoSA.Q,p.Tx_AoSA.Q,num_subcarries);
Alpha_LoS = cell(p.Rx_AoSA.Q,p.Tx_AoSA.Q,num_subcarries);
Alpha_NLoS = cell(p.Rx_AoSA.Q,p.Tx_AoSA.Q,num_subcarries,p.L);
Gamma = cell(p.Rx_AoSA.Q,p.Tx_AoSA.Q,num_subcarries,p.L);

% Main Loop
% Loop over subcarries
for indx_subc = 1:num_subcarries
    % Loops over SA MIMO arrays    
    for mr = 1:p.Rx_AoSA.Qdim(1)
        for nr = 1:p.Rx_AoSA.Qdim(2)
            for mt = 1:p.Tx_AoSA.Qdim(1)
                for nt = 1:p.Tx_AoSA.Qdim(2)
                    
                    Tx_indx = (mt-1)*p.Tx_AoSA.Qdim(2)+nt;
                    Rx_indx = mr*p.Rx_AoSA.Qdim(2)-nr+1;
                   
                    [d_SAAEs, AoD_SV_LOS_loc, AoA_SV_LOS_loc] = get_Distance_Angle_LoS(p,mr,nr,mt,nt);
                    [AlphaLoS,AlphaLoSGain] = get_PathLoss(p, indx_subc, d_SAAEs, K_abs);
                    
                    AoD_LoS_loc{Rx_indx,Tx_indx,indx_subc} = AoD_SV_LOS_loc;
                    AoA_LoS_loc{Rx_indx,Tx_indx,indx_subc} = AoA_SV_LOS_loc;
                    Dist_LoS{Rx_indx,Tx_indx,indx_subc} = d_SAAEs;
                    Alpha_LoS{Rx_indx,Tx_indx,indx_subc} = AlphaLoS;
                    Alpha_NLoS{Rx_indx,Tx_indx,indx_subc,1} = AlphaLoS;
                    
                    Gamma{Rx_indx,Tx_indx,indx_subc,1} = 1;
                    Gammatmp = zeros(p.L,1);Gammatmp(1) = 1;
                    AlphaNLoS = zeros([size(AlphaLoSGain),p.L]);AlphaNLoS(:,:,1) = AlphaLoSGain;
                    HLoStmp = zeros(size(AlphaLoS));HNLoStmp = zeros(size(AlphaLoS));
                    
                    if strcmp(p.channelType,'LoS') || strcmp(p.channelType,'Multipath+LoS')
                        switch p.AntennaModel
                            % Frequency Independent Model
                            case 'Isotropic'
                                if strcmp(p.WaveModelSA,'Sphere') && strcmp(p.WaveModelAE,'Sphere')
                                    Gt_LOS = p.Tx.Gain*ones(p.Rx_AoSA.Qbar,p.Tx_AoSA.Qbar);
                                    Gr_LOS = p.Rx.Gain*ones(p.Rx_AoSA.Qbar,p.Tx_AoSA.Qbar);
                                else
                                    Gt_LOS = p.Tx.Gain;
                                    Gr_LOS = p.Rx.Gain;
                                    Gt_NLOS = Gt_LOS;
                                    Gr_NLOS = Gr_LOS;
                                end
                            case 'Directional'
                                % todo: add get_SectorGain with check sizes
                                % Gt_LOS = get_SectorGain(p, AoD_SV_LOS, 'LoS');
                                % Gr_LOS = get_SectorGain(p, AoA_SV_LOS, 'LoS');
                            otherwise
                                error('This Antenna Model isn''t implemented !!');
                        end
                        if strcmp(p.WaveModelAE,'Plane')
                            a_sv_t_LoS = get_ArrayResponse(p,indx_subc,p.Tx_AoSA.Qbardim(1),p.Tx_AoSA.Qbardim(2),AoD_SV_LOS_loc,'SV','T');
                            a_sv_t_LoS_tmp = reshape(a_sv_t_LoS,p.Tx_AoSA.Qbar,1,num_freqs_per_subcarr);
                            a_sv_r_LoS = get_ArrayResponse(p,indx_subc,p.Rx_AoSA.Qbardim(1),p.Rx_AoSA.Qbardim(2),AoA_SV_LOS_loc,'SV','R');
                            a_sv_r_LoS_tmp = reshape(a_sv_r_LoS,p.Rx_AoSA.Qbar,1,num_freqs_per_subcarr);
                            H_LoS_tmp = zeros(p.Rx_AoSA.Qbar,p.Tx_AoSA.Qbar, num_freqs_per_subcarr);
                            for indx_numfreqpersubc = 1:num_freqs_per_subcarr
                                H_LoS_tmp(:,:,indx_numfreqpersubc) = Gr_LOS*Gt_LOS*AlphaLoS(:,:,indx_numfreqpersubc)*a_sv_r_LoS_tmp(:,:,indx_numfreqpersubc)*a_sv_t_LoS_tmp(:,:,indx_numfreqpersubc).';
                            end
                            % LoS Channel, Frequency Domain Implementation
                            H_LoS{Rx_indx,Tx_indx,indx_subc} = sqrt(p.Rx_AoSA.Qbar*p.Tx_AoSA.Qbar)*H_LoS_tmp;
                            B_T_tot{Rx_indx,Tx_indx,indx_subc}(:,1) = a_sv_t_LoS_tmp;
                            B_R_tot{Rx_indx,Tx_indx,indx_subc}(:,1) = a_sv_r_LoS_tmp;
                        elseif strcmp(p.WaveModelAE,'Sphere')
                            HLoStmp = repmat(Gt_LOS,1,1,num_freqs_per_subcarr).*AlphaLoS.*repmat(Gr_LOS,1,1,num_freqs_per_subcarr);
                            H_LoS{Rx_indx,Tx_indx,indx_subc} =HLoStmp;
                            Dist_SWM_Ch = Dist_LoS{Rx_indx,Tx_indx,indx_subc};
                            % This is just an approximation for getting the near-field array response vector valid for
                            % some communication distances (not always accurate: example for short distances and UM-MIMO at both ends)
                            bt_LoS = get_Nearfield_manifold(p.Tx_AoSA.Qbar, Dist_SWM_Ch(floor(p.Rx_AoSA.Qbar/2),:).',p.d_tx_rx,p.lambda_kc(indx_subc));
                            br_LoS = get_Nearfield_manifold(p.Rx_AoSA.Qbar,Dist_SWM_Ch(:,floor(p.Tx_AoSA.Qbar)/2),p.d_tx_rx,p.lambda_kc(indx_subc));
                            B_T_tot{Rx_indx,Tx_indx,indx_subc}(:,1) = bt_LoS;
                            B_R_tot{Rx_indx,Tx_indx,indx_subc}(:,1) = br_LoS;
                        else
                            % AE SWM need to add SV and/or BF
                            error('This wavemodel for AE level isn''t implemented. Supported Options are: SWM/PWM');
                        end
                        if p.L > 1
                            for ell = 2:p.L
                                gamma = (cos(p.varphi_in(ell)) - p.kappa_T*cos(p.varphi_ref(ell))) / (cos(p.varphi_in(ell)) + p.kappa_T*cos(p.varphi_ref(ell)));
                                rho = exp(-8*(pi^2)*(p.freq(indx_subc,1)^2)*(p.sigma_rough^2)*(cos(p.varphi_in(ell))^2) / (p.c^2));
                                Gammatmp(ell) = gamma*rho;   % Eq (17) in the manuscript
                                Gamma{Rx_indx,Tx_indx,indx_subc,ell} = Gammatmp(ell);
                                d_NLoS  = d_SAAEs + (p.d_ell(ell)-p.d_tx_rx);
                                [~,AlphaNLoSGain] = get_PathLoss(p, indx_subc, d_NLoS, K_abs);
                                AlphaNLoS(:,:,ell) = abs(Gammatmp(ell))*AlphaNLoSGain*exp(1j*2*pi*rand(1));  % Eq (18) in the manuscript
                                Alpha_NLoS{Rx_indx,Tx_indx,indx_subc,ell} = AlphaNLoS(:,:,ell);
                                % For NLoS, L>1, the approximation of using two near-field array response vector is valid 
                                bt_NLoS = get_ArrayResponse_Spherical(p,indx_subc,mt,nt,p.d_ell_Tx(ell),p.theta_TxNLoS(ell),'SV','T');
                                br_NLoS = get_ArrayResponse_Spherical(p,indx_subc,mr,nr,p.d_ell_Rx(ell),p.phi_RxNLoS(ell),'SV','R');
                                B_T_tot{Rx_indx,Tx_indx,indx_subc}(:,ell) = bt_NLoS;
                                B_R_tot{Rx_indx,Tx_indx,indx_subc}(:,ell) = br_NLoS;
                                b_tot = br_NLoS*bt_NLoS';
                                HNLoStmp = HNLoStmp + AlphaNLoS(:,:,ell).*b_tot*exp(-1j*2*pi*p.freq(indx_subc,:)*p.tau_LoS_NLoS(ell));
                            end
                            H_NLoS{Rx_indx,Tx_indx,indx_subc} = HNLoStmp;
                        end
                        H{Rx_indx,Tx_indx,indx_subc} = (HLoStmp + HNLoStmp); % Eq (10) in the manuscript but the LoS is accurate here based on Eq (8)
                    end
                end
            end
        end
    end
end
% Channel Response
switch p.channelType
    case 'LoS'
        CH_Response.H = H;
        CH_Response.H_LoS = H_LoS;
        CH_Response.H_NLoS = H_NLoS;
        CH_Response.AT_tot = B_T_tot;
        CH_Response.AR_tot = B_R_tot;
        CH_Response.AoD_LoS_loc = AoD_LoS_loc;
        CH_Response.AoA_LoS_loc = AoA_LoS_loc;
        CH_Response.Dist_LoS = Dist_LoS;
        CH_Response.Alpha_LoS = Alpha_LoS;
        CH_Response.Gamma = Gamma;
        CH_Response.Alpha_NLoS = Alpha_NLoS;
    otherwise
        error('This channel configuration is still not implemented in this version');
end
end