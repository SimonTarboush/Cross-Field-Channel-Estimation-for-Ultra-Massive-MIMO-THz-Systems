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
% This function computes the distances and local AoD/AoA. This function supports the models: 
%         1) Hybrid Spherical Planar Wave Model (HSPWM): i.e. Spherical Wave Model (SWM)/Plane Wave Model (PWM) (SWM/PWM) on the level of SA/AE.
%         2) Spherical Wave Model (SWM): i.e. Spherical Wave Model (SWM) on the level of SA and AE.
% The other models are implemented in the original code of TeraMIMO
% This is a modified version of the original code:
% https://github.com/Simon-Tarboush/TeraMIMO || https://github.com/hasarieddeen/TeraMIMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% p: Channel struct that contains the channel parameters
% mr: Index of the number of the receiver SAs (rows)
% nr:  Index of the number of the receiver SAs (columns)
% mt: Index of the number of the transmitter SAs (rows)
% nt:  Index of the number of the transmitter SAs (columns)
% Output Arguments:
% Dist_SAAEs: Distance between Tx-Rx SAs/AEs:
%                     1) HSPWM (SWM-SA/PWM-AE): we need the distance between SAs
%                     2) SWM: Distance Matrix between all Tx-Rx AEs pair inside a SA of size (Qbar_R, Qbar_T) 
% AoD_SV_LOS_loc: Local angle-of-departure for LoS beamsteering vector; size is based on the used model for the SA/AE:
%                     1) HSPWM (SWM-SA/PWM-AE): vector of size (2, 1), (Azimuth; Elevation)
%                     2) SWM: 3D-Array of size (2, Qbar_R, Qbar_T), (Azimuth; Elevation) for every Tx-Rx AEs pair inside a SA
% AoA_SV_LOS_loc: Local angle-of-arrival for LoS beamsteering vector; size is based on the used model for the SA/AE: 
%                     1) HSPWM (SWM-SA/PWM-AE): vector of size(2, 1), (Azimuth; Elevation)
%                     2) SWM: 3D-Array of size (2, Qbar_R, Qbar_T), (Azimuth; Elevation) for every Tx-Rx AEs pair inside a SA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dist_SAAEs, AoD_SV_LOS_loc, AoA_SV_LOS_loc] = get_Distance_Angle_LoS(p,mr,nr,mt,nt)

positionLocalSARx = [0; (nr-1-(p.Rx_AoSA.Qdim(2)-1)/2)*p.Rx_AoSA.Delta(2); (mr-1-(p.Rx_AoSA.Qdim(1)-1)/2)*p.Rx_AoSA.Delta(1)];
positionLocalSATx = [0; (nt-1-(p.Tx_AoSA.Qdim(2)-1)/2)*p.Tx_AoSA.Delta(2); (mt-1-(p.Tx_AoSA.Qdim(1)-1)/2)*p.Tx_AoSA.Delta(1)];

if strcmp(p.WaveModelSA,'Sphere') && strcmp(p.WaveModelAE,'Plane') 
    p.RxSATxSA.Rot_mat = p.Rx.Rot_mat;
    p.TxSARxSA.Rot_mat = p.Tx.Rot_mat;
    p.RxSATxSA.Pos = p.RxSATxSA.Rot_mat*positionLocalSARx + p.Rx.Pos;
    p.TxSARxSA.Pos = p.TxSARxSA.Rot_mat*positionLocalSATx + p.Tx.Pos;
    p.RxSATxSAt = get_unitdirvec_glob_loc(p.RxSATxSA,p.TxSARxSA);
    p.TxSARxSAt = get_unitdirvec_glob_loc(p.TxSARxSA,p.RxSATxSA);
    % Physical AoD/AoA
    AoA_SV_LOS_loc = get_anglevec_from_unitdirvec(p.RxSATxSAt.t_Loc);
    AoD_SV_LOS_loc = get_anglevec_from_unitdirvec(p.TxSARxSAt.t_Loc);
    Dist_SAAEs = norm(p.RxSATxSA.Pos-p.TxSARxSA.Pos);
elseif strcmp(p.WaveModelSA,'Plane') && strcmp(p.WaveModelAE,'Plane')
    % See TeraMIMO codes
elseif strcmp(p.WaveModelSA,'Sphere') && strcmp(p.WaveModelAE,'Sphere')
    Dist_SAAEs = zeros(p.Rx_AoSA.Qbar,p.Tx_AoSA.Qbar);
    AoA_SV_LOS_loc = zeros(2,p.Rx_AoSA.Qbar,p.Tx_AoSA.Qbar);
    AoD_SV_LOS_loc = zeros(2,p.Rx_AoSA.Qbar,p.Tx_AoSA.Qbar);
    p.RxAETxAE.Rot_mat = p.Rx.Rot_mat;
    p.TxAERxAE.Rot_mat = p.Tx.Rot_mat;
    for aenr = 1:p.Rx_AoSA.Qbardim(2)
        for aemr = 1:p.Rx_AoSA.Qbardim(1)
            for aent = 1:p.Tx_AoSA.Qbardim(2)
                for aemt = 1:p.Tx_AoSA.Qbardim(1)
                    positionLocalAERx = [0; (aenr-1-(p.Rx_AoSA.Qbardim(2)-1)/2)*p.Rx_AoSA.delta(2); (aemr-1-(p.Rx_AoSA.Qbardim(1)-1)/2)*p.Rx_AoSA.delta(1)] + positionLocalSARx;
                    p.RxAETxAE.Pos = p.RxAETxAE.Rot_mat*positionLocalAERx + p.Rx.Pos;
                    positionLocalAETx = [0; (aent-1-(p.Tx_AoSA.Qbardim(2)-1)/2)*p.Tx_AoSA.delta(2); (aemt-1-(p.Tx_AoSA.Qbardim(1)-1)/2)*p.Tx_AoSA.delta(1)] + positionLocalSATx;
                    p.TxAERxAE.Pos = p.TxAERxAE.Rot_mat*positionLocalAETx + p.Tx.Pos;
                    p.RxAETxAEt = get_unitdirvec_glob_loc(p.RxAETxAE,p.TxAERxAE);
                    p.TxAERxAEt = get_unitdirvec_glob_loc(p.TxAERxAE,p.RxAETxAE);
                    % Physical AoD/AoA
                    AoA_SV_LOS_loc(:,(aenr-1)*p.Rx_AoSA.Qbardim(1)+aemr,(aent-1)*p.Tx_AoSA.Qbardim(1)+aemt) = get_anglevec_from_unitdirvec(p.RxAETxAEt.t_Loc);
                    AoD_SV_LOS_loc(:,(aenr-1)*p.Rx_AoSA.Qbardim(1)+aemr,(aent-1)*p.Tx_AoSA.Qbardim(1)+aemt) = get_anglevec_from_unitdirvec(p.TxAERxAEt.t_Loc);
                    % Distances Matrix
                    Dist_SAAEs((aenr-1)*p.Rx_AoSA.Qbardim(1)+aemr,(aent-1)*p.Tx_AoSA.Qbardim(1)+aemt) = norm(p.RxAETxAE.Pos-p.TxAERxAE.Pos);
                end
            end
        end
    end
else
    error('This wavemodel on SA/AE level isn''t implemented. Supported options are: SWM/SWM, SWM/PWM, PWM/PWM');
end