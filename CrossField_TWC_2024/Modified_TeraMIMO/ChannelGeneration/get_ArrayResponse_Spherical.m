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
% This function computes the array response vector (AE level) for the spherical wave model based on the choosed antenna structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% p: Channel struct that contains the channel parameters
% indx_subc: Index of the subcarrier
% mSA: index of SA (rows) based on TRx 
% nSA: index of SA (columns) based on TRx
% Sca_Dist: Scaterrer Distance
% AngleIn: AoD/AoA for NLoS
% Psi_type: Define whether to calculate a beamsteering or beamforming vector 'SV', 'BF', as follows:
%                 1) @ Tx SV ---> AoD, @ Rx SV ---> AoA
%                 2) @ Tx BF ---> AoD_BF, @ Rx BF ---> AoA_BF
%                 In this version we only use SV
% TRx: Defines the direction of the link, 'T' at Tx side, 'R' at Rx side
% Output Arguments:
% b: Spherical wave model array response, a 3D-Array of size(M, N, num_freqs_per_subcarr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b = get_ArrayResponse_Spherical(p, indx_subc, mSA,nSA, Sca_Dist, AngleIn, Psi_type, TRx)
num_freqs_per_subcarr = p.nFreq(2);
% Initialize Output
if strcmp(TRx, 'T')
    delta1 = p.Tx_AoSA.delta(1);delta2 = p.Tx_AoSA.delta(2);
    Delta1 = p.Tx_AoSA.Delta(1);Delta2 = p.Tx_AoSA.Delta(2);
    Rot_mat = p.Tx.Rot_mat;
    SrcPos = p.Tx.Pos;
    M = p.Tx_AoSA.Qbardim(1);N = p.Tx_AoSA.Qbardim(2);
    M_SAs = p.Tx_AoSA.Qdim(1);N_SAs = p.Tx_AoSA.Qdim(2);
elseif strcmp(TRx, 'R')
    delta1 = p.Rx_AoSA.delta(1); delta2 = p.Rx_AoSA.delta(2);
    Delta1 = p.Rx_AoSA.Delta(1); Delta2 = p.Rx_AoSA.Delta(2);
    Rot_mat = p.Rx.Rot_mat;
    SrcPos = p.Rx.Pos;
    M = p.Rx_AoSA.Qbardim(1);N = p.Rx_AoSA.Qbardim(2);
    M_SAs = p.Rx_AoSA.Qdim(1);N_SAs = p.Rx_AoSA.Qdim(2);
else
    error('TRx has only two options: T/R');
end
b = zeros(M,N,num_freqs_per_subcarr);
positionLocalSA = [0; (nSA-1-(N_SAs-1)/2)*Delta2; (mSA-1-(M_SAs-1)/2)*Delta1];
if strcmp(Psi_type,'SV')
    lambda = p.lambda_subb(indx_subc,:);
end
for m = 1:M
    for n = 1:N
        positionLocalAE = [0; (n-1-(N-1)/2)*delta2; (m-1-(M-1)/2)*delta1] + positionLocalSA;
        TRxAE_Pos = Rot_mat*positionLocalAE + SrcPos;
        UnitDirVec = get_unitdirvec_from_anglevec([0; AngleIn]);
        ScaPos = Sca_Dist*UnitDirVec;
        D = sqrt((ScaPos(1)-TRxAE_Pos(1))^2 + (ScaPos(2)-TRxAE_Pos(2))^2 + (ScaPos(3)-TRxAE_Pos(3))^2);
        b(m,n,:) = exp(-1j*2*pi*D./lambda);
    end
end
end