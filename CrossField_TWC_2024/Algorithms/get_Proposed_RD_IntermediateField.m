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
% This function implements Algorithm 2 in the manuscript to construct the proposed Reduced Dictionary (RD) method by selecting the most correlated
% columns from an oversampled dictionary. More details in Sec. III-C3 HSPWM-based channel estimation with RD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% Abar_OVS: The oversampled far-field DFT dictionary (defined in Step 1 of Algorithm 2)
% Abar_Est: The selected columns of the Polar-domain dictionary based on the first "reference" T-R SA channel estimation
% G_RD: The size of the reduced dictionary
% Output Arguments:
% Abar_RD: The RD far-field DFT dictionary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Abar_RD = get_Proposed_RD_IntermediateField(Abar_OVS, Abar_Est, G_RD)
U_corr = Abar_OVS'*Abar_Est;
U_corr_norm = sqrt(diag(U_corr*U_corr'));
[~ , sel_ind] = sort(U_corr_norm,'descend');
% Read more in Sec. III-C3 HSPWM-based channel estimation with RD
% Detect different clusters centered around the beams corresponding to each path component in the beamspace
sel_ind = sort(sel_ind(1:size(Abar_Est,2)));
IndexSet = zeros(2*size(Abar_Est,2),1);
IndexSet(1:2:end) = sel_ind(1:size(Abar_Est,2));
IndexSet(2:2:end) = sel_ind(1:size(Abar_Est,2))+1;
N_Clu = numel(find(abs([0; diff(IndexSet)])>1)) +1;
loc_each_cluster = unique([1; find(abs([0; diff(IndexSet)])>1)]);
cluster_members = cell(N_Clu,1);
% Detect each set of consecutive indices and group them into a cluster
for indx_num_clu = 1:N_Clu
    if indx_num_clu == N_Clu
        cluster_members{indx_num_clu} = IndexSet(loc_each_cluster(indx_num_clu):end);
    else
        cluster_members{indx_num_clu} = IndexSet(loc_each_cluster(indx_num_clu):loc_each_cluster(indx_num_clu+1)-1);
    end
end
% The center and the number of indices for each cluster
% The number of columns for this cluster is determined by floor(G_RD*L_n_clu/L_dom)
center_each_cluster = floor(cellfun(@mean,cluster_members));
weight_each_cluster = floor(G_RD*(cellfun(@numel,cluster_members))/length(IndexSet));
index_each_cluster = cell(N_Clu,1);
for indx_num_clu = 1:N_Clu
    tmp_indx =  center_each_cluster(indx_num_clu)- floor(weight_each_cluster(indx_num_clu)/2): center_each_cluster(indx_num_clu)+ ceil(weight_each_cluster(indx_num_clu)/2)-1;
    tmp_indx = tmp_indx(tmp_indx>0);
    tmp_indx = tmp_indx(tmp_indx<=size(Abar_OVS,2));
    index_each_cluster{indx_num_clu} = tmp_indx;
end
Indx_RD_each_clu = [];
for indx_num_clu = 1:N_Clu
    Indx_RD_each_clu = [Indx_RD_each_clu index_each_cluster{indx_num_clu}];
end
Indx_RD_each_clu = unique(Indx_RD_each_clu(:));
Abar_RD = Abar_OVS(:,Indx_RD_each_clu);