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
% This function generates Figs. 9b (or 8b), 9c (8c), and 9d (8d) in the manuscript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add Paths
path(pathdef); addpath(pwd);
cd Algorithms; addpath(genpath(pwd)); cd ..;
cd Modified_TeraMIMO; addpath(genpath(pwd)); cd ..;
cd Molecular_Absorption;addpath(genpath(pwd)); cd ..;
cd Utilities; addpath(genpath(pwd)); cd ..;
%%
% Note that the simulation is too heavy and requires strong computing capabilities 
% and long simulation time due to the used size of arrays and dictionaries
clear;clc;
Dist_Vec = [1 20 50 1000]; % To reduce the simulation time I used these values
% To get Fig. 9, you need to adjust the Tx/Rx SA separations and you can use the following vector for evaluation:
% [1, 1.5, 2, 4, 5, 7, 10, 15, 30, 50, 80, 100, 120, 250, 450, 650, 950, 1050, 1250, 1600, 1900, 2050, 2500, 2800, 3000, 3200];
% The values used in the manuscript Scenario 1 (Fig. 8) are:
% [1, 1.5, 2, 3, 4, 5, 7, 10, 15, 20, 35, 40, 50, 60, 70, 80, 90, 100, 220, 320, 450, 650, 800, 900, 1000, 1200, 1400, 1500];
DistVec_Len = length(Dist_Vec);
RotAng_Vec = [-30, -22 -15, -7, 0, 7, 15, 22, 30]; % You can use less number of rotations to reduce the simulation time
% RotAng_Vec = [-30, -22 -15, -7, 0, 7, 15, 22, 30]; % this is the vector that was used in the manuscript
R = length(RotAng_Vec);
for indx_dist = 1 :DistVec_Len
    for indx_rot = 1:R
        get_Proposed_CrossField_RD_Parallel(Dist_Vec(indx_dist), RotAng_Vec(indx_rot));
    end
end
% Wait until the simulation finishes and all files are generated.
%% Load required data
% Note that the simulation is too heavy and requires strong computing capabilities 
% and long simulation time due to the used size of arrays and dictionaries
clear;clc;
cd Res_Online
% Use the same values you put in get_Proposed_CrossField_RD_Parallel while generating the dataset
% and check "p_ch" to know the value of many of these parameters
Dist_Vec = [1 20 50 1000]; % To reduce the simulation time I used these values
% Dist_Vec = [1 5 10 20 50 100 600 1000]; % To reduce the simulation time I used these values
% The values used in the manuscript Scenario 1 (Fig. 8) are:
% [1, 1.5, 2, 3, 4, 5, 7, 10, 15, 20, 35, 40, 50, 60, 70, 80, 90, 100, 220, 320, 450, 650, 800, 900, 1000, 1200, 1400, 1500];
% To get Fig. 9, you need to adjust the Tx/Rx SA separations and you can use the following vector for evaluation:
% [1, 1.5, 2, 4, 5, 7, 10, 15, 30, 50, 80, 100, 120, 250, 450, 650, 950, 1050, 1250, 1600, 1900, 2050, 2500, 2800, 3000, 3200];
DistVec_Len = length(Dist_Vec);
RotAng_Vec = [-30, -22 -15, -7, 0, 7, 15, 22, 30];
R = length(RotAng_Vec);
RequiredSNR = 6; 
E = 10; % Change based on the generated scenario
p_ch.L = 3;
K = 16;
Q_T_v = 4;Q_T_h = 1;
Q_R_v = 4;Q_R_h = 1;
Qbar_T_v = 256;Qbar_T_h = 1;
Qbar_R_v = 16;Qbar_R_h = 1;
p_ch.Sep_factTx = 2;p_ch.Sep_factRx = 10;
Qbar_T = Qbar_T_v*Qbar_T_h;
Qbar_R = Qbar_R_v*Qbar_R_h;
M_T = Qbar_T/2; 
M_R = Qbar_R;
% Initialize array for results:
% NMSE
NMSE_SOMP_NF_Rot_AllTrial = zeros(DistVec_Len,R,E,K);
NMSE_SOMP_FF_Rot_AllTrial = zeros(DistVec_Len,R,E,K);
NMSE_SOMP_Hybrid_NF_FF_Rot_AllTrial = zeros(DistVec_Len,R,E,K);
NMSE_SOMP_Hybrid_FF_NF_Rot_AllTrial = zeros(DistVec_Len,R,E,K);
NMSE_SOMP_Cross_Rot_AllTrial = zeros(DistVec_Len,R,E,K);
% AR
AR_PCSI_Rot_AllTrial = zeros(DistVec_Len,R,E,K);
AR_SOMP_NF_Rot_AllTrial = zeros(DistVec_Len,R,E,K);
AR_SOMP_FF_Rot_AllTrial = zeros(DistVec_Len,R,E,K);
AR_SOMP_Hybrid_NF_FF_Rot_AllTrial = zeros(DistVec_Len,R,E,K);
AR_SOMP_Hybrid_FF_NF_Rot_AllTrial = zeros(DistVec_Len,R,E,K);
AR_SOMP_Cross_Rot_AllTrial = zeros(DistVec_Len,R,E,K);
% Complexity
Complexity_NF_Rot_AllTrial = zeros(DistVec_Len,R,E);
Complexity_FF_Rot_AllTrial = zeros(DistVec_Len,R,E);
Complexity_Cross_Rot_AllTrial = zeros(DistVec_Len,R,E);
SelectedModel_Rot_AllTrial = strings(DistVec_Len,R,E);
for indx_dist = 1:DistVec_Len
    for indx_rot = 1:R
        % define the .mat file name
        simulationname = strcat('OnEst',num2str(Dist_Vec(indx_dist)),'m',num2str(RotAng_Vec(indx_rot)),'RotAng',num2str(E),'Iter',...
        'L',num2str(p_ch.L),'Paths','SNR',num2str(RequiredSNR),'dB',num2str(Q_T_v),'by',num2str(Q_T_h),'x',num2str(Q_R_v),'by',num2str(Q_R_h),'SAs',...
        num2str(Qbar_T_v),'by',num2str(Qbar_T_h),'x',num2str(Qbar_R_v),'by',num2str(Qbar_R_h),'AEs',...
        num2str(M_T),'by',num2str(M_R),'T_R_meas',num2str(p_ch.Sep_factTx),'TxSep', num2str(p_ch.Sep_factRx),'RxSep');
        newfilename=strcat(simulationname,'.mat');
        % load the file
        load(newfilename);
        % Keep all simulation trials to plot the average value
        % NMSE
        NMSE_SOMP_NF_Rot_AllTrial(indx_dist,indx_rot,:,:) = NMSE_SOMP_NF;
        NMSE_SOMP_FF_Rot_AllTrial(indx_dist,indx_rot,:,:) = NMSE_SOMP_FF;
        NMSE_SOMP_Hybrid_NF_FF_Rot_AllTrial(indx_dist,indx_rot,:,:) = NMSE_SOMP_Hybrid_NF_FF;
        NMSE_SOMP_Hybrid_FF_NF_Rot_AllTrial(indx_dist,indx_rot,:,:) = NMSE_SOMP_Hybrid_FF_NF;
        NMSE_SOMP_Cross_Rot_AllTrial(indx_dist,indx_rot,:,:) = NMSE_SOMP_Cross;
        % AR
        AR_PCSI_Rot_AllTrial(indx_dist,indx_rot,:,:) = AR_PCSI;
        AR_SOMP_NF_Rot_AllTrial(indx_dist,indx_rot,:,:) = AR_SOMP_NF;
        AR_SOMP_FF_Rot_AllTrial(indx_dist,indx_rot,:,:) = AR_SOMP_FF;
        AR_SOMP_Hybrid_NF_FF_Rot_AllTrial(indx_dist,indx_rot,:,:) = AR_SOMP_Hybrid_NF_FF;
        AR_SOMP_Hybrid_FF_NF_Rot_AllTrial(indx_dist,indx_rot,:,:) = AR_SOMP_Hybrid_FF_NF;
        AR_SOMP_Cross_Rot_AllTrial(indx_dist,indx_rot,:,:) = AR_SOMP_Cross;
        % Complexity
        Complexity_NF_Rot_AllTrial(indx_dist,indx_rot,:) = Complexity_NF;
        Complexity_FF_Rot_AllTrial(indx_dist,indx_rot,:) = Complexity_FF;
        Complexity_Cross_Rot_AllTrial(indx_dist,indx_rot,:) = Complexity_Cross;
        SelectedModel_Rot_AllTrial(indx_dist,indx_rot,:) = SelectedModel;
    end
end
% to generate Fig. 9a (8a), please check the codes for generating Fig. 7c and 7d (7a and 7b) in "Fig7_OfflineTraining.m"
%% Generate Fig. 9b (8b)
%% NMSE versus distance
d_gamma_S_H = 35;
d_gamma_H_P = 400;
figure('color',[1,1,1]);
semilogx(Dist_Vec,10*log10(mean(NMSE_SOMP_NF_Rot_AllTrial,[2 3 4])),'-diamondb','linewidth',2,'MarkerSize',10);hold on;
semilogx(Dist_Vec,10*log10(mean(NMSE_SOMP_FF_Rot_AllTrial,[2 3 4])),'-square','Color',[0.85,0.33,0.10],'linewidth',2,'MarkerSize',10);
semilogx(Dist_Vec,10*log10(mean(NMSE_SOMP_Hybrid_NF_FF_Rot_AllTrial,[2 3 4])),'-xg','linewidth',2,'MarkerSize',10);
semilogx(Dist_Vec,10*log10(mean(NMSE_SOMP_Hybrid_FF_NF_Rot_AllTrial,[2 3 4])),'-+','Color',[0.30,0.75,0.93],'linewidth',2,'MarkerSize',10);
semilogx(Dist_Vec,10*log10(mean(NMSE_SOMP_Cross_Rot_AllTrial,[2 3 4])),'-ok','linewidth',2,'MarkerSize',10);
dmin_all = min(10*log10(mean(NMSE_SOMP_Cross_Rot_AllTrial,[2 3 4])));
dmax_all = max(10*log10(mean(NMSE_SOMP_FF_Rot_AllTrial,[2 3 4])));
plot([d_gamma_S_H d_gamma_S_H], [dmin_all dmax_all], '--','Color',[0.47,0.67,0.19],'Linewidth', 2);hold on;
plot([d_gamma_H_P d_gamma_H_P], [dmin_all dmax_all], '-','Color',[0.47,0.67,0.19],'Linewidth', 2);hold on;
xlabel('Distance (m)','FontSize',14,'Interpreter','latex');
ylabel('NMSE (dB)','FontSize',14,'Interpreter','latex');
hplt = zeros(7,1);
hplt(1) = plot(NaN,NaN,'-diamondb','linewidth',2,'MarkerSize',7);
hplt(2) = plot(NaN,NaN,'-square','Color',[0.85,0.33,0.10],'linewidth',2,'MarkerSize',7);
hplt(3) = plot(NaN,NaN,'-xg','linewidth',2,'MarkerSize',7);
hplt(4) = plot(NaN,NaN,'-+','Color',[0.30,0.75,0.93],'linewidth',2,'MarkerSize',7);
hplt(5) = plot(NaN,NaN,'-ok','linewidth',2,'MarkerSize',7);
hplt(6) = plot(NaN,NaN,'--','Color',[0.47,0.67,0.19],'linewidth',2,'MarkerSize',7);
hplt(7) = plot(NaN,NaN,'-','Color',[0.47,0.67,0.19],'linewidth',2,'MarkerSize',7);
L_lngd=legend(hplt,'SOMP NF-polar [22]','SOMP FF-angular [33]','Hybrid-field NF-FF [26]', 'Hybrid-field FF-NF [26]', 'Proposed cross-field RD', ...
    '$d_{\gamma_{\mathrm{S}\--\mathrm{H}}}$','$d_{\gamma_{\mathrm{H}\--\mathrm{P}}}$', 'NumColumns',2);
L_lngd.Interpreter='latex';L_lngd.Location='southwest';L_lngd.FontSize=7;
ax=gca;ax.TickLabelInterpreter='latex';ax.FontSize=12;
box off;
grid on;
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
if length(Dist_Vec)>1
    xlim([min(Dist_Vec) max(Dist_Vec)]);
end
axis tight
axlims=axis;
x_range=axlims(2)-axlims(1);
y_range=axlims(4)-axlims(3);
loseness=5;
axis([axlims(1)-Dist_Vec(1)/5 axlims(2)+loseness/1e2*x_range...
    axlims(3)-loseness/1e2*y_range axlims(4)+loseness/1e2*y_range])
savefig('Fig9b.fig');
%% AR versus distance
figure('color',[1,1,1]);
semilogx(Dist_Vec,mean(AR_PCSI_Rot_AllTrial,[2 3 4]),'-pentagram','Color',[0 0 0.5],'linewidth',2,'MarkerSize',10);hold on;
semilogx(Dist_Vec,mean(AR_SOMP_NF_Rot_AllTrial,[2 3 4]),'-diamondb','linewidth',2,'MarkerSize',10);
semilogx(Dist_Vec,mean(AR_SOMP_FF_Rot_AllTrial,[2 3 4]),'-square','Color',[0.85,0.33,0.10],'linewidth',2,'MarkerSize',10);
semilogx(Dist_Vec,mean(AR_SOMP_Hybrid_NF_FF_Rot_AllTrial,[2 3 4]),'-xg','linewidth',2,'MarkerSize',10);
semilogx(Dist_Vec,mean(AR_SOMP_Hybrid_FF_NF_Rot_AllTrial,[2 3 4]),'-+','Color',[0.30,0.75,0.93],'linewidth',2,'MarkerSize',10);
semilogx(Dist_Vec,mean(AR_SOMP_Cross_Rot_AllTrial,[2 3 4]),'-ok','linewidth',2,'MarkerSize',10);
xlabel('Distance (m)','FontSize',14,'Interpreter','latex');
ylabel('Achievable Rate (bits/sec/Hz)','FontSize',14,'Interpreter','latex');
hplt = zeros(6,1);
hplt(1) = plot(NaN,NaN,'-pentagram','Color',[0 0 0.5],'linewidth',2,'MarkerSize',7);
hplt(2) = plot(NaN,NaN,'-diamondb','linewidth',2,'MarkerSize',7);
hplt(3) = plot(NaN,NaN,'-square','Color',[0.85,0.33,0.10],'linewidth',2,'MarkerSize',7);
hplt(4) = plot(NaN,NaN,'-xg','linewidth',2,'MarkerSize',7);
hplt(5) = plot(NaN,NaN,'-+','Color',[0.30,0.75,0.93],'linewidth',2,'MarkerSize',7);
hplt(6) = plot(NaN,NaN,'-ok','linewidth',2,'MarkerSize',7);
L_lngd=legend(hplt, 'PCSI','SOMP NF-polar [22]','SOMP FF-angular [33]','Hybrid-field NF-FF [26]','Hybrid-field FF-NF [26]', 'Proposed cross-field RD');
L_lngd.Interpreter='latex';L_lngd.Location='northeast';L_lngd.FontSize=10;
ax=gca;ax.TickLabelInterpreter='latex';ax.FontSize=14;
box off;
grid on;
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
if length(Dist_Vec)>1
    xlim([min(Dist_Vec) max(Dist_Vec)]);
end
axis tight
axlims=axis;
x_range=axlims(2)-axlims(1);
y_range=axlims(4)-axlims(3);
loseness=5;
axis([axlims(1)-Dist_Vec(1)/5 axlims(2)+loseness/1e2*x_range...
    axlims(3)-loseness/1e2*y_range axlims(4)+loseness/1e2*y_range])
savefig('Fig9c.fig');
 %% Computational Complexity versus distance
SWM_SelectedModel_comp = strcmp(SelectedModel_Rot_AllTrial,"SWM");
HSPWM_SelectedModel_comp = strcmp(SelectedModel_Rot_AllTrial,"HSPWM");
PWM_SelectedModel_comp = strcmp(SelectedModel_Rot_AllTrial,"PWM");
Lbar = 10;Ldom = 9;
G_OVS = 2;
G_T_Polar = 2077;G_R_Polar = 16;
G_T = 256;G_R = 16;
Q_T = 4;Q_R = 4;
d_gamma_S_H = 35;
d_gamma_H_P = 400;
Cross_SWM_RD = G_OVS*G_T_Polar*Lbar *(Q_R-1)*Q_T+G_OVS*G_R_Polar*Lbar*(Q_R-1)*Q_T;
Cross_HSPWM_RD = G_OVS*G_T*Ldom *(Q_R-1)*Q_T+G_OVS*G_R*Ldom*(Q_R-1)*Q_T;
Cross_SWM_RD_Comp = Cross_SWM_RD*SWM_SelectedModel_comp;
Cross_HSPWM_RD_Comp = Cross_HSPWM_RD*HSPWM_SelectedModel_comp;
AVG_Complexity_NF = mean(Complexity_NF_Rot_AllTrial,[2 3]);
AVG_Complexity_FF = mean(Complexity_FF_Rot_AllTrial,[2 3]);
AVG_Complexity_Hybrid_NF_FF = (AVG_Complexity_NF+AVG_Complexity_FF)/2;
AVG_Complexity_Hybrid_FF_NF = (AVG_Complexity_NF+AVG_Complexity_FF)/2;
AVG_Complexity_Cross = mean(Complexity_Cross_Rot_AllTrial+Cross_SWM_RD_Comp+Cross_HSPWM_RD_Comp,[2 3]);
figure('color',[1,1,1]);
semilogx(Dist_Vec,AVG_Complexity_NF,'-diamondb','linewidth',2,'MarkerSize',10);hold on;
semilogx(Dist_Vec,AVG_Complexity_FF,'-square','Color',[0.85,0.33,0.10],'linewidth',2,'MarkerSize',10);
semilogx(Dist_Vec,AVG_Complexity_Hybrid_NF_FF,'-xg','linewidth',2,'MarkerSize',10);
semilogx(Dist_Vec,AVG_Complexity_Hybrid_FF_NF,'-+','Color',[0.30,0.75,0.93],'linewidth',2,'MarkerSize',10);
semilogx(Dist_Vec,AVG_Complexity_Cross,'-ok','linewidth',2,'MarkerSize',10);
dmin_all = min(AVG_Complexity_Cross);
dmax_all = max(AVG_Complexity_NF);
plot([d_gamma_S_H d_gamma_S_H], [dmin_all dmax_all], '--','Color',[0.47,0.67,0.19],'Linewidth', 2);hold on;
plot([d_gamma_H_P d_gamma_H_P], [dmin_all dmax_all], '-','Color',[0.47,0.67,0.19],'Linewidth', 2);hold on;
xlabel('Distance (m)','FontSize',14,'Interpreter','latex');
ylabel('Computational Complexity','FontSize',14,'Interpreter','latex');
hplt = zeros(7,1);
hplt(1) = plot(NaN,NaN,'-diamondb','linewidth',2,'MarkerSize',7);
hplt(2) = plot(NaN,NaN,'-square','Color',[0.85,0.33,0.10],'linewidth',2,'MarkerSize',7);
hplt(3) = plot(NaN,NaN,'-xg','linewidth',2,'MarkerSize',7);
hplt(4) = plot(NaN,NaN,'-+','Color',[0.30,0.75,0.93],'linewidth',2,'MarkerSize',7);
hplt(5) = plot(NaN,NaN,'-ok','linewidth',2,'MarkerSize',7);
hplt(6) = plot(NaN,NaN,'--','Color',[0.47,0.67,0.19],'linewidth',2,'MarkerSize',7);
hplt(7) = plot(NaN,NaN,'-','Color',[0.47,0.67,0.19],'linewidth',2,'MarkerSize',7);
L_lngd=legend(hplt, 'SOMP NF-polar [22]','SOMP FF-angular [33]','Hybrid-field NF-FF [26]','Hybrid-field FF-NF [26]','Proposed cross-field RD', ...
    '$d_{\gamma_{\mathrm{S}\--\mathrm{H}}}$','$d_{\gamma_{\mathrm{H}\--\mathrm{P}}}$','NumColumns',2);
L_lngd.Interpreter='latex';L_lngd.Location='southwest';L_lngd.FontSize=8;
ax=gca;ax.TickLabelInterpreter='latex';ax.FontSize=12;
box off;
grid on;
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
if length(Dist_Vec)>1
    xlim([min(Dist_Vec) max(Dist_Vec)]);
end
axis tight
axlims=axis;
x_range=axlims(2)-axlims(1);
y_range=axlims(4)-axlims(3);
loseness=5;
axis([axlims(1)-Dist_Vec(1)/5 axlims(2)+loseness/1e2*x_range...
    axlims(3)-loseness/1e2*y_range axlims(4)+loseness/1e2*y_range])
savefig('Fig9d.fig');
cd ..