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
% This function generates Figs. 7c and 7d (or Figs. 7a and 7b, depending on the simulation parameters) in the manuscript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add Paths
path(pathdef); addpath(pwd);
cd Algorithms; addpath(genpath(pwd)); cd ..;
cd Modified_TeraMIMO; addpath(genpath(pwd)); cd ..;
cd Molecular_Absorption;addpath(genpath(pwd)); cd ..;
cd Utilities; addpath(genpath(pwd)); cd ..;
%%
clear;clc;
Dist_Vec = [1 4 7 20 100 500 1000]; % To reduce the simulation time I used these values
% The values used in the manuscript Scenario 1 (Figs. 7a and 7b) are:
% [1, 1.5, 2, 3, 4, 5, 7, 10, 20, 30, 40, 60, 80, 100, 200, 300, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1300, 1500];
% The values used in the manuscript Scenario 2 (Figs. 7c and 7d) are:
% [1, 1.5, 2, 3, 4, 5, 7, 10, 20, 40, 60, 80, 100, 200, 300, 500, 700, 900,1000, 1200, 1600, 1900, 2100, 2400, 2800, 3000, 3300];
DistVec_Len = length(Dist_Vec);
RotAng_Vec = -30:2:30; 
R_offline = length(RotAng_Vec);
% You can use randi([-30,30],1,31) but make sure that the generated vector is unique to avoid re-simulating the same value
for indx_dist = 1 :DistVec_Len
    for indx_rot = 1:R_offline
        get_ProposedMetric_OfflineTraining_Parallel(Dist_Vec(indx_dist), RotAng_Vec(indx_rot));
    end
end
% Wait until the simulation finishes and all files are generated.
%% Load required data
clear;clc;
cd Res_Offline
% Use the same values you put in get_ProposedMetric_OfflineTraining_Parallel while generating the dataset
Dist_Vec = [1 4 7 20 100 500 1000]; 
% The values used in the manuscript Scenario 1 (Figs. 7a and 7b) are:
% [1, 1.5, 2, 3, 4, 5, 7, 10, 20, 30, 40, 60, 80, 100, 200, 300, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1300, 1500];
% The values used in the manuscript Scenario 2 (Figs. 7c and 7d) are:
% [1, 1.5, 2, 3, 4, 5, 7, 10, 20, 40, 60, 80, 100, 200, 300, 500, 700, 900,1000, 1200, 1600, 1900, 2100, 2400, 2800, 3000, 3300];
DistVec_Len = length(Dist_Vec);
RotAng_Vec = -30:2:30; 
R_offline = length(RotAng_Vec);
RequiredSNR = [0 6 8 14 16]; 
LenSNR = length(RequiredSNR);
E_offline = 10; 
% Change based on the generated scenario and parameters
p_ch.L=3;
Q_T_v = 4;Q_T_h = 1;
Q_R_v = 4;Q_R_h = 1;
Qbar_T_v = 256;Qbar_T_h = 1;
Qbar_R_v = 16;Qbar_R_h = 1;
p_ch.Sep_factTx = 1;p_ch.Sep_factRx = 1;
Qbar_T = Qbar_T_v*Qbar_T_h;
M_T = [Qbar_T/2 Qbar_T/4 Qbar_T/8]; 
MeasTx_vec = length(M_T);
% Initialize array for results:
EtaMax_Dist_Rot_AvgTrial = zeros(LenSNR,MeasTx_vec,DistVec_Len,R_offline);
EtaMax_Dist_Rot_AllTrial = zeros(LenSNR,MeasTx_vec,DistVec_Len,R_offline,E_offline);
for indx_dist = 1:DistVec_Len
    for indx_rot = 1:R_offline
        % define the .mat file name
        simulationname = strcat('OffTr',num2str(Dist_Vec(indx_dist)),'m',num2str(RotAng_Vec(indx_rot)),'RotAng',num2str(E_offline),'Iter',...
        'L',num2str(p_ch.L),'Paths',num2str(Q_T_v),'by',num2str(Q_T_h),'x',num2str(Q_R_v),'by',num2str(Q_R_h),'SAs',...
        num2str(Qbar_T_v),'by',num2str(Qbar_T_h),'x',num2str(Qbar_R_v),'by',num2str(Qbar_R_h),'AEs',...
        num2str(p_ch.Sep_factTx),'TxSep', num2str(p_ch.Sep_factRx),'RxSep');
        newfilename=strcat(simulationname,'.mat');
        % load the file
        load(newfilename);
        % averaging over all simulation trials "E_offline" to plot the CDF
        EtaMax_Dist_Rot_AvgTrial(:,:,indx_dist,indx_rot) = mean(Eta_MAT_Max,3);
        % Keep all simulation trials "E_offline" to plot the CDF
        EtaMax_Dist_Rot_AllTrial(:,:,indx_dist,indx_rot,:) = Eta_MAT_Max;
    end
end
%% Generate Fig. 7c (or 7a)
% Proposed metric versus distance
% Since we generate the files only for one value for the number of multi-paths (e.x. L=1 or L=3), you must also perform the simulation again for another L.
EtaMax_Dist_AvgRot_AvgTrial = mean(EtaMax_Dist_Rot_AvgTrial,4);
% Plot for one L value
indx_meas =1; % Select the number of used beams/measurements to plot the proposed metric, this should not be greater than MeasTx_vecfigure('color',[1,1,1]);
ColorsCell = {[0.85,0.33,0.10], 'b','r','m','g'};
Line_Style_Marker = '-*'; % '--o', ':diamond'
for indx_SNR =1:LenSNR
        semilogx(Dist_Vec,10*log10(squeeze(EtaMax_Dist_AvgRot_AvgTrial(indx_SNR,indx_meas,:))), ...
            Line_Style_Marker,'Color',ColorsCell{indx_SNR},'linewidth',2,'MarkerSize',8);hold on;
        xlabel('Distance (m)','FontSize',14,'Interpreter','latex');ylabel('$\eta$ (dB)','FontSize',14,'Interpreter','latex');
        ax=gca;ax.TickLabelInterpreter='latex';ax.FontSize=14;
        box off;grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
        if length(Dist_Vec)>1
            xlim([min(Dist_Vec) max(Dist_Vec)]);
        end
        axis tight;axlims=axis;
        x_range=axlims(2)-axlims(1);
        y_range=axlims(4)-axlims(3);
        loseness=5;
        axis([axlims(1)-Dist_Vec(1)/5 axlims(2)+loseness/1e2*x_range axlims(3)-loseness/1e2*y_range axlims(4)+loseness/1e2*y_range]);
end
Dist_gamma_H_P = 35;
Dist_gamma_S_H = 1000;
hold on;
xpwm = [0.8 Dist_gamma_H_P Dist_gamma_H_P 0.8];
ypwm = [-17,-17,-1,-1];
p_pwm = patch(xpwm,ypwm,[0.96,0.74,0.91]);
p_pwm.EdgeColor = 'none';
p_pwm.FaceAlpha = 0.4;
xhspwm = [Dist_gamma_H_P Dist_gamma_S_H Dist_gamma_S_H Dist_gamma_H_P];
yhspwm = [-17,-17,-1,-1];
p_hspwm = patch(xhspwm,yhspwm,[0.77,0.96,0.74]);
p_hspwm.EdgeColor = 'none';
p_hspwm.FaceAlpha = 0.4;
xswm = [Dist_gamma_S_H 3359.95 3359.95 Dist_gamma_S_H];
yswm = [-17,-17,-1,-1];
p_swm = patch(xswm,yswm,[0.97,0.63,0.63]);
p_swm.EdgeColor = 'none';
p_swm.FaceAlpha = 0.4;
% I provide some codes for generating the legend, you should before generate all possible values for L
hplt = zeros(11, 1);
hplt(1) = plot(NaN,NaN,'-*k','linewidth',2,'MarkerSize',8);
hplt(2) = plot(NaN,NaN,'--ok','linewidth',2,'MarkerSize',8);
hplt(3) = plot(NaN,NaN,':diamondk','linewidth',2,'MarkerSize',8);
hplt(4) = plot(NaN,NaN,'Color',[0.85,0.33,0.10],'linewidth',2,'MarkerSize',8);
hplt(5) = plot(NaN,NaN,'b','linewidth',2,'MarkerSize',8);
hplt(6) = plot(NaN,NaN,'r','linewidth',2,'MarkerSize',8);
hplt(7) = plot(NaN,NaN,'m','linewidth',2,'MarkerSize',8);
hplt(8) = plot(NaN,NaN,'g','linewidth',2,'MarkerSize',8);
hplt(9) = p_pwm;
hplt(10) = p_hspwm;
hplt(11) = p_swm;
L_lngd=legend(hplt,'L=1','L=3','L=5','SNR=0 dB','SNR=6 dB','SNR=8 dB','SNR=14 dB','SNR=16 dB','SWM','HSPWM','PWM');
L_lngd.Interpreter='latex';L_lngd.Location='southwest';L_lngd.FontSize=11;
savefig('Fig7c.fig');
%% Generate Fig. 7d (or 7b)
indx_meas =1; 
indx_SNR = 2; % we defined RequiredSNR = [0 6 8 14 16], then indx_SNR = 2 corresponds to SNR = 6 dB
EtaMax_CDFs = squeeze(EtaMax_Dist_Rot_AllTrial(indx_SNR,indx_meas,:,:,:));
EtaMax_CDFs = reshape(EtaMax_CDFs,size(EtaMax_CDFs,1),size(EtaMax_CDFs,2)*size(EtaMax_CDFs,3));
Threshold_SWM_HSPM = 0.4;
Threshold_HSPM_PWM = 0.12; 
cmap = jet(DistVec_Len) ;
h_prnt = figure();
hold on;
xpwm = [0 Threshold_HSPM_PWM Threshold_HSPM_PWM 0];
ypwm = [0    0 1 1];
p_pwm = patch(xpwm,ypwm,[0.97,0.63,0.63]);
p_pwm.EdgeColor = 'none';
p_pwm.FaceAlpha = 0.4;
xhspwm = [Threshold_HSPM_PWM Threshold_SWM_HSPM Threshold_SWM_HSPM Threshold_HSPM_PWM];
yhspwm = [0    0 1 1];
p_hspwm = patch(xhspwm,yhspwm,[0.77,0.96,0.74]);
p_hspwm.EdgeColor = 'none';
p_hspwm.FaceAlpha = 0.4;
xswm = [Threshold_SWM_HSPM 0.78 0.78 Threshold_SWM_HSPM];
yswm = [0    0 1 1];
p_swm = patch(xswm,yswm,[0.96,0.74,0.91]);
p_swm.EdgeColor = 'none';
p_swm.FaceAlpha = 0.4;
for indx_dist = 1:DistVec_Len
    h_plt(indx_dist) = cdfplot(EtaMax_CDFs(indx_dist,:));
    set(h_plt(indx_dist),'Color',cmap(indx_dist,:))
end
xlabel('X ($\eta$)','FontSize',14,'Interpreter','latex');
ylabel('F(X)','FontSize',14,'Interpreter','latex');
title('');
Xarr1 = [0.192142857142857,0.221428571428571];
Yarr1 = [0.420952380952388,0.527619047619048];
arr1 = annotation('textarrow',Xarr1,Yarr1);
arr1.String='$\gamma_{\mathrm{H}\--\mathrm{P}}$';
arr1.Interpreter = 'latex';
arr1.FontSize=14;
Xarr2 = [0.433571428571428,0.468928571428571];
Yarr2 = [0.705714285714289,0.795714285714286];
arr2 = annotation('textarrow',Xarr2,Yarr2);
arr2.String='$\gamma_{\mathrm{S}\--\mathrm{H}}$';
arr2.Interpreter = 'latex';
arr2.FontSize=14;
cm = colormap(cmap);
cb = colorbar;
caxis([log10(Dist_Vec(1)) log10(Dist_Vec(end))]);
cb.Label.String = 'log(distance)';
cb.FontSize=14;
cb.TickLabelInterpreter='latex';
cb.Label.Interpreter='latex';
hplt = zeros(3, 1);
hplt(1) = p_pwm;
hplt(2) = p_hspwm;
hplt(3) = p_swm;
le=legend(hplt,'PWM','HSPWM','SWM');
le.Interpreter='latex';le.Location='northwest';le.FontSize=14;
ax=gca;ax.TickLabelInterpreter='latex';ax.FontSize=14;
box off;grid on;
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
axis tight;
set(h_prnt,'Units','Inches');
pos = get(h_prnt,'Position');
set(h_prnt,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(h_prnt,'Fig7d','-dpdf','-r0');
savefig('Fig7d.fig');
cd ..