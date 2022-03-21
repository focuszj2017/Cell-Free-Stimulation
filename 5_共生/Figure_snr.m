clear all;
clc;

%%Define simulation setup

%Number of setups with random UE locations
nbrOfSetups = 1;

%Number of channel realizations per setup
nbrOfRealizations = 10^4;

%Number of APs 
M = 10;

%Number of BDs in the network
K = 4;

% number of receiver antenna
N = 4; 

%Number of transmit antenna
PT_tx = 1;

%Transmit power: W
p = 10^(-1);

%BD power reflection coefficient
alpha = 1;

%ratio between the symbol duration of the BD symbols and that of the PT symbols
J = 100;

%power of noise(W)
noise_var = 10^(-96/10)/1000;

%List of SNR
SNR_list = 0 : 2 : 20;

%Prepare to save simulation results
rate_noBD_avg = zeros(1,length(SNR_list));
rate_direct_avg = zeros(1,length(SNR_list));
rate_BD_MRC_avg = zeros(1,length(SNR_list));
rate_BD_MMSE_avg = zeros(1,length(SNR_list));
SE_MRC_tot = zeros(K,length(SNR_list));
SE_MMSE_tot = zeros(K,length(SNR_list));

%Go through all SNR
for snr = 1 : length(SNR_list)

    %Display simulation progress
    %disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);

    %Generate one setup with UEs at random locations
    [beta_PT_AP,beta_PT_BD,beta_BD_AP] = generateSetup(M,K,N);
    
    %Generate channel realizations
    %channel gain for MMSE
    [H_PT_AP,H_PT_BD,H_BD_AP] = functionComputeChannelGain(nbrOfRealizations,beta_PT_AP,beta_PT_BD, beta_BD_AP, alpha, M, K, N);

    %Compute SE for the Cell-free mMIMO system with Monte Carlo simulations
    [SE_noBD,SE_direct,SE_BD_MRC, SE_BD_MMSE] = functionComputeSE(nbrOfRealizations,beta_BD_AP,H_PT_AP,H_PT_BD,H_BD_AP,M,K,N,alpha,noise_var,p,J);

    %Save SE values
    rate_noBD_avg(snr) = SE_noBD;
    rate_direct_avg(snr) = SE_direct;
    rate_BD_MRC_avg(snr) = sum(SE_BD_MRC)/J;
    rate_BD_MMSE_avg(snr) = sum(SE_BD_MMSE)/J;

end

%% Plot simulation results

figure(1)
plot(SNR_list,rate_noBD_avg,'k-','Linewidth',2);hold on;
plot(SNR_list,rate_direct_avg,'b','Linewidth',2);
plot(SNR_list,rate_direct_avg+rate_BD_MRC_avg,'r','Linewidth',2);
xlim([10,20]);
grid on;
legend('单一主动传输系统速率','共生传输系统-主动速率','共生传输系统-整体速率');



