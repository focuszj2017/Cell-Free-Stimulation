function [SE_CF_MMSE_tot, SE_CF_MMSE_sum, SE_CF_MRC_tot, SE_CF_MRC_sum] = function_BD_alpha(K, alpha_f, nbrOfSetups)

%% Define simulation setup

%Number of setups with random UE locations
% nbrOfSetups = 200;

%Number of channel realizations per setup
nbrOfRealizations = 1000;

% Number of APs in the cell-free network
M = 16;

% Number of UEs in the network
% K = 4;

% Number of antennas per AP
N = 4;

%Uplink transmit power emmit by each Tag (W)
p = 0.1;

%Reflection coefficient
% alpha_f = 1;

%Carrier frequency(MHz)
f = 900;

%power of noise (W)
segma = 10^(-96/10)/1000;

%Wavelength (m)
lambda = (3*1e8)/(f*1e6);

%Distance between antennas (m)
l = lambda * 0.5;

%Prepare to save simulation results
SE_CF_MRC_tot = zeros(K,nbrOfSetups);
SE_CF_MMSE_tot = zeros(K,nbrOfSetups);

%Go through all setups
for n = 1:nbrOfSetups

    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);

    %% Cell-Free mMIMO
   %Generate one setup with UEs at random locations
    [Beta,dist,theta] = generateSetup(M,K,N);

    %Generate channel realizations
    %channel gain for MMSE
    [Hhat] = functionComputeChannelGain(nbrOfRealizations,Beta,dist,theta,alpha_f,M,K,N,l,lambda);

    %Compute SE for the Cell-free mMIMO system with Monte Carlo simulations
    [SE_CF_MRC, SE_CF_MMSE] = functionComputeSE(nbrOfRealizations,Beta,dist,theta,Hhat,M,K,N,alpha_f,segma,p,l,lambda);
    
    %Save SE values
    SE_CF_MRC_tot(:,n) = SE_CF_MRC;
    SE_CF_MMSE_tot(:,n) = SE_CF_MMSE;
    
    %Remove large matrices at the end of analyzing this setup
    %clear Hhat;

end

SE_CF_MRC_sum = sum(SE_CF_MRC_tot,1);
SE_CF_MMSE_sum = sum(SE_CF_MMSE_tot,1);

end


