close all;
clear;

%% Define simulation setup

%Number of setups with random UE locations
nbrOfSetups = 1;

%Number of channel realizations per setup
nbrOfRealizations = 200;

%Number of APs in the cell-free network
M = 100;

%Number of UEs in the network
K = 40;

%Number of antennas per AP
N = 4;

%Uplink transmit power emmit by each Tag (W)
p = 1000;

%power control coefficients (W)
alpha_f = 1;

%power of noise (W)
segma = 1e-14;

%Prepare to save simulation results
SE_CF_MRC_tot = zeros(K,nbrOfSetups);
SE_CF_MMSE_tot = zeros(K,nbrOfSetups);
SE_CL_MRC_tot = zeros(K,nbrOfSetups);
SE_CL_MMSE_tot = zeros(K,nbrOfSetups);


%Go through all setups
for n = 1:nbrOfSetups

    %Display simulation progress
    %disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);

    %Generate one setup with UEs at random locations
    [pathLoss] = generateSetup(M,K,N);

    %Generate channel realizations
    %channel gain for MMSE
    [Hhat] = functionComputeChannelGain(nbrOfRealizations, pathLoss, alpha_f, M, K, N);

    %Compute SE for the Cell-free mMIMO system with Monte Carlo simulations
    [SE_CF_MRC, SE_CF_MMSE] = functionComputeSE(nbrOfRealizations,pathLoss,Hhat,M,K,N,alpha_f,segma,p);

    %Save SE values
    SE_CF_MRC_tot(:,n) = SE_CF_MRC;
    SE_CF_MMSE_tot(:,n) = SE_CF_MMSE;

    %Remove large matrices at the end of analyzing this setup
    %clear Hhat;

end

%% Plot simulation results

figure;
hold on; box on;
plot(sort(reshape(SE_CF_MRC_tot,[K*nbrOfSetups,1])), linspace(0,1,K*nbrOfSetups),'b-','LineWidth',2);
plot(sort(reshape(SE_CF_MMSE_tot,[K*nbrOfSetups,1])), linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'CellFree (MRC)','CellFree (MMSE)'},'Interpreter','Latex','Location','NorthWest');
%xlim([0 10]);