function [H_PT_AP,H_PT_BD,H_BD_AP] = functionComputeChannelGain(nbrOfRealizations,beta_PT_AP,beta_PT_BD, beta_BD_AP, alpha, M, K, N)
    %Generate the channel realizations and compute the channel gain for all channels.The channels are modeled 
    %as uncorrelated Rayleigh fading and the MMSE estimator is used.
    %
    %INPUT:
    %nbrOfRealizations = Number of channel realizations
    %beta              = Matrix with dimesion M x K where element (m,k) is the 
    %                    pathloss between AP m and UE k
    %alpha             = BD power reflection coefficient
    %M                 = Number of APs
    %K                 = Number of UEs in the network
    %N                 = Number of antennas per AP
    %
    %OUTPUT:
    %H_PT_AP           = Matrix with dimension M*N x nbrOfRealizations x 1 where
    %                    (:,n,1) is the estimated collective channel from PT to AP at
    %                    channel realization n. 
    %H_PT_BD           = Matrix with dimension K x nbrOfRealizations x 1 where
    %                    (:,n,k) is the estimated collective channel frin PT to UE k at
    %                    channel realization n. 
    %H_BD_AP           = Matrix with dimension M*N x nbrOfRealizations x K where
    %                    (k,n,1) is the estimated collective channel to UE k at
    %                    channel realization n. 
    %

    %% Generate channel gain from PT to AP
    %Rayleigh cofficients
    rayleigh_PT_AP = (randn(1,nbrOfRealizations,M*N)+1i*randn(1,nbrOfRealizations,M*N))/sqrt(2);
    %channel gain
    Gain = zeros(M*N, 1);
    %Beta in each AP replicates the number of atennas
    for m = 1:M
        Gain((m-1)*N+1 : m*N, :) = repmat(sqrt(beta_PT_AP(m,:)), N, 1);
    end
    %reshape Gain matrix to [M*N nbrOfRealizations 1]
    Gain = reshape(Gain, [M*N, 1, 1]);
    Gain = repmat(Gain, [1 nbrOfRealizations 1]);    
    %Channel gain from PT to AP
    H_PT_AP = zeros(M*N, nbrOfRealizations, 1)
    for  ii = 1 : M*N
        for n = 1 : nbrOfRealizations
            H_PT_AP(ii,n,:) = Gain(ii,n,1) * rayleigh_PT_AP(1,n,ii);
        end
    end
    
    %% Generate channel gain from PT to BD   
    %Rayleigh cofficients
    rayleigh_PT_BD = (randn(K,nbrOfRealizations,1)+1i*randn(K,nbrOfRealizations,1))/sqrt(2);
    %reshape Gain matrix to [1 nbrOfRealizations K]
    Gain = reshape(sqrt(beta_PT_BD),[K 1 1]);
    Gain = repmat(Gain, [1 nbrOfRealizations 1]);
    %Channel gain from PT to BD
    H_PT_BD = rayleigh_PT_BD.* Gain;
    
    %% Generate channel gain from BD to AP
    %Rayleigh cofficients
    rayleigh_BD_AP = (randn(M*N,nbrOfRealizations,K)+1i*randn(M*N,nbrOfRealizations,K))/sqrt(2);
    %channel gain
    Gain = zeros(M*N, K);
    %Beta in each AP replicates the number of atennas
    for m = 1:M
        Gain((m-1)*N+1 : m*N, :) = repmat(sqrt(beta_BD_AP(m,:)), N, 1);
    end
    %reshape matrix to [M*N nbrOfRealizations K]
    Gain = reshape(Gain, [M*N, 1, K]);
    Gain = repmat(Gain, [1 nbrOfRealizations 1]);
    %channel gain form BD to AP
    H_BD_AP = Gain.*rayleigh_BD_AP;

    