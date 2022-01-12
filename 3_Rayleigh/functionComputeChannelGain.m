function [Hhat] = functionComputeChannelGain(nbrOfRealizations, Beta, alpha_f, M, K, N)
    %Generate the channel realizations and compute the channel gain for all channels.The channels are modeled 
    %as uncorrelated Rayleigh fading and the MMSE estimator is used.
    %
    %INPUT:
    %nbrOfRealizations = Number of channel realizations
    %Beta              = Matrix with dimesion M x K where element (m,k) is the 
    %                    pathloss between AP m and UE k
    %alpha_f           = Power control cofficients in the forward channels
    %                    (same for everyone)
    %M                 = Number of APs
    %K                 = Number of UEs in the network
    %N                 = Number of antennas per AP
    %
    %OUTPUT:
    %Hhat              = Matrix with dimension M*N x nbrOfRealizations x K where
    %                    (:,n,k) is the estimated collective channel to UE k at
    %                    channel realization n. (for MMSE estimate) 
    %

    %% Generate channel gain
    
    %forward channel 
    H_f = (randn(1,nbrOfRealizations,K)+1i*randn(1,nbrOfRealizations,K))/sqrt(2);
    a_ch = alpha_f * H_f;

    %backward channel
    %Rayleigh cofficients
    H_b = (randn(M*N,nbrOfRealizations,K)+1i*randn(M*N,nbrOfRealizations,K))/sqrt(2);
    %channel gain
    Gain = zeros(M*N, K);
    %Beta in each AP replicates the number of atennas
    for m = 1:M
        Gain((m-1)*N+1 : m*N, :) = repmat(Beta(m,:), N, 1);
    end
    %reshape matrix to [M*N nbrOfRealizations K]
    Gain = reshape(Gain, [M*N, 1, K]);
    Gain = repmat(Gain, [1 nbrOfRealizations 1]);
    %backward channel gain
    Gain = Gain.*H_b;
    
    %compute the total channel gain
    Hhat = zeros(M*N, nbrOfRealizations, K);
    %hk = ak * gk
    for k = 1:K
        for n = 1:nbrOfRealizations
            Hhat(:,n,k) = a_ch(1,n,k) * Gain(:,n,k);
        end
    end
    