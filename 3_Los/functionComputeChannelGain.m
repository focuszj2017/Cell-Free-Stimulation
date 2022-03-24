function [Hhat] = functionComputeChannelGain(nbrOfRealizations, Beta, theta, alpha_f, M, K, N)
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
    a_ch = alpha_f

    %backward channel
    %channel gain
    Gain = zeros(M*N, K);
    %Beta in each AP replicates the number of atennas
    for m = 1:M
        Gain((m-1)*N+1 : m*N, :) = repmat(sqrt(Beta(m,:)), N, 1);
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
    
    
    Hhat = zeros(M*N, K);
    Gain = zeros(M,K,N);
    %% ToDO
    for idx = 1:N
        Gain(:,:,idx) = BETAA.^(1/2).*exp(-2i*pi*( dist + (idx-1)*l*cos(theta) ) / lambda);
    end
    
    for k = 1:K
        tmp = [];
        for m = 1:M
            gN = reshape(Gain(m,k,:), [N,1]);
            tmp = [tmp;gN];
        end
        Hhat(:,k) = tmp;
    end
    
    