function [Hhat] = functionComputeChannelGain(nbrOfRealizations,Beta,dist,theta,alpha_f,M,K,N,l,lambda)
    %Generate the channel realizations and compute the channel gain for all channels.The channels are modeled 
    %as uncorrelated Line of Sight and the MMSE estimator is used.
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
    Hhat = zeros(M*N,K);
    Gain = zeros(M,K,N);
    for idx = 1:N
        Gain(:,:,idx) = Beta.^(1/2).*exp(-2i*pi*( dist + (idx-1)*l*cos(theta) ) / lambda);
    end
    
    for k = 1:K
        tmp = [];
        for m = 1:M
            gN = reshape(Gain(m,k,:), [N,1]);
            tmp = [tmp;gN];
        end
        Hhat(:,k) = tmp;
    end
    
    