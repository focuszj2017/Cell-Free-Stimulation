function [SE_MRC, SE_MMSE] = functionComputeSE_parfor(nbrOfRealizations,Beta,Hhat,M,K,N,alpha_f,segma,p)
    %Compute uplink SE for Cell-Free mMIMO network in Rayleigh channel.
    %Combining methods including MRC and MMSE combining are used.
    %
    %INPUT:
    %nbrOfRealizations = Number of channel realizations
    %Beta              = Matrix with dimesion M x K where element (m,k) is the 
    %                    pathloss between AP m and UE k (for MRC estimate)
    %Hhat              = Matrix with dimension M*N x nbrOfRealizations x K where
    %                    (:,n,k) is the estimated collective channel to UE k at
    %                    channel realization n.(for MMSE estimate) 
    %alpha_f           = Power control cofficients in forward channels (same for everyone)
    %segma             = Power of Guassian noise (mW)
    %M                 = Number of APs
    %K                 = Number of UEs in the network
    %N                 = Number of antennas per AP
    %p                 = Uplink transmit power per UE (same for everyone)
    %
    %OUTPUT:
    %SE_MR     = K x 1 matrix where the element (k,1) is the uplink SE of 
    %            UE k achieved with MRC combining 
    %SE_MMMSE  = Same as SE_MR but with MMSE or MMSE combining    
    %

    %If only one transmit power is provided, use the same for all the UEs
    if length(p) == 1
        p = p*ones(K,1);
    end
    
    %Store identity matrices of different sizes
    eyeMN = eye(M*N);

    %Prepare to sotre simulation results
    SE_MRC_tot = zeros(K,nbrOfRealizations);
    SE_MMSE_tot = zeros(K,nbrOfRealizations);
    
%     %% MRC combining (Mathmatics methods)
%     for k = 1:K
%         numerator = p(k) * N * alpha_f^2 * sum(Beta(:,k))^2;
%         denominator = p(k) * alpha_f^2 * sum(  (Beta(:,k)' * Beta)  ) + segma * sum(Beta(:,k));
%         SE_MRC(k) = log2(1 + numerator / denominator);
%     end
    
    %% MMSE combining

    %Diagonal matrix with transmit powers and its square root
    Dp = diag(p);
    Dp12 = diag(sqrt(p));

    %Go through all channel realizations
    parfor n = 1:nbrOfRealizations
    
        %Display simulation progress
        n

        %Extract channel estimate realizations from all UEs to all APs
        Hhatallj = reshape(Hhat(:,n,:),[M*N, K]);

        %Compute MMSE combing
        V_MMSE = ((Hhatallj*Dp*Hhatallj')+segma*eyeMN) \ (Hhatallj*Dp);

        %Go through all UEs
        for k = 1:K

            %%MRC combining
            v = Hhatallj(:,k); %Extract combining vector
        
            %Compute numerator and denominator of instantaneous SINR at Level 4
            numerator = p(k)*abs(v'*Hhatallj(:,k))^2;
            denominator = norm(v'*Hhatallj*Dp12)^2 + v'*(segma*eyeMN)*v - numerator;
        
            %Compute instantaneous SE for one channel realization
            SE_MRC_tot(k,n) =  real(log2(1+numerator/denominator))/nbrOfRealizations;
            
            %Combining vector
            v = V_MMSE(:,k);

            %Compute numerator and denominator of instantaneous SINR
            numerator = p(k)*abs(v'*Hhatallj(:,k))^2;
            denominator = norm(v'*Hhatallj*Dp12)^2 + v'*segma*eyeMN*v - numerator;
            
            %Compute instantaneous SE for one channel realization
            SE_MMSE_tot(k,n) =  real(log2(1 + numerator / denominator)) / nbrOfRealizations;

        end
    end
    
    SE_MRC = sum(SE_MRC_tot,2);
    SE_MMSE = sum(SE_MMSE_tot,2);

    