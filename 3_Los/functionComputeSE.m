function [SE_MRC, SE_MMSE] = functionComputeSE(nbrOfRealizations,Beta,dist,theta,Hhat,M,K,N,alpha_f,segma,p,l,lambda)
    %Compute uplink SE for Cell-Free mMIMO network in LoS channel.
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
    %SE_MMSE  = Same as SE_MR but with MMSE or MMSE combining    
    %
    
    
    %If only one transmit power is provided, use the same for all the UEs
    if length(p) == 1
        p = p*ones(K,1);
    end
    
    %Store identity matrices of different sizes
    p_tag = alpha_f^2 * p;
    diagP = diag(p_tag);
    eyeMN = segma * eye(M*N);
    l1 = l / lambda;

    %Prepare to sotre simulation results
    SE_MRC_tot = zeros(K,nbrOfRealizations);
    SE_MMSE_tot = zeros(K,nbrOfRealizations);
    
    %Go through all realizations
    parfor n = 1:nbrOfRealizations
        %% MRC combining
        
        %%Inter-symbol interference
        PC = zeros(K,K);
    
        %compute Partial sum
        for k = 1:K
            for ii = 1:K
                deno = 0;
                for v = 1:M
                    q = 0;
                    if cos(theta(v,k)) == cos(theta(v,ii))
                        q = N;
                    else
                        q = (1 - exp(2i * pi * l1 * N * (cos(theta(v,k)) - cos(theta(v,ii))))) / (1 - exp(2i * pi * l1 * (cos(theta(v,k)) - cos(theta(v,ii)))));
                    end
                
                    deno = deno + ( Beta(v,k) * Beta(v,ii) )^(1/2) * exp(2i * pi * (dist(v,k) - dist(v,ii)) / lambda) * q;
                end
                PC(ii,k) = deno;
            end
        end
        PC1 = (abs(PC)).^2;
    
        %Ssave temp results
        SE_MRC_tot_tmp = zeros(K,1);
        for k = 1:K
            deno = 0;
            for ii = 1:K
                deno = deno + p(k) * alpha_f^2 * PC1(ii,k) / (N * sum(Beta(:,k)));
            end
            
            %Compute numerator and denominator of instantaneous SINR
            numerator = p(k) * alpha_f^2 * N * sum(Beta(:,k));
            denominator = deno - p(k) * alpha_f^2 * N * sum(Beta(:,k)) + segma;
            
            %Compute instantaneous SE for one channel realization
            SE_MRC_tot_tmp(k) =  real(log2(1+numerator/denominator))/nbrOfRealizations;
        end
        SE_MRC_tot(:,n) = SE_MRC_tot_tmp;
        
        %% MMSE combining
        %Ssave temp results
        SE_MMSE_tot_tmp = zeros(K,1);
        
        %Combning vector
        V_MMSE = ((Hhat*diagP*Hhat')+eyeMN)\(Hhat*diagP); %Matrix MN * K
        
        %Go through all UEs
        for k = 1:K
            v = V_MMSE(:,k);
            deno = 0;
            for ii = 1:K
                deno = deno + p(ii) * alpha_f^2 * abs(v' * Hhat(:,ii))^2;
            end
            
            %Compute numerator and denominator of instantaneous SINR            
            numerator = p(k) * alpha_f^2 * abs(v' * Hhat(:,k))^2;
            denominator = deno - numerator + segma * norm(v)^2;
            
            %Compute instantaneous SE for one channel realization
            SE_MMSE_tot_tmp(k) =  real(log2(1+numerator/denominator))/nbrOfRealizations;
        end
        SE_MMSE_tot(:,n) = SE_MMSE_tot_tmp;
        
    end

    SE_MRC = sum(SE_MRC_tot,2);
    SE_MMSE = sum(SE_MMSE_tot,2);