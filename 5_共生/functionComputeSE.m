function [SE_noBD,SE_direct,SE_BD_MRC, SE_BD_MMSE] = functionComputeSE(nbrOfRealizations,beta_BD_AP,H_PT_AP,H_PT_BD,H_BD_AP,M,K,N,alpha,segma,p,J)
    %Compute uplink SE for Cell-Free mMIMO network.
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

    for n = 1 : nbrOfRealizations
        
        % Display simulation progress
        n
        
        % generate transmit symbols of PTs and backscatters
        PT_signals = (randn(J,1) + 1i*randn(J,1))/sqrt(2); % CN(0,1)
        BD_signals = (randn(1,K) + 1i*randn(1,K))/sqrt(2); % CN(0,1)
        
        % generate noise
        noise = sqrt(noise_var)*(randn(M*N,J) + 1i*randn(M*N,J))/sqrt(2);
        
       %% NoBD-MRC
%         V_MRC = H_PT_AP(:,n,:)/norm(H_PT_AP(:,n,:));
%         gamma_no_BD = PT_power*abs(V_MRC'*H_PT_AP(:,n,:))^2/noise_var;
%         rate_noBD(i) = log2(1+gamma_no_BD);        
    end


