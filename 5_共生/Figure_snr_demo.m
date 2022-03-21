clear all;
clc;

PT_tx = 1; % number of transmit antenna
PT_power = 10^(-1); % transmit power: W
BD_num = 4; % number of BD devices
AP_num = 100; %number of AP receiver
AP_rx = 4; % number of receiver antenna
alpha = 1; % BD power reflection coefficient
J = 100; % ratio between the symbol duration of the BD symbols and that of the PT symbols

SNR_list = 0 : 2 : 20;
ite_num = 10^4;

rate_noBD_avg = zeros(1,length(SNR_list));
rate_direct_avg = zeros(1,length(SNR_list));
rate_BD_avg = zeros(1,length(SNR_list));

for snr = 1 : length(SNR_list)
    SNR = SNR_list(snr);
    noise_var = PT_power/10^(SNR/10); % power of the noise
    rate_direct = zeros(1,ite_num);
    rate_BD = zeros(1,ite_num);
    rate_noBD = zeros(1,ite_num);
    rate_direct_MMSE = zeros(1,ite_num);
    rate_BD_MMSE = zeros(1,ite_num);
    rate_noBD_MMSE = zeros(1,ite_num);    
    parfor i = 1 : ite_num
        disp([' SNR ' num2str(SNR) ' round ' num2str(i)]);
        LSF_PT_AP = -120; % large-scale channel gain from PT to AP in dB
        LSF_PT_BD = -110; % large-scale channel gain from PT to BD in dB
        LSF_BD_AP = -20; % large-scale channel gain from BD to AP in dB
%         var_noise = -110; % noise variance in dBm
        
        % generate transmit symbols of PTs and backscatters
        PT_signals = (randn(J,1) + 1i*randn(J,1))/sqrt(2); % CN(0,1)
        BD_signals = (randn(1,BD_num) + 1i*randn(1,BD_num))/sqrt(2); % CN(0,1)
        
        % generate channel gains
        H_PT_AP = sqrt(10^(LSF_PT_AP/10))*(randn(AP_num*AP_rx,PT_tx) + 1i*randn(AP_num*AP_rx,PT_tx))/sqrt(2);
        H_PT_BD = sqrt(10^(LSF_PT_BD/10))*(randn(BD_num,PT_tx) + 1i*randn(BD_num,PT_tx))/sqrt(2);
        H_BD_AP = sqrt(10^(LSF_BD_AP/10))*(randn(AP_num*AP_rx,BD_num) + 1i*randn(AP_num*AP_rx,BD_num))/sqrt(2);
%         H_PT_AP = (randn(AP_num*AP_rx,PT_tx) + 1i*randn(AP_num*AP_rx,PT_tx))/sqrt(2);
%         H_PT_BD = (randn(BD_num,PT_tx) + 1i*randn(BD_num,PT_tx))/sqrt(2);
%         H_BD_AP = sqrt(1/10)*(randn(AP_num*AP_rx,BD_num) + 1i*randn(AP_num*AP_rx,BD_num))/sqrt(2);

        % generate noise
        noise = sqrt(noise_var)*(randn(AP_num*AP_rx,J) + 1i*randn(AP_num*AP_rx,J))/sqrt(2);
        
       %% NoBD-MRC
        v_MRC = H_PT_AP/norm(H_PT_AP); 
        gamma_no_BD = PT_power*abs(v_MRC'*H_PT_AP)^2/noise_var;
        rate_noBD(i) = log2(1+gamma_no_BD);
        
       %% CSR-MRC
        % compute the rate of active link
        h_equal = H_PT_AP;
        for ii = 1 : BD_num
            h_equal = h_equal + sqrt(alpha)*H_BD_AP(:,ii)*H_PT_BD(ii,:)*BD_signals(1,ii);
        end
        v_MRC = h_equal/norm(h_equal); % assuming MRC combining is used
        gamma_direct = PT_power*abs(v_MRC'*h_equal)^2/ noise_var;
        rate_direct(i) = log2(1 + gamma_direct);
                
        % compute the rate of passive link
        % Monte carlo
        H_eq = zeros(AP_num*AP_rx,BD_num);
        gamma_BD = zeros(1,BD_num);
        for k = 1 : BD_num
            H_eq(:,k) = sqrt(alpha*J)*H_PT_BD(k,:)*H_BD_AP(:,k);
            v_MRC = H_eq(:,k)/norm(H_eq(:,k));
            deno = 0;
            for j = (k + 1) : BD_num
                deno = deno + PT_power*J*alpha*abs(v_MRC'*H_PT_BD(j,:)*H_BD_AP(:,j)*BD_signals(:,j))^2;
            end
            gamma_BD(:,k) = PT_power*abs(v_MRC'*H_eq(:,k)*BD_signals(1,k))^2/(deno + noise_var);
        end
        rate_BD(i) = sum(log2(1+gamma_BD))/J; % sum of passive link rate
    end
    
    rate_noBD_avg(snr) = mean(rate_noBD);
    rate_direct_avg(snr) = mean(rate_direct);
    rate_BD_avg(snr) = mean(rate_BD);
end

figure(1)
plot(SNR_list,rate_noBD_avg,'k-','Linewidth',2);hold on;
plot(SNR_list,rate_direct_avg,'b','Linewidth',2);
plot(SNR_list,rate_direct_avg+rate_BD_avg,'r','Linewidth',2);
xlim([10,20]);
%ylim([3,6]);
grid on;
legend('单一主动传输系统速率','共生传输系统-主动速率','共生传输系统-整体速率');

figure(2)
grid on;
plot(SNR_list,rate_BD_avg,'m','Linewidth',2);
xlim([10,20]);
%ylim([0,3]);
legend('共生传输系统-被动速率');
