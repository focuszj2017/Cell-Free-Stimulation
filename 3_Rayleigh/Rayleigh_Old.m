clc 
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Uplink
%Consider a square are of DxD m^2
%M distributed APs serves K terminals, they all randomly located in the area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Inital parameters

M = 100; %number of access points in cell free
K = 40; %number of terminals
N = 4; %number of antennas per AP in cell free
M_cl = 1;%number of APs in collocated massive MIMO
N_cl = 400;%number of antennas per AP in collocated

D=1; %in kilometer

Hb = 15; % Base station height in m
Hm = 1.65; % Mobile height in m
f = 900; % Frequency in MHz
aL = (1.1*log10(f)-0.7)*Hm-(1.56*log10(f)-0.8);
L = 46.3+33.9*log10(f)-13.82*log10(Hb)-aL;

d0=0.01;%km
d1=0.05;%km

nbrOfRealizations=200;%number of loops

power_tag = ones(1,K);% power of signals emmit by tag k
half_wavelengh = 1 / 2;
segma = 1e-14; %power of noise: 10^-11 mW
c = 3 * 1e8;%speed of light
lambda = c / (f * 1e6);%wavelength
l = lambda * half_wavelengh;%Distance between antennas scale in m
l1 = l / lambda;

alpha_f = 0.01 * ones(1,K);%power control coefficients
H_f = (randn(1,K)+1i*randn(1,K))/sqrt(2);%the gain of rayleigh forward channel,matrix with 1 * K
a_ch = alpha_f.*H_f;% channel gain forward

%to store the result

R_cf_MR_avg = zeros(1,nbrOfRealizations);%min rate, cell-free
R_cf_MR_sum = zeros(1,nbrOfRealizations);%capacity
R_cf_MMSE_avg = zeros(1,nbrOfRealizations);
R_cf_MMSE_sum = zeros(1,nbrOfRealizations);

R_cl_MR_avg = zeros(1,nbrOfRealizations);%collocated massive MIMO
R_cl_MR_sum = zeros(1,nbrOfRealizations);%capacity
R_cl_MMSE_avg = zeros(1,nbrOfRealizations);
R_cl_MMSE_sum = zeros(1,nbrOfRealizations);

%% Scattering points

%%%%%% APs in Cell Free %%%%%

% Randomly locations of M APs
AP=zeros(M,2,9);
%Grid fixed
APpositions = functionLocalScattering(M, D);
for m = 1 : M
    AP(m,1,1) = real(APpositions(m));
    AP(m,2,1) = imag(APpositions(m));
end
%Wrapped around (8 neighbor cells)
D1=zeros(M,2);
D1(:,1)=D1(:,1)+ D*ones(M,1);
AP(:,:,2)=AP(:,:,1)+D1;

D2=zeros(M,2);
D2(:,2)=D2(:,2)+ D*ones(M,1);
AP(:,:,3)=AP(:,:,1)+D2;

D3=zeros(M,2);
D3(:,1)=D3(:,1)- D*ones(M,1);
AP(:,:,4)=AP(:,:,1)+D3;

D4=zeros(M,2);
D4(:,2)=D4(:,2)- D*ones(M,1);
AP(:,:,5)=AP(:,:,1)+D4;

D5=zeros(M,2);
D5(:,1)=D5(:,1)+ D*ones(M,1);
D5(:,2)=D5(:,2)- D*ones(M,1);
AP(:,:,6)=AP(:,:,1)+D5;

D6=zeros(M,2);
D6(:,1)=D6(:,1)- D*ones(M,1);
D6(:,2)=D6(:,2)+ D*ones(M,1);
AP(:,:,7)=AP(:,:,1)+D6;

D7=zeros(M,2);
D7=D7+ D*ones(M,2);
AP(:,:,8)=AP(:,:,1)+D7;

D8=zeros(M,2);
D8=D8- D*ones(M,2);
AP(:,:,9)=AP(:,:,1)+D8;


%%%%%% APs in Collocated Massive MIMO %%%%%

AP_cl=zeros(M_cl,2,9);
AP_cl(1,:,1) = [0 0];

%Wrapped around (8 neighbor cells)
D1=zeros(M_cl,2);
D1(:,1)=D1(:,1)+ D*ones(M_cl,1);
AP_cl(:,:,2)=AP_cl(:,:,1)+D1;

D2=zeros(M_cl,2);
D2(:,2)=D2(:,2)+ D*ones(M_cl,1);
AP_cl(:,:,3)=AP_cl(:,:,1)+D2;

D3=zeros(M_cl,2);
D3(:,1)=D3(:,1)- D*ones(M_cl,1);
AP_cl(:,:,4)=AP_cl(:,:,1)+D3;

D4=zeros(M_cl,2);
D4(:,2)=D4(:,2)- D*ones(M_cl,1);
AP_cl(:,:,5)=AP_cl(:,:,1)+D4;

D5=zeros(M_cl,2);
D5(:,1)=D5(:,1)+ D*ones(M_cl,1);
D5(:,2)=D5(:,2)- D*ones(M_cl,1);
AP_cl(:,:,6)=AP_cl(:,:,1)+D5;

D6=zeros(M_cl,2);
D6(:,1)=D6(:,1)- D*ones(M_cl,1);
D6(:,2)=D6(:,2)+ D*ones(M_cl,1);
AP_cl(:,:,7)=AP_cl(:,:,1)+D6;

D7=zeros(M_cl,2);
D7=D7+ D*ones(M_cl,2);
AP_cl(:,:,8)=AP_cl(:,:,1)+D7;

D8=zeros(M_cl,2);
D8=D8- D*ones(M_cl,2);
AP_cl(:,:,9)=AP_cl(:,:,1)+D8;


%% Loops of realizations

parfor n = 1:nbrOfRealizations
    n

%Randomly locations of K terminals:
Ter=zeros(K,2,9);
Ter(:,:,1)=unifrnd(-D/2,D/2,K,2);

%Wrapped around (8 neighbor cells)
D1=zeros(K,2);
D1(:,1)=D1(:,1)+ D*ones(K,1);
Ter(:,:,2)=Ter(:,:,1)+D1;

D2=zeros(K,2);
D2(:,2)=D2(:,2)+ D*ones(K,1);
Ter(:,:,3)=Ter(:,:,1)+D2;

D3=zeros(K,2);
D3(:,1)=D3(:,1)- D*ones(K,1);
Ter(:,:,4)=Ter(:,:,1)+D3;

D4=zeros(K,2);
D4(:,2)=D4(:,2)- D*ones(K,1);
Ter(:,:,5)=Ter(:,:,1)+D4;

D5=zeros(K,2);
D5(:,1)=D5(:,1)+ D*ones(K,1);
D5(:,2)=D5(:,2)- D*ones(K,1);
Ter(:,:,6)=Ter(:,:,1)+D5;

D6=zeros(K,2);
D6(:,1)=D6(:,1)- D*ones(K,1);
D6(:,2)=D6(:,2)+ D*ones(K,1);
Ter(:,:,7)=Ter(:,:,1)+D6;

D7=zeros(K,2);
D7=D7+ D*ones(K,2);
Ter(:,:,8)=Ter(:,:,1)+D7;

D8=zeros(K,2);
D8=D8- D*ones(K,2);
Ter(:,:,9)=Ter(:,:,1)+D8;


%% Cell-Free

%%%%%%%%%%%%%%%%%%%%%%%%%     MRC BEGIN     %%%%%%%%%%%%%%%%%%%%%%%%%

%%calculation of Beta

%Create an MxK large-scale coefficients beta_mk
BETAA = zeros(M, K);
dist = zeros(M, K);
theta = zeros(M, K);%the angle from user k to AP m

for m = 1:M  
    for k = 1:K
    [dist(m,k),index] = min([norm(AP(m,:,1)-Ter(k,:,1)), norm(AP(m,:,2)-Ter(k,:,1)),norm(AP(m,:,3)-Ter(k,:,1)),norm(AP(m,:,4)-Ter(k,:,1)),norm(AP(m,:,5)-Ter(k,:,1)),norm(AP(m,:,6)-Ter(k,:,1)),norm(AP(m,:,7)-Ter(k,:,1)),norm(AP(m,:,8)-Ter(k,:,1)),norm(AP(m,:,9)-Ter(k,:,1)) ]); %distance between Terminal k and AP m
    if dist(m,k)<d0
        betadB = -L - 35*log10(d1) + 20*log10(d1) - 20*log10(d0);
    elseif ((dist(m,k)>=d0) && (dist(m,k)<=d1))
        betadB = -L - 35*log10(d1) + 20*log10(d1) - 20*log10(dist(m,k));
    else
        betadB = -L - 35*log10(dist(m,k)); %large-scale in dB
    end
    
    BETAA(m,k) = 10^(betadB / 10); 

    theta(m,k) = atan(abs(Ter(k,2,1) - AP(m,2,index)) / abs(Ter(k,1,1) - AP(m,1,index)));
    end

end

%calculation of SINR
SINR = zeros(1,K);
R_cf = zeros(1,K);

%%â€”â?â€”â?â€”â?â€”â?â€”â?â€”â?â€”â?â€”â?â€”â?â€”â?â€”â?TODO

R_cf_MR_avg(n) = min(R_cf(1,:));
R_cf_MR_sum(n) = sum(R_cf(1,:));
%%%%%%%%%%%%%%%%%%%%%%%%%     MRC END     %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
%compute channel gain
Gain = zeros(M*N, K);
%Beta in each AP replicates the number of atennas
for m = 1:M
    Gain((m-1)*N+1 : m*N, :) = repmat(BETAA(m,:), N, 1);
end
%Rayleigh cofficients
H_b = (randn(M*N,K)+1i*randn(M*N,K))/sqrt(2);
Gain = Gain.*H_b;

%hk = ak * gk
Hhatallj = zeros(M*N, K);
for idx = 1:K
    Hhatallj(:,k) = a_ch(k) * Gain(:,k);
end

%power cofficients
p = power_tag(1) * ones(K,1);
Dp = diag(p);
Dp12 = diag(sqrt(p));
eyeMNa = segama * eye(M*N);

%compute MMSE combining
V_MMSE = ((Hhatallj*Dp*Hhatallj')+eyeMNa) \ (Hhatallj*Dp);

%Go through all APs
for k = 1:K

    %MMSE combing
    v = V_MMSE(:,k);
    
    %Compute numerator and denominator of instantaneous SINR
    numerator = p(k)*abs(v'*Hhatallj(:,k))^2;
    denominator = norm(v'*Hhatallj*Dp12)^2 + v'*eyeMNa*v - numerator;

    %Compute instantaneous SE for one channel realization
    R_cf_MMSE_avg(k) = R_cf_MMSE_avg(k) + real(log2(1+numerator/denominator))/nbrOfRealizations;                        

end
%%%%%%%%%%%%%%%%%%%%%%%%%     MMSE END     %%%%%%%%%%%%%%%%%%%%%%%%%



%% Collocated Massive MIMO

%%%%%%%%%%%%%%%%%%%%%%%%%     MRC BEGIN     %%%%%%%%%%%%%%%%%%%%%%%%%   

%calculation of Beta
%Create an MxK large-scale coefficients beta_mk
BETAA = zeros(M_cl,K);
dist = zeros(M_cl,K);
theta = zeros(M_cl,K);%the angle from user k to AP m
for m = 1:M_cl  
    for k = 1:K
    [dist(m,k),index] = min([norm(AP_cl(m,:,1)-Ter(k,:,1)), norm(AP_cl(m,:,2)-Ter(k,:,1)),norm(AP_cl(m,:,3)-Ter(k,:,1)),norm(AP_cl(m,:,4)-Ter(k,:,1)),norm(AP_cl(m,:,5)-Ter(k,:,1)),norm(AP_cl(m,:,6)-Ter(k,:,1)),norm(AP_cl(m,:,7)-Ter(k,:,1)),norm(AP_cl(m,:,8)-Ter(k,:,1)),norm(AP_cl(m,:,9)-Ter(k,:,1)) ]); %distance between Terminal k and AP m
    if dist(m,k)<d0
        betadB = -L - 35*log10(d1) + 20*log10(d1) - 20*log10(d0);
    elseif ((dist(m,k)>=d0) && (dist(m,k)<=d1))
        betadB = -L - 35*log10(d1) + 20*log10(d1) - 20*log10(dist(m,k));
    else
        betadB = -L - 35*log10(dist(m,k)); %large-scale in dB
    end
    
    BETAA(m,k) = 10^(betadB / 10); 

    theta(m,k) = atan(abs(Ter(k,2,1) - AP_cl(m,2,index)) / abs(Ter(k,1,1) - AP_cl(m,1,index)));
    end

end

%calculation of SINR
SINR = zeros(1,K);
R_cl = zeros(1,K);
%___________TODO

%%%%%%%%%%%%%%%%%%%%%%%%%     MRC END     %%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%     MMSE BEGIN     %%%%%%%%%%%%%%%%%%%%%%%%%
%compute channel gain
Gain = zeros(M_cl*N_cl, K);
for m = 1:M_cl
    Gain((m-1)*N_cl+1 : m*N_cl, :) = repmat(BETAA(m,:), N_cl, 1);
end
H_b = (randn(M_cl*N_cl,K)+1i*randn(M_cl*N_cl,K))/sqrt(2);%Rayleigh cofficients
Gain = Gain.*H_b;

%hk = ak * gk
Hhatallj = zeros(M_cl*N_cl, K);
for idx = 1:K
    Hhatallj(:,k) = a_ch(k) * Gain(:,k);
end


p = power_tag(1) * ones(K,1);
Dp = diag(p);
Dp12 = diag(sqrt(p));
eyeMNa = segama * eye(M_cl*N_cl);

%compute MMSE combining
V_MMSE = ((Hhatallj*Dp*Hhatallj')+eyeMNa) \ (Hhatallj*Dp);

%Go through all APs
for k = 1:K

    %MMSE combing
    v = V_MMSE(:,k);
    
    %Compute numerator and denominator of instantaneous SINR
    numerator = p(k)*abs(v'*Hhatallj(:,k))^2;
    denominator = norm(v'*Hhatallj*Dp12)^2 + v'*eyeMNa*v - numerator;

    %Compute instantaneous SE for one channel realization
    R_cl_MMSE_avg(k) = R_cl_MMSE_avg(k) + real(log2(1+numerator/denominator))/nbrOfRealizations;                        

end
%%%%%%%%%%%%%%%%%%%%%%%%%     MMSE END     %%%%%%%%%%%%%%%%%%%%%%%%%


end

%% plot
Y=linspace(0,1,K);

figure(1)
hold on;
plot(R_cf_MR_avg, Y(:), 'r-', 'LineWidth',1.25);
plot(R_cl_MR_avg, Y(:), 'r--', 'LineWidth',1.25);
plot(R_cf_MMSE_avg, Y(:), 'b-', 'LineWidth',1.25);
plot(R_cl_MMSE_avg, Y(:), 'b--', 'LineWidth',1.25);
xlabel('Spectral efficiency [bit/s/Hz]','FontSize',13);
ylabel('CDF','FontSize',13);
legend('MRC(Cell-Free mMIMO)','MRC(Collocated mMIMO)','MMSE(Cell-Free mMIMO)','MMSE(Collocated mMIMO)');



% figure(1)
% hold on;
% plot(sort(R_cf_MR_avg), Y(:), 'r-', 'LineWidth',1.25);
% plot(sort(R_cl_MR_avg), Y(:), 'r--', 'LineWidth',1.25);
% plot(sort(R_cf_MMSE_avg), Y(:), 'b-', 'LineWidth',1.25);
% plot(sort(R_cl_MMSE_avg), Y(:), 'b--', 'LineWidth',1.25);
% xlabel('Spectral efficiency [bit/s/Hz]','FontSize',13);
% ylabel('CDF','FontSize',13);
% legend('MRC(Cell-Free mMIMO)','MRC(Collocated mMIMO)','MMSE(Cell-Free mMIMO)','MMSE(Collocated mMIMO)');


% figure(2)
% hold on;
% plot(sort(R_cf_MR_sum), Y(:), 'r-', 'LineWidth',1.25);
% plot(sort(R_cl_MR_sum), Y(:), 'r--', 'LineWidth',1.25);
% plot(sort(R_cf_MMSE_sum), Y(:), 'b-', 'LineWidth',1.25);
% plot(sort(R_cl_MMSE_sum), Y(:), 'b--', 'LineWidth',1.25);
% xlabel('Sum spectral efficiency [bit/s/Hz]','FontSize',13);
% ylabel('CDF','FontSize',13);
% legend('MRC(Cell-Free mMIMO)','MRC(Collocated mMIMO)','MMSE(Cell-Free mMIMO)','MMSE(Collocated mMIMO)');

