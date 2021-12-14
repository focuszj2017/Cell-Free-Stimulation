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
Na = 4; %number of antennas per AP in cell free
M_cl = 1;%number of APs in collocated massive MIMO
Na_cl = 400;%number of antennas per AP in collocated


D=1; %in kilometer

%B=20; %Mhz
Hb = 15; % Base station height in m
Hm = 1.65; % Mobile height in m
f = 900; % Frequency in MHz
aL = (1.1*log10(f)-0.7)*Hm-(1.56*log10(f)-0.8);
L = 46.3+33.9*log10(f)-13.82*log10(Hb)-aL;

d0=0.01;%km
d1=0.05;%km

N=200;%number of loops

%Some adjustable parameters
power_tag = zeros(1,K);% power of signals emmit by tag k
a_ch = zeros(1,K);% channel gain forward
half_wavelengh = 1 / 2;
segma = 1e-14; %power of noise: 10^-11 mW
c = 3 * 1e8;%speed of light
lambda = c / (f * 1e6);%wavelength
l = lambda * half_wavelengh;%Distance between antennas scale in m
l1 = l / lambda;

%initalize
for k=1:K
    power_tag(k) = 1;
    a_ch(k) = 0.01;
end

R_cf_MR_min=zeros(1,N);%min rate, cell-free
R_cf_MR_sum=zeros(1,N);%capacity
R_cf_MMSE_min=zeros(1,N);
R_cf_MMSE_sum=zeros(1,N);


R_cl_MR_min=zeros(1,N);%collocated massive MIMO
R_cl_MR_sum=zeros(1,N);%capacity
R_cl_MMSE_min=zeros(1,N);
R_cl_MMSE_sum=zeros(1,N);

%%%%%% APs in Cell Free%%%%%
%%% Randomly locations of M APs%%%
AP=zeros(M,2,9);
%Grid fixed
APpositions = FixedAPSetup(M, D);
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


%%%%%% APs in Collocated Massive MIMO%%%%%
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


parfor n=1:N
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                       CELL FREE                     %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%calculation of SINR
SINR = zeros(1,K);
R_cf = zeros(1,K);

%%Inter-symbol interference
PC = zeros(K,K);

%compute Partial sum
for k = 1:K
    for ii = 1:K  
        deno = 0;
        for v = 1:M
            q = 0;
            if cos(theta(v,k)) == cos(theta(v,ii))
                q = Na;
            else
                q = (1 - exp(2i * pi * l1 * Na * (cos(theta(v,k)) - cos(theta(v,ii))))) / (1 - exp(2i * pi * l1 * (cos(theta(v,k)) - cos(theta(v,ii)))));
            end
            
            deno = deno + ( BETAA(v,k) * BETAA(v,ii) )^(1/2) * exp(2i * pi * (dist(v,k) - dist(v,ii)) / lambda) * q;
        end
        PC(ii,k) = deno;
    end
end
PC1 = (abs(PC)).^2;
   
for k = 1:K
    deno = 0;
    for ii = 1:K
        deno = deno + power_tag(ii) * a_ch(ii)^2 * PC1(ii,k) / (Na * sum(BETAA(:,k)));
    end
    SINR(k) = power_tag(k) * (a_ch(k))^2 * Na * sum(BETAA(:,k)) / (deno - power_tag(k) * a_ch(k)^2 * Na * sum(BETAA(:,k)) + segma);
    %Rate
    R_cf(k) = log2(1 + SINR(k));
end

R_cf_MR_min(n) = min(R_cf(1,:));
R_cf_MR_sum(n) = sum(R_cf(1,:));
%%%%%%%%%%%%%%%%%%%%%%%%%     MRC END     %%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%     MMSE BEGIN     %%%%%%%%%%%%%%%%%%%%%%%%%
%compute channel gain
gain = zeros(M*Na, K);
gMK = zeros(M,K,Na);
for idx = 1:Na
    gMK(:,:,idx) = BETAA.^(1/2).*exp(-2i*pi*( dist + (idx-1)*l*cos(theta) ) / lambda);
end

for k = 1:K
    tmp = [];
    for m = 1:M
        gNa = reshape(gMK(m,k,:), [Na,1]);
        tmp = [tmp;gNa];
    end
    gain(:,k) = tmp;
end

%some parameters
p = power_tag(1) * a_ch(1)^2 * ones(K,1);
diagP = diag(p);
eyeMNa = segma * eye(M*Na);
%compute MMSE combining
V_MMSE = ((gain*diagP*gain')+eyeMNa)\(gain*diagP); %Matrix MNa * K
%compute SINR
SINR = zeros(1,K);
R_cf = zeros(1,K);
for k = 1:K
    v = V_MMSE(:,k);
    deno = 0;
    for ii = 1:K
        deno = deno + power_tag(ii) * a_ch(ii)^2 * abs(v' * gain(:,ii))^2;
    end
    partialSum = power_tag(k) * a_ch(k)^2 * abs(v' * gain(:,k))^2;
    SINR(k) = partialSum / (deno - partialSum + segma * norm(v)^2);
    R_cf(k) = log2(1 + SINR(k));
end
R_cf_MMSE_min(n) = min(R_cf(1,:));
R_cf_MMSE_sum(n) = sum(R_cf(1,:));
%%%%%%%%%%%%%%%%%%%%%%%%%     MMSE END     %%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                 COLLOCATED MASSIVE MIMO             %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%     MRC BEGIN     %%%%%%%%%%%%%%%%%%%%%%%%%

%%calculation of Beta
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

%%%calculation of SINR
SINR = zeros(1,K);
R_cl = zeros(1,K);

%%Inter-symbol interference
PC = zeros(K,K);

%compute Partial sum
for k = 1:K
    for ii = 1:K  
        deno = 0;
        for v = 1:M_cl
            q = 0;
            if cos(theta(v,k)) == cos(theta(v,ii))
                q = Na_cl;
            else
                q = (1 - exp(2i * pi * l1 * Na_cl * (cos(theta(v,k)) - cos(theta(v,ii))))) / (1 - exp(2i * pi * l1 * (cos(theta(v,k)) - cos(theta(v,ii)))));
            end
            
            deno = deno + ( BETAA(v,k) * BETAA(v,ii) )^(1/2) * exp(2i * pi * (dist(v,k) - dist(v,ii)) / lambda) * q;
        end
        PC(ii,k) = deno;
    end
end
PC1 = (abs(PC)).^2;
   
for k = 1:K
    deno = 0;
    for ii = 1:K
        deno = deno + power_tag(ii) * a_ch(ii)^2 * PC1(ii,k) / (Na_cl * sum(BETAA(:,k)));
    end
    SINR(k) = power_tag(k) * (a_ch(k))^2 * Na_cl * sum(BETAA(:,k)) / (deno - power_tag(k) * a_ch(k)^2 * Na_cl * sum(BETAA(:,k)) + segma);
    %Rate
    R_cl(k) = log2(1 + SINR(k));
end

R_cl_MR_min(n) = min(R_cl(1,:));
R_cl_MR_sum(n) = sum(R_cl(1,:));
%%%%%%%%%%%%%%%%%%%%%%%%%     MRC END     %%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%     MMSE BEGIN     %%%%%%%%%%%%%%%%%%%%%%%%%
%compute channel gain
gain = zeros(M_cl*Na_cl, K);
gMK = zeros(M_cl,K,Na_cl);
for idx = 1:Na_cl
    gMK(:,:,idx) = BETAA.^(1/2).*exp(-2i*pi*( dist + (idx-1)*l*cos(theta) ) / lambda);
end

for k = 1:K
    tmp = [];
    for m = 1:M_cl
        gNa = reshape(gMK(m,k,:), [Na_cl,1]);
        tmp = [tmp;gNa];
    end
    gain(:,k) = tmp;
end
%some parameters
p = power_tag(1) * a_ch(1)^2 * ones(K,1);
diagP = diag(p);
eyeMNa = segma * eye(M_cl*Na_cl);
%compute MMSE combining
V_MMSE = ((gain*diagP*gain')+eyeMNa)\(gain*diagP); %Matrix M_cl*Na_cl ✖️ K
%compute SINR
SINR = zeros(1,K);
R_cl = zeros(1,K);
for k = 1:K
    v = V_MMSE(:,k);
    deno = 0;
    for ii = 1:K
        deno = deno + power_tag(ii) * a_ch(ii)^2 * abs(v' * gain(:,ii))^2;
    end
    partialSum = power_tag(k) * a_ch(k)^2 * abs(v' * gain(:,k))^2;
    SINR(k) = partialSum / (deno - partialSum + segma * norm(v)^2);
    R_cl(k) = log2(1 + SINR(k));
end
R_cl_MMSE_min(n) = min(R_cl(1,:));
R_cl_MMSE_sum(n) = sum(R_cl(1,:));

%%%%%%%%%%%%%%%%%%%%%%%%%     MMSE END     %%%%%%%%%%%%%%%%%%%%%%%%%
end

Y=linspace(0,1,N);


figure(1)
hold on;
plot(sort(R_cf_MR_min), Y(:), 'r-', 'LineWidth',1.25);
plot(sort(R_cl_MR_min), Y(:), 'r--', 'LineWidth',1.25);
plot(sort(R_cf_MMSE_min), Y(:), 'b-', 'LineWidth',1.25);
plot(sort(R_cl_MMSE_min), Y(:), 'b--', 'LineWidth',1.25);
xlabel('Spectral efficiency [bit/s/Hz]','FontSize',13);
ylabel('CDF','FontSize',13);
legend('MRC(Cell-Free mMIMO)','MRC(Collocated mMIMO)','MMSE(Cell-Free mMIMO)','MMSE(Collocated mMIMO)');


figure(2)
hold on;
plot(sort(R_cf_MR_sum), Y(:), 'r-', 'LineWidth',1.25);
plot(sort(R_cl_MR_sum), Y(:), 'r--', 'LineWidth',1.25);
plot(sort(R_cf_MMSE_sum), Y(:), 'b-', 'LineWidth',1.25);
plot(sort(R_cl_MMSE_sum), Y(:), 'b--', 'LineWidth',1.25);
xlabel('Sum spectral efficiency [bit/s/Hz]','FontSize',13);
ylabel('CDF','FontSize',13);
legend('MRC(Cell-Free mMIMO)','MRC(Collocated mMIMO)','MMSE(Cell-Free mMIMO)','MMSE(Collocated mMIMO)');