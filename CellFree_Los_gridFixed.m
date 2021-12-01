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

R_cf_min=zeros(1,N);%min rate, cell-free

R_cl_min=zeros(1,N);%collocated massive MIMO


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


%%%%%% APs in Collocative Massive MIMO%%%%%
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

%some parameters
power_tag = zeros(1,K);
a_ch = zeros(1,K);
half_wavelengh = 1 / 2;
segma = 1e-14; %power of noise: 10^-11 mW
c = 3 * 10e8;%speed of light
lambda = c / (f * 10e6);%wavelength
l = lambda * half_wavelengh;%Distance between antennas scale in m
l1 = l / lambda;

%initalize
for k=1:K
    power_tag(k) = 1;
    a_ch(k) = 0.01;
end
%%Inter-symbol interference
PC = zeros(K,K);

%compute Partial sum
for ii = 1:K
    for k = 1:K  
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

R_cf_min(n) = min(R_cf(1,:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                 COLLOCATED MASSIVE MIMO             %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%some parameters
% power_tag = zeros(1,K);
% a_ch = zeros(1,K);
% half_wavelengh = 1 / 2;
% segma = 1e-14; %power of noise: 10^-11 mW
% c = 3 * 10e8;%speed of light
% lambda = c / (f * 10e6);%wavelength
% l = lambda * half_wavelengh;%Distance between antennas scale in m
% l1 = l / lambda;

% %initalize
% for k=1:K
%     power_tag(k) = 1;
%     a_ch(k) = 0.01;
% end

%%Inter-symbol interference
PC = zeros(K,K);

%compute Partial sum
for ii = 1:K
    for k = 1:K  
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

R_cl_min(n) = min(R_cl(1,:));

end

Y=linspace(0,1,N);

hold on 
plot(sort(R_cf_min),Y(:),'r');
plot(sort(R_cl_min),Y(:),'b');
%plot(sort(R_sc_opt_min),Y(:),'b')