clc
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Uplink
%Consider a square are of DxD m^2
%M distributed APs serves K terminals, they all randomly located in the area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Inital parameters
M=100; %number of access points
K=40; %number of terminals
Na=4; %number of antennas per AP

D=1; %in kilometer
tau=20;%training length
[U,S,V]=svd(randn(tau,tau));%U includes tau orthogonal sequences 

B=20; %Mhz
Hb = 15; % Base station height in m
Hm = 1.65; % Mobile height in m
f = 1900; % Frequency in MHz
aL = (1.1*log10(f)-0.7)*Hm-(1.56*log10(f)-0.8);
L = 46.3+33.9*log10(f)-13.82*log10(Hb)-aL;

power_f=0.1; %uplink power: 100 mW
noise_p = 10^((-203.975+10*log10(20*10^6)+9)/10); %noise power
Pu = power_f/noise_p;%nomalized receive SNR
Pp=Pu;%pilot power

d0=0.01;%km
d1=0.05;%km

N=200;
R_cf_min=zeros(1,N);%min rate, cell-free, without power allocation
R_cf_opt_min=zeros(1,N); %min rate, cell-free, with power allocation

R_cl_min=zeros(1,N);%collacted massive MIMO
R_cl_opt_min=zeros(1,N);

R_cf_user=zeros(N,K);
R_cl_user=zeros(N,K);

for n=1:N
    n
%%%%%Randomly locations of M APs%%%%
AP=zeros(M,2,9);
AP(:,:,1)=unifrnd(-D/2,D/2,M,2);

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

%%calculation of Beta
%Create an MxK large-scale coefficients beta_mk
BETAA = zeros(M,K);
dist=zeros(M,K);
for m=1:M  
    for k=1:K
    [dist(m,k),index] = min([norm(AP(m,:,1)-Ter(k,:,1)), norm(AP(m,:,2)-Ter(k,:,1)),norm(AP(m,:,3)-Ter(k,:,1)),norm(AP(m,:,4)-Ter(k,:,1)),norm(AP(m,:,5)-Ter(k,:,1)),norm(AP(m,:,6)-Ter(k,:,1)),norm(AP(m,:,7)-Ter(k,:,1)),norm(AP(m,:,8)-Ter(k,:,1)),norm(AP(m,:,9)-Ter(k,:,1)) ]); %distance between Terminal k and AP m
    if dist(m,k)<d0
        betadB=-L - 35*log10(d1) + 20*log10(d1) - 20*log10(d0);
    elseif ((dist(m,k)>=d0) && (dist(m,k)<=d1))
        betadB= -L - 35*log10(d1) + 20*log10(d1) - 20*log10(dist(m,k));
    else
        betadB = -L - 35*log10(dist(m,k)); %large-scale in dB
    end
    
    BETAA(m,k)=10^(betadB/10); 
    end

end

%%%calculation of SINR
SINR = zeros(1,K);
R_cf = zeros(1,K);

%some parameters
power_tag = ones(1,K);
a_ch = ones(1,K);
sign = ones(1,K);
segma = 1; %variance
l = 0.1;%Distance between antennas scale in m
c = 3 * 10e8;%speed of light
lambda = c / f;
l1 = l / lambda;

%%Inter-symbol interference
PC = zeros(K,K);
q_theta = zeros(K,K);
theta = zeros(1,K);
%i denotes k'
for ii=1:K
    for k=1:K  
        PC(ii,k) = sum(  (BETAA(:,k) .*BETAA(:,ii)).^(1/2) .* exp(2i * pi * (dist(:,k) - dist(:,ii)) / lambda)  );
    end
end
PC1 = (abs(PC)).^2;
   
for k=1:K
    deno = 0;
    for ii=1:K
        if cos(theta(k)) == cos(theta(ii))
            q_theta(ii,k) = Na;
        else
            q_theta(ii,k) = (sin(pi * l1 * Na * ( cos(theta(k)) - cos(theta(ii)) )))^2 / (Na * sin(pi * l1 * (  cos(theta(k)) - cos(theta(ii)) ))^2 );
        end
        deno = deno + sign(ii) * a_ch(ii) * PC1(ii,k) * q_theta(ii,k) / (sum(BETAA(:,k)));
    end
    SINR(k) = power_tag(k) * (a_ch(k))^2 * Na * sum(BETAA(:,k)) / (deno - sign(k) * a_ch(k) * Na * sum(BETAA(:,k)) + segma);
    %Rate
    R_cf = log2(1 + SINR(k));
end


end