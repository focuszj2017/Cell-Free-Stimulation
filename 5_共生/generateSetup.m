function[beta_PT_AP,beta_PT_BD,beta_BD_AP, APpositions, UEpositions] = generateSetup(M,K,N)

%INPUT:
%M                  = Number of APs for the Cell-free system
%K                  = Number of UEs in the network
%N                  = Number of antennas per APpositions
%
%OUTPUT:
%beta               = Matrix with dimension M x K where
%                     element (m,k) is the pathLoss between APpositions m and UE k 
%UEpositions        = Matrix with dimension K x 2 x 9 with UE locations, where the element (k,:,1) real
%                     part is the horizontal position and the imaginary
%                     part is the vertical position
%APpositions        = Matrix with dimension M x 2 x 9 with APpositions locations, measured in
%                     the same way as UEpositions

%% Define simulation setup

%Size of the coverage area (as a square with wrap-around)
D = 1; %km

%Pathloss parameters
Hb = 15; %Base station height in m
Hm = 1.65; % Mobile height in m
f = 900; %Frequency in MHz
aL = (1.1*log10(f)-0.7)*Hm-(1.56*log10(f)-0.8);
L = 46.3+33.9*log10(f)-13.82*log10(Hb)-aL;
d0=0.01; %km
d1=0.05; %km

%% BD to AP

%Deploy APs on the grid
APpositions = zeros(M,2,9);
positions = functionLocalScattering(M,D);
for m = 1:M
    APpositions(m,1,1) = real(positions(m));
    APpositions(m,2,1) = imag(positions(m));
end

%Wrapped around (8 neighbor cells)
D1=zeros(M,2);
D1(:,1)=D1(:,1)+ D*ones(M,1);
APpositions(:,:,2)=APpositions(:,:,1)+D1;

D2=zeros(M,2);
D2(:,2)=D2(:,2)+ D*ones(M,1);
APpositions(:,:,3)=APpositions(:,:,1)+D2;

D3=zeros(M,2);
D3(:,1)=D3(:,1)- D*ones(M,1);
APpositions(:,:,4)=APpositions(:,:,1)+D3;

D4=zeros(M,2);
D4(:,2)=D4(:,2)- D*ones(M,1);
APpositions(:,:,5)=APpositions(:,:,1)+D4;

D5=zeros(M,2);
D5(:,1)=D5(:,1)+ D*ones(M,1);
D5(:,2)=D5(:,2)- D*ones(M,1);
APpositions(:,:,6)=APpositions(:,:,1)+D5;

D6=zeros(M,2);
D6(:,1)=D6(:,1)- D*ones(M,1);
D6(:,2)=D6(:,2)+ D*ones(M,1);
APpositions(:,:,7)=APpositions(:,:,1)+D6;

D7=zeros(M,2);
D7=D7+ D*ones(M,2);
APpositions(:,:,8)=APpositions(:,:,1)+D7;

D8=zeros(M,2);
D8=D8- D*ones(M,2);
APpositions(:,:,9)=APpositions(:,:,1)+D8;

%Deploy UE in the grid
UEpositions=zeros(K,2,9);
UEpositions(:,:,1)=unifrnd(-D/2,D/2,K,2);

%Wrapped around (8 neighbor cells)
D1=zeros(K,2);
D1(:,1)=D1(:,1)+ D*ones(K,1);
UEpositions(:,:,2)=UEpositions(:,:,1)+D1;

D2=zeros(K,2);
D2(:,2)=D2(:,2)+ D*ones(K,1);
UEpositions(:,:,3)=UEpositions(:,:,1)+D2;

D3=zeros(K,2);
D3(:,1)=D3(:,1)- D*ones(K,1);
UEpositions(:,:,4)=UEpositions(:,:,1)+D3;

D4=zeros(K,2);
D4(:,2)=D4(:,2)- D*ones(K,1);
UEpositions(:,:,5)=UEpositions(:,:,1)+D4;

D5=zeros(K,2);
D5(:,1)=D5(:,1)+ D*ones(K,1);
D5(:,2)=D5(:,2)- D*ones(K,1);
UEpositions(:,:,6)=UEpositions(:,:,1)+D5;

D6=zeros(K,2);
D6(:,1)=D6(:,1)- D*ones(K,1);
D6(:,2)=D6(:,2)+ D*ones(K,1);
UEpositions(:,:,7)=UEpositions(:,:,1)+D6;

D7=zeros(K,2);
D7=D7+ D*ones(K,2);
UEpositions(:,:,8)=UEpositions(:,:,1)+D7;

D8=zeros(K,2);
D8=D8- D*ones(K,2);
UEpositions(:,:,9)=UEpositions(:,:,1)+D8;

%compute the pathLoss

%Store the results
beta_BD_AP = zeros(M,K);
dist = zeros(M,K);

for m = 1:M  
    for k = 1:K
    [dist(m,k),~] = min([norm(APpositions(m,:,1)-UEpositions(k,:,1)), norm(APpositions(m,:,2)-UEpositions(k,:,1)),norm(APpositions(m,:,3)-UEpositions(k,:,1)),norm(APpositions(m,:,4)-UEpositions(k,:,1)),norm(APpositions(m,:,5)-UEpositions(k,:,1)),norm(APpositions(m,:,6)-UEpositions(k,:,1)),norm(APpositions(m,:,7)-UEpositions(k,:,1)),norm(APpositions(m,:,8)-UEpositions(k,:,1)),norm(APpositions(m,:,9)-UEpositions(k,:,1)) ]); %distance between Terminal k and APpositions m
    if dist(m,k)<d0
        betadB = -L - 35*log10(d1) + 20*log10(d1) - 20*log10(d0);
    elseif ((dist(m,k)>=d0) && (dist(m,k)<=d1))
        betadB = -L - 35*log10(d1) + 20*log10(d1) - 20*log10(dist(m,k));
    else
        betadB = -L - 35*log10(dist(m,k)); %large-scale in dB
    end
    
    beta_BD_AP(m,k) = 10^(betadB / 10); 

    end

end

%% PT to AP
PT_positions = zeros(1,2);

%Store the results
beta_PT_AP = zeros(M,1);
dist_PT_AP = zeros(M,1);

for m = 1 : M
    [dist_PT_AP(m),~] = min([norm(APpositions(m,:,1)-PT_positions), norm(APpositions(m,:,2)-PT_positions),norm(APpositions(m,:,3)-PT_positions),norm(APpositions(m,:,4)-PT_positions),norm(APpositions(m,:,5)-PT_positions),norm(APpositions(m,:,6)-PT_positions),norm(APpositions(m,:,7)-PT_positions),norm(APpositions(m,:,8)-PT_positions),norm(APpositions(m,:,9)-PT_positions)]); %distance between PT and APpositions m
    if dist_PT_AP(m)<d0
        betadB = -L - 35*log10(d1) + 20*log10(d1) - 20*log10(d0);
    elseif ((dist_PT_AP(m)>=d0) && (dist_PT_AP(m)<=d1))
        betadB = -L - 35*log10(d1) + 20*log10(d1) - 20*log10(dist_PT_AP(m));
    else
        betadB = -L - 35*log10(dist_PT_AP(m)); %large-scale in dB
    end
    beta_PT_AP = 10^(betadB / 10);
end

%% PT to BD
beta_PT_BD = zeros(K,1);

% Go through all UEs
for k = 1 : K
    dist_PT_BD = norm(PT_positions - UEpositions(k,:,1));
    if dist_PT_BD<d0
        betadB = -L - 35*log10(d1) + 20*log10(d1) - 20*log10(d0);
    elseif ((dist_PT_BD>=d0) && (dist_PT_BD<=d1))
        betadB = -L - 35*log10(d1) + 20*log10(d1) - 20*log10(dist_PT_BD);
    else
        betadB = -L - 35*log10(dist_PT_BD); %large-scale in dB
    end
    beta_PT_BD = 10^(betadB / 10);    
end


