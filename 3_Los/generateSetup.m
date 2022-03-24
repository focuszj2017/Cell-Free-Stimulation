function[pathLoss, theta, APpositions, UEpositions] = generateSetup(M,K,N)

%INPUT:
%M                  = Number of APs for the Cell-free system
%K                  = Number of UEs in the network
%N                  = Number of antennas per APpositions
%
%OUTPUT:
%pathLoss           = Matrix with dimension M x K where
%                     element (m,k) is the pathLoss between APpositions m and UE k 
%theta              = Matrix with dimension M x K where
%                     element (m,k) is the angle between APpositions m and UE k 
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
pathLoss = zeros(M,K);
dist = zeros(M,K);
theta = zeros(M, K);%the angle from user k to AP m

for m = 1:M  
    for k = 1:K
    [dist(m,k),index] = min([norm(APpositions(m,:,1)-UEpositions(k,:,1)), norm(APpositions(m,:,2)-UEpositions(k,:,1)),norm(APpositions(m,:,3)-UEpositions(k,:,1)),norm(APpositions(m,:,4)-UEpositions(k,:,1)),norm(APpositions(m,:,5)-UEpositions(k,:,1)),norm(APpositions(m,:,6)-UEpositions(k,:,1)),norm(APpositions(m,:,7)-UEpositions(k,:,1)),norm(APpositions(m,:,8)-UEpositions(k,:,1)),norm(APpositions(m,:,9)-UEpositions(k,:,1)) ]); %distance between Terminal k and APpositions m
    if dist(m,k) < 0.01
        dist(m,k) = 0.01;
    end
    betadB = -30.5-36.7*log10(dist(m,k)*1000);
    pathLoss(m,k) = 10^(betadB / 10); 
    theta(m,k) = atan(abs(Ter(k,2,1) - AP(m,2,index)) / abs(Ter(k,1,1) - AP(m,1,index)));
    end

end

