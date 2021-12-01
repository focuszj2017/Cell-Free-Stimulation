function APpositions = FixedAPSetup(L, D)
% L = 64;% number of APs
% D=2;% in kilometer

nbrAPsPerDim = sqrt(L);
%Distance between APs in vertical/horizontal direction
interAPDistance = D/nbrAPsPerDim;
locationsGridHorizontal = repmat(interAPDistance/2:interAPDistance:D-interAPDistance/2,[nbrAPsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
APpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);
APpositions = APpositions - D/2 - 1i*D/2;
% scatter(real(APpositions),imag(APpositions),'rs'),hold on;

end