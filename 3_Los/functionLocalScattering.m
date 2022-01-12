function positions = functionLocalScattering(L, D)
% L = 64;% number of APs
% D=2;% in kilometer

nbrAPsPerDim = sqrt(L);
%Distance between APs in vertical/horizontal direction
interAPDistance = D/nbrAPsPerDim;
locationsGridHorizontal = repmat(interAPDistance/2:interAPDistance:D-interAPDistance/2,[nbrAPsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
positions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);
positions = positions - D/2 - 1i*D/2;

%figure(3)
%scatter(real(positions),imag(positions),'rs'),hold on;

end