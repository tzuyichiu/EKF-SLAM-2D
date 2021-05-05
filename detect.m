function xra = detect(scan)
% scan = LASER(idx, :)

%Mask13 = uint16(2^13 - 1);
%RR = double(bitand(Mask13, scan));
xra = detectTreesI16(double(scan)/100);
xra(2, :) = xra(2, :) - pi/2;
xra = xra(:, xra(1, :) < 75);
