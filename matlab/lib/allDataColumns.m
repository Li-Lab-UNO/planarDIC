function C = allDataColumns(F, dFx, dFy, wx, wy, H)
% The kinematic and stress columns of CellScale's AllData export, from deformation gradients F (2 x 2 x n,
% relative to the first frame of the cycle; x = tester X = axial, y = tester Y = circumferential) and the
% forces tared at that frame (dFx, dFy, mN, n x 1). The X force is divided by the width wx across it
% (mm), the Y force by wy, and both by the thickness H (mm):
%   lx = F11, ly = F22, kx = F12, ky = F21
%   txx = dFx / (wx H) * lx, tyy = dFy / (wy H) * ly          (Cauchy stress, kPa)
%   txy = dFy / (wy H) * kx, tyx = dFx / (wx H) * ky
%   ShearAngleDeg: angle between the deformed x and y axes (90 undeformed)
% CellScale's own export follows these definitions with its 4-point F, the rake separations at the start
% of the test as widths and a typed thickness (lx ... tyy, kx ... tyx exact; the angle within 0.02 deg).

F11 = squeeze(F(1, 1, :)); F12 = squeeze(F(1, 2, :)); F21 = squeeze(F(2, 1, :)); F22 = squeeze(F(2, 2, :));
Pxx = dFx(:) / (wx * H); Pyy = dFy(:) / (wy * H);
C11 = F11 .^ 2 + F21 .^ 2; C22 = F12 .^ 2 + F22 .^ 2; C12 = F11 .* F12 + F21 .* F22;
C = table(F11, Pxx .* F11, F22, Pyy .* F22, acosd(C12 ./ sqrt(C11 .* C22)), F12, Pyy .* F12, F21, Pxx .* F21, ...
          'VariableNames', {'lx', 'txx', 'ly', 'tyy', 'ShearAngleDeg', 'kx', 'txy', 'ky', 'tyx'});
end
