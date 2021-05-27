% Ocean acoustic normal modes.

% Copyright (C) 2021 Houwang Tu
% -------------------------------------------------------------------------
% This program is free software: you can redistribute it and/or modify it |
% under the terms of the GNU General Public License as published by the   |
% Free Software Foundation, either version 3 of the License, or (at your  |
% option) any later version.                                              |
%                                                                         |
% This code is distributed in the hope that it will be useful, but without|
% any warranty; without even the implied warranty of merchantability or   |
% fitness for a particular purpose. See the GNU General Public License for|
% more details.                                                           |
%                                                                         |
% You should have received a copy of the GNU General Public License along |
% with this program. If not, see <http://www.gnu.org/licenses/>.          |
%                                                                         |
% Originally developed as part of the author's article (Y.Wang, H.Tu, W.  |
% Liu et al., Application of a Chebyshev Collocation Method to Solve a    |
% Parabolic Equation Model of Underwater Acoustic Propagation, Acoustics  |
% Australia, https://doi.org/10.1007/s40857-021-00218-5) under the        |
% supervision of Prof. Yongxian Wang, National University of Defense      |
% Technology, China.                                                      |
%																		  |
% This Matlab/Scilab style code computes the range-independent acoustic   |
% field using the Chebyshev collocation spectral method based on wide-    |
% angle PE model.                                                         |
% -------------------------------------------------------------------------
clc;
clear;
close all;
tic;
% edit 'input_SMPE.txt';

[casename, N, np, f, zs, zr, rmax, dr, H, dz, tlmin, tlmax, dep, ...
 c, rho, alpha] = ReadEnvParameter('input_SMPE.txt');
    
c0  = 1500;
ns  = 1;
r   = dr : dr : rmax;
nr  = length(r);
w   = 2 * pi * f;
k0  = w / c0;

[D, x] = DifferenceMatrix(N);
z      = (1.0 - x) * H / 2;          
cs     = interp1(dep, c,     z, 'linear');
rho    = interp1(dep, rho,   z, 'linear');
alpha  = interp1(dep, alpha, z, 'linear');
n      = (c0 ./ cs .* (1.0 + 1i * alpha / (40.0 * pi * log10( exp(1.0) ) ) ) ) .^ 2 - 1.0;
X      = 4.0 / H ^ 2 / k0 ^ 2 * diag(rho) * D * diag(1.0 ./ rho) * D + diag(n);

%*********calculated the initial field*************
zd = 0 : 0.1 * dz : H;
cw = interp1(dep, c, zd, 'linear');
[~, ~, ~, ~, ~, starter] = selfstarter(zs, 0.1 * dz, k0, w, cw', np, ns, c0, dr, length(zd));

psi = zeros(N + 1, nr);
psi(:, 1) = interp1(zd, starter, z, 'linear');
[pade1, pade2] = epade(np, ns, 1, k0, dr);

%*****************split-step interation******************  
A = zeros(N - 1, N - 1);
B = zeros(N - 1, N - 1);
T = eye(N - 1);
for ip = 1 : np
     A = eye(N - 1) + pade1(ip) * X(2 : N, 2 : N);   
     B = eye(N - 1) + pade2(ip) * X(2 : N, 2 : N);
     T = B \ A * T;
end

for ir = 2 : nr
    psi(2 : N, ir) = T * psi(2 : N, ir - 1);
end

u = exp(1i * k0 * dr) .* psi * diag( 1 ./ sqrt(r) );
    
%********************plot the results**************************
tl = - 20 * log10( abs( u ) );   
zl = 0 : dz : H;
tl = interp1(z, tl, zl, 'linear');
tl_zr = interp1(zl, tl, zr, 'linear'); 
ShowSoundField(r, zl, tl, tlmin, tlmax, casename);
ShowTLcurve(r, zr, tl_zr);    
toc;
