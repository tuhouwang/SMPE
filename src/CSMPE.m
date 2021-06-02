% Ocean acoustics.

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
% Originally developed as part of the author's article (H.Tu, Y.Wang, X.  |
% Ma et al., Applying the Chebyshev-Tau spectral method to solve the      |
% parabolic equation model of wide-angle rational approximation in ocean  |
% acoustics, Journal of Theoretical and Computational Acoustics, 2021,    |
% https://doi.org/10.1142/S2591728521500134) under the supervision of     |
% Prof. Yongxian Wang, National University of Defense Technology, China.  |
% This Matlab/Scilab style code computes the range-independent acoustic   |
% field using the Chebyshev-Tau spectral method based on wide-angle       |
% PE model.                                                               |
% -------------------------------------------------------------------------
clc
close all
clear 
tic;
% edit 'input_SMPE.txt';

[casename, N, np, f, zs, zr, rmax, dr, H, dz, tlmin, tlmax, dep, ...
c, rho, alpha] = ReadEnvParameter('input_SMPE.txt');

c0 = 1500;
ns = 1;
r  = dr : dr : rmax;
nr = length(r);
w  = 2 * pi * f;
k0 = w / c0;

x  = cos( (0 : N) * pi / N )';  
z  = (1.0 - x) * H / 2;         
cs     = interp1(dep, c,     z, 'linear');
rho    = interp1(dep, rho,   z, 'linear');
alpha  = interp1(dep, alpha, z, 'linear');
n      = (c0 ./ cs .* (1.0 + 1i * alpha / (40.0 * pi * ...
         log10( exp(1.0) ) ) ) ) .^ 2 - 1.0;

C  = ConvolutionMatrix( ChebTransFFT(N, n) );
D  = DerivationMatrix ( N + 1);
X  = 4.0 / H ^ 2 / k0 ^ 2 * ConvolutionMatrix( ChebTransFFT(N, rho) ) ...  
     * D  * ConvolutionMatrix( ChebTransFFT(N, 1.0 ./ rho) ) * D + C;

%*********calculated the initial field*************
zd = 0 : 0.1 * dz : H;

cw = interp1(dep, c, zd, 'linear');
[~, ~, ~, ~, ~, starter] = selfstarter(zs, 0.1 * dz, k0, ...
                           w, cw', np, ns, c0, dr, length(zd));

starter        = interp1(zd, starter, z, 'linear');
psi            = zeros(N + 1, nr);
psi(:, 1)      = ChebTransFFT(N, starter);
[pade1, pade2] = epade(np, ns, 1, k0, dr);

%*****************split-step interation******************  
B = zeros(N + 1, N + 1);
A = zeros(N + 1, N + 1);
T = eye(N + 1);

for ip = 1 : np
    A = eye(N + 1) + pade1(ip) * X;    
    B = eye(N + 1) + pade2(ip) * X;
    B(N  ,       :) =  1.0;
    B(N+1, 1:2:N+1) =  1.0; 
    B(N+1, 2:2:N+1) = -1.0; 
    T               =  B \ A * T;
end
    T(N  ,       :) =  1.0;
    T(N+1, 1:2:N+1) =  1.0; 
    T(N+1, 2:2:N+1) = -1.0; 

for ir = 2 : nr
    psi(:, ir) = T * [psi(1 : N - 1, ir - 1); 0; 0];
end

psi = psi .* exp(1i * k0 * dr);
zl  = 0 : dz : H;
xl  = 1 - 2 ./ H * zl' ;
u   = InvChebTrans(psi, xl);
u   = u * diag( 1 ./ sqrt(r) ); 

%********************plot the results**************************

tl    = - 20 * log10( abs( u ));
tl_zr = interp1(zl,  tl,  zr,  'linear');
ShowSoundField(r,  zl,  tl,  tlmin,  tlmax,  casename);
ShowTLcurve(r,  zr,  tl_zr);   
toc;
