% Ocean acoustics.

% Copyright (C) 2023 Houwang Tu
% -------------------------------------------------------------------------
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the   
% Free Software Foundation, either version 3 of the License, or (at your  
% option) any later version.                                              
%                                                                         
% This code is distributed in the hope that it will be useful, but without
% any warranty; without even the implied warranty of merchantability or   
% fitness for a particular purpose. See the GNU General Public License for
% more details.                                                           
%                                                                         
% You should have received a copy of the GNU General Public License along 
% with this program. If not, see <http://www.gnu.org/licenses/>.          
%                                                                         
%                                                                         
% Originally developed as part of the author's article (H.Tu, Y. Wang, X. 
% Ma et al., Applying the Chebyshev-Tau spectral method to solve the      
% parabolic equation model of wide-angle rational approximation in ocean  
% acoustics, Journal of Theoretical and Computational Acoustics, 2021,    
% https://doi.org/10.1142/S2591728521500134) under the supervision of     
% Prof. Yongxian Wang, National University of Defense Technology, China.  
% This Matlab/Scilab style code computes the range-independent acoustic   
% field using the Chebyshev-Tau spectral method based on wide-angle       
% PE model.                                                               
%                                                                         
%                                                                         
% In 2023, we extended this model to be able to solve inhomogeneous oceans
% with any number of layers. Please refer to (H.Tu,Y.Wang,Y.Zhang et al.,A
% spectrally discretized wide-angle parabolic equation model for simulating
% acoustic propagation in laterally inhomogeneous oceans, The Journal of 
% Acoustical Society of America, 2023, https://doi.org//10.1121/10.0019748)
% for details.
% -------------------------------------------------------------------------
clc
% close all
clear 
tic;
% edit 'input_SMPE.txt';

[casename, Layers, np, ns, c0, freq, zs, zr, rmax, dr, depth, ... 
 dz, tlmin, tlmax, N, range, Coll, dep, c, rho, alpha, ...
 Lowerboundary] = ReadEnvParameter('input_SMPE.txt');

if(Lowerboundary == 'A')
    % Acoustic half-space.
    
    % Add another layer of Perfectly Matched Layer with thickness of two 
    % wavelengths under the original layers.
    lambda = 2 * c0 / freq;
    depth  = [depth; depth(end) + lambda];
    Layers = Layers + 1;
    Coll   = [Coll; max(ceil(Coll(end) * 0.1), 20)];
    
    % The sound speed, density and attenuating coefficient in the PML are 
    % consistent with the seabed (the first SSP).
    PMLd = {depth(end-1:end)'};
    PMLc = {[    c{end,1}(end),    c{end,1}(end)]};
    PMLr = {[  rho{end,1}(end),  rho{end,1}(end)]};
    PMLa = {[alpha{end,1}(end),alpha{end,1}(end)]};
    
    PMLd = repmat(PMLd, 1, N);
    PMLc = repmat(PMLc, 1, N);
    PMLr = repmat(PMLr, 1, N);
    PMLa = repmat(PMLa, 1, N);
    
    dep   = [dep;   PMLd];
    c     = [c;     PMLc];
    rho   = [rho;   PMLr];
    alpha = [alpha; PMLa];
end

[w, k0, z, dep, c, rho, alpha] = Initialization(freq, c0, dz, depth, ... 
 Layers, N, Coll, dep, c, rho, alpha);

if(Lowerboundary == 'A')
    % PML uses complex coordinate transformation technique to absorb 
    % downward radiation energy. 
    tau = (dep{end} - depth(end-1)) / lambda;
    tau = 100 * tau .^ 3 ./ (1 + tau .^ 2);
    % The factor g is used to modify the 
    % Depth Operator of the last layer.
    g = 1.0 ./ (1 + tau + 2i * tau );
else
    g = 1;
end
%***************************obtain the initial field***********************
zd = (0 : 0.1 * dz : depth(end))';
[zl, cl] = Column(Layers, dep(:,1), c(:,1));
cw = interp1(zl, cl, zd, 'linear');
starter = selfstarter(zs, 0.1 * dz, w, cw, np, ns, c0, dr);

phi = {};
r   = dr;
for m = 1 : Layers
    % The self-starter is projected onto the spectral space.
    phi(m,1) = {interp1(zd, starter, dep{m,1}, 'linear')};
    phi(m,1) = {ChebTransFFT(Coll(m,1), phi{m,1})};
end

[pade1, pade2] = epade(np, ns, 1, k0, dr);

range = [dr, range(2:end), rmax];
for j = 1 : length(range)-1
    X = DepthOperator(Layers, k0, c0, Coll, dep(:,j), ...
                      c(:,j), rho(:,j), alpha(:,j), g);

    if( range(j+1) - range(j) > dr )
        % Long flat segment.
        [phit, rt] = FlatStep(X, Layers, Coll, pade1, pade2, np, ... 
                              range(j), range(j+1), dr, dep(:,j), ...
                              rho(:,j), phi(:,end), Lowerboundary);

         phi = [phi, phit(:,2:end)];
         r   = [r,     rt(  2:end)];
    else
         % Long flat section, multiple steps can be calculated without 
         % updating the Depth Operator.     
         phit = OneStep(X, Layers, Coll, pade1, pade2, np, ...
                        dep(:,j), rho(:,j), phi(:,end), Lowerboundary);               
         
         phi = [phi,   phit];
         r   = [r, range(j)];
    end
end

% Inversion of spectral coefficients into physical space.
phi2 = cell(Layers,1);
for i = 1 : Layers
    phi2(i) = {cell2mat(phi(i,:))};
end
psi = KernelFunc(phi2, dz, dep(:,1), Layers);

psi = psi .* exp(1i * k0 * dr);
u   = psi ./ sqrt(r);

% Cut off the sound field of PML.
if(Lowerboundary == 'A')
    ind = find( z <= depth(end-1) );
    z   = z(ind);
    u   = u(ind,:);
end

%*********************calculate and plot the results***********************
tl    = - 20 * log10(abs(u));
% tl_zr = interp1(z, tl, zr, 'linear');
ShowSoundField(r, z, tl, tlmin, tlmax, casename);
% ShowTLcurve(r, zr, tl_zr);   
toc;
