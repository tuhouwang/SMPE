% A simple version RAM developed based on ram1.5.f
% Only suitable for range-independent situations.
% Only one layer of water is considered,  so the 
% default density is 1 g/cm^3,  no attenuation.
% Developed by Houwang Tu from NUDT,  Feb 15, 2020.
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

z   = 0 : dz : H;
nz  = length(z);
cw  = interp1(dep,   c, z, 'linear');
rho = interp1(dep, rho, z, 'linear');

[u, alpw, alpb, ksqw, ksqb, psi1] = selfstarter(zs, dz, k0, w, cw, np, ns, c0, dr, nz);

psi = complex(zeros(nz, nr));
psi(:, 1) = psi1;

[pdu, pdl] = epade(np, ns, 1, k0, dr);
iz = nz - 2;
for ir= 2 : nr
    [r1, r2, r3, s1, s2, s3, f3] = matrc(k0, dz, iz, rho, alpw, alpb, ... 
                                                   ksqw, ksqb, pdu, pdl);
    u = solvetri(iz, u, r1, r2, r3, s1, s2, s3);
    psi(:, ir) = u .* f3;
end

p  = psi * diag( 1.0 ./ sqrt(r) );
tl = -20 * log10( abs(p) );
ShowSoundField(r, z, tl, tlmin, tlmax, casename);
toc



 