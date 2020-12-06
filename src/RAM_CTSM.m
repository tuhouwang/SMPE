
clc
close all
clear 
tic;
% edit 'input_SMPE.txt';

[casename,N,np,f,zs,zr,rmax,dr,H,dz,tlmin,tlmax,dep,...
c,rho,alpha] = ReadEnvParameter('input_SMPE.txt');

c0 = 1500;
ns = 1;
r  = dr:dr:rmax;
nr = length(r);
w  = 2 * pi * f;
k0 = w / c0;

x  = cos( (0 : N) * pi / N )';  
z  = (1.0 - x) * H / 2;         
cs = interp1(dep, c, z, 'linear');
n  = (c0 ./ cs) .^ 2 - 1.0;

C  = ConvolutionMatrix( ChebTransFFT(N, n) );
D  = DerivationMatrix ( N + 1);
X  = 4.0 / H ^2 / k0^2 * D * D + C;

%*********calculated the initial field*************
zd = 0 : 0.1 * dz : H;

cw = interp1(dep, c, zd, 'linear');
[~,~,~,~,~,starter]=selfstarter(zs,0.1 * dz,k0,w,cw',np,ns,c0,dr,length(zd));

starter  = interp1(zd, starter, z, 'linear');
psi      = zeros(N+1, nr);
psi(:,1) = ChebTransFFT(N, starter);
[pade1, pade2] = epade(np, ns, 1, k0, dr);

%*****************split-step interation******************  
B = zeros(np, N+1, N+1);
A = zeros(np, N+1, N+1);
L = zeros(np, N+1, N+1);
U = zeros(np, N+1, N+1); 
for ip = 1 : np
    B(ip, :, :) = eye(N + 1) + pade2(ip) * X;
    B(ip, N, :)         = 1.0;
    B(ip, N+1, 1:2:N+1) = 1.0; 
    B(ip, N+1, 2:2:N+1) = -1.0; 
    [L(ip, :, :), U(ip, :, :)] = lu( shiftdim( B(ip, :, :) ) );
    A(ip, :, :) = eye(N + 1) + pade1(ip) * X;  
end
    q  = psi(:,1);
for ir = 2 : nr
    for ip = 1 : np
        R  = shiftdim( A(ip, :, :) ) * q;
        R(N : N+1) = 0; 
        y = shiftdim( L(ip, :, :) ) \ R; 
        q = shiftdim( U(ip, :, :) ) \ y;           
    end
    psi(:, ir) = q;
end

psi = psi .* exp(1i * k0 * dr);
zl  = 0 : dz : H;
xl  = 1 - 2 ./ H * zl' ;
u   = InvChebTrans(psi, xl);
u   = u * diag( 1 ./ sqrt(r) ); 

%********************plot the results**************************

tl    = - 20 * log10( abs( u ));
ShowSoundField(r, zl, tl, tlmin, tlmax, casename);
tl_zr = interp1(zl, tl, zr, 'linear');
ShowTLcurve(r,zr,tl_zr);   
toc;
