%Compute the single layer of range-independent modal acoustic 
%field, using the Chebyshev-Collation spectral method based on
%PE model. This code was written by Tu Houwang at 08/06/2020 
% from NUDT. psi=sqrt(r)*pressure
clc;
clear;
close all;
tic;
% edit 'input_SMPE.txt';

[casename,N,np,f,zs,zr,rmax,dr,H,dz,tlmin,tlmax,...
        dep,c,rho,alpha] = ReadEnvParameter('input_SMPE.txt');
    
c0             = 1500;
ns             = 1;
r              = dr:dr:rmax;
nr             = length(r);
w              = 2 * pi * f;
k0             = w / c0;

[D,x] = DifferenceMatrix(N);
z  = (1.0 - x) * H / 2;          % GL节点的深度
cs = interp1(dep,c,z,'linear');  % GL节点处的声速
n  = (c0 ./ cs) .^ 2 - 1.0;
X  = 4.0 / H^2 / k0 ^2 * D * D + diag(n);

%*********calculated the initial field*************
zd = 0 : dz : H;
cw = interp1(dep,c,zd,'linear'); % 等距节点处的声速
[~,~,~,~,~,starter] = selfstarter(zs,dz,k0,w,cw',np,ns,c0,dr,length(zd));

psi = zeros(N+1,nr);
psi(:,1) = interp1(zd,starter,z,'linear','extrap');
[pade1, pade2] = epade(np, ns, 1, k0, dr);

%*****************split-step interation******************  
B=zeros(np,N-1,N-1);
A=zeros(np,N-1,N-1);
L=zeros(np,N-1,N-1);
U=zeros(np,N-1,N-1); 
for ip = 1 : np
  B(ip,:,:)=eye(N-1)+pade2(ip)*X(2:N,2:N);
  [L(ip,:,:),U(ip,:,:)]=lu( shiftdim(B(ip,:,:)));
  A(ip,:,:)=eye(N-1)+pade1(ip)*X(2:N,2:N);  
end

for ir = 2 : nr
    q=psi(2:N,ir-1);
    for ip=1:np
        R=shiftdim(A(ip,:,:))*q;
        y=shiftdim(L(ip,:,:))\R; 
        q=shiftdim(U(ip,:,:))\y;
    end
    psi(2:N, ir)=exp(1i*k0*dr)*q;
end

    u=exp(1i*k0*dr).*psi*diag(1./sqrt(r));
    
%********************plot the results**************************

    tl = -20 * log10( abs( u ));
    
    zl = 0 : dz : H;
    tl = interp1(z,tl,zl,'linear'); % GL节点处的声速
    ShowSoundField(r,zl,tl,tlmin,tlmax,casename);
    
    toc;