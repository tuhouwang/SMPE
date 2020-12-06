% Self starter，psi1是初始距离处的psi乘sqrt(dr)声压
%当dz过小时，本程序计算的结果与RAM不同，但应以本程序为准，因为RAM输入声速剖面时
%是插值的，本程序是用Munk声速剖面算的，要多精确有多精确。
function [u,alpw,alpb,ksqw,ksqb,psi1]=selfstarter(zs,dz,k0,w,cw,np,ns,c0,dr,nz)
alpw=sqrt(cw/c0);
alpb=ones(1,nz);
ksqw=(w./cw).^2-k0^2;
ksqb=zeros(1,nz);
ksqb=complex(ksqb);
u=zeros(nz,1);
rhob=ones(1,nz);
iz=nz-2;

si=1+zs/dz;
is=fix(si);
dis=si-double(is);
u(is)=(1-dis)*sqrt(2*pi/k0)/(dz*alpw(is));
u(is+1)=dis*sqrt(2*pi/k0)/(dz*alpw(is));
clear is si dis
% Divide the delta function by (1-X)**2 to get a smooth rhs.
pd1=complex(0); pd2=complex(-1);
[r1,r2,r3,s1,s2,s3,~]=matrc(k0,dz,iz,rhob,alpw,alpb,ksqw,ksqb,pd1,pd2);
u=complex(u);
u=solvetri(iz,u,r1,r2,r3,s1,s2,s3);
u=solvetri(iz,u,r1,r2,r3,s1,s2,s3);
%Apply the operator (1-X)**2*(1+X)**(-1/2)*exp(ci*k0*dr*sqrt(1+X)).
[pd1, pd2]=epade( np, ns, 0, k0, dr);
[r1,r2,r3,s1,s2,s3,f3]=matrc(k0,dz,iz,rhob,alpw,alpb,ksqw,ksqb,pd1,pd2);
u=solvetri(iz,u,r1,r2,r3,s1,s2,s3);
psi1=u.*f3;
end