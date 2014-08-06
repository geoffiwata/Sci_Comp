function [ ] = final_prob3(  )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
close all
clear all
%First lets define some variables
a = 0.01;   %grid spacing (0.01)
L = 2;      %length of well
n = L/a+1;  %number of grid points
V = 16;    %potential (30)
m = 1;      %mass
W = 3;      %width of barrier
tau = 0.01;  %timestep (0.01)
N=n+2;
time = 310; %number of timesteps (forward)

position = zeros(1,n+2);

for m=1:n+2
    position(1,m) = -1+(m-2)*a;
end

psi = zeros(N,1);
psi(1) = 0;
psi(n+2) = 0;

VV=zeros(N,1);
for g = floor(3*N/7):floor(4*N/7)
    if g < (floor(3*N/7)+W) || g > (floor(4*N/7)-W)
        VV(g) = V;
    else
        VV(g) = -V;
    end    
end

%Initial wavefunction
k = 108;
sigma = .07;

for p=2:n+1
    psi(p) = exp(+1i*k*position(1,p)) * exp(-(position(1,p)+0.5)^2/(2*sigma^2));
end

%Hamiltonian things
front = floor(3*N/7);
back = floor(4*N/7);
potup = -2*m*a^2*V; %the upper potential term
potlow = +2*m*a^2*V; %the lower potential term
midd=-2+4*1i*m*a^2/(tau);
ff = 8*1i*m*a^2/(tau);
%1-d laplacian
A = sparse(1:n,1:n,midd,N,N);
B = sparse(2:n,1:n-1,1,N,N);
C = sparse(front:front+W,front:front+W,potup,N,N);
D = sparse(front+W:back-W,front+W:back-W,potlow,N,N);
E = sparse(back-W:back,back-W:back,potup,N,N);

Ham = (A+C+D+E+B+B');

chi = zeros(N,1);

psiold=psi;

for t=1:N
psi_plot(t,1) = (psi(t))*conj(psi(t));
end

m=1;
w=1;

while m<=time
   m     
chi = zeros(N,1);
   %Calculate Chi(m)
    r = ff*psiold - Ham*chi; %residual
    p = r;
    rr = r'*r;
    for g=1:100
        %Set up variables
        Hamp=Ham*p;
        alpha = rr / (p'*Hamp);
        T=alpha*p;
        %Update value of chi
        for g=3:N-3
            chi(g) = chi(g) + T(g);
        end
        r =  r - alpha*Hamp; %update r
        rrnew=r'*r;
%         %If r is too small, quit.
%             if sqrt(rrnew)<1e5
%                 i
%                 break;
%             end
        p = r + rrnew*p./(rr); %update p
        rr=rrnew;
    end
%     if (psi)'*conj(psi) > 10
%         break
%     end
    psi = chi - psiold; %update psi
    if w==2
        q = 1+m/w;
        for t=1:N
            psi_plot(t,q) = (psi(t))*conj(psi(t));
            w=0;
        end        
        plot(position,psi_plot(:,q), position,VV);
        axis([-1.1,1.1,-0.25,1.25]);
        M(q) = getframe;
        last = q;
    end
    psiold = psi;
    m=m+1;
    w=w+1;
end
w=0;


plot(position,psi_plot(:,1));
axis([-1.1,1.1,-0.25,1.25]);
M(1) = getframe;

movie(M,1)
title('Initial wavepacket');
figure
plot(position,psi_plot(:,1),position, VV, position, psi_plot(:,q),position, psi_plot(:,85))
axis([-1.1,1.1,-0.25,1.25])


end
