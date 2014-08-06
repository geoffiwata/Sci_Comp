function [] = QM_box()
%HW 5-2. The function will calculate the wavefunction.


clear all
a = 0.02; %spacing of grid
length = 1; %region length
n = (length ./ a); %row length
N = n^2; %total number of grid points
m = .1; %particle mass

%position of source
xs = 0.51;
ys = 0.26;
%position of square
xsq = 0.625;
ysq = 0.75;
% 3.75<jsq<8.75


%x,y positions of all the grid points
for j = 1:n
    for k= 1:n
        g = (k-1)*n +j;
        gpos(g,1) = a*j;
        gpos(g,2) = a*k;
    end
end


% Make the laplacian operator
D = sparse(1:N,1:N,4,N,N);
B = sparse(2:N,1:N-1,-1,N,N);
C = sparse(1:N-n,n+1:N,-1,N,N);
A = (1/(2*m))*(D+B+B'+C+C')./(a^2);
%Lap = full(Laps);



psi = zeros(N,1); %the potential we are solving for

rho = zeros(N,1); %the sources

gdo = zeros(N,1); %whether we are at a boundary point or not(not 0 means we are)

%distance from source to each point
for p=1:N
    dist(p,1)=sqrt((gpos(p,1)-xs)^2 + (gpos(p,2)-ys)^2);
    if dist(p,1)<=(a*.75)
        gdo(p,1)=3;
        %Make the source
    end  
end



for j=1:n
    for k=1:n
        g = (k-1)*n +j;        
        if k==1||j==1||j==n||k==n
            gdo(g,1)=1;
            psi(g) = 0;
            A(:,g) = 0;
        end
%        Make the box - infinite potential means psi is zero
        if j >= (0.5*n) && j <= (0.75*n) && k >= (n*0.625) && k <= (.875*n)
            gdo(g,1)=2;
            psi(g,1)=0;
            A(:,g) = 0;
        end 
    end
end

mu = 170;
v = zeros(N,1);
v(1,1) =1;
Ar = (A-mu*eye(N,N));
f=0;
while f<100
    i=0; %iteration
    t=50; %tolerance
    W = ones(N,1);
    r = -v + Ar*W; %residual
    p = r;
    rr = r'*r;
    while i < t
        %Set up variables
        Arp=Ar*p;
        alpha = rr / (p'*Arp);
        T=alpha*p;
        W = W + T;
        r = Ar*W - alpha*Arp; %update r
        rrnew=r'*r;
        %If r is too small, quit.
        if sqrt(rrnew)<1e-5
            i
            break;
        end
        p = r + rrnew*p./(rr); %update p
        rr=rrnew;
        i=i+1;
    end
    v=W./norm(W);
    d = (A*v)./v;
    eigenvalue = d(4);
    
    f=f+1;  
end


[V,D] = eig(full(A));
[a,b] = size(V);

energies = [10000,10000,10000];
which = [0,0,0];
for p=1:b
    if D(p,p) ~= 0
        if D(p,p) < energies(1)
            energies(2) = energies(1);
            energies(1) = D(p,p);
            which(1) = p;
        elseif D(p,p) < energies(2)
            energies(3) = energies(2);
            energies(2) = D(p,p);            
            which(2) = p;
        elseif D(p,p) < energies(3)
            energies(3) = D(p,p);
            which(3) = p;
        end
    end
end

energies
which
        
        
for j=1:n
    for k=1:n
        g = (k-1)*n +j;
        if gdo(g) == 1 || gdo(g) == 2
            V(g,:) = 0;
            v(g) = 0;
        end
        psi_plot1(j,k) = (V(g,which(1)))^2;
        psi_plot2(j,k) = (V(g,which(2)))^2;
        psi_plot3(j,k) = (V(g,which(3)))^2;
        psi_plot4(j,k) = (v(g))^2;
    end
end
figure
surf(psi_plot1)
figure
surf(psi_plot2)
figure
surf(psi_plot3)
figure
surf(psi_plot4)



end
