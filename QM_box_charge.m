function [] = QM_box_charge()
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
    
        gdo(p,1)=3;
        %Make the source (-143 for 90%)
        U(p,1)= (-143)./dist(p,1);
        A(p,p) = A(p,p)+U(p);
     
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

[V,D] = eig(full(A));
[a,b] = size(V);

energies = [10000,10000,10000];
which = [0,0,0];
%find the three lowest energy eigenvalues
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
        
prob=zeros(N,1); %probability estimator for lower half   
        
for j=1:n
    for k=1:n
        g = (k-1)*n +j;
        if gdo(g) == 1 || gdo(g) == 2
            V(g,:) = 0;
        end
        psi_plot1(j,k) = (V(g,which(1)))^2;
%         psi_plot2(j,k) = (V(g,which(2)))^2;
%         psi_plot3(j,k) = (V(g,which(3)))^2;
        if j<=0.5*n
            prob(g) = (V(g,which(1)))^2;
        end
    end
end
%what is probability that particle is found in lower half?
totalprob = sum((V(:,which(1)))'*(V(:,which(1))))
plow = sum(prob)/totalprob 


figure
surf(psi_plot1)




end
