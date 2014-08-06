function [ p ] = argon_md_MC()
close all

NP = 864;      %  Number of particles
tau = 0.032;    %  Time step
tauhalf = tau/2;
tau2half = tau * tau / 2;
Tsteps = 1200;       %  Number of time steps
T = Tsteps + 1;
den = 0.75;     %  density
L = (NP/den)^(1/3);
Li = 1/L;
Lo2i = 2/L;
UpdatePlotInterval = 200;
NPlot = 5;      %  Number of particles to track
NPlotIndex = floor( [1 : NPlot ] * 864/(NPlot + 1));    %  Indicies for particle(s) to track in 3d
InitFigures = 1;
mu = 150; %stabalized iteration step

p.NP = NP;
p.tau = tau;
p.T = T;
p.den = den;
p.L = L;

fprintf('Number of particles = %d\n', NP);
fprintf('Density = %5.3f\n', den);
fprintf('Linear size of box = %6.3f\n', L);

fprintf('Index of particles whose location is being tracked in 3d\n')

for i=1:NPlot
    fprintf('%d\n',NPlotIndex(i));
end

%  A structure array for the particles
%  p.p contains the particle's position at a given time
%  p.v contains the particle's velocity at a given time
%  p.vh contains the particle's velocity at a given time - 1/2
%  p.f contains the force on a particle at a given time
%  p.virial contains the virial at a given time

p.p = zeros(3,NP,T);
p.v = zeros(3,NP,T);
p.vh = zeros(3,NP,T);
p.f = zeros(3,NP,T);
p.virial = zeros(T,1);
p.ke = zeros(T,1);
p.pe = zeros(T,1);
p.e = zeros(T,1);

%  Put the particles roughly uniformly spaced
%  Find integer number of particles in each dimensions (NL)
%  that ensures NP particles in box.

NL = floor(NP^(1/3)) + 1;
sep = L/NL;
rsep = 0.1;
n = 1;
for i=1:NL
    for j=1:NL
        for k=1:NL
            if (n <= NP )
            p.p(1,n,1) = ((i - 0.5) + rsep * (rand(1) - 0.5) ) * sep;
            p.p(2,n,1) = ((j - 0.5) + rsep * (rand(1) - 0.5) ) * sep;
            p.p(3,n,1) = ((k - 0.5) + rsep * (rand(1) - 0.5) ) * sep;
            
            p.v(1,n,1) = 0.105*(rand(1) - 0.5);
            p.v(2,n,1) = 0.105*(rand(1) - 0.5);
            p.v(3,n,1) = 0.105*(rand(1) - 0.5);
                        
            end
            n = n + 1;
        end
    end
end
KE(1);
Temp(1) = p.ke(1)/NP;
tic, PEaF(1), toc;
p.e(1) = p.ke(1) + p.pe(1);

PrintE(1);
for t=2:T
    %update all positions
    p.p(:,:,t) = p.p(:,:,t-1);
    p.v(:,:,t) = p.v(:,:,t-1);
    %save the positions/velocities of all the particles:
    p.po(:,:,t) = p.p(:,:,t-1);
    p.vo(:,:,t) = p.v(:,:,t-1);
     
    %make some random numbers
    which1 =ceil((NP).*rand());
    which2 =ceil((NP).*rand());
    dispx = 0.1*(rand()-0.5);
    dispy = 0.1*(rand()-0.5);
    dispz = 0.1*(rand()-0.5);
    vdispx = .5*(rand()-0.5);
    vdispy = .5*(rand()-0.5);
    vdispz = .5*(rand()-0.5);
    
    %make new positions for this particle
    p.p(1,which1,t) = p.p(1,which1,t-1) + dispx;
    p.p(2,which1,t) = p.p(2,which1,t-1) + dispy;
    p.p(3,which1,t) = p.p(3,which1,t-1) + dispz;
    %  Make particles have coordinates between 0 and L    
    p.p(:,which1,t) = p.p(:,which1,t) - floor(p.p(:,which1,t) * Li ) * L;        


    %make new velocities for this particle
    p.v(1,which2,t) = p.v(1,which2,t) + vdispx;
    p.v(2,which2,t) = p.v(2,which2,t) + vdispy;
    p.v(3,which2,t) = p.v(3,which2,t) + vdispz;

    %  Calculate new Potential energy and force    
    PEaF(t);
    %   Calculate new kinetic energy
    KE(t);
        Temp(t) = 16*p.ke(t)/NP;
    p.e(t) = p.ke(t) + p.pe(t);
    
    R = rand();
    deltV =p.e(t)-p.e(t);
    cond = exp(-deltV./Temp(t));
    if deltV >= 0 &&  cond>=R
        p.p(1,which1,t) = p.p(1,which1,t-1);
        p.p(2,which1,t) = p.p(2,which1,t-1);
        p.p(3,which1,t) = p.p(3,which1,t-1);
      p.v(:,which2,t)= p.vo(:,which2,t); 
    end
    
        
    PrintE(t)
    
end

EvolutionPlots(t)
k=2;
St=Tsteps-mu-k; %steps after stabilized point
tempavg = mean(Temp)
PEavg = mean(p.pe)
virialavg = mean(p.virial)

for f = 1:Tsteps-k
    Tacf(f) = (Temp(f)-tempavg)*(Temp(f+k)-tempavg);
    Tnorm = (Temp(f)-tempavg)^2;
    PEacf(f)= (p.pe(f)-PEavg)*(p.pe(f+k)-PEavg);
    PEnorm = (p.pe(f)-PEavg)^2;
    virialacf(f) = (p.virial(f)-virialavg)*(p.virial(f+k)-virialavg);
    VIRnorm = (p.virial(f)-virialavg)^2;
    
    GammaT(t) = Tacf(f)/Tnorm;
    GammaPE(t) = PEacf(f)/PEnorm;
    GammaVIR(t) = virialacf(f)/VIRnorm;
    
end

tauT = tau*trapz(GammaT)
tauPE = tau*trapz(GammaPE);
tauVIR = tau*trapz(GammaVIR);
alphT = sqrt(2*tauT);
alphPE = sqrt(2*tauPE);
alphVIR = sqrt(abs(2*tauVIR));

varT = std(Temp);
varPE = std(p.pe);
varVIR = std(p.virial);
Temperr = varT*alphT
PEerr = varPE*alphPE
VIRerr = varVIR*alphVIR


    function [] = PEaF(tt)
        
        for n=1:NP-1
            for m=n+1:NP
                
                r2 = 0.0;
                
                for i=1:3
                    
                    dr(i) = p.p(i,n,tt) - p.p(i,m,tt);
                    dr(i) = dr(i) - floor(abs(dr(i) * Lo2i)) * L * sign(dr(i));
                    r2 = r2 + dr(i) * dr(i);
                    
                end
                
                r2i = 1.0/r2;
                r6i = r2i * r2i * r2i;
                
                %  Accumulate potential energy
                
                p.pe(tt) = p.pe(tt) + r6i * r6i - r6i;
                
                %  Find force
                
                ft = r6i * r2i * ( r6i - 0.5);
                
                for i=1:3
                    p.f(i,n,tt) = p.f(i,n,tt) + ft * dr(i);
                    p.f(i,m,tt) = p.f(i,m,tt) - ft * dr(i);
                end
                
            end
        end
      
        p.pe(tt) = 4.0 * p.pe(tt);
        
        p.virial(tt) = sum( p.p(:,tt) .* p.f(:,tt));
        
    end

    function [] = KE(tt)
        
        for n=1:NP
            for i=1:3
                p.ke(tt) = p.ke(tt) + p.v(i,n,tt) * p.v(i,n,tt);
            end
        end
        p.ke(tt) = 24.0 * p.ke(tt);
        
    end

    function [] = PrintE(tt)
        
        fprintf('Time step %d\tKE = %13.7e\tPE = %13.7e\tE = %13.7e\n Temp = %6.3e\n', ...
            tt-1, p.ke(tt), p.pe(tt), p.e(tt),Temp(tt));
        
    end

    function [] = EvolutionPlots(tt)
        
        %  Save figure handles between calls, to close and reopen figures.
        persistent fh0 fh1 fh2;
        
        mdtime = repmat([0:tt-1]', 1, 3 );
        EtoPlot = [p.e(1:tt,1), p.ke(1:tt,1), p.pe(1:tt,1)];
        
        if (InitFigures == 0)
            close(fh0);
            close(fh1);
            close(fh2);
        else
            InitFigures = 0;
        end
        
        fh0 = figure('Units','inches','Position', [1 7 7 6]);
        plot(mdtime,EtoPlot);
        title('Evolution of Energies with MD time');
        xlabel('MD time');
        ylabel('Energy');
        legend('Total Energy','Kinetic Energy','Potential Energy','Location','NorthWest');
        drawnow; pause(1);
        
        x3d = squeeze(p.p(1,NPlotIndex,1:tt))';
        y3d = squeeze(p.p(2,NPlotIndex,1:tt))';
        z3d = squeeze(p.p(3,NPlotIndex,1:tt))';
        
        fh1 = figure('Units','inches','Position', [9 7 7 6]);
        plot3(x3d,y3d,z3d);
        title('Evolution of sample particles with MD time');
        drawnow; pause(1);
        
        fh2 = figure('Units','inches','Position', [17 7 7 6]);
        plot(0:tt-1,p.virial(1:tt));
        title('Evolution of virial with MD time');
        xlabel('MD time');
        ylabel('Virial');
        drawnow; pause(1);
        
    end
        

end