%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script computes a time step for the voltage souce term according
% to the Aliev-Panfilov model (1996).
%
% Implemented using a second-order implicit scheme (Crank Nicolson)
%
% Created by Emilio Garcia-Blanco, 06/06/2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%--------------------------------------------------------------------------
% Constants of the model
%--------------------------------------------------------------------------
k           =  8;
alpha       =  0.15; 
beta        =  0.20;
gamma       =  0.30;
StimAm      =  2;

%--------------------------------------------------------------------------
% Initialising arrays and initial conditions 
%--------------------------------------------------------------------------
dt          =  0.01;
t_max       =  1000;
steps       =  round(t_max/dt)+1;
scale       =  10;

U           =  zeros(steps,1);
V           =  zeros(steps,1);
PHI         =  zeros(steps,1);
T           =  linspace(0,t_max,steps);

u           =  0;
u1          =  0;
v           =  0;
v1          =  0;
In1         =  0;

U(1)        =  u;
V(1)        =  v;
PHI(1)      =  u*110-86.4;

%--------------------------------------------------------------------------
% Evolution of the voltage aand internal variables up to t_max
%--------------------------------------------------------------------------

for i=2:steps
    
    iter   =  0;
    Res    =  1e6;

    Istim  = 0;
    if (T(i-1)<2); Istim  = -StimAm; end;
    
    while abs(Res)>1e-10 && iter < 14  
        %------------------------------------------------------------------
        % Computing parameters
        %------------------------------------------------------------------
        epsil   = 0.002;
        %epsil  = 1; if (u1>0.05); epsil = 0.1; end;
        
        aux         =  (1 + dt/scale/2*(epsil + beta*k*u1*(u1-alpha-1)/(u1+gamma)));
        v1          =  (-aux + sqrt(aux^2 - dt/scale*beta/(u1+gamma)*(dt/scale*epsil*k*u1*(u1-alpha-1)- 2*v + dt/scale*(epsil + beta*v/(u+gamma))*(v + k*u*(u-alpha-1))) ))/(dt/scale*beta/(u1+gamma));
        dvdu        =  beta*v1/(u1+gamma)*(v1 + k*u1*(u1-alpha-1)) - k*(2*u1-alpha-1)*(epsil*(u1+gamma)+beta*v1)/(beta*(v1+k*u1*(u1-alpha-1)) + (2/dt/scale+epsil)*(u1+gamma));
        
        Itotal      =  k*u1*(u1-alpha)*(u1-1) + u1*v1;
        dIt         =  k*(u1*(3*u1-2) - alpha*(2*u1-1)) + v1 + u1*dvdu;
    
        Res         =  (u1 - u)/dt + 0.5*(Itotal + Istim + In1);
        u1          =  u1 - Res/(1/dt + 0.5*dIt);
        iter        =  iter + 1;
    
    end

    %----------------------------------------------------------------------
    %  Save in memory
    %----------------------------------------------------------------------
    U(i)         =  u1;
    V(i)         =  v1;
    PHI(i)       =  u1*110-86.4;

    %----------------------------------------------------------------------
    % Rename variables for the next time-step
    %----------------------------------------------------------------------
    v            =  v1;
    u            =  u1;
    In1          =  Itotal + Istim;
    
end


%--------------------------------------------------------------------------
% Postprocessing
%--------------------------------------------------------------------------
set(groot, 'defaulttextinterpreter',         'latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter',       'latex');  
set(groot, 'defaulttextinterpreter',         'latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter',       'latex');  
set(gca,   'FontSize', 8)


figure(1)
subplot(1,2,1)
plot(T,U,'-k','LineWidth',1.5)
set(gcf,'color','w');
xlim([0,t_max])
ylim([0,1])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(1,2,1);
xlabel('Time (ms)','FontSize',12) 
ylabel('Voltage $$u$$','FontSize',12) 

subplot(1,2,2)
plot(T,V,'-k','LineWidth',1.5)
set(gcf,'color','w');
xlim([0,t_max])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(1,2,2);
xlabel('Time (ms)','FontSize',12) 
ylabel('Refractoriness variable $$v$$','FontSize',12) 

hFig  =  figure(1);
set(hFig,'Position', [0 0 800 420])


figure(2)
plot(T,PHI,'-k','LineWidth',1.5)
set(gcf,'color','w');
xlim([0,t_max])
ylim([-90,60])
grid on
set(gca,'GridLineStyle','--')
hFig  =  figure(2);  
xlabel('Time (ms)','FontSize',12) 
ylabel('Voltage $$\phi$$ (mV)','FontSize',12) 

