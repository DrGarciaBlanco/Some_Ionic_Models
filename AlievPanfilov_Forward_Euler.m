%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script computes a time step for the voltage souce term according
% to the Aliev-Panfilov model (1996).
%
% Implemented using a first-order explicit scheme (Forward Euler)
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
v           =  0;

U(1)        =  u;
V(1)        =  v;
PHI(1)      =  u*110-86.4;

%--------------------------------------------------------------------------
% Evolution of the voltage aand internal variables up to t_max
%--------------------------------------------------------------------------

for i=2:steps
    %----------------------------------------------------------------------
    % Computing parameters
    %----------------------------------------------------------------------
    Istim  = 0; if (mod(T(i-1),1000)<10); Istim  = -StimAm; end;
    
    Itotal       =  k*u*(u - alpha)*(u - 1) + u*v;
    
    epsil   = 0.002;
    %epsil  = 1; if (u>0.05); epsil = 0.1; end;
    v1           =  v - dt/scale*(epsil + beta*v/(u+gamma))*(v+k*u*(u-alpha-1));
    u1           =  u - dt*(Itotal + Istim);

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
set(gca,   'FontSize', 10)


figure(1)
subplot(1,2,1)
plot(T,U,'-k','LineWidth',1.5)
set(gcf,'color','w');
xlim([0,t_max])
ylim([-0.4,1])
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
ylim([-120,60])
grid on
set(gca,'GridLineStyle','--')
hFig  =  figure(2);  
xlabel('Time (ms)','FontSize',12) 
ylabel('Voltage $$\phi$$ (mV)','FontSize',12) 

