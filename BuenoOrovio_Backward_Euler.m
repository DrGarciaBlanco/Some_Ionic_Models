%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script computes a time step for the voltage souce term according
% to the BuenoOrovio-Cherry-Fenton model (2008).
%
% Using a first-order implicit scheme (Backward Euler)
%
% Created by Emilio Garcia-Blanco, 06/06/2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%--------------------------------------------------------------------------
% Constants of the model
%--------------------------------------------------------------------------
uo          =  0;
uu          =  1.55;
thev        =  0.3;
thew        =  0.13;
thevmin     =  0.006;
theo        =  0.006;
tauv1min    =  60;
tauv2min    =  1150;
tauvplus    =  1.4506;
tauw1min    =  60;
tauw2min    =  15;
kwmin       =  65;
uwmin       =  0.03;
tauwplus    =  200;
taufi       =  0.11;
tauo1       =  400;
tauo2       =  6;
tauso1      =  30.0181;
tauso2      =  0.9957;
kso         =  2.0458;
uso         =  0.65;
taus1       =  2.7342;
taus2       =  16;
ks          =  2.0994;
us          =  0.9087;
tausi       =  1.8875;
tauwinf     =  0.07;
winfast     =  0.94;
StimAm      =  1;


%--------------------------------------------------------------------------
% Initialising arrays and initial conditions 
%--------------------------------------------------------------------------
dt          =  0.001;
t_max       =  1000;
steps       =  round(t_max/dt)+1;

U           =  zeros(steps,1);
V           =  zeros(steps,1);
W           =  zeros(steps,1);
S           =  zeros(steps,1);
PHI         =  zeros(steps,1);
ISI         =  zeros(steps,1);
ISO         =  zeros(steps,1);
IFI         =  zeros(steps,1);
T           =  linspace(0,t_max,steps);

u           =  0;
u1          =  0;
v           =  1;
w           =  1;  
s           =  0; 

U(1)        =  u;
V(1)        =  v;
W(1)        =  w;
S(1)        =  s;
PHI(1)      =  u*85.7-84;
ISI(1)      =  0;
ISO(1)      =  0;
IFI(1)      =  0;


%--------------------------------------------------------------------------
% Evolution of the voltage aand internal variables up to t_max
%--------------------------------------------------------------------------

for i=2:steps
    %----------------------------------------------------------------------
    % Computing gating variables and parameters
    %----------------------------------------------------------------------
    
    iter   =  0;
    Res    =  1e6;

    Istim  = 0;
    if (T(i-1)<2); Istim  = -StimAm; end;
    
    while abs(Res)>1e-7 && iter < 14
        
        tauwmin = tauw1min + (tauw2min-tauw1min)*(1+tanh(kwmin*(u-uwmin)))*0.5;
        tauso   = tauso1   + (tauso2-tauso1)*(1+tanh(kso*(u-uso)))*0.5;
       
        if (u1>thev)
            
            v1      =  tauvplus*v/(dt+tauvplus);
            w1      =  tauwplus*w/(dt+tauwplus);
            s1      =  (2*s*taus2+(tanh(ks*(u1-us))+1)*dt)/(2*(taus2+dt)); 
            Isi     =  - w1*s1/tausi;
            Iso     =  1/tauso;
            Ifi     =  (thev-u1)*(uu-u1)*v1/taufi;
            dIt     =  v1*(2*u1-uu-thev)/taufi + (2*kso*(tauso1-tauso2)*sech(kso*(u1-uso))^2)/(2*tauso)^2 - w1/tausi*(ks*(sech(ks.*(u1-us))).^2*dt)./(2*(taus2+dt));

        elseif (u1>thew)
            
            v1      =  tauv2min*v/(dt+tauv2min);
            if (u1<thevmin); v1 = (dt+v*tauv1min)/(dt+tauv1min); end;
            w1      =  tauwplus*w/(dt+tauwplus);
            s1      =  (2*s*taus2+(tanh(ks*(u1-us))+1)*dt)/(2*(taus2+dt)); 
            Isi     =  - w1*s1/tausi;
            Iso     =  1/tauso; 
            Ifi     =  0;
            dIt     =  (2*kso*(tauso1-tauso2)*sech(kso*(u1-uso))^2)/(2*tauso)^2 - w1/tausi*(ks*(sech(ks.*(u1-us))).^2*dt)./(2*(taus2+dt));

        elseif (u1>theo)
            
            v1      =  tauv2min*v/(dt+tauv2min);
            if (u1<thevmin); v1 = (dt+v*tauv1min)/(dt+tauv1min); end;
            w1      =  (tauwmin*w+dt*winfast)/(dt+tauwmin);
            s1      =  (2*s*taus1+(tanh(ks*(u1-us))+1)*dt)/(2*(taus1+dt)); 
            Isi     =  0;
            Iso     =  (u1-uo)/tauo2;
            Ifi     =  0;
            dIt     =  1/tauo2;

        else
            
            v1      =  (dt+v*tauv1min)/(dt+tauv1min); 
            w1      =  (tauwmin*w+dt*(-(u1/tauwinf)+1))/(dt+tauwmin);
            s1      =  (2*s*taus1+(tanh(ks*(u1-us))+1)*dt)/(2*(taus1+dt));
            Isi     =  0;
            Iso     =  (u-uo)/tauo1;
            Ifi     =  0;
            dIt     =  1/tauo1;

        end
        
        Res         =  (u1 - u)/dt + (Isi + Iso + Ifi + Istim);
        u1          =  u1 - Res/(1/dt + dIt);
        iter        =  iter + 1;
    
    end

    %----------------------------------------------------------------------
    %  Save in memory
    %----------------------------------------------------------------------
    U(i)         =  u1;
    V(i)         =  v1;
    W(i)         =  w1;
    S(i)         =  s1;
    PHI(i)       =  u1*85.7-84;
    ISI(i)       =  Isi;
    ISO(i)       =  Iso;
    IFI(i)       =  Ifi;

    %----------------------------------------------------------------------
    % Rename variables for the next time-step
    %----------------------------------------------------------------------
    v            =  v1;
    w            =  w1;
    s            =  s1;
    u            =  u1;
    
end


%--------------------------------------------------------------------------
% Postprocessing
%--------------------------------------------------------------------------
set(groot, 'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');  
set(groot, 'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');  
set(gca, 'FontSize', 10)

figure(1)
subplot(2,2,1)
plot(T,U,'-k','LineWidth',1.5)
set(gcf,'color','w');
xlim([0,t_max])
ylim([0,1.5])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(2,2,1);
xlabel('Time (ms)','FontSize',12) 
ylabel('Voltage $$u$$','FontSize',12) 


subplot(2,2,2)
plot(T,V,'-k','LineWidth',1.5)
set(gcf,'color','w');
xlim([0,t_max])
ylim([0,1])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(2,2,2);
xlabel('Time (ms)','FontSize',12) 
ylabel('Gate variable $$v$$','FontSize',12) 


subplot(2,2,3)
plot(T,W,'-k','LineWidth',1.5)
set(gcf,'color','w');
xlim([0,t_max])
ylim([0,1])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(2,2,3);
xlabel('Time (ms)','FontSize',12) 
ylabel('Gate variable $$w$$','FontSize',12) 


subplot(2,2,4)
plot(T,S,'-k','LineWidth',1.5)
set(gcf,'color','w');
xlim([0,t_max])
ylim([0,1])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(2,2,4);
xlabel('Time (ms)','FontSize',12) 
ylabel('Gate variable $$s$$','FontSize',12) 

hFig  =  figure(1);
set(hFig,'Position', [0 0 800 720])




figure(2)
subplot(2,2,1)
plot(T,PHI,'-k','LineWidth',1.5)
set(gcf,'color','w');
xlim([0,t_max])
ylim([-90,60])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(2,2,1);  
xlabel('Time (ms)','FontSize',12) 
ylabel('Voltage $$\phi$$ (mV)','FontSize',12) 


subplot(2,2,2)
plot(T,ISI,'-k','LineWidth',1.5)
set(gcf,'color','w');
xlim([0,t_max])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(2,2,2);
xlabel('Time (ms)','FontSize',12) 
ylabel('Intensity $$I_{si}$$ (mA)','FontSize',12) 


subplot(2,2,3)
plot(T,ISO,'-k','LineWidth',1.5)
set(gcf,'color','w');
xlim([0,t_max])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(2,2,3);
xlabel('Time (ms)','FontSize',12) 
ylabel('Intensity $$I_{so}$$ (mA)','FontSize',12) 


subplot(2,2,4)
plot(T,IFI,'-k','LineWidth',1.5)
set(gcf,'color','w');
xlim([0,10])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(2,2,4);
xlabel('Time (ms)','FontSize',12)
ylabel('Intensity $$I_{fi}$$ (mA)','FontSize',12) 


hFig  =  figure(2);
set(hFig,'Position', [0 0 800 720])