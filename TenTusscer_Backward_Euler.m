%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script computes a time step for the voltage souce term according
% to the Ten Tusscher model (2004).
%
% Implemented using a first-order implicit scheme (Backward Euler). For
% more information about its implementation, the reader is referred to the
% Appendix D of my PhD thesis.
%
% Created by Emilio Garcia-Blanco, 06/06/2019.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%--------------------------------------------------------------------------
% Constants of the model (values for epicardium)
%--------------------------------------------------------------------------
R          =  8.3143;           %(mJ K mmol-1)
T          =  310;              %(K)
F          =  96.4867;          %(C/mmol)
Cm         =  0.000185;         %(muF cm-2)
S          =  0.2;              %(mum-1)
pho        =  162;              %(Ohm cm)
VC         =  0.016404;         %(cm3)
VSR        =  0.001094;         %(cm3)
KO         =  5.4;              %(mM) 
NaO        =  140;              %(mM) 
CaO        =  2;                %(mM) 
GNa        =  14.838;           %(nS/pF) 
GK1        =  5.405;            %(nS/pF) 
Gto        =  0.294;            %(nS/pF)    (0.073 for endocardium)
GKr        =  0.096;            %(nS/pF) 
GKs        =  0.245;            %(nS/pF)  (0.245 for epi and endocardium)
pKNa       =  0.03;             %(-)
GCaL       =  1.75e-1;          %(mm3 muF-1 s-1) 
kNaCa      =  1000;             %(pA/pF)
gamma      =  0.35;             %(-)
KmCa       =  1.38;             %(mM)
KmNai      =  87.5;             %(mM)
ksat       =  0.1;              %(-)
alpha      =  2.5;              %(-)
PNaK       =  1.362;            %(pA/pF)
KmK        =  1;                %(mM) 
KmNa       =  40;               %(mM) 
GpK        =  0.0146;           %(nS/pF) 
GpCa       =  0.825;            %(nS/pF) 
KpCa       =  0.0005;           %(mM) 
GbNa       =  0.00029;          %(nS/pF) 
GbCa       =  0.000592;         %(nS/pF) 
Vmaxup     =  0.000425;         %(mM/ms) 
Kup        =  0.00025;          %(mM)
arel       =  0.016464;         %(mM/ms) 
brel       =  0.25;             %(mM) 
crel       =  0.008232;         %(mM/ms) 
Vleak      =  0.00008;          %(ms-1)  
Bufc       =  0.15;             %(mM)
Kbufc      =  0.001;            %(mM) 
Bufsr      =  10;               %(mM)
Kbufsr     =  0.3;              %(mM) 
StimAm     =  72;        
StimDu     =  2; 
StimSt     =  0.3;

%--------------------------------------------------------------------------
% Initialising arrays and initial conditions 
%--------------------------------------------------------------------------
dt               =  0.01;
t_max            =  1000;
steps            =  round(t_max/dt)+1;

Volt             =  zeros(steps,1);

gate_m           =  zeros(steps,1);
gate_g           =  zeros(steps,1);
gate_d           =  zeros(steps,1);
gate_h           =  zeros(steps,1);
gate_f           =  zeros(steps,1);
gate_s           =  zeros(steps,1);
gate_r           =  zeros(steps,1);
gate_j           =  zeros(steps,1);
gate_fCa         =  zeros(steps,1);
gate_xs          =  zeros(steps,1);
gate_xr1         =  zeros(steps,1);
gate_xr2         =  zeros(steps,1);
gate_xK1inf      =  zeros(steps,1);

ion_CaSR         =  zeros(steps,1);
ion_Ca           =  zeros(steps,1);
ion_K            =  zeros(steps,1);
ion_Na           =  zeros(steps,1);

inten_Ina        =  zeros(steps,1);   
inten_ICaL       =  zeros(steps,1);   
inten_Ito        =  zeros(steps,1);  
inten_IKs        =  zeros(steps,1);  
inten_IKr        =  zeros(steps,1);  
inten_IK1        =  zeros(steps,1);  
inten_INaCa      =  zeros(steps,1);
inten_INaK       =  zeros(steps,1); 
inten_IpCa       =  zeros(steps,1);
inten_IpK        =  zeros(steps,1); 
inten_IbNa       =  zeros(steps,1); 
inten_IbCa       =  zeros(steps,1); 
inten_Ileak      =  zeros(steps,1); 
inten_Iup        =  zeros(steps,1); 
inten_Irel       =  zeros(steps,1); 

Time             =  linspace(0,t_max,steps);


Iax              =  0;
V                =  -86.2;
V1               =  -86.2;       
m1               =  0.000;
h1               =  0.750;
j1               =  0.750;
xr11             =  0.000;
xr21             =  0.000;
xs1              =  0.000;
r1               =  0.000;
s1               =  1.000;
d1               =  0.000;
f1               =  1.000;
xK1inf1          =  0.052;
fCa1             =  1.000;
g1               =  1.000;
Nai1             =  11.60;
Ki1              =  138.3;
Cai1             =  8e-5; 
Casr1            =  0.560;
m                =  0.000;
h                =  0.750;
j                =  0.750;
xr1              =  0.000;
xr2              =  0.000;
xs               =  0.000;
r                =  0.000;
s                =  1.000;
d                =  0.000;
f                =  1.000;
xK1inf           =  0.052;
fCa              =  1.000;
g                =  1.000;
Nai              =  11.60;
Ki               =  138.3;
Cai              =  8e-5; 
Casr             =  0.560;   


gate_m(1)        =  m1;
gate_g(1)        =  g1;
gate_d(1)        =  d1;
gate_h(1)        =  h1;
gate_f(1)        =  f1;
gate_s(1)        =  s1;
gate_r(1)        =  r1;
gate_j(1)        =  j1;
gate_fCa(1)      =  fCa1;
gate_xs(1)       =  xs1;
gate_xr1(1)      =  xr11;
gate_xr2(1)      =  xr21;
gate_xK1inf(1)   =  xK1inf1;

ion_CaSR(1)      =  Casr1;
ion_Ca(1)        =  Cai1;
ion_K(1)         =  Ki1;
ion_Na(1)        =  Nai1;


%--------------------------------------------------------------------------
% Evolution of the voltage aand internal variables up to t_max
%--------------------------------------------------------------------------

for i=2:steps

%--------------------------------------------------------------------------
% Newton-Raphson algorithm for gating variables and ion concentrations.                                                                                                                                                                                                                  
%-------------------------------------------------------------------------- 
Var       =  [m;h;j;xr1;xr2;xs;r;s;d;f;xK1inf;fCa;g;Nai;Ki;Cai;Casr];
Res       =  1e6;
iter      =  0;

while Res>1e-11
    
%--------------------------------------------------------------------------
% Computing Nerst potenials                                                                                                                                                                                                                   
%-------------------------------------------------------------------------- 
ENa     = R*T/F*log(NaO/Nai);
EK      = R*T/F*log(KO/Ki);
ECa     = R*T/2/F*log(CaO/Cai);
EKs     = R*T/F*log((KO + pKNa*NaO)/(Ki + pKNa*Nai));

%--------------------------------------------------------------------------
% Computing gating variables                                                                                                                                                                                                                   
%-------------------------------------------------------------------------- 
minf    = 1/(1 + exp((-56.86-V)/9.03))^2;
alpham  = 1/(1 + exp((-60-V)/5));
betam   = 0.1/(1 + exp((V+35)/5)) + 0.1/(1 + exp((V-50)/200));
taum    = alpham*betam;

hinf    = 1/(1 + exp((V+71.55)/7.43))^2;
tauh    = (0.057*exp(-(V+80)/6.8)+2.7*exp(0.079*V)+3.1e5*exp(0.3485*V)).^(-1); if (V>=-40) tauh = (0.77./0.13./(1 + exp(-(V+10.66)/11.1))).^(-1); end;

jinf    = 1/(1 + exp((V+71.55)/7.43))^2;
alphaj  = (-2.5428e4*exp(0.2444*V) - 6.948e-6*exp(-0.04391*V))*(V + 37.78)/(1 + exp(0.311*(V+79.23))); if (V>=-40) alphaj = 0; end;
betaj   = 0.6*exp(0.057*V)/(1 + exp(-0.1*(V+32))); if (V<-40) betaj = 0.02424*exp(-0.01052*V)/(1 + exp(-0.1378*(V+40.14))); end;
tauj    = 1/(alphaj + betaj);

dinf    = 1/(1 + exp(-(V+5)/7.5));
alphad  = 1.4/(1 + exp(-(35+V)/13)) + 0.25;
betad   = 1.4/(1 + exp((V+5)/5));
gammad  = 1/(1 + exp((50-V)/20));
taud    = alphad*betad + gammad;

finf    = 1/(1 + exp((V+20)/7));
tauf    = 1125*exp(-(V+27)^2/240) + 165/(1 + exp((25-V)/10)) + 80;

rinf    = 1/(1 + exp((20-V)/6));
taur    = 9.5*exp(-(V+40)^2/1800) + 0.8;

sinf    = 1/(1 + exp((V+20)/5));
taus    = 85*exp(-(V+45)^2/320) + 5/(1 + exp((V-20)/5)) + 3;

xsinf   = 1/(1 + exp(-(V+5)/14));
alphaxs = 1100/sqrt(1 + exp(-(V+10)/6)); 
betaxs  = 1/(1 + exp((V-60)/20));
tauxs   = alphaxs*betaxs;

xr1inf  = 1/(1 + exp(-(26+V)/7));
alxr1   = 450/(1 + exp(-(45+V)/10));
bexr1   = 6/(1 + exp((V+30)/11.5));
tauxr1  = alxr1*bexr1;

xr2inf  = 1/(1 + exp((V+88)/24));
alxr2   = 3/(1 + exp(-(V+60)/20));
bexr2   = 1.12/(1 + exp((V-60)/20));
tauxr2  = alxr2*bexr2;

alpfca  = 1/(1 + (Cai/0.000325)^8);
betfca  = 0.1/(1 + exp((Cai-0.0005)/0.0001));
gamfca  = 0.2/(1 + exp((Cai-0.00075)/0.0008));
fcainf  = (alpfca + betfca + gamfca + 0.23)/1.46;
taufca  = 2; kfca       = 1; if (fcainf>fCa && V>-60) kfca = 0; end;

alK1    = 0.1/(1 + exp(0.06*(V-EK-200)));
beK1    = (3*exp(0.0002*(V-EK+100)) + exp(0.1*(V-EK-10)))/(1 + exp(-0.5*(V-EK)));

ginf    = 1/(1 + Cai^6/0.00035^6); if (Cai>=0.00035) ginf = 1/(1 + Cai^16/0.00035^16); end;
taug    = 2; kg       = 1; if (ginf>g && V>-60) kg = 0; end;

%--------------------------------------------------------------------------
% Computing intensities                                                                                                                                                                                                                   
%-------------------------------------------------------------------------- 
INa     = GNa*m^3*h*j*(V-ENa);
ICaL    = GCaL*d*f*fCa*4*V*F^2/R/T*(Cai*exp(2*V*F/R/T) - 0.341*CaO)/(exp(2*V*F/R/T) - 1);
Ito     = Gto*r*s*(V - EK);
IKs     = GKs*xs^2*(V - EKs);
IKr     = GKr*sqrt(KO/5.4)*xr1*xr2*(V - EK);
IK1     = GK1*sqrt(KO/5.4)*xK1inf*(V - EK);
INaCa   = kNaCa*(exp(gamma*V*F/R/T)*Nai^3*CaO -exp((gamma-1)*V*F/R/T)*NaO^3*Cai*alpha)/(KmNai^3 + NaO^3)/(KmCa + CaO)/(1 + ksat*exp((gamma-1)*V*F/R/T));
INaK    = PNaK*KO*Nai/(KO + KmK)/(Nai + KmNa)/(1 + 0.1245*exp(-0.1*V*F/R/T) + 0.0353*exp(-V*F/R/T));
IpCa    = GpCa*Cai/(KpCa + Cai);
IpK     = GpK*(V - EK)/(1 + exp((25-V)/5.98));
IbNa    = GbNa*(V - ENa);
IbCa    = GbCa*(V - ECa);
Ileak   = Vleak*(Casr - Cai);
Iup     = Vmaxup/(1 + Kup^2/Cai^2);
Irel    = (arel*Casr^2/(brel^2 + Casr^2) + crel)*d*g;
Istim   = 0;
if (Time(i-1)<1); Istim  = -StimAm; end;


%--------------------------------------------------------------------------
% Additional variables for ion concentrations                                                                                                                                                                                                                   
%-------------------------------------------------------------------------- 
Caibufc      = 1/(1 + (Bufc*Kbufc)/(Cai+Kbufc)^2);
Casrbufsr    = 1/(1 + (Bufsr*Kbufsr)/(Casr+Kbufsr)^2);

%--------------------------------------------------------------------------
% Computing derivatives respect to the internal variables (gates and concentrations)                                                                                                                                                                                                                 
%-------------------------------------------------------------------------- 
dEKdKi       =  -R*T/F/Ki;
dECadCai     =  -R*T/2/F/Cai;
dENadNai     =  -R*T/F/Nai;
dEKsdKi      =  -R*T/F/(Ki + pKNa*Nai);
dEKsdNai     =  -R*T*pKNa/F/(Ki + pKNa*Nai);

dINadNai     =  (-1)*GNa*m^3*h*j*dENadNai;
dIbNadNai    =  (-1)*GbNa*dENadNai;
dINaKdNai    =  PNaK*KO*KmNa*exp(1.1*V*F/R/T)/(KO + KmK)/(Nai + KmNa)^2/(exp(1.1*V*F/R/T) + 0.1245*exp(V*F/R/T) + 0.0353*exp(0.1*V*F/R/T));
dINaCadCai   =  (-1)*kNaCa*(exp((gamma-1)*V*F/R/T)*NaO^3*alpha)/(KmNai^3 + NaO^3)/(KmCa + CaO)/(1 + ksat*exp((gamma-1)*V*F/R/T));
dINaCadNai   =  kNaCa*(3*exp(gamma*V*F/R/T)*Nai^2*CaO)/(KmNai^3 + NaO^3)/(KmCa + CaO)/(1 + ksat*exp((gamma-1)*V*F/R/T));
dIK1dKi      =  (-1)*GK1*sqrt(KO/5.4)*xK1inf*dEKdKi;
dIKrdKi      =  (-1)*GKr*sqrt(KO/5.4)*xr1*xr2*dEKdKi;
dIKsdKi      =  (-1)*GKs*xs^2*dEKsdKi; 
dIKsdNai     =  (-1)*GKs*xs^2*dEKsdNai;
dIpKdKi      =  (-1)*GpK*dEKdKi/(1 + exp((25-V)/5.98));
dItodKi      =  (-1)*Gto*r*s*dEKdKi;
dICaLdCai    =  GCaL*d*f*V*F^2*(4*fCa*exp(2*V*F/R/T))/R/T/(exp(2*V*F/R/T) - 1); 
dIbCadCai    =  (-1)*GbCa*dECadCai; 
dIpCadCai    =  GpCa*KpCa/(KpCa + Cai)^2;
dIleakdCai   =  (-1)*Vleak;
dIleakdCasri =  Vleak;
dIupdCai     =  Vmaxup*2*Kup^2*Cai/(Kup^2 + Cai^2)^2;
dIreldCasri  =  2*arel*brel^2*Casr/(brel^2 + Casr^2)^2*d*g;

dINadm       =  3*GNa*m^2*h*j*(V-ENa);
dINadh       =  GNa*m^3*j*(V-ENa);
dINadj       =  GNa*m^3*h*(V-ENa);
dIK1dK1      =  GK1*sqrt(KO/5.4)*(V - EK);
dIKrdxr1     =  GKr*sqrt(KO/5.4)*xr2*(V - EK);
dIKrdxr2     =  GKr*sqrt(KO/5.4)*xr1*(V - EK);
dIKsdxs      =  2*GKs*xs*(V - EKs);
dItodr       =  Gto*s*(V - EK);
dItods       =  Gto*r*(V - EK);
dICaLdd      =  GCaL*f*fCa*4*V*F^2/R/T*(Cai*exp(2*V*F/R/T) - 0.341*CaO)/(exp(2*V*F/R/T) - 1);
dICaLdf      =  GCaL*d*fCa*4*V*F^2/R/T*(Cai*exp(2*V*F/R/T) - 0.341*CaO)/(exp(2*V*F/R/T) - 1);
dICaLdfCa    =  GCaL*d*f*4*V*F^2/R/T*(Cai*exp(2*V*F/R/T) - 0.341*CaO)/(exp(2*V*F/R/T) - 1);
dIreldd      =  (arel*Casr^2/(brel^2 + Casr^2) + crel)*g; 
dIreldg      =  (arel*Casr^2/(brel^2 + Casr^2) + crel)*d; 

dK1infdKi    =  (exp(4*EK/25+13)*(3000*exp(14*EK/25+3*V/5+12) + 3300*exp(EK/2+33*V/50) + 500*exp(3*EK/50+11*V/10+12) +...
                7503*exp((3299*EK+2501*V+65100)/5000) + 8403*exp((2999*EK+2801*V+5100)/5000) + 3*exp((799*EK+5001*V+65100)/5000) +...
                903*exp((499*EK+5301*V+5100)/5000) + 800*exp(29*V/25)))/...
                500/(exp(4*EK/25+V/2+13) + 10*exp(3*EK/50+3*V/5+12) + 30*exp((799*EK+2501*V+65100)/5000) +...
                30*exp((499*EK+2801*V+5100)/5000) + exp(33*EK/50+13) + 10*exp(33*V/50))^2*dEKdKi;
dfCainfdCai  =  -8*Cai^7*0.000325^8/(Cai^8 + 0.000325^8)^2/1.46 - 67.0588*exp(1250*Cai)/(0.391606*exp(1250*Cai) + 1)^2 - 4.61503*exp(10000*Cai)/(0.00673795*exp(10000*Cai) + 1)^2;
dginfdCai    =  -6*Cai^5*0.00035^6/(Cai^6 + 0.00035^6)^2; if (Cai>=0.00035) dginfdCai    = -16*Cai^15*0.00035^16/(Cai^16 + 0.00035^16)^2; end;

dCaibufcdCai      =  2*Bufc*Kbufc*(Cai+Kbufc)/(Bufc*Kbufc+(Cai+Kbufc)^2)^2;
dCasrbufsrdCasri  =  2*Bufsr*Kbufsr*(Casr+Kbufsr)/(Bufsr*Kbufsr+(Casr+Kbufsr)^2)^2;

Stiffness    =  [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
                 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
                 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
                 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;...
                 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;...
                 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;...
                 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;...
                 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;...
                 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;...
                 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;...
                 0 0 0 0 0 0 0 0 0 0 1 0 0 0 dK1infdKi 0 0;...
                 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 -dt*dfCainfdCai/(taufca+dt) 0;...
                 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 -dt*dginfdCai/(taug+dt) 0;...
                 dt*Cm/VC/F*dINadm dt*Cm/VC/F*dINadh dt*Cm/VC/F*dINadj 0 0 0 0 0 0 0 0 0 0 1+dt*Cm/VC/F*(dINadNai+dIbNadNai+3*dINaKdNai+3*dINaCadNai) 0 dt*Cm/VC/F*(3*dINaCadCai) 0;...
                 0 0 0 dt*Cm/VC/F*dIKrdxr1 dt*Cm/VC/F*dIKrdxr2 dt*Cm/VC/F*dIKsdxs dt*Cm/VC/F*dItodr dt*Cm/VC/F*dItods 0 0 dt*Cm/VC/F*dIK1dK1 0 0 dt*Cm/VC/F*(dIKsdNai-2*dINaKdNai) 1+dt*Cm/VC/F*(dIK1dKi+dIKrdKi+dIKsdKi+dIpKdKi+dItodKi) 0 0;...
                 0 0 0 0 0 0 0 0 -dt*Caibufc*(-Cm/2/VC/F*(dICaLdd)+dIreldd) -dt*Caibufc*(-Cm/2/VC/F*dICaLdf) 0 -dt*Caibufc*(-Cm/2/VC/F*dICaLdfCa) -dt*Caibufc*(dIreldg) -dt*Caibufc*2*dINaCadNai*Cm/2/VC/F 0 1-dt*Caibufc*(-(dICaLdCai + dIbCadCai + dIpCadCai - 2*dINaCadCai)*Cm/2/VC/F + dIleakdCai - dIupdCai)-dt*dCaibufcdCai*(-(ICaL + IbCa + IpCa - 2*INaCa)*Cm/2/VC/F + Ileak - Iup + Irel) -dt*Caibufc*(dIleakdCasri+dIreldCasri);...
                 0 0 0 0 0 0 0 0 -dt*Casrbufsr*VC/VSR*(-dIreldd) 0 0 0 -dt*Casrbufsr*VC/VSR*(-dIreldg) 0 0 -dt*Casrbufsr*VC/VSR*(-dIleakdCai+dIupdCai) 1-dt*Casrbufsr*VC/VSR*(-dIleakdCasri-dIreldCasri)-dt*dCasrbufsrdCasri*VC/VSR*(-Ileak+Iup-Irel)];
                    
Residual     =  [m - (m1*taum + minf*dt)/(taum + dt);...
                 h - (h1*tauh + hinf*dt)/(tauh + dt);...
                 j - (j1*tauj + jinf*dt)/(tauj + dt);...
                 xr1 - (xr11*tauxr1 + xr1inf*dt)/(tauxr1 + dt);...
                 xr2 - (xr21*tauxr2 + xr2inf*dt)/(tauxr2 + dt);...
                 xs - (xs1*tauxs + xsinf*dt)/(tauxs + dt);...
                 r - (r1*taur + rinf*dt)/(taur + dt);...
                 s - (s1*taus + sinf*dt)/(taus + dt);...
                 d - (d1*taud + dinf*dt)/(taud + dt);...
                 f - (f1*tauf + finf*dt)/(tauf + dt);...
                 xK1inf - alK1/(alK1 + beK1);...
                 fCa - (fCa1*taufca + fcainf*kfca*dt)/(taufca + kfca*dt);...
                 g - (g1*taug + ginf*kg*dt)/(taug + kg*dt);...
                 Nai - Nai1 + dt*(INa + IbNa + 3*INaK + 3*INaCa)*Cm/VC/F;...
                 Ki - Ki1 + dt*(IK1 + Ito + IKr + IKs - 2*INaK + IpK + Istim - Iax)*Cm/VC/F;...
                 Cai - Cai1 - dt*Caibufc*( -(ICaL + IbCa + IpCa - 2*INaCa)*Cm/2/VC/F + Ileak - Iup + Irel);...
                 Casr - Casr1 - dt*Casrbufsr*VC/VSR*(-Ileak + Iup - Irel)];
             

Res         =  norm(Residual);

Var         =  Var - Stiffness\Residual; 

m = Var(1);    h = Var(2);  j  = Var(3);  xr1 = Var(4);     xr2 = Var(5);   xs = Var(6);
r = Var(7);    s = Var(8);  d  = Var(9);    f = Var(10); xK1inf = Var(11); fCa = Var(12);
g = Var(13); Nai = Var(14); Ki = Var(15); Cai = Var(16);   Casr = Var(17);


end

%--------------------------------------------------------------------------
% Calculate the source term and its derivative for the action potential
% evolution. Derivatives respect to electric potential are required
%-------------------------------------------------------------------------- 
source     =  INa + IbNa + INaK + INaCa + IK1   + IKr   +...
              IKs + IpK  + Ito  + ICaL  + IpCa  + IbCa  + Istim;

dENadphi   =  0;
dEKdphi    =  0;
dEKsdphi   =  0;
dECadphi   =  0;

dINadphi   =  GNa*m^3*h*j*(1-dENadphi);
dIbNadphi  =  GbNa*(1-dENadphi);
dINaKdphi  =  (-1)*PNaK*KO*Nai*F*(-0.0353*exp(-V*F/T/R)-0.01245*exp(-0.1*V*F/R/T))/R/T/(KO + KmK)/(Nai + KmNa)/(1 + 0.1245*exp(-0.1*V*F/R/T) + 0.0353*exp(-V*F/R/T))^2;
dINaCadphi =  kNaCa*F*exp((gamma-1)*V*F/R/T)*(Nai^3*CaO*(ksat*exp(gamma*V*F/R/T)+gamma*exp(V*F/R/T))-(gamma-1)*Cai*NaO^3*alpha)/(KmNai^3 + NaO^3)/(KmCa + CaO)/(1 + ksat*exp((gamma-1)*V*F/R/T))^2/R/T;
dIK1dphi   =  GK1*xK1inf*sqrt(KO/5.4)*(1 - dEKdphi);
dIKrdphi   =  GKr*sqrt(KO/5.4)*xr1*xr2*(1 - dEKdphi);
dIKsdphi   =  GKs*xs^2*(1 - dEKsdphi);
dIpKdphi   =  GpK*exp(0.167224*V)*(65.4052+exp(0.167224*V)-10.9373*EK+10.9373*V)/(65.4052+exp(0.167224*V))^2;
dItodphi   =  Gto*s*r*(1 - dEKdphi);
dICaLdphi  =  GCaL*f*d*fCa*4*F^2/R^2/T^2*(-exp(2*V*F/R/T)*(Cai*(2*F*V+R*T) + 0.341*CaO*(R*T-2*F*V))+R*T*Cai*exp(4*V*F/R/T)+0.341*R*T*CaO)/(exp(2*V*F/R/T) - 1)^2;
dIbCadphi  =  GbCa*(1-dECadphi);

dminfdphi   =  0.00040802028375750473*exp(0.221483942414175*V)/(exp(0.1107419712070875*V)+0.0018422115811651339)^3;
dtaumdphi   =  (exp(-12-V/5)/10*(1/(1+exp(7+V/5))+1/(1+exp((V-50)/200))) + 5*(1+exp(-12-V/5))*(-exp(7+V/5)/50/(1+exp(7+V/5))^2 - exp((50+V)/200)/2000/(exp(1/4)+exp(V/200))^2))/5/(1+exp(-12-V/5))^2;
dhinfdphi   =  -2*exp((V+71.55)/7.43)/7.43/(1+exp((V+71.55)/7.43))^3;
dtauhdphi   =  -0.00582171*exp(-0.0900901*V); if (V<=-40); dtauhdphi   =  (8.93905416279455e-9*exp(0.147059*V)-0.02925925925925926*exp(0.37311764705882355*V)-14819.615912208501*exp(0.6426176470588235*V))/(exp(0.22605882352941176*V)+114814.8148148148*exp(0.4955588235294117*V)+1.6412103442890793e-7)^2; end;
djinfdphi   =  -2*exp((V+71.55)/7.43)/7.43/(1+exp((V+71.55)/7.43))^3;
dtaujdphi   =  -0.01067*exp(-0.157*V) - 0.095*exp(-0.057*V); if (V<=-40); dtaujdphi   =  -(-0.00025505*exp(0.12728*V)/(0.00396087+exp(0.1378*V)) + 0.0000132*exp(0.12728*V)/(0.00396087+exp(0.1378*V))^2 + (-6.948e-6*exp(-0.04391*V) + exp(0.2444*V) - 25428)/(5.02627e10*exp(0.311*V) + 1) - 1.56317e10*exp(0.311*V)*(V+37.78)*(-6.948e-6*exp(-0.04391*V) + exp(0.2444*V) - 25428)/(5.02627e10*exp(0.311*V) + 1)^2 + (V+37.78)*(3.05087e-7*exp(-0.04391*V) + 0.2444*exp(0.2444*V))/(5.02627e10*exp(0.311*V) + 1))/(0.02424*exp(0.12728*V)/(0.00396087+exp(0.1378*V)) + (V+37.78)*(-6.948e-6*exp(-0.04391*V) + exp(0.2444*V) - 25428)/(5.02627e10*exp(0.311*V) + 1))^2; end;
ddinfdphi   =  exp((V+5)/7.5)/7.5/(1+exp((V+5)/7.5))^2;
dtauddphi   =  -0.07*exp((V+5)/5)/(1+exp((V+5)/5))^2 - 0.392*exp((V+5)/5)/(1+exp(-(35+V)/13))/(1+exp((V+5)/5))^2 + exp((50-V)/20)/20/(1+exp((50-V)/20))^2 + 0.150769*exp(-(35+V)/13)/(1+exp(-(35+V)/13))^2/(1+exp((V+5)/5));
dfinfdphi   =  -exp((V+20)/7)/7/(1+exp((V+20)/7))^2;
dtaufdphi   =  3/8*(44*exp((25+V)/10)/(exp(V/10)+exp(5/2))^2 -25*exp(-(V+27)^2/240)*(V+27));
drinfdphi   =  exp((20-V)/6)/6/(1+exp((20-V)/6))^2;
dtaurdphi   =  -19*exp(-(V+40)^2/1800)*(V+40)/1800;
dsinfdphi   =  -exp((V+20)/5)/5/(1+exp((V+20)/5))^2;
dtausdphi   =  -17/32*exp(-(V+45)^2/320)*(V+45)-exp((V-20)/5)/(1+exp((V-20)/5))^2;
dxsinfdphi  =  exp(-(V+5)/14)/14/(1+exp(-(V+5)/14))^2;
dtauxsdphi  =  275*exp(-(V+10)/6)/3/((1+exp(-(V+10)/6))^(3/2))/(1+exp((V-60)/20)) - 55*exp((V-60)/20)/sqrt(1+exp(-(V+10)/6))/(1+exp((V-60)/20))^2;
dxr1infdphi =  exp((V+26)/7)/7/(1+exp((V+26)/7))^2;
dtauxr1dphi =  270*exp((V+45)/10)*(exp((V+30)/11.5)*(3-20*exp((V+45)/10))+23)/23/(exp((V+30)/11.5)+1)^2/(1+exp((V+45)/10))^2;
dxr2infdphi =  -exp((V+88)/24)/24/(1+exp((V+88)/24))^2;
dtauxr2dphi =  -21*exp((V+180)/20)*(exp(V/10)-1)/125/(1+exp((V+60)/20))^2/(exp(V/20)+exp(3))^2;
dK1infdphi  =  -dK1infdKi/dEKdKi;

%--------------------------------------------------------------------------
% Compute the derivatives dCdphi and dWdphi by solving a system of equation                                                                                                                                                                                                                  
%-------------------------------------------------------------------------- 
dFdphi      = -[-((taum+dt)*dminfdphi     + (m1 - minf)*dtaumdphi)*dt/(taum + dt)^2;...
                -((tauh+dt)*dhinfdphi     + (h1 - hinf)*dtauhdphi)*dt/(tauh + dt)^2;...
                -((tauj+dt)*djinfdphi     + (j1 - jinf)*dtaujdphi)*dt/(tauj + dt)^2;...
                -((tauxr1+dt)*dxr1infdphi + (xr11 - xr1inf)*dtauxr1dphi)*dt/(tauxr1 + dt)^2;...
                -((tauxr2+dt)*dxr2infdphi + (xr21 - xr2inf)*dtauxr2dphi)*dt/(tauxr2 + dt)^2;...
                -((tauxs+dt)*dxsinfdphi   + (xs1 - xsinf)*dtauxsdphi)*dt/(tauxs + dt)^2;...
                -((taur+dt)*drinfdphi     + (r1 - rinf)*dtaurdphi)*dt/(taur + dt)^2;...
                -((taus+dt)*dsinfdphi     + (s1 - sinf)*dtausdphi)*dt/(taus + dt)^2;...
                -((taud+dt)*ddinfdphi     + (d1 - dinf)*dtauddphi)*dt/(taud + dt)^2;...
                -((tauf+dt)*dfinfdphi     + (f1 - finf)*dtaufdphi)*dt/(tauf + dt)^2;...
                -dK1infdphi;...
                0;...
                0];
            
dGdphi     = -[ dt*(dINadphi + dIbNadphi + 3*dINaKdphi + 3*dINaCadphi)*Cm/VC/F;...
                dt*(dIK1dphi + dItodphi + dIKrdphi + dIKsdphi - 2*dINaKdphi + dIpKdphi)*Cm/VC/F;...
               -dt*Caibufc*( -(dICaLdphi + dIbCadphi - 2*dINaCadphi)*Cm/2/VC/F);...
                0];

Neoresidual = [dFdphi;dGdphi];
            
solution   =  Stiffness\Neoresidual;
dwdphi     =  solution(1:13)';
dcdphi     =  solution(14:17)';

dInadC    =   [dINadNai;0;0;0];
dIbNadC   =   [dIbNadNai;0;0;0];
dINaKdC   =   [dINaKdNai;0;0;0];
dINaCadC  =   [dINaCadNai;0;dINaCadCai;0];
dIK1dC    =   [0;dIK1dKi;0;0];
dIKrdC    =   [0;dIKrdKi;0;0];
dIKsdC    =   [dIKsdNai;dIKsdKi;0;0];
dIpKdC    =   [0;dIpKdKi;0;0];
dItodC    =   [0;dItodKi;0;0];
dICaLdC   =   [0;0;dICaLdCai;0];
dIbCadC   =   [0;0;dIbCadCai;0];
dIpCadC   =   [0;0;dIpCadCai;0];
dINadW    =   [dINadm;dINadh;dINadj;0;0;0;0;0;0;0;0;0;0];
dIKrdW    =   [0;0;0;dIKrdxr1;dIKrdxr2;0;0;0;0;0;0;0;0];
dIKsdW    =   [0;0;0;0;0;dIKsdxs;0;0;0;0;0;0;0];
dItodW    =   [0;0;0;0;0;0;dItodr;dItods;0;0;0;0;0];
dICaLdW   =   [0;0;0;0;0;0;0;0;dICaLdd;dICaLdf;0;dICaLdfCa;0];
dIK1dW    =   [0;0;0;0;0;0;0;0;0;0;dIK1dK1;0;0];
dsource   =   dINadphi + dIbNadphi + dINaKdphi + dINaCadphi + dIK1dphi + dIKrdphi + dIKsdphi + dIpKdphi + dItodphi + dICaLdphi + dIbCadphi + ...
              dcdphi*dInadC + dcdphi*dIbNadC + dcdphi*dINaKdC + dcdphi*dINaCadC + dcdphi*dIK1dC + dcdphi*dIKrdC + dcdphi*dIKsdC + dcdphi*dIpKdC + ...
              dcdphi*dItodC + dcdphi*dICaLdC + dcdphi*dIbCadC + dcdphi*dIpCadC + dwdphi*dINadW + dwdphi*dIKrdW + dwdphi*dIKsdW + dwdphi*dItodW + ...
              dwdphi*dICaLdW + dwdphi*dIK1dW;
          
%--------------------------------------------------------------------------
% Calculate voltage                                                                                                                                                                                                                
%-------------------------------------------------------------------------- 
Res        =  (V1 - V)/dt + source;
V1         =  V1 - Res/(1/dt + dsource);
iter       =  iter + 1;

%--------------------------------------------------------------------------
%  Save in memory
%--------------------------------------------------------------------------
Volt(i)          =  V1;

gate_m(i)        =  m;
gate_g(i)        =  g;
gate_d(i)        =  d;
gate_h(i)        =  h;
gate_f(i)        =  f;
gate_s(i)        =  s;
gate_r(i)        =  r;
gate_j(i)        =  j;
gate_fCa(i)      =  fCa;
gate_xs(i)       =  xs;
gate_xr1(i)      =  xr1;
gate_xr2(i)      =  xr2;
gate_xK1inf(i)   =  xK1inf; 

ion_Ca(i)        =  Cai;
ion_CaSR(i)      =  Casr;
ion_Na(i)        =  Nai;
ion_K(i)         =  Ki; 

inten_Ina(i)     =  INa;   
inten_ICaL(i)    =  ICaL;   
inten_Ito(i)     =  Ito;  
inten_IKs(i)     =  IKs;  
inten_IKr(i)     =  IKr;  
inten_IK1(i)     =  IK1;  
inten_INaCa(i)   =  INaCa;
inten_INaK(i)    =  INaK; 
inten_IpCa(i)    =  IpCa; 
inten_IpK(i)     =  IpK;  
inten_IbNa(i)    =  IbNa; 
inten_IbCa(i)    =  IbCa; 
inten_Ileak(i)   =  Ileak;
inten_Iup(i)     =  Iup; 
inten_Irel(i)    =  Irel; 

%----------------------------------------------------------------------
% Rename variables for the next time-step
%----------------------------------------------------------------------
m1               =  m;
g1               =  g;
d1               =  d;
h1               =  h;
f1               =  f;
s1               =  s;
r1               =  r;
j1               =  j;
fCa1             =  fCa;
xs1              =  xs;
xr11             =  xr1;
xr21             =  xr2; 
xK1inf1          =  xK1inf;
Cai1             =  Cai; 
Casr1            =  Casr; 
Nai1             =  Nai;
Ki1              =  Ki;
V                =  V1;

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
subplot(5,3,1)
plot(Time,gate_m,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([0,1])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,1);
xlabel('Time (ms)','FontSize',12) 
ylabel('Gate variable $$m$$','FontSize',12) 

subplot(5,3,2)
plot(Time,gate_g,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([0,1])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,2);
xlabel('Time (ms)','FontSize',12) 
ylabel('Gate variable $$g$$','FontSize',12) 

subplot(5,3,3)
plot(Time,gate_d,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([0,1])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,3);
xlabel('Time (ms)','FontSize',12) 
ylabel('Gate variable $$d$$','FontSize',12) 

subplot(5,3,4)
plot(Time,gate_h,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([0,0.8])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,4);
xlabel('Time (ms)','FontSize',12) 
ylabel('Gate variable $$h$$','FontSize',12) 

subplot(5,3,5)
plot(Time,gate_f,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([0.2,1])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,5); 
xlabel('Time (ms)','FontSize',12) 
ylabel('Gate variable $$f$$','FontSize',12) 

subplot(5,3,6)
plot(Time,gate_s,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([0,1])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,6); 
xlabel('Time (ms)','FontSize',12) 
ylabel('Gate variable $$s$$','FontSize',12) 

subplot(5,3,7)
plot(Time,gate_r,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([0,0.6])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,7); 
xlabel('Time (ms)','FontSize',12) 
ylabel('Gate variable $$r$$','FontSize',12) 

subplot(5,3,8)
plot(Time,gate_j,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([0,0.8])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,8); 
xlabel('Time (ms)','FontSize',12) 
ylabel('Gate variable $$j$$','FontSize',12) 

subplot(5,3,9)
plot(Time,gate_fCa,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([0.2,1])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,9);
xlabel('Time (ms)','FontSize',12) 
ylabel('Gate variable $$f_{Ca}$$','FontSize',12) 

subplot(5,3,10)
plot(Time,gate_xs,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])  
ylim([0,0.16])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,10);
xlabel('Time (ms)','FontSize',12) 
ylabel('Gate variable $$x_{s}$$','FontSize',12) 

subplot(5,3,11)
plot(Time,gate_xr1,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([0,1])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,11);
xlabel('Time (ms)','FontSize',12)
ylabel('Gate variable $$x_{r1}$$','FontSize',12) 

subplot(5,3,12)
plot(Time,gate_xr2,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([0,0.5])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,12);
xlabel('Time (ms)','FontSize',12) 
ylabel('Gate variable $$x_{r2}$$','FontSize',12) 

subplot(5,3,13)
plot(Time,gate_xK1inf,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([0,0.2])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,13);
xlabel('Time (ms)','FontSize',12) 
ylabel('Gate variable $$K_1^{\infty}$$','FontSize',12)


figure(2)
subplot(2,2,1)
plot(Time,ion_Ca,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([0,0.001])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(2,2,1);
xlabel('Time (ms)','FontSize',12) 
ylabel('Concentration $$Ca$$ (mM)','FontSize',12) 

subplot(2,2,2)
plot(Time,ion_CaSR,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([0.4,0.6])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(2,2,2);
xlabel('Time (ms)','FontSize',12) 
ylabel('Concentration $$Ca^{sr}$$ (mM)','FontSize',12) 

subplot(2,2,3)
plot(Time,ion_Na,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])  
ylim([11.58,11.61])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(2,2,3);
xlabel('Time (ms)','FontSize',12) 
ylabel('Concentration $$Na$$ (mM)','FontSize',12)

subplot(2,2,4)
plot(Time,ion_K,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])  
ylim([138.288,138.306])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(2,2,4);
xlabel('Time (ms)','FontSize',12) 
ylabel('Concentration $$K$$ (mM)','FontSize',12) 


figure(3)
subplot(5,3,1)
plot(Time,inten_Ina,'-k','LineWidth',2) 
set(gcf,'color','w');
xlim([0,2])
ylim([-300,0])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,1);
xlabel('Time (ms)','FontSize',12)
ylabel('Intensity $$I_{Na}$$ (pA/pF)','FontSize',12)

subplot(5,3,2)
plot(Time,inten_ICaL,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([-7,0])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,2); 
xlabel('Time (ms)','FontSize',12)
ylabel('Intensity $$I_{CaL}$$ (pA/pF)','FontSize',12) 

subplot(5,3,3)
plot(Time,inten_Ito,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([0,12])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,3); 
xlabel('Time (ms)','FontSize',12) 
ylabel('Intensity $$I_{to}$$ (pA/pF)','FontSize',12)

subplot(5,3,4)
plot(Time,inten_IKs,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000]) 
ylim([0,0.35])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,4);
xlabel('Time (ms)','FontSize',12) 
ylabel('Intensity $$I_{Ks}$$ (pA/pF)','FontSize',12) 

subplot(5,3,5)
plot(Time,inten_IKr,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])  
ylim([0,0.6])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,5); 
xlabel('Time (ms)','FontSize',12) 
ylabel('Intensity $$I_{Kr}$$ (pA/pF)','FontSize',12)

subplot(5,3,6)
plot(Time,inten_IK1,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([0,2])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,6); 
xlabel('Time (ms)','FontSize',12) 
ylabel('Intensity $$I_{K1}$$ (pA/pF)','FontSize',12) 

subplot(5,3,7)
plot(Time,inten_INaCa,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([-0.4,0.4])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,7); 
xlabel('Time (ms)','FontSize',12) 
ylabel('Intensity $$I_{NaCa}$$ (pA/pF)','FontSize',12) 

subplot(5,3,8)
plot(Time,inten_INaK,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([0.12,0.24])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,8);
xlabel('Time (ms)','FontSize',12)
ylabel('Intensity $$I_{NaK}$$ (pA/pF)','FontSize',12) 

subplot(5,3,9)
plot(Time,inten_IpCa,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([0.1,0.55])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,9);
xlabel('Time (ms)','FontSize',12) 
ylabel('Intensity $$I_{pCa}$$ (pA/pF)','FontSize',12) 

subplot(5,3,10)
plot(Time,inten_IpK,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([0,1.4])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,10); 
xlabel('Time (ms)','FontSize',12)
ylabel('Intensity $$I_{pK}$$ (pA/pF)','FontSize',12) 

subplot(5,3,11)
plot(Time,inten_IbNa,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([-0.045,-0.005])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,11);
xlabel('Time (ms)','FontSize',12)
ylabel('Intensity $$I_{bNa}$$ (pA/pF)','FontSize',12) 

subplot(5,3,12)
plot(Time,inten_IbCa,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([-0.14,-0.04])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,12);
xlabel('Time (ms)','FontSize',12)
ylabel('Intensity $$I_{bCa}$$ (pA/pF)','FontSize',12) 

subplot(5,3,13)
plot(Time,inten_Ileak,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([0.000032,0.000048])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,13);
xlabel('Time (ms)','FontSize',12) 
ylabel('Intensity $$I_{leak}$$ (mM/ms)','FontSize',12) 

subplot(5,3,14)
plot(Time,inten_Iup,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([0,0.0004])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,14);
xlabel('Time (ms)','FontSize',12) 
ylabel('Intensity $$I_{up}$$ (mM/ms)','FontSize',12) 

subplot(5,3,15)
plot(Time,inten_Irel,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,100])
ylim([0,0.025])
grid on
set(gca,'GridLineStyle','--')
hFig  =  subplot(5,3,15); 
xlabel('Time (ms)','FontSize',12) 
ylabel('Intensity $$I_{rel}$$ (mM/ms)','FontSize',12) 

figure(4)
plot(Time,Volt,'-k','LineWidth',2)
set(gcf,'color','w');
xlim([0,1000])
ylim([-90,60])
grid on
set(gca,'GridLineStyle','--')
hFig  =  figure(4);
xlabel('Time (ms)','FontSize',12) 
ylabel('Voltage $$\phi$$ (mV)','FontSize',12) 





% ENa     = R*T/F*log(NaO/Nai1);
% EK      = R*T/F*log(KO/Ki1);
% ECa     = R*T/2/F*log(CaO/Cai1);
% EKs     = R*T/F*log((KO + pKNa*NaO)/(Ki1 + pKNa*Nai1));
% alK1    = 0.1/(1 + exp(0.06*(V-EK-200)));
% beK1    = (3*exp(0.0002*(V-EK+100)) + exp(0.1*(V-EK-10)))/(1 + exp(-0.5*(V-EK)));
% xK1inf  = alK1/(alK1 + beK1);
% 
% 
% INa     = GNa*m^3*h*j*(V-ENa);
% ICaL    = GCaL*d*f*fCa*4*V*F^2/R/T*(Cai1*exp(2*V*F/R/T) - 0.341*CaO)/(exp(2*V*F/R/T) - 1);
% Ito     = Gto*r*s*(V - EK);
% IKs     = GKs*xs^2*(V - EKs);
% IKr     = GKr*sqrt(KO/5.4)*xr1*xr2*(V - EK);
% IK1     = GK1*sqrt(KO/5.4)*xK1inf*(V - EK);
% INaCa   = kNaCa*(exp(gamma*V*F/R/T)*Nai1^3*CaO -exp((gamma-1)*V*F/R/T)*NaO^3*Cai1*alpha)/(KmNai^3 + NaO^3)/(KmCa + CaO)/(1 + ksat*exp((gamma-1)*V*F/R/T));
% INaK    = PNaK*KO*Nai1/(KO + KmK)/(Nai1 + KmNa)/(1 + 0.1245*exp(-0.1*V*F/R/T) + 0.0353*exp(-V*F/R/T));
% IpCa    = GpCa*Cai1/(KpCa + Cai1);
% IpK     = GpK*(V - EK)/(1 + exp((25-V)/5.98));
% IbNa    = GbNa*(V - ENa);
% IbCa    = GbCa*(V - ECa);
% Ileak   = Vleak*(Casr1 - Cai1);
% Iup     = Vmaxup/(1 + Kup^2/Cai1^2);
% Irel    = (arel*Casr1^2/(brel^2 + Casr1^2) + crel)*d*g;
% 
% StimAm  = 72; 
% Istim   = -StimAm;
%     
% source     =  INa + IbNa + INaK + INaCa + IK1 + IKr + IKs + IpK + Ito + ICaL  + IpCa  + IbCa  + Istim;
% 
% Caibufc      = 1/(1 + (Bufc*Kbufc)/(Cai1+Kbufc)^2);
% Casrbufsr    = 1/(1 + (Bufsr*Kbufsr)/(Casr1+Kbufsr)^2);
% 
% 
% Nai1    = Nai1 - dt*(INa + IbNa + 3*INaK + 3*INaCa)*Cm/VC/F;
% Ki1     = Ki1 - dt*(IK1 + Ito + IKr + IKs - 2*INaK + IpK + Istim - Iax)*Cm/VC/F;
% Cai1    = Cai1 + dt*Caibufc*( -(ICaL + IbCa + IpCa - 2*INaCa)*Cm/2/VC/F + Ileak - Iup + Irel);
% Casr1   = Casr1 + dt*Casrbufsr*VC/VSR*(-Ileak + Iup - Irel);
% 
% V         = V - dt*source;

% Cai_step   = Cai  + dt*Caibufc.*( -(ICaL + IbCa + IpCa - 2*INaCa).*(Cm/2/VC/F) + Ileak - Iup + Irel);
% Casr_step  = Casr + dt*Casrbufsr*VC/VSR.*(-Ileak + Iup - Irel);
% Nai_step   = Nai  - dt*(INa + IbNa + 3*INaK + 3*INaCa)*Cm/VC/F;
% Ki_step    = Ki   - dt*(IK1 + Ito + IKr + IKs - 2*INaK + IpK + Istim - Iax)*Cm/VC/F;
