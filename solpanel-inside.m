clear all
clc

%Modelling parameters
timestep=120;                      %seconds
Nd=365;                    %Numbers of days to run
N=1;                       %Numbers of years to run
endtime=N*Nd*24*3600;     %Numbers of seconds to run
steps=endtime/timestep;

%Fundamental constants
h=6.63E-34;   %Planc
c=3.00E8;     %Light speed
k=1.38E-23;   %Boltzmann
sigma=5.67E-8; % W/(m^2*K^4) Stefan-Boltzmann constant

%Orbital paramters
wj=7.272205217E-5;		%Angular frequency earth spin	(rad/s)
wjs=7.291666667E-5; 	%Sidereal angular frequency earth spin	(rad/s)
ws=1.991021278E-7;		%Angular frequency earth orbit	(rad/s)
axistilt=-0.40909;		%rad

%Atmospheric paramteters
IoSun=1361;			    %W/m^2
abscoeff=0.357;			%no unit
sunchance=0.3;      %Average percentage of sunshine
%Tair=5+273.15;      %K Air temperture

%Geografical parameters
la=55.607443*pi/180;		%Latitude (rad)
lo=12.362914*pi/180;		%Longitude (rad)
%la=-5*pi/180;		%Latitude (rad)
%lo=3*pi/180;		%Longitude (rad)
timezone=15*pi/180; 		%The longitude where the Sun would be exact south at 12. (rad)

%Solar panel parameters
etacell=0.136;			  %no unit
area=2.85;			      %4.65 % m^2
area1=1.80;           %
alpha=90*pi/180;			%rad (inclination of panel (angle between horizontal and solar panel))
beta=-46.6*pi/180;		%rad (angle between sunpanel and south (negative = east))
alpha1=90*pi/180;
beta1=-43.4*pi/180;
lacel=la-alpha; 			%rad (Sunpanel vertical angle)
locel=lo-beta; 			  %rad (Sunpanel horizontal angle)
lacel1=la-alpha1; 			%rad (Sunpanel vertical angle)
locel1=lo-beta1; 			  %rad (Sunpanel horizontal angle)
absbp=0.95;           % No unit. Absorption coefficient of black paint
emibp=0.90;           % No unit. Emissivity of black paint
areaemit=2*(area+area1);      % m^2
conductivityPE=0.42;  % W/(m*K) Thermal conductivity of PE
thicknessPE=0.005;    % m Thickness of panel lamels

%Greenhouse parameters
Agreenhouse=33.81;    % m^2. Glass surface, not floor
Uglass=6;             % W/(m^2*K)
Pelectricradiator=1700; % W
transmissivitygreenhousevis=0.80; % Estimated value

%Water tank parameters
Vw=1;                 % m^3 Volume of water tank
Aw=6;                 % m^2 Surface area of water tank
dinsulate=0.01;        % m Thickness of insulation around water tank
mw=1000*Vw;           % kg Mass of water in water tank
cw=4180;              % J/(kg*K) Specific heat capacity of water
conductivityrockwool=0.033;      % W/(m*K) Thermal conductivity of rockwool



%Load temperature data
timedata=(textread('finaltime.dat', '%f',"headerlines",5))';
outtempdata=(textread('finalouttemp.dat', '%f',"headerlines",5))';
greentempdata=(textread('finalgreentemp.dat', '%f',"headerlines",5))';
arduinotimedata=(textread('arduinotime.dat', '%f',"headerlines",5))';
arduinogreentempdata=(textread('arduinogreentemp.dat', '%f',"headerlines",5))';
arduinoradiatorstatus=(textread('arduinoradiatorstatus.dat', '%f',"headerlines",5))';


%%%%%%%%%%

% Preassign variables:
%Tin=zeros(1,steps+1);
%tite=zeros(1,steps+1);
%Tair=zeros(1,steps+1);
%inclination=zeros(1,steps+1);
%orientation=zeros(1,steps+1);
%ISun=zeros(1,steps+1);
%Pcell=zeros(1,steps+1);
%Pemit=zeros(1,steps+1);
%Pconduct=zeros(1,steps+1);
%Twater=zeros(1,steps+1);
%Plosstank=zeros(1,steps+1);
%Qtot=zeros(1,steps+1);
%time=zeros(1,steps+1);
%dailytvar=zeros(1,steps+1);

% Initial parameters
Twaterstart=7+273.15; % K
Tinstart=0+273.15;    % K
t=0;
orientationsum=0;
Twater(1)=Twaterstart;
Tin(1)=Tinstart;
Qradiatortotal(1)=0;
totalsunhours=0;
Qradiatortotalmeasured(1)=0;
totalpumptime=0;
Qradiatorless(1)=0;

% Time loop
for j=1:steps
  time(j)=t;
  i=axistilt*cos(ws*(t+9*24*60*60));
  fi=-wj*t+timezone+pi;
  Tin(j)=greentempdata(j)+273.15;    % K Temperature of greenhouse   %273.15+15.0+15*sin((2*pi/(365*24*3600))*t-pi/2)+5*sin(fi+pi/2);
  Tair(j)=outtempdata(j)+273.15;     % K Outside air temperature     %273.15+10.0+10*sin((2*pi/(365*24*3600))*t-pi/2)+5*sin(fi+pi/2);
  %dailytvar(j)=5*sin(fi+pi/2);
  inclination(j)=90-(180/pi)*acos(cos(la)*cos(i)*(cos(lo)*cos(fi)+sin(lo)*sin(fi))+sin(la)*sin(i));               %Elevetion of Sun over the horizon
  orientation(j)=90-(180/pi)*acos(cos(lacel)*cos(i)*(cos(locel)*cos(fi)+sin(locel)*sin(fi))+sin(lacel)*sin(i));   %Angle between Sun and solar panels
  orientation1(j)=90-(180/pi)*acos(cos(lacel1)*cos(i)*(cos(locel1)*cos(fi)+sin(locel1)*sin(fi))+sin(lacel1)*sin(i));   %Angle between Sun and solar panels
  
  if Tair(j) > 2+273.15
    if Tin(j) - Tair(j) > 5
      sunshine(j)=1;
      totalsunhours=totalsunhours+timestep/3600;  % h Total hours of sunshine
    else
      sunshine(j)=0;
    endif
  else
    sunshine(j)=0;
  endif
  
  if inclination(j) < 0
    orientation(j)=0.0;
    airmass=0.0;
    ISun(j)=0.0;
  else
    airmass=1/(sin((inclination(j)+244.0/(165.0+47.0*(inclination(j)**1.1)))*pi/180));
    ISun(j)=IoSun*exp(-abscoeff*airmass);
  endif

  if orientation(j) > 0
    Pcell(j)=transmissivitygreenhousevis*ISun(j)*sunshine(j)*absbp*(area*sin(orientation(j)*pi/180)+area1*sin(orientation1(j)*pi/180));
  else
    orientation(j)=0.0;
    Pcell(j)=0;
    Pemit(j)=0.0;
    Pconduct(j)=0.0;
  endif

  Pemit(j)=areaemit*emibp*sigma*(Twater(j)^4-Tin(j)^4);                % W Emitted radiation of panel  
  Pconduct(j)=conductivityPE*areaemit*(Twater(j)-Tin(j))/thicknessPE;  % W Emitted radiation of panel
  Plosstank(j)=(conductivityrockwool*Aw/dinsulate)*(Twater(j)-Tin(j));  % W Heat loss from water tank
  
  if Tair(j) < 2+273.15
    Pradiator(j)=(Uglass*Agreenhouse)*(Tin(j)-Tair(j));
    Pradiatoreffective(j)=conductivityPE*areaemit*(Twater(j)-Tin(j))/thicknessPE;  % W Emitted radiation of radiator
    totalpumptime=totalpumptime + timestep;
  elseif Pcell(j) > (Pemit(j)+Pconduct(j))
    Pradiator(j)=0;
    Pradiatoreffective(j)=0;
    totalpumptime=totalpumptime + timestep;
  else
    Pradiatoreffective(j)=0;
    Pradiator(j)=0;
  endif
  
  if Pcell(j) > (Pemit(j)+Pconduct(j))
    Qtot(j)=(Pcell(j)-Pemit(j)-Pconduct(j)-Plosstank(j)-Pradiator(j))*timestep;
  else
    Qtot(j)=(-Plosstank(j)-Pradiator(j))*timestep;
  endif

  if isnan(greentempdata(j)) == 1
    if isnan(outtempdata(j)) == 1
      Qtot(j)=0;
    endif
  endif

  if isnan(arduinoradiatorstatus(j)) == 0
    Qradiatortotalmeasured(j+1)=Qradiatortotalmeasured(j) + arduinoradiatorstatus(j)*Pelectricradiator*timestep/(3.6E6);
  else
    Qradiatortotalmeasured(j+1)=Qradiatortotalmeasured(j);
  endif
  
  if Twater(j) < 2 + 273.15
    Qradiatorless(j+1) = Qradiatorless(j) + Pradiator(j)*timestep/(3.6E6);
  else
    Qradiatorless(j+1) = Qradiatorless(j);
  endif
  
  Twater(j+1)=Twater(j)+Qtot(j)/(mw*cw);
  
  
  Qradiatortotal(j+1)=Qradiatortotal(j)+Pradiator(j)*timestep/(3.6E6); % kWh Total amount of heating energy required   
  t=t+timestep;
endfor

%    Pemit(j)=areaemit*emibp*sigma*(Tw^4-Tair^4); % W Emitted radiation of panel  
%    Eabs=(Pcell(j)-Pemit(j))*timestep;

%  orientationsum=orientationsum+orientation(j);
%  Tmaxtheory(j)=(Pcell(j)/(areaemit*sigma*emibp*sunchance))^(1/4);  


Tin(steps+1)=Tin(steps);
%time(steps+1)=time(steps);
Tair(steps+1)=Tair(steps);
inclination(steps+1)=inclination(steps);
orientation(steps+1)=orientation(steps);
ISun(steps+1)=ISun(steps);
Pcell(steps+1)=Pcell(steps);
Pemit(steps+1)=Pemit(steps);
Pconduct(steps+1)=Pconduct(steps);
Twater(steps+1)=Twater(steps);
Plosstank(steps+1)=Plosstank(steps);
Qtot(steps+1)=Qtot(steps);
time(steps+1)=time(steps);
sunshine(steps+1)=sunshine(steps);
Pradiator(steps+1)=Pradiator(steps);
Qradiatorless(steps+1)=Qradiatorless(steps);
Pradiatoreffective(steps+1)=Pradiatoreffective(steps);
%dailytvar(steps+1)=dailytvar(steps);
%Tmaxtheory(steps+1)=Tmaxtheory(steps);

%plotyy(t,Pcell,t,inclination)
plot(time,Tair)
%axis ([700 25000 0 1E7])
grid on

%Qradiatortotal(j)
totalpumptime/(3600*24)
Qradiatorless(j)