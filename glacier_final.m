% This code generates the growth and evolution of an eroding glacier, and
% plots the volume and ELA as a function of time

% ELA oscillates

% written by AGT 2/18/2016

clear all
figure(1)
clf
figure(2)
clf
figure(3)
clf

%% initialize

%constants
rhoi = 917; %this is in kg/m3
g = 9.81; % gravity m/s2
usl = 0.05; %sliding speed, m/yr 
a = 2.1e-16; %flow law parameter Pa-3 yr-1
gamma = 0.01; %gradient in net accumulation rate m/yr per m altitude (aka 1/yr)
E = 0.001; %erosion rate m/yr
amp = 100; %amplitude of oscillation
p = 2000; %period of oscillation, years

% create distance array
xmax = 20000; %m 
dx = 100; %m
x = (dx/2):dx:xmax-(dx/2); %so that the x value is in the middle of each 'box'

% set up time array
tmax = 2000; %years
dt = 0.005; 
t = 0:dt:tmax;

% define bedrock profile
zbmax = 1500;
zbmin = 0;
xstar=12000;
zb = zbmin+((zbmax-zbmin)*exp(-x/xstar)); %equation for bedrock profile
width = 1000; %constant width of glacier (for calculating ice volume)

%set up initial glacier parameters
h = zeros(size(x)); %initial thickness of glacier
z = zb+h; %height of glacier on top of bedrock
Vi = zeros(size(t)); %empty array for volume
ela = zeros(size(t)); %equilibrium line altitude m

imax = length(t);
nplots = 50;
tplot = tmax/nplots;

%% run

for i = 1:imax
    
    ela(i) = (zbmax-300) + amp*sin((2*pi*t(i))/p); %oscillating ela
    
    %calculate growth of glacier
    slope = diff(z)/dx; %slope of glacier, n-1 elements
    hedge = h(1:end-1) + (diff(h)/2); %this yields the new edge h's
    q = (usl.*hedge)+((a.*((rhoi*g.*abs(slope)).^3)).*((hedge.^5)/5)); %flux of ice, n-1 elements
    q = [0 q 0]; %pads q on either side, n+1 elements
    dqdx = diff(q)/dx; %change in flux, n elements
    b = gamma.*(z-ela(i)); %dynamic local mass balance, n elements
    dhdt = b-dqdx; %change in height, n elements
    h = h+(dhdt*dt); %new height of glacier, n elements
    h = max(0,h); %keeps h from becoming negative
    
    %erode bedrock
    xglacier = find(h>0); %find x values where the glacier exists
    zb(xglacier) = zb(xglacier)-E; %erode the bedrock under the glacier
    
    %calculate new height and volume of glacier
    z = zb+h; %new height of glacier on top of bedrock, n elements
    vi(i) = width*mean(h)*(dx*length(xglacier)); %volume of ice

    if(rem(t(i),tplot)==0)
        figure(1)
        plot(x/1000,zb,'k','linewidth',1.5)
        hold on
        plot(x/1000,z,'b','linewidth',1.5)
        axis([0 xmax/1000 min(zb) zbmax+500])
        xlabel('Distance (km)','fontname','arial','fontsize', 18)
        ylabel('Height (m)', 'fontname', 'arial', 'fontsize', 18)
        set(gca, 'fontsize', 14, 'fontname', 'arial') 
        time=num2str(t(i)); %convert time of each plot to 'letters'
        timetext=strcat(time,' years'); %add years to the time
        text(2,400,timetext,'fontsize',14) %shows time on each plot
        hold off
        pause(0.1)
    end

end

%% finalize

%equilibrium time scale 
% teq = 600; %time it takes to get to equilibrium, need to find a way to calculate this
% charts = teq/(1-(1/exp(1))); %characteristic time scale
% eqts = charts*3; %equilibrium time scale

% plot volume as a function of time
vimax = max(vi);
figure(2)
plot(t,vi/1e9)
axis([0 tmax 0 ((vimax/1e9)+0.1)])
xlabel('Time (years)','fontname','arial','fontsize',18)
ylabel('Volume (km^3)','fontname','arial','fontsize',18)
set(gca,'fontsize',14,'fontname','arial')

% plot ela as a function of time
figure(3)
plot(t,ela)
xlabel('Time (years)','fontname','arial','fontsize',18)
ylabel('ELA (km)','fontname','arial','fontsize',18)
set(gca,'fontsize',14,'fontname','arial')



