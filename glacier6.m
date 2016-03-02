% This code generates the growth and evolution of an eroding glacier, and
% plots the volume and ELA as a function of time

% ELA is determined by d18O values for past 2000k years

% written by AGT 2/25/2016

clear all
figure(1)
clf
figure(2)
clf

%% initialize

%constants
rhoi = 917; %this is in kg/m3
g = 9.81; % gravity m/s2
usl = 0.05; %sliding speed, m/yr 
a = 2.1e-16; %flow law parameter Pa-3 yr-1
gamma = 0.01; %gradient in net accumulation rate m/yr per m altitude (aka 1/yr)
E = 0.001; %erosion rate m/yr

% create distance array
xmax = 20000; %m 
dx = 100; %m
x = (dx/2):dx:xmax-(dx/2); %so that the x value is in the middle of each 'box'
N=length(x);

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
vi = zeros(size(t)); %empty array for volume
ela = zeros(size(t)); %equilibrium line altitude m
q = zeros(1,(N+1));

%load data file
load pleist_del18O_2000yrs.txt
age_data = pleist_del18O_2000yrs(:,1); %age column
d18O_data = pleist_del18O_2000yrs(:,2); %d18O column

%convert d18O to ela
age0 = age_data; %age in kyrs
age=flipud(age0); %flip ages so it goes from past to present
rangeO = range(d18O_data); %range of isotope data 
scaled = (d18O_data./rangeO) - min(d18O_data./rangeO); %scales d18O data from 0-1
ela0 = (zbmax-500)+(scaled.*(200)); %scales data to 200m of ela change
ela = interp1(age,ela0,t); %interpolate ela for each age being calculated

imax = length(t);
nplots = 200;
tplot = tmax/nplots;
nframe = 0;

%% run

for i = 1:imax
    
    elal = ela(i)*ones(size(x)); %line of ela for plotting
    
    %calculate growth of glacier
    slope = diff(z)/dx; %slope of glacier, n-1 elements
    hedge = h(1:end-1) + (diff(h)/2); %this yields the new edge h's
    q(2:N) = (usl.*hedge)+((a.*((rhoi*g.*abs(slope)).^3)).*((hedge.^5)/5)); %flux of ice, n-1 elements
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
        nframe = nframe+1;
        figure(1)
        subplot(2,1,1)%subplot (nrows,ncols,plot_number)
        plot(x/1000,z,'b','linewidth',1.5)
        hold on
        plot(x/1000,zb,'k','linewidth',1.5)
        plot(x/1000,elal,'r','linewidth',1)
        title('Glacier growth')
        axis([0 xmax/1000 min(zb) zbmax+500])
        xlabel('Distance (km)','fontname','arial','fontsize', 18)
        ylabel('Height (m)', 'fontname', 'arial', 'fontsize', 18)
        set(gca, 'fontsize', 14, 'fontname', 'arial') 
        time=num2str(t(i)); %convert time of each plot to 'letters'
        timetext=strcat(time,' kyrs'); %add years to the time
        text(2,400,timetext,'fontsize',14) %shows time on each plot
        legend('Glacier height','Bedrock height','ELA')
        hold off
        
        subplot(2,1,2)
        pause(0.1)
        plot(x/1000,h,'b','linewidth',1.5)
        axis([0 xmax/1000 0 350])
        title('Ice thickness')
        xlabel('Distance (km)','fontname','arial','fontsize', 18)
        ylabel('Height (m)', 'fontname', 'arial', 'fontsize', 18)
        set(gca, 'fontsize', 14, 'fontname', 'arial') 
        time=num2str(t(i)); %convert time of each plot to 'letters'
        timetext=strcat(time,' kyrs'); %add years to the time
        text(15,300,timetext,'fontsize',14) %shows time on each plot
        M(:,nframe) = getframe(gcf);
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
subplot(2,1,1)
plot(t,vi/1e9)
axis([0 tmax 0 ((vimax/1e9)+0.1)])
title('Volume')
xlabel('Time (kyrs)','fontname','arial','fontsize',18)
ylabel('Volume (km^3)','fontname','arial','fontsize',18)
set(gca,'fontsize',14,'fontname','arial')

% plot ela as a function of time
subplot(2,1,2)
plot(t,ela)
title('Equilibrium line altitude')
xlabel('Time (kyrs)','fontname','arial','fontsize',18)
ylabel('ELA (km)','fontname','arial','fontsize',18)
set(gca,'fontsize',14,'fontname','arial')


movie2avi(M,'Glacier evolution','fps',10) %


