%% AEM 4602W Fluids lab
% Due Nov 22
% Daniel Kim
clear
close all
clc
%% Constants and stuff

% NACA 0016
% max camber 0% of chord
% location of max camber: 0% of chord from LE
% thickness: 16% of chord, symm => @ 30% chord from LE

% Load Data Files
Data = load('Group6Data.mat');

% Uncertainties
Unc.Force = 0.1;            % +/- 0.1 N
Unc.Speed = 0.4;            % +/- 0.4 m/s
Unc.Density = 0.02;         % +/- 2%
Unc.Viscosity = 0.01;       % +/- 1%
Unc.Dimensions = 0.001;     % +/- 1 mm
Unc.AoA = 0.2;              % +/- 0.2 degrees



%% Forces
forcesData = Data.stingData;
AoA = zeros(1,length(forcesData));
forceComp = zeros(3,length(forcesData));
windSpeed = zeros(1,length(forcesData));
for i = 1:length(forcesData)
    AoA(1,i) = forcesData(i).aoa;
    forceComp(:,i) = forcesData(i).forces;
    windSpeed(:,i) = forcesData(i).windSpeed;
end
forceNormal = forceComp(1,:);
forceTransverse = forceComp(2,:);
forceAxial = forceComp(3,:);

results = Forces(forceNormal,forceAxial,AoA,windSpeed);

cLAlpha = figure;
hold on
grid on
axis([-10 25 -0.5 1.5])
plot(AoA,results.cL,'-o')
title('c_L vs Alpha')
xlabel('Angle of Attack (degrees)')
ylabel('c_L')
xl = xline(16,'--',{'Alpha Stall','16 deg'});
xl.LabelVerticalAlignment = 'bottom';
snapnow
hold off
saveas(cLAlpha,'cLAlpha.png')

cDAlpha = figure;
hold on
grid on
axis([-10 25 0 0.5])
plot(AoA,results.cD,'-o')
title('c_D vs Alpha')
xlabel('Angle of Attack (degrees)')
ylabel('c_D')
xl = xline(16,'--',{'Alpha Stall','16 deg'});
xl.LabelVerticalAlignment = 'bottom';
snapnow
hold off
saveas(cDAlpha,'cDAlpha.png')

%% Hot Wire Calibration
calibrationData = Data.hotWireCalData;
voltage = zeros(length(calibrationData(1).recDat),length(calibrationData));
meanVoltage = zeros(1,length(calibrationData));
airSpeed = zeros(1,length(calibrationData));
for i = 1:length(calibrationData)
    voltage(:,i) = calibrationData(i).recDat;
    airSpeed(1,i) = calibrationData(i).windSpeed;
end

x = voltage.^2;
for i = 1:11
    meanVoltage(i) = mean(x(:,i));
end

HWCal = figure;
hold on
grid on
plot(airSpeed.^0.5.*20,meanVoltage,'-o')
plot(airSpeed.^0.75.*9,meanVoltage,'-o')
plot(airSpeed.^0.9.*6,meanVoltage,'-o')
plot(airSpeed.*4,meanVoltage,'-o')
title('Hot Wire Calibration')
xlabel('U^n')
ylabel('E_0^2')
legend('E_0^2 = 20*U^{0.5}','E_0^2 = 9*U^{0.75}','E_0^2 = 6*U^{0.9}','E_0^2 = 4*U^1')
legend('Location','northwest')
snapnow
hold off
saveas(HWCal,'HWCal.png')

A = 2;
B = 4;
n = 1;

%% Wake Data
wakeNeg5 = WakeProcess(Data.wakeDatan5,A,B,n);
wake5 = WakeProcess(Data.wakeData5,A,B,n);
wake10 = WakeProcess(Data.wakeData10,A,B,n);
wake15 = WakeProcess(Data.wakeData15,A,B,n);
wake17 = WakeProcess(Data.wakeData17,A,B,n);

L0 = 6;

Wake = figure;
hold on
grid on
plot(wakeNeg5.speed/wakeNeg5.windSpeed,wakeNeg5.pos/L0,'-o')
plot(wake5.speed/wake5.windSpeed,wake5.pos/L0,'-+')
plot(wake10.speed/wake10.windSpeed,wake10.pos/L0,'-s')
plot(wake15.speed/wake15.windSpeed,wake15.pos/L0,'-d')
title('Wake Profiles')
xlabel('U_x/U_{inf}')
ylabel('y/L_0')
legend('-5 deg AoA','5 deg AoA','10 deg AoA','15 deg AoA')
legend('Location','northwest')
snapnow
hold off
saveas(Wake,'Wake.png')

AoA17 = figure;
hold on
grid on
plot(wake17.speed/wake17.windSpeed,wake17.pos/L0,'-o')
title('Wake Profile at 17 degrees AoA')
xlabel('U_x/U_{inf}')
ylabel('y/L_0')
snapnow
hold off
saveas(AoA17,'AoA17.png')

TurbInt = figure;
hold on
grid on
plot(wakeNeg5.dev,wakeNeg5.pos/L0,'-o')
plot(wake5.dev,wake5.pos/L0,'-+')
plot(wake10.dev,wake10.pos/L0,'-s')
plot(wake15.dev,wake15.pos/L0,'-d')
plot(wake17.dev,wake17.pos/L0,'-*')
title('Turbulence Intensity')
xlabel('std(U_x)')
ylabel('y/L_0')
legend('-5 deg AoA','5 deg AoA','10 deg AoA','15 deg AoA','17 deg AoA')
legend('Location','northeast')
snapnow
hold off
saveas(TurbInt,'TurbInt.png')

%% Wake Data processing
function out = WakeProcess(data,A,B,n)

pos = zeros(length(data),1);
voltage = zeros(length(data(1).recDat),length(data));
for i = 1:length(data)
    pos(i) = data(i).pos(3);
    voltage(:,i) = data(i).recDat;
end
speed = ((voltage.^2-A)./B).^(1/n);
for i = 1:length(data)
    realSpeed(i) = mean(speed(:,i));
    dev(i) = std(speed(:,i));
end
windSpeed = (realSpeed(1)+realSpeed(end))/2;

[posSorted, order] = sort(pos);
voltSorted = voltage(order);
speedSorted = realSpeed(order);

out.pos = posSorted;
out.voltage = voltSorted;
out.windSpeed = windSpeed;
out.speed = speedSorted;
out.dev = dev;

% out.pos = pos;
% out.voltage = voltage;
% out.windSpeed = windSpeed;
% out.speed = speed;
end

%% Converts forces
function out = Forces(forceNormal,forceAxial,AoA,U)
span = 0.13;
chord = 0.19;
rho = 1.17274;

lift = forceNormal.*cosd(AoA)-forceAxial.*sind(AoA);
drag = forceNormal.*sind(AoA)+forceAxial.*cosd(AoA);
q = rho*U.^2;
S = span*chord;
cL = lift./q./S;
cD = drag./q./S;

out.Lift = lift;
out.Drag = drag;
out.q = q;
out.S = S;
out.cL = cL;
out.cD = cD;
end