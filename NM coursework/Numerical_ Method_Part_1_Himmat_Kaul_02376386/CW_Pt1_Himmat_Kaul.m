
%clears comand window
clc
%clears veh structure
clear veh
%closes all files
close all
%formats long format for numerical data types
format long

set(groot,'defaultLineLineWidth',2)  %sets graph line width as 2
set(groot,'defaultAxesFontSize',12)  %sets graph axes font size as 12
set(groot,'defaulttextfontsize',12)  %sets graph text font size as 12
set(groot,'defaultLineMarkerSize',10) %sets line marker size as 10
set(groot,'defaultAxesXGrid','on')   %sets X axis grid on 
set(groot,'defaultAxesYGrid','on')   %sets Y axis grid on
set(groot, 'DefaultAxesBox', 'on')   %sets Axes boxes on

%loading planetdata.mat files as structured arrays
load("planetdata.mat");
%loading MH_data.mat file as structured array
load("MH_data.mat");

%reads cl and cd data from imported .csv files into matricies under veh
%structure
veh.cd = readmatrix('MH_Cddata.csv');
%clears the NaN values in the initial 2 rows
veh.cd([1,2],:) = [];
%reads cl and cd data from imported .csv files into matricies under veh
%structure
veh.cl = readmatrix('MH_Cldata.csv');
%clears the NaN values in the initial 2 rows
veh.cl([1,2],:) = [];

%reads chord and twist distridutions for the blade and writes to the
%variables MH_chord and MH_twist 
MH_chord = readmatrix('MH_chord.csv');
MH_twist = readmatrix("MH_twist.csv");

%Defines the number of descritisation points as N
N = 1000;
%defines the radial station values along the blade
r = linspace(0.09, 1, N);
%stores radial station values into structure veh
veh.r = r;
%interpolates and extrapolates the MH_chord and MH_twist data for the Radial stations
%along the blade
veh.chorddistribution = interp1(MH_chord(:,1),MH_chord(:,2),r,"linear","extrap");
veh.twistdistribution = interp1(MH_twist(:,1),MH_twist(:,2),r,"linear","extrap");

%defines the table of given Bilinear chord distribution 
bilinearchorddistribution = [0.09,0.34,1.0;0.05,0.2,0.07];
%defines the table of given bilinear twist distribution
bilineartwistdistribution = [0.09,0.2,1;16,18,0];

%interpolates Bilinear chord distribution for ploting
veh.bilineartwist = interp1(bilineartwistdistribution(1,:),bilineartwistdistribution(2,:),r);
%interpolates Bilinear twist distribution for ploting
veh.bilinearchord = interp1(bilinearchorddistribution(1,:),bilinearchorddistribution(2,:),r);

%evaluates the BET function for the MH data
MH_actualdata = BET(env,veh);


%Chord plots for MH data, Interpolated distribution from MH data and Simplified bilinear distribution  
ChordDistributionFigure = figure;
set(ChordDistributionFigure,'WindowState','maximized')
plot(MH_chord(:,1),MH_chord(:,2),"x");
hold on
plot(r,veh.chorddistribution, LineStyle=":")
hold on
plot(r,veh.bilinearchord, LineStyle="--")

title('Mars Helicopter Chord Distributions')
xlabel("Non-Dimentional radial distance, (r/R)")
ylabel("Normalised chord length, (c/R)")
legend("Mars Helicopter Distribution",'Interpolated Distribution','Simplified Bilinear Distribution',Location = 'best')
grid on
hold off

saveas(ChordDistributionFigure,'Mars Helicopter Chord Distributions.png')

%Twist plots for MH data, Interpolated distribution from MH data and Simplified bilinear distribution  
TwistDistributionFigure = figure;
set(TwistDistributionFigure,'WindowState','maximized')
plot(MH_twist(:,1),MH_twist(:,2),"x");
hold on
plot(r,veh.twistdistribution, LineStyle=":")
hold on
plot(r,veh.bilineartwist, LineStyle="--")
title('Mars Helicopter Twist Distributions')
xlabel("Non-Dimentional radial distance, (r/R)")
ylabel("twist angle, (°)")
legend("Mars Helicopter Distribution",'Interpolated Distribution','Simplified Bilinear Distribution',Location='best')
grid on
hold off

saveas(TwistDistributionFigure,'Mars Helicopter Twist Distributions.png')

%Cl plots for each mach number from the veh.cl array
Clplot = figure;
set(Clplot,'WindowState','maximized')
plot(veh.cl(:,1),veh.cl(:,2),LineStyle="-");
hold on
plot(veh.cl(:,3),veh.cl(:,4),LineStyle="--");
hold on
plot(veh.cl(:,5),veh.cl(:,6),LineStyle=":");
hold on
plot(veh.cl(:,7),veh.cl(:,8),LineStyle="-.");
hold on
plot(veh.cl(:,9),veh.cl(:,10),'-p');
hold on

title("Cl vs α, for different Mach Numbers")
xlabel("Angle of Attack, α (°) ")
ylabel("Lift Coefficient, Cl")
legend('M = 0.2','M = 0.4','M = 0.6','M = 0.8','M = 0.9',Location='best')
grid on
hold off

saveas(Clplot,'Cl vs alpha, for different Mach Numbers.png')

%Cl plots for each mach number from the veh.cl array
Cdplot = figure;
set(Cdplot,'WindowState','maximized')
plot(veh.cd(:,1),veh.cd(:,2),LineStyle="-")
hold on
plot(veh.cd(:,3),veh.cd(:,4),LineStyle="--")
hold on
plot(veh.cd(:,5),veh.cd(:,6),LineStyle=":")
hold on
plot(veh.cd(:,7),veh.cd(:,8),LineStyle="-.")
hold on
plot(veh.cd(:,9),veh.cd(:,10),"-p")
hold on

title("Cd vs α, for different Mach Numbers")
xlabel("Angle of Attack, α (°) ")
ylabel("Drag Coefficient, Cd")
legend('M = 0.2','M = 0.4','M = 0.6','M = 0.8','M = 0.9',Location='best')
grid on
hold off

saveas(Cdplot,'Cd vs alpha, for different Mach Numbers.png')

%SENSETIVITY STUDY%
%initialises arrays for storing BET function outputs for all studies 
study1 = []; %Varying mid Point radial location in the bilinear chord distribution
study2 = []; %Varying the mid point Non dimentional chord length in the bilinear chord distribution 
study3 = []; %Varying mid point radial location of the Bilinear twist Distributiom
study4 = []; %Varying the mid point twist in the bilinears twist distribution

%Evaluates the BET funtion for Study 1 and stores outputs into array study1
for i = 1:21
    bilinearchorddistribution = [0.09,0.34*(0.9+(i-1)*0.01),1.0;0.05,0.2,0.07];
    bilineartwistdistribution = [0.09,0.2,1;16,18,0];
    veh.twistdistribution = interp1(bilineartwistdistribution(1,:),bilineartwistdistribution(2,:),r);
    veh.chorddistribution = interp1(bilinearchorddistribution(1,:),bilinearchorddistribution(2,:),r);
    study1(i,:) = BET(env,veh);
end

%Evaluates the BET funtion for Study 2 and stores outputs into array study2
for i = 1:21
    bilinearchorddistribution = [0.09,0.34,1.0;0.05,0.2*(0.9+(i-1)*0.01),0.07];
    bilineartwistdistribution = [0.09,0.2,1;16,18,0];
    veh.twistdistribution = interp1(bilineartwistdistribution(1,:),bilineartwistdistribution(2,:),r);
    veh.chorddistribution = interp1(bilinearchorddistribution(1,:),bilinearchorddistribution(2,:),r);
    study2(i,:) = BET(env,veh);
end

%Evaluates the BET funtion for Study 3 and stores outputs into array study3
for i = 1:21
    bilinearchorddistribution = [0.09,0.34,1.0;0.05,0.2,0.07];
    bilineartwistdistribution = [0.09,0.2*(0.9+(i-1)*0.01),1;16,18,0];
    veh.twistdistribution = interp1(bilineartwistdistribution(1,:),bilineartwistdistribution(2,:),r);
    veh.chorddistribution = interp1(bilinearchorddistribution(1,:),bilinearchorddistribution(2,:),r);
    study3(i,:) = BET(env,veh);
end

%Evaluates the BET funtion for Study 4 and stores outputs into array study4
for i = 1:21
    bilinearchorddistribution = [0.09,0.34,1.0;0.05,0.2,0.07];
    bilineartwistdistribution = [0.09,0.2,1;16,18*(0.9+(i-1)*0.01),0];
    veh.twistdistribution = interp1(bilineartwistdistribution(1,:),bilineartwistdistribution(2,:),r);
    veh.chorddistribution = interp1(bilinearchorddistribution(1,:),bilinearchorddistribution(2,:),r);
    study4(i,:) = BET(env,veh);
end

%plotting Study 1

BilinearChordDist = figure;
set(BilinearChordDist,'WindowState','maximized')
subplot(2,2,1)
variation = 0.34*(0.9+0.01.*(1:21));
plot(variation([1:9,11:21]), study1([1:9,11:21],1),linestyle = "none",marker="+")
hold on
plot(variation(10),MH_actualdata(1),"o",'color','r')
hold on
plot(variation(10),study1(10,1),'pentagram','color','g')

title("Sensetivity study of Radial location of Station B1")
xlabel("Station Radial Location (y/R)")
ylabel("Figure of Merit (FM)")
legend('Sensetivity Study (+/-10%)','MH FM','Bilinear Distribution',Location = 'best')
grid on
hold off
%%%

subplot(2,2,2);
variation = 0.34*(0.9+0.01.*(1:21));
plot(variation([1:9,11:21]), study1([1:9,11:21],3),linestyle = "none",marker="+")
hold on
plot(variation(10),MH_actualdata(3),"o",'color','r')
hold on
plot(variation(10),study1(10,3),'pentagram','color','g')

title("Sensetivity study of Radial location of Station B1")
xlabel("Station Radial Location (y/R)")
ylabel("Power Coefficient, (CP)")
legend('Sensetivity Study (+/-10%)','MH CP','Bilinear Distribution',Location = 'best')
grid on
hold off

%plotting Study 2

subplot(2,2,3)
variation = 0.34*(0.9+0.01.*(1:21));
plot(variation([1:9,11:21]), study2([1:9,11:21],1),linestyle = "none",marker="+")
hold on
plot(variation(10),MH_actualdata(1),"o",'color','r')
hold on
plot(variation(10),study2(10,1),'pentagram','color','g')

title("Sensetivity study of Chord value of Station B1")
xlabel("Station Radial Location (y/R)")
ylabel("Figure of Merit (FM)")
legend('Sensetivity Study (+/-10%)','MH FM','Bilinear Distribution',Location='best')
grid on
hold off
%%%

subplot(2,2,4);
variation = 0.34*(0.9+0.01.*(1:21));
plot(variation([1:9,11:21]), study2([1:9,11:21],3),linestyle = "none",marker="+")
hold on
plot(variation(10),MH_actualdata(3),"o",'color','r')
hold on
plot(variation(10),study2(10,3),'pentagram','color','g')

title("Sensetivity study of Chord value of Station B1")
xlabel("Station Radial Location (y/R)")
ylabel("Power Coefficient, (CP)")
legend('Sensetivity Study (+/-10%)','MH CP','Bilinear Distribution',Location='best')
grid on
hold off

saveas(BilinearChordDist,'Sensetivity study of Chord Distribution.png')

%plotting Study 3

BilinearTwistDist = figure;
set(BilinearTwistDist,'WindowState','maximized')
subplot(2,2,1)
variation = 0.34*(0.9+0.01.*(1:21));
plot(variation([1:9,11:21]), study3([1:9,11:21],1),linestyle = "none",marker="+")
hold on
plot(variation(10),MH_actualdata(1),"o",'color','r')
hold on
plot(variation(10),study3(10,1),'pentagram','color','g')

title("Sensetivity study of Radial location of Station B2")
xlabel("Station Radial Location (y/R)")
ylabel("Figure of Merit (FM)")
legend('Sensetivity Study (+/-10%)','MH FM','Bilinear Distribution',Location='best')
grid on
hold off
%%%

subplot(2,2,2);
variation = 0.34*(0.9+0.01.*(1:21));
plot(variation([1:9,11:21]), study3([1:9,11:21],3),linestyle = "none",marker="+")
hold on
plot(variation(10),MH_actualdata(3),"o",'color','r')
hold on
plot(variation(10),study3(10,3),'pentagram','color','g')

title("Sensetivity study of Radial location of Station B2")
xlabel("Station Radial Location (y/R)")
ylabel("Power Coefficient, (CP)")
legend('Sensetivity Study (+/-10%)','MH CP','Bilinear Distribution',Location='best')
grid on
hold off

%plotting Study 4 

subplot(2,2,3)
variation = 0.34*(0.9+0.01.*(1:21));
plot(variation([1:9,11:21]), study4([1:9,11:21],1),linestyle = "none",marker="+")
hold on
plot(variation(10),MH_actualdata(1),"o",'color','r')
hold on
plot(variation(10),study4(10,1),'pentagram','color','g')

title("Sensetivity study of Twist value of Station B2")
xlabel("Station Radial Location (y/R)")
ylabel("Figure of Merit (FM)")
legend('Sensetivity Study (+/-10%)','MH FM','Bilinear Distribution',Location='best')
grid on
hold off
%%%

subplot(2,2,4);
variation = 0.34*(0.9+0.01.*(1:21));
plot(variation([1:9,11:21]), study4([1:9,11:21],3),linestyle = "none",marker="+")
hold on
plot(variation(10),MH_actualdata(3),"o",'color','r')
hold on
plot(variation(10),study4(10,3),'pentagram','color','g')

title("Sensetivity study of Twist value of Station B2")
xlabel("Station Radial Location (y/R)")
ylabel("Power Coefficient, (CP)")
legend('Sensetivity Study (+/-10%)','MH CP','Bilinear Distribution',Location='best')
grid on
hold off

saveas(BilinearTwistDist,'Sensetivity study of Twist Distribution.png')

