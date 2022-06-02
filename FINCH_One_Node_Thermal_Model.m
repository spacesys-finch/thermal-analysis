% FINCH One Node Thermal Model

close all;
clear all;
clc;

%% Overview

% This program is a transient orbital thermal simulation that models the FINCH satellite as a single
% node. The purpose of this program is to identify the on orbit temperature ranges that the FINCH 
% satellite will encounter. Two cases are analyzed, a "hot" case corresponding to winter in the
% northern hemisphere, and a "cold" case corresponding to summer in the northern hemisphere. To
% modify the program, adjust values in the "Constants" section.

% The design of this program is based on the paper "Preliminary Thermal Analysis of Small Satellites" 
% by Casper Versteeg and David L. Cotten

%% Assumptions & Approximations

% FINCH is a sphere made of 6061-T6 aluminum with the same surface area as a 3U CubeSat
% This sphere is a lumped mass (i.e., all points within the sphere are at the same temperature at each time step)
% FINCH's mass is constant
% FINCH's specific heat capacity is constant
% FINCH's absorptivity and emissivity are the same across its entire surface and do not change over time
% FINCH generates the same amount of power each second for the entire duration of this analysis
% FINCH is deployed at the halfway point of the sunlit portion of its orbit
% FINCH is orbiting at a low altitude (i.e., altitude << Earth's radius) and in a circular orbit
% The Earth's albedo is constant over its entire surface
% The Earth has a constant infrared flux emission over its entire surface
% The solar heat flux is constant for the durations of the hot and cold cases
% Initial temperatures for hot and cold cases are the same
% The orbital altitude is high enough such that the effects of free molecular heating are negligible
% The effects of charged-particle heating are negligible


%% Constants

% Orbital Parameters

altitude = 550; % in km
orbital_period = 5708; % in seconds
beta_angle = 22.5; % in degrees (beta angle is constant because FINCH is in a sun-synchronous orbit)
no_orbits = 10; % number of orbits to be completed in this analysis

% FINCH Parameters

mass = 4; % in kg (maximum possible mass)
specific_heat = 896; % specific heat capacity of 6061-T6 aluminum in J/(kg-K)
abosrptivity = 0.7; % unitless
emissivity = 0.52; % unitless
FINCH_radius = 0.1079; % radius of spherical FINCH in m

% Other Constants

sb = 5.67*(10^-8); % Stefan-Boltzmann Constant in W/(m^2-K^4)
earth_albedo = 0.14; % unitless
earth_IR_heat_flux = 228; % in W/m^2
solar_flux_hot = 1414; % solar flux for hot case in W/m^2
solar_flux_cold = 1322; % solar flux for cold case in W/m^2
space_temp = -270.45; % temperature of space in degrees Celsius
Qdot_gen_Wh = 9.28; % internally generated heat load in Watt-hours over a single orbit (based on power budget)
earth_radius = 6371; % in km
delta_t = 1; % define time step in seconds

%% Analysis

% Calculate the internally generated heat load in Watts (need to convert from Wh to Watts)

Qdot_gen = Qdot_gen_Wh/(orbital_period/3600);

% Calculate the critical beta angle and eclipse fraction

critical_beta_angle = asind(earth_radius/(earth_radius+altitude));

if beta_angle < critical_beta_angle
    eclipse_fraction = (1/180)*acosd((sqrt(altitude^2+2*earth_radius*altitude))/((earth_radius+altitude)*cosd(beta_angle)));
else
    eclipse_fraction = 0;
end

% Calculate the area of FINCH exposed to radiation sources in m^2 (a cirlce in this case)
area = pi*(FINCH_radius)^2;

% Initialize the index variable
i = 1;

% Initialize the spacecraft's temperature in degrees Celsius for both cases
T_hot(i) = 20;
T_cold(i) = 20;

% Initialize time vector at halfway point of sunlit portion of the orbit (units: seconds)
t(i) = 0;

% Calculate the heat transfer rates at time t=0
Qdot_earth_IR = earth_IR_heat_flux*area;
% Hot Case
Qdot_solar_hot(i) = solar_flux_hot*area*abosrptivity;
Qdot_albedo_hot(i) = earth_albedo*Qdot_solar_hot(i);
% Cold Case
Qdot_solar_cold(i) = solar_flux_cold*area*abosrptivity;
Qdot_albedo_cold(i) = earth_albedo*Qdot_solar_cold(i);

% Calculate the spacecraft temperatures over the desired number of orbits for both cases

for j = 0 : 1 : (no_orbits-1)

    while t(i) <= (orbital_period*(j+1) - 1)

        % On to the next time step
        i = i+1;
        t(i) = t(i-1)+delta_t;
        
        % Extra "if statement" in case we choose an increment such that an additional temperature point is calculated beyond the specified number of orbits
        if t <= orbital_period*no_orbits

            % Decide whether the solar radiation is turned on or off
            % Solar radiation is "turned off" when FINCH enters eclipse and "turned on" when it exits eclipse
            % This "if statement" outputs a factor that, when applied to the solar heat term, either turns solar heating on or off

            if t(i) < orbital_period*0.5*(1-eclipse_fraction)+(orbital_period*j) || t(i) > orbital_period*0.5*(1+eclipse_fraction)+(orbital_period*j)
                solar_switch = 1;
            else
                solar_switch = 0;
            end

            % Calculate the net heat transfer rate for each time step of the hot case
            Qdot_solar_hot(i) = solar_flux_hot*area*solar_switch*abosrptivity;
            Qdot_albedo_hot(i) = earth_albedo*Qdot_solar_hot(i);
            Qdot_radiation_hot = 4*area*sb*emissivity*(((T_hot(i-1)+273.15)^4)-((space_temp+273.15)^4));
            Qdot_hot = Qdot_earth_IR + Qdot_solar_hot(i) + Qdot_albedo_hot(i) + Qdot_gen - Qdot_radiation_hot;

            % Calculate the temperature for each time step of the hot case
            T_hot(i) = T_hot(i-1) + (delta_t*Qdot_hot)/(specific_heat*mass);

            % Calculate the net heat transfer rate for each time step of the cold case
            Qdot_solar_cold(i) = solar_flux_cold*area*solar_switch*abosrptivity;
            Qdot_albedo_cold(i) = earth_albedo*Qdot_solar_cold(i);
            Qdot_radiation_cold = 4*area*sb*emissivity*(((T_cold(i-1)+273.15)^4)-((space_temp+273.15)^4));
            Qdot_cold = Qdot_earth_IR + Qdot_solar_cold(i) + Qdot_albedo_cold(i) + Qdot_gen - Qdot_radiation_cold;

            % Calculate the temperature for each time step of the cold case
            T_cold(i) = T_cold(i-1) + (delta_t*Qdot_cold)/(specific_heat*mass);

        else
            % If we accidentally step beyond our desired number of orbits, the final time step is removed from the time vector
            t(i) = [];
            break
        end
        
    end

end

%% Plotting

% Plot the heat loads over one orbit for the hot case
figure(1);
set(1,'WindowStyle','Docked');
grid on
title('FINCH One Node Hot Case Heat Loads over a Single Orbit','fontsize', 15)
xlabel('Time from Deployment (s)', 'fontsize', 15);
ylabel('Heat Load (W)', 'fontsize',15);
hold on;
plot([0,orbital_period],[Qdot_gen,Qdot_gen],'g','LineWidth',2);
hold on;
plot([0,orbital_period],[Qdot_earth_IR,Qdot_earth_IR],'LineWidth',2);
hold on;
plot(t(1:orbital_period),Qdot_solar_hot(1:orbital_period),'LineWidth',2);
hold on;
plot(t(1:orbital_period),Qdot_albedo_hot(1:orbital_period),'LineWidth',2); 
legend('Generated Heat','Earth IR','Direct Solar','Reflected Solar','Location','southeast');

% Plot the heat loads over one orbit for the cold case
figure(2);
set(2,'WindowStyle','Docked');
grid on
title('FINCH One Node Cold Case Heat Loads over a Single Orbit','fontsize', 15)
xlabel('Time from Deployment (s)', 'fontsize', 15);
ylabel('Heat Load (W)', 'fontsize',15);
hold on;
plot([0,orbital_period],[Qdot_gen,Qdot_gen],'g','LineWidth',2);
hold on;
plot([0,orbital_period],[Qdot_earth_IR,Qdot_earth_IR],'LineWidth',2);
hold on;
plot(t(1:orbital_period),Qdot_solar_cold(1:orbital_period),'LineWidth',2);
hold on;
plot(t(1:orbital_period),Qdot_albedo_cold(1:orbital_period),'LineWidth',2);
legend('Generated Heat','Earth IR','Direct Solar','Reflected Solar','Location','southeast');

% Plot the temperatures for both cases
figure(3);
set(3,'WindowStyle','Docked');
grid on
title('FINCH One Node Temperatures over Ten Orbits','fontsize', 15);
xlabel('Time from Deployment (s)', 'fontsize', 15);
ylabel('Temperature (°C)', 'fontsize',15);
hold on;
plot(t,T_hot,'r');
hold on
plot(t,T_cold,'b');
legend('Hot Case','Cold Case')