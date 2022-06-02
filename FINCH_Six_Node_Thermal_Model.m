% FINCH Six Node Thermal Model

close all;
clear all;
clc;

%% Overview

% This program is a transient orbital thermal simulation that models the FINCH satellite as six
% thermally coupled nodes. The purpose of this program is to identify the on orbit temperature ranges 
% that the FINCH satellite will encounter. Two cases are analyzed, a "hot" case corresponding to winter 
% in the northern hemisphere, and a "cold" case corresponding to summer in the northern hemisphere. To
% modify the program, adjust values in the "Constants" section.

% The design of this program is based on the paper "Preliminary Thermal Analysis of Small Satellites" 
% by Casper Versteeg and David L. Cotten
% The view factors for Earth IR and reflected solar radiation are from the paper "A Simplified, Closed-Form 
% Method for Screening Spacecraft Orbital Heating Variations" by S. L. Rickman

%% Assumptions & Approximations

% FINCH consists of six thermally coupled nodes
% Each of these nodes is a lumped mass (i.e., all points within each node are at the same temperature at each time step)
% The mass of each node is constant
% Nodes are thermally coupled along spacecraft edges (i.e., each node includes a single external spacecraft panel)
% The total mass of the satellite is distributed among the nodes such that the mass-area ratio remains constant
% The specific heat capacity (at constant pressure) of each node is that of aluminum 6061-T6 and is constant
% Each node's absorptivity and emissivity are the same across its entire external surface and do not change over time
% FINCH's internal power generation is constant
% This power is distributed among the nodes such that the power-area ratio remains constant
% All internal radiation is ignored
% FINCH is deployed at the halfway point of the sunlit portion of its orbit
% FINCH is orbiting at a low altitude (i.e., altitude << Earth's radius) and in a circular orbit
% FINCH is oriented in a fixed attitude with the -Z node always facing Earth, the +Z node always facing away 
% from the Earth, the +Y node acting as the "leading face" (i.e., its surface normal points in the direction of motion), 
% and the -X node never facing the Sun
% The Earth's albedo is constant over its entire surface
% The Earth has a constant infrared flux emission over its entire surface
% The solar heat flux is constant for the durations of the hot and cold cases
% Initial temperatures for the hot and cold cases for all nodes are the same
% The orbital altitude is high enough such that the effects of free molecular heating are negligible
% The effects of charged-particle heating are negligible


%% Constants

% Orbital Parameters

altitude = 550; % in km
orbital_period = 5708; % in seconds
beta_angle = 22.5; % in degrees (beta angle is constant because FINCH is in a sun-synchronous orbit)
no_orbits = 10; % number of orbits to be completed in this analysis

% FINCH Parameters

specific_heat = 896; % specific heat capacity of 6061-T6 aluminum in J/(kg-K)
large_node_mass = 0.8675; % in kg
small_node_mass = 0.2649; % in kg
large_node_bolts = 5; % number of bolts joining adjacent large nodes
small_node_bolts = 3; % number of bolts joining small nodes to adjacent large nodes
% Note: "large" and "small" refer to the relative sizes of the areas of the panels contained within each node (i.e., ±X and ±Y panels are 
% contained within "large" nodes, while the ±Z panels are contained within the "small" nodes)

% +X Node
plusX_absorptivity = 0.7507; % unitless
plusX_emissivity = 0.619; % unitless

% -X Node (includes radiator panel, coated in Z93 paint - property values from Gilmore Thermal Control Handbook Vol. 1)
minusX_absorptivity = 0.17; % unitless
minusX_emissivity = 0.92; % unitless

% +Y Node
plusY_absorptivity = 0.7507; % unitless
plusY_emissivity = 0.619; % unitless

% -Y Node
minusY_absorptivity = 0.7507; % unitless
minusY_emissivity = 0.619; % unitless

% +Z Node
plusZ_absorptivity = 0.379; % unitless
plusZ_emissivity = 0.08; % unitless

% -Z Node
minusZ_absorptivity = 0.379; % unitless
minusZ_emissivity = 0.08; % unitless

small_node_area = 0.01; % area in m^2 of ±Z panels contained within the "small" ±Z nodes
large_node_area = 0.03275; % area in m^2 of ±X and ±Y panels contained within the "large" ±X and ±Y nodes
bolt_conductance = 0.26; % thermal conductance of bolted connections between node in W/K (also from Gilmore's book)

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

% Distribute the internally generated heat load among the nodes

small_node_power = Qdot_gen/15.1;
large_node_power = 3.275*small_node_power;

% Calculate the critical beta angle and eclipse fraction

critical_beta_angle = asind(earth_radius/(earth_radius+altitude));

if beta_angle < critical_beta_angle
    eclipse_fraction = (1/180)*acosd((sqrt(altitude^2+2*earth_radius*altitude))/((earth_radius+altitude)*cosd(beta_angle)));
else
    eclipse_fraction = 0;
end

% Calculate the view factors for Earth IR and reflected solar radiation (from Rickman paper)
% Since the +Z (zenith) node is never oriented towards the Earth, its view factor is always 0
VF_nadir =(earth_radius/(earth_radius+altitude))^2; % Earth Facing
a=sqrt(1-VF_nadir);
b=2*asin(a);
c=sin(b);
VF_walls = (pi-b-c)/(2*pi); % Perpendicular to nadir

% Initialize the index variable
i = 1;

% Initialize the spacecraft's temperature in degrees Celsius for all nodes for both cases
T_hot_plusX(i) = 20;
T_cold_plusX(i) = 20;
T_hot_minusX(i) = 20;
T_cold_minusX(i) = 20;
T_hot_plusY(i) = 20;
T_cold_plusY(i) = 20;
T_hot_minusY(i) = 20;
T_cold_minusY(i) = 20;
T_hot_plusZ(i) = 20;
T_cold_plusZ(i) = 20;
T_hot_minusZ(i) = 20;
T_cold_minusZ(i) = 20;

% Initialize time vector at halfway point of sunlit portion of the orbit (units: seconds)
t(i) = 0;

% Calculate the heat transfer rates for the +X node at time t=0
Qdot_earth_IR_plusX = VF_walls*earth_IR_heat_flux*large_node_area;
% Hot Case
Qdot_solar_hot_plusX(i) = sind(beta_angle)*solar_flux_hot*large_node_area*plusX_absorptivity;
Qdot_albedo_hot_plusX(i) = VF_walls*earth_albedo*solar_flux_hot*plusX_absorptivity*large_node_area*cosd(beta_angle);
% Cold Case
Qdot_solar_cold_plusX(i) = sind(beta_angle)*solar_flux_cold*large_node_area*plusX_absorptivity;
Qdot_albedo_cold_plusX(i) = VF_walls*earth_albedo*solar_flux_cold*plusX_absorptivity*large_node_area*cosd(beta_angle);

% Calculate the heat transfer rates for the -X node at time t=0
Qdot_earth_IR_minusX = VF_walls*earth_IR_heat_flux*large_node_area;
% Hot Case
Qdot_solar_hot_minusX(i) = 0;
Qdot_albedo_hot_minusX(i) = VF_walls*earth_albedo*solar_flux_hot*minusX_absorptivity*large_node_area*cosd(beta_angle);
% Cold Case
Qdot_solar_cold_minusX(i) = 0;
Qdot_albedo_cold_minusX(i) = VF_walls*earth_albedo*solar_flux_cold*minusX_absorptivity*large_node_area*cosd(beta_angle);

% Calculate the heat transfer rates for the +Y node at time t=0
Qdot_earth_IR_plusY = VF_walls*earth_IR_heat_flux*large_node_area;
% Hot Case
Qdot_solar_hot_plusY(i) = 0;
Qdot_albedo_hot_plusY(i) = VF_walls*earth_albedo*solar_flux_hot*plusY_absorptivity*large_node_area*cosd(beta_angle);
% Cold Case
Qdot_solar_cold_plusY(i) = 0;
Qdot_albedo_cold_plusY(i) = VF_walls*earth_albedo*solar_flux_cold*plusY_absorptivity*large_node_area*cosd(beta_angle);

% Calculate the heat transfer rates for the -Y node at time t=0
Qdot_earth_IR_minusY = VF_walls*earth_IR_heat_flux*large_node_area;
% Hot Case
Qdot_solar_hot_minusY(i) = 0;
Qdot_albedo_hot_minusY(i) = VF_walls*earth_albedo*solar_flux_hot*minusY_absorptivity*large_node_area*cosd(beta_angle);
% Cold Case
Qdot_solar_cold_minusY(i) = 0;
Qdot_albedo_cold_minusY(i) = VF_walls*earth_albedo*solar_flux_cold*minusY_absorptivity*large_node_area*cosd(beta_angle);

% Calculate the heat transfer rates for the +Z node at time t=0
% The +Z node never faces the Earth, so Earth IR and reflected solar radiation heat transfer to this node is 0
% Hot Case
Qdot_solar_hot_plusZ(i) = cos((2*pi*t(i))/orbital_period)*cosd(beta_angle)*solar_flux_hot*small_node_area*plusZ_absorptivity;
% Cold Case
Qdot_solar_cold_plusZ(i) = cos((2*pi*t(i))/orbital_period)*cosd(beta_angle)*solar_flux_cold*small_node_area*plusZ_absorptivity;

% Calculate the heat transfer rates for the -Z node at time t=0
Qdot_earth_IR_minusZ = VF_nadir*earth_IR_heat_flux*small_node_area;
% Hot Case
Qdot_solar_hot_minusZ(i) = 0;
Qdot_albedo_hot_minusZ(i) = VF_nadir*earth_albedo*solar_flux_hot*minusZ_absorptivity*small_node_area*cosd(beta_angle);
% Cold Case
Qdot_solar_cold_minusZ(i) = 0;
Qdot_albedo_cold_minusZ(i) = VF_nadir*earth_albedo*solar_flux_cold*minusZ_absorptivity*small_node_area*cosd(beta_angle);

% Calculate the total heat loads for both cases

total_heat_load_plusX_hot(i) = Qdot_earth_IR_plusX + Qdot_solar_hot_plusX(i) + Qdot_albedo_hot_plusX(i) + large_node_power;
total_heat_load_minusX_hot(i) = Qdot_earth_IR_minusX + Qdot_albedo_hot_minusX(i) + large_node_power;
total_heat_load_plusY_hot(i) = Qdot_earth_IR_plusY + Qdot_solar_hot_plusY(i) + Qdot_albedo_hot_plusY(i) + large_node_power;
total_heat_load_minusY_hot(i) = Qdot_earth_IR_minusY + Qdot_solar_hot_minusY(i) + Qdot_albedo_hot_minusY(i) + large_node_power;
total_heat_load_plusZ_hot(i) = Qdot_solar_hot_plusZ(i) + small_node_power;
total_heat_load_minusZ_hot(i) = Qdot_earth_IR_minusZ + Qdot_solar_hot_minusZ(i) + Qdot_albedo_hot_minusZ(i) + small_node_power;

total_heat_load_plusX_cold(i) = Qdot_earth_IR_plusX + Qdot_solar_cold_plusX(i) + Qdot_albedo_cold_plusX(i) + large_node_power;
total_heat_load_minusX_cold(i) = Qdot_earth_IR_minusX + Qdot_albedo_cold_minusX(i) + large_node_power;
total_heat_load_plusY_cold(i) = Qdot_earth_IR_plusY + Qdot_solar_cold_plusY(i) + Qdot_albedo_cold_plusY(i) + large_node_power;
total_heat_load_minusY_cold(i) = Qdot_earth_IR_minusY + Qdot_solar_cold_minusY(i) + Qdot_albedo_cold_minusY(i) + large_node_power;
total_heat_load_plusZ_cold(i) = Qdot_solar_cold_plusZ(i) + small_node_power;
total_heat_load_minusZ_cold(i) = Qdot_earth_IR_minusZ + Qdot_solar_cold_minusZ(i) + Qdot_albedo_cold_minusZ(i) + small_node_power;

% Calculate the spacecraft temperatures over the desired number of orbits for both cases

for j = 0 : 1 : (no_orbits-1)

    while t(i) <= (orbital_period*(j+1) - 1)

        % On to the next time step
        i = i+1;
        t(i) = t(i-1)+delta_t;
        
        % Extra "if statement" in case we choose an increment such that an additional temperature point is calculated beyond the specified number of orbits
        if t <= orbital_period*no_orbits
             
            %% Calculate the direct solar radiation view factor for each node
            % The -X node never faces the Sun, so its view factor is always 0)
            
            if t(i) < orbital_period*0.5*(1-eclipse_fraction)+(orbital_period*j) || t(i) > orbital_period*0.5*(1+eclipse_fraction)+(orbital_period*j)
                plusX_VF = sind(beta_angle);
            else
                plusX_VF = 0;
            end

            if t(i) > orbital_period*0.5*(1+eclipse_fraction)+(orbital_period*j)
                plusY_VF = -sin((2*pi*t(i))/orbital_period)*cosd(beta_angle);
            else
                plusY_VF = 0;
            end

            if t(i) < orbital_period*0.5*(1-eclipse_fraction)+(orbital_period*j)
                minusY_VF = sin((2*pi*t(i))/orbital_period)*cosd(beta_angle);
            else
                minusY_VF = 0;
            end

            if t(i) < 0.25*orbital_period+(orbital_period*j) || t(i) > 0.75*orbital_period+(orbital_period*j)
                plusZ_VF = cos((2*pi*t(i))/orbital_period)*cosd(beta_angle);
            else
                plusZ_VF = 0;
            end

            if (t(i) > 0.25*orbital_period+(orbital_period*j) && t(i) < 0.5*orbital_period*(1-eclipse_fraction)+(orbital_period*j)) || (t(i) > 0.5*orbital_period*(1+eclipse_fraction)+(orbital_period*j) && t(i) < 0.75*orbital_period+(orbital_period*j))
                minusZ_VF = -cos((2*pi*t(i))/orbital_period)*cosd(beta_angle);
            else
                minusZ_VF = 0;
            end
            %% Calculate the albedo cosine factor
            % This factor either activates or deactivates the reflected solar heating in orbit

            if t(i) < 0.25*orbital_period+(orbital_period*j) || t(i) > 0.75*orbital_period+(orbital_period*j)
                albedo_switch = 1;
            else
                albedo_switch = 0;
            end
            %% Hot Case

            % Calculate the net heat transfer rate for the +X node for each time step of the hot case
            Qdot_solar_hot_plusX(i) = plusX_VF*solar_flux_hot*large_node_area*plusX_absorptivity;
            Qdot_albedo_hot_plusX(i) = albedo_switch*VF_walls*earth_albedo*solar_flux_hot*plusX_absorptivity*large_node_area*cosd(beta_angle)*cos((2*pi*t(i))/orbital_period);
            Qdot_radiation_hot_plusX = large_node_area*sb*plusX_emissivity*(((T_hot_plusX(i-1)+273.15)^4)-((space_temp+273.15)^4));
            Qdot_bolts_hot_plusX = bolt_conductance*(large_node_bolts*(T_hot_plusY(i-1)+T_hot_minusY(i-1)-2*T_hot_plusX(i-1))+small_node_bolts*(T_hot_plusZ(i-1)+T_hot_minusZ(i-1)-2*T_hot_plusX(i-1)));
            Qdot_hot_plusX = Qdot_earth_IR_plusX + Qdot_solar_hot_plusX(i) + Qdot_albedo_hot_plusX(i) + large_node_power - Qdot_radiation_hot_plusX + Qdot_bolts_hot_plusX;
            
            % Calculate the +X node temperature for each time step of the hot case
            T_hot_plusX(i) = T_hot_plusX(i-1) + (delta_t*Qdot_hot_plusX)/(specific_heat*large_node_mass);

            % Calculate the net heat transfer rate for the -X node for each time step of the hot case
            Qdot_albedo_hot_minusX(i) = albedo_switch*VF_walls*earth_albedo*solar_flux_hot*minusX_absorptivity*large_node_area*cosd(beta_angle)*cos((2*pi*t(i))/orbital_period);
            Qdot_radiation_hot_minusX = large_node_area*sb*minusX_emissivity*(((T_hot_minusX(i-1)+273.15)^4)-((space_temp+273.15)^4));
            Qdot_bolts_hot_minusX = bolt_conductance*(large_node_bolts*(T_hot_plusY(i-1)+T_hot_minusY(i-1)-2*T_hot_minusX(i-1))+small_node_bolts*(T_hot_plusZ(i-1)+T_hot_minusZ(i-1)-2*T_hot_minusX(i-1)));
            Qdot_hot_minusX = Qdot_earth_IR_minusX + Qdot_albedo_hot_minusX(i) + large_node_power - Qdot_radiation_hot_minusX + Qdot_bolts_hot_minusX;
            
            % Calculate the -X node temperature for each time step of the hot case
            T_hot_minusX(i) = T_hot_minusX(i-1) + (delta_t*Qdot_hot_minusX)/(specific_heat*large_node_mass);

            % Calculate the net heat transfer rate for the +Y node for each time step of the hot case
            Qdot_solar_hot_plusY(i) = plusY_VF*solar_flux_hot*large_node_area*plusY_absorptivity;
            Qdot_albedo_hot_plusY(i) = albedo_switch*VF_walls*earth_albedo*solar_flux_hot*plusY_absorptivity*large_node_area*cosd(beta_angle)*cos((2*pi*t(i))/orbital_period);
            Qdot_radiation_hot_plusY = large_node_area*sb*plusY_emissivity*(((T_hot_plusY(i-1)+273.15)^4)-((space_temp+273.15)^4));
            Qdot_bolts_hot_plusY = bolt_conductance*(large_node_bolts*(T_hot_plusX(i-1)+T_hot_minusX(i-1)-2*T_hot_plusY(i-1))+small_node_bolts*(T_hot_plusZ(i-1)+T_hot_minusZ(i-1)-2*T_hot_plusY(i-1)));
            Qdot_hot_plusY = Qdot_earth_IR_plusY + Qdot_solar_hot_plusY(i) + Qdot_albedo_hot_plusY(i) + large_node_power - Qdot_radiation_hot_plusY + Qdot_bolts_hot_plusY;
            
            % Calculate the +Y node temperature for each time step of the hot case
            T_hot_plusY(i) = T_hot_plusY(i-1) + (delta_t*Qdot_hot_plusY)/(specific_heat*large_node_mass);

            % Calculate the net heat transfer rate for the -Y node for each time step of the hot case
            Qdot_solar_hot_minusY(i) = minusY_VF*solar_flux_hot*large_node_area*minusY_absorptivity;
            Qdot_albedo_hot_minusY(i) = albedo_switch*VF_walls*earth_albedo*solar_flux_hot*minusY_absorptivity*large_node_area*cosd(beta_angle)*cos((2*pi*t(i))/orbital_period);
            Qdot_radiation_hot_minusY = large_node_area*sb*minusY_emissivity*(((T_hot_minusY(i-1)+273.15)^4)-((space_temp+273.15)^4));
            Qdot_bolts_hot_minusY = bolt_conductance*(large_node_bolts*(T_hot_plusX(i-1)+T_hot_minusX(i-1)-2*T_hot_minusY(i-1))+small_node_bolts*(T_hot_plusZ(i-1)+T_hot_minusZ(i-1)-2*T_hot_minusY(i-1)));
            Qdot_hot_minusY = Qdot_earth_IR_minusY + Qdot_solar_hot_minusY(i) + Qdot_albedo_hot_minusY(i) + large_node_power - Qdot_radiation_hot_minusY + Qdot_bolts_hot_minusY;

            % Calculate the -Y node temperature for each time step of the hot case
            T_hot_minusY(i) = T_hot_minusY(i-1) + (delta_t*Qdot_hot_minusY)/(specific_heat*large_node_mass);

            % Calculate the net heat transfer rate for the +Z node for each time step of the hot case
            Qdot_solar_hot_plusZ(i) = plusZ_VF*solar_flux_hot*small_node_area*plusZ_absorptivity;
            Qdot_radiation_hot_plusZ = small_node_area*sb*plusZ_emissivity*(((T_hot_plusZ(i-1)+273.15)^4)-((space_temp+273.15)^4));
            Qdot_bolts_hot_plusZ = bolt_conductance*large_node_bolts*(T_hot_plusX(i-1)+T_hot_minusX(i-1)+T_hot_plusY(i-1)+T_hot_minusY(i-1)-4*T_hot_plusZ(i-1));
            Qdot_hot_plusZ = Qdot_solar_hot_plusZ(i) + small_node_power - Qdot_radiation_hot_plusZ + Qdot_bolts_hot_plusZ;
            
            % Calculate the +Z node temperature for each time step of the hot case
            T_hot_plusZ(i) = T_hot_plusZ(i-1) + (delta_t*Qdot_hot_plusZ)/(specific_heat*small_node_mass);

            % Calculate the net heat transfer rate for the -Z node for each time step of the hot case
            Qdot_solar_hot_minusZ(i) = minusZ_VF*solar_flux_hot*small_node_area*minusZ_absorptivity;
            Qdot_albedo_hot_minusZ(i) = albedo_switch*VF_nadir*earth_albedo*solar_flux_hot*minusZ_absorptivity*small_node_area*cosd(beta_angle)*cos((2*pi*t(i))/orbital_period);
            Qdot_radiation_hot_minusZ = small_node_area*sb*minusZ_emissivity*(((T_hot_minusZ(i-1)+273.15)^4)-((space_temp+273.15)^4));
            Qdot_bolts_hot_minusZ = bolt_conductance*large_node_bolts*(T_hot_plusX(i-1)+T_hot_minusX(i-1)+T_hot_plusY(i-1)+T_hot_minusY(i-1)-4*T_hot_minusZ(i-1));
            Qdot_hot_minusZ = Qdot_solar_hot_minusZ(i) + Qdot_albedo_hot_minusZ(i) + small_node_power - Qdot_radiation_hot_minusZ + Qdot_bolts_hot_minusZ;

            % Calculate the -Z node temperature for each time step of the hot case
            T_hot_minusZ(i) = T_hot_minusZ(i-1) + (delta_t*Qdot_hot_minusZ)/(specific_heat*small_node_mass);

            % Calculate the total heat loads for each node
            total_heat_load_plusX_hot(i) = Qdot_earth_IR_plusX + Qdot_solar_hot_plusX(i) + Qdot_albedo_hot_plusX(i) + large_node_power;
            total_heat_load_minusX_hot(i) = Qdot_earth_IR_minusX + Qdot_albedo_hot_minusX(i) + large_node_power;
            total_heat_load_plusY_hot(i) = Qdot_earth_IR_plusY + Qdot_solar_hot_plusY(i) + Qdot_albedo_hot_plusY(i) + large_node_power;
            total_heat_load_minusY_hot(i) = Qdot_earth_IR_minusY + Qdot_solar_hot_minusY(i) + Qdot_albedo_hot_minusY(i) + large_node_power;
            total_heat_load_plusZ_hot(i) = Qdot_solar_hot_plusZ(i) + small_node_power;
            total_heat_load_minusZ_hot(i) = Qdot_earth_IR_minusZ + Qdot_solar_hot_minusZ(i) + Qdot_albedo_hot_minusZ(i) + small_node_power;
            %% Cold Case

            % Calculate the net heat transfer rate for the +X node for each time step of the cold case
            Qdot_solar_cold_plusX(i) = plusX_VF*solar_flux_cold*large_node_area*plusX_absorptivity;
            Qdot_albedo_cold_plusX(i) = albedo_switch*VF_walls*earth_albedo*solar_flux_cold*plusX_absorptivity*large_node_area*cosd(beta_angle)*cos((2*pi*t(i))/orbital_period);
            Qdot_radiation_cold_plusX = large_node_area*sb*plusX_emissivity*(((T_cold_plusX(i-1)+273.15)^4)-((space_temp+273.15)^4));
            Qdot_bolts_cold_plusX = bolt_conductance*(large_node_bolts*(T_cold_plusY(i-1)+T_cold_minusY(i-1)-2*T_cold_plusX(i-1))+small_node_bolts*(T_cold_plusZ(i-1)+T_cold_minusZ(i-1)-2*T_cold_plusX(i-1)));
            Qdot_cold_plusX = Qdot_earth_IR_plusX + Qdot_solar_cold_plusX(i) + Qdot_albedo_cold_plusX(i) + large_node_power - Qdot_radiation_cold_plusX + Qdot_bolts_cold_plusX;
            
            % Calculate the +X node temperature for each time step of the cold case
            T_cold_plusX(i) = T_cold_plusX(i-1) + (delta_t*Qdot_cold_plusX)/(specific_heat*large_node_mass);

            % Calculate the net heat transfer rate for the -X node for each time step of the cold case
            Qdot_albedo_cold_minusX(i) = albedo_switch*VF_walls*earth_albedo*solar_flux_cold*minusX_absorptivity*large_node_area*cosd(beta_angle)*cos((2*pi*t(i))/orbital_period);
            Qdot_radiation_cold_minusX = large_node_area*sb*minusX_emissivity*(((T_cold_minusX(i-1)+273.15)^4)-((space_temp+273.15)^4));
            Qdot_bolts_cold_minusX = bolt_conductance*(large_node_bolts*(T_cold_plusY(i-1)+T_cold_minusY(i-1)-2*T_cold_minusX(i-1))+small_node_bolts*(T_cold_plusZ(i-1)+T_cold_minusZ(i-1)-2*T_cold_minusX(i-1)));
            Qdot_cold_minusX = Qdot_earth_IR_minusX + Qdot_albedo_cold_minusX(i) + large_node_power - Qdot_radiation_cold_minusX + Qdot_bolts_cold_minusX;
            
            % Calculate the -X node temperature for each time step of the cold case
            T_cold_minusX(i) = T_cold_minusX(i-1) + (delta_t*Qdot_cold_minusX)/(specific_heat*large_node_mass);

            % Calculate the net heat transfer rate for the +Y node for each time step of the cold case
            Qdot_solar_cold_plusY(i) = plusY_VF*solar_flux_cold*large_node_area*plusY_absorptivity;
            Qdot_albedo_cold_plusY(i) = albedo_switch*VF_walls*earth_albedo*solar_flux_cold*plusY_absorptivity*large_node_area*cosd(beta_angle)*cos((2*pi*t(i))/orbital_period);
            Qdot_radiation_cold_plusY = large_node_area*sb*plusY_emissivity*(((T_cold_plusY(i-1)+273.15)^4)-((space_temp+273.15)^4));
            Qdot_bolts_cold_plusY = bolt_conductance*(large_node_bolts*(T_cold_plusX(i-1)+T_cold_minusX(i-1)-2*T_cold_plusY(i-1))+small_node_bolts*(T_cold_plusZ(i-1)+T_cold_minusZ(i-1)-2*T_cold_plusY(i-1)));
            Qdot_cold_plusY = Qdot_earth_IR_plusY + Qdot_solar_cold_plusY(i) + Qdot_albedo_cold_plusY(i) + large_node_power - Qdot_radiation_cold_plusY + Qdot_bolts_cold_plusY;
            
            % Calculate the +Y node temperature for each time step of the cold case
            T_cold_plusY(i) = T_cold_plusY(i-1) + (delta_t*Qdot_cold_plusY)/(specific_heat*large_node_mass);

            % Calculate the net heat transfer rate for the -Y node for each time step of the cold case
            Qdot_solar_cold_minusY(i) = minusY_VF*solar_flux_cold*large_node_area*minusY_absorptivity;
            Qdot_albedo_cold_minusY(i) = albedo_switch*VF_walls*earth_albedo*solar_flux_cold*minusY_absorptivity*large_node_area*cosd(beta_angle)*cos((2*pi*t(i))/orbital_period);
            Qdot_radiation_cold_minusY = large_node_area*sb*minusY_emissivity*(((T_cold_minusY(i-1)+273.15)^4)-((space_temp+273.15)^4));
            Qdot_bolts_cold_minusY = bolt_conductance*(large_node_bolts*(T_cold_plusX(i-1)+T_cold_minusX(i-1)-2*T_cold_minusY(i-1))+small_node_bolts*(T_cold_plusZ(i-1)+T_cold_minusZ(i-1)-2*T_cold_minusY(i-1)));
            Qdot_cold_minusY = Qdot_earth_IR_minusY + Qdot_solar_cold_minusY(i) + Qdot_albedo_cold_minusY(i) + large_node_power - Qdot_radiation_cold_minusY + Qdot_bolts_cold_minusY;

            % Calculate the -Y node temperature for each time step of the cold case
            T_cold_minusY(i) = T_cold_minusY(i-1) + (delta_t*Qdot_cold_minusY)/(specific_heat*large_node_mass);

            % Calculate the net heat transfer rate for the +Z node for each time step of the cold case
            Qdot_solar_cold_plusZ(i) = plusZ_VF*solar_flux_cold*small_node_area*plusZ_absorptivity;
            Qdot_radiation_cold_plusZ = small_node_area*sb*plusZ_emissivity*(((T_cold_plusZ(i-1)+273.15)^4)-((space_temp+273.15)^4));
            Qdot_bolts_cold_plusZ = bolt_conductance*large_node_bolts*(T_cold_plusX(i-1)+T_cold_minusX(i-1)+T_cold_plusY(i-1)+T_cold_minusY(i-1)-4*T_cold_plusZ(i-1));
            Qdot_cold_plusZ = Qdot_solar_cold_plusZ(i) + small_node_power - Qdot_radiation_cold_plusZ + Qdot_bolts_cold_plusZ;
            
            % Calculate the +Z node temperature for each time step of the cold case
            T_cold_plusZ(i) = T_cold_plusZ(i-1) + (delta_t*Qdot_cold_plusZ)/(specific_heat*small_node_mass);

            % Calculate the net heat transfer rate for the -Z node for each time step of the cold case
            Qdot_solar_cold_minusZ(i) = minusZ_VF*solar_flux_cold*small_node_area*minusZ_absorptivity;
            Qdot_albedo_cold_minusZ(i) = albedo_switch*VF_nadir*earth_albedo*solar_flux_cold*minusZ_absorptivity*small_node_area*cosd(beta_angle)*cos((2*pi*t(i))/orbital_period);
            Qdot_radiation_cold_minusZ = small_node_area*sb*minusZ_emissivity*(((T_cold_minusZ(i-1)+273.15)^4)-((space_temp+273.15)^4));
            Qdot_bolts_cold_minusZ = bolt_conductance*large_node_bolts*(T_cold_plusX(i-1)+T_cold_minusX(i-1)+T_cold_plusY(i-1)+T_cold_minusY(i-1)-4*T_cold_minusZ(i-1));
            Qdot_cold_minusZ = Qdot_solar_cold_minusZ(i) + Qdot_albedo_cold_minusZ(i) + small_node_power - Qdot_radiation_cold_minusZ + Qdot_bolts_cold_minusZ;

            % Calculate the -Z node temperature for each time step of the cold case
            T_cold_minusZ(i) = T_cold_minusZ(i-1) + (delta_t*Qdot_cold_minusZ)/(specific_heat*small_node_mass);

            % Calculate the total heat loads for each node
            total_heat_load_plusX_cold(i) = Qdot_earth_IR_plusX + Qdot_solar_cold_plusX(i) + Qdot_albedo_cold_plusX(i) + large_node_power;
            total_heat_load_minusX_cold(i) = Qdot_earth_IR_minusX + Qdot_albedo_cold_minusX(i) + large_node_power;
            total_heat_load_plusY_cold(i) = Qdot_earth_IR_plusY + Qdot_solar_cold_plusY(i) + Qdot_albedo_cold_plusY(i) + large_node_power;
            total_heat_load_minusY_cold(i) = Qdot_earth_IR_minusY + Qdot_solar_cold_minusY(i) + Qdot_albedo_cold_minusY(i) + large_node_power;
            total_heat_load_plusZ_cold(i) = Qdot_solar_cold_plusZ(i) + small_node_power;
            total_heat_load_minusZ_cold(i) = Qdot_earth_IR_minusZ + Qdot_solar_cold_minusZ(i) + Qdot_albedo_cold_minusZ(i) + small_node_power;
        else
            % If we accidentally step beyond our desired number of orbits, the final time step is removed from the time vector
            t(i) = [];
            break
        end
        
    end

end

%% Min and Max Nodal Temperatures

% Hot Case

T_hot_plusX_min = min(T_hot_plusX);
T_hot_plusX_max = max(T_hot_plusX);
T_hot_minusX_min = min(T_hot_minusX);
T_hot_minusX_max = max(T_hot_minusX);

T_hot_plusY_min = min(T_hot_plusY);
T_hot_plusY_max = max(T_hot_plusY);
T_hot_minusY_min = min(T_hot_minusY);
T_hot_minusY_max = max(T_hot_minusY);

T_hot_plusZ_min = min(T_hot_plusZ);
T_hot_plusZ_max = max(T_hot_plusZ);
T_hot_minusZ_min = min(T_hot_minusZ);
T_hot_minusZ_max = max(T_hot_minusZ);

% Cold Case

T_cold_plusX_min = min(T_cold_plusX);
T_cold_plusX_max = max(T_cold_plusX);
T_cold_minusX_min = min(T_cold_minusX);
T_cold_minusX_max = max(T_cold_minusX);

T_cold_plusY_min = min(T_cold_plusY);
T_cold_plusY_max = max(T_cold_plusY);
T_cold_minusY_min = min(T_cold_minusY);
T_cold_minusY_max = max(T_cold_minusY);

T_cold_plusZ_min = min(T_cold_plusZ);
T_cold_plusZ_max = max(T_cold_plusZ);
T_cold_minusZ_min = min(T_cold_minusZ);
T_cold_minusZ_max = max(T_cold_minusZ);

%% Plotting

% Plot the direct solar radiation heat loads over one orbit for the hot case
figure(1);
set(1,'WindowStyle','Docked');
grid on
title('Hot Case Direct Solar Radiation Heat Loads over a Single Orbit','fontsize', 15)
xlabel('Time from Deployment (s)', 'fontsize', 15);
ylabel('Heat Load (W)', 'fontsize',15);
hold on;
plot(t(1:orbital_period),Qdot_solar_hot_plusX(1:orbital_period));
hold on;
plot([0,orbital_period],[0,0]);
hold on;
plot(t(1:orbital_period),Qdot_solar_hot_plusY(1:orbital_period));
hold on;
plot(t(1:orbital_period),Qdot_solar_hot_minusY(1:orbital_period));
hold on;
plot(t(1:orbital_period),Qdot_solar_hot_plusZ(1:orbital_period));
hold on;
plot(t(1:orbital_period),Qdot_solar_hot_minusZ(1:orbital_period));
legend('+X Node','-X Node','+Y Node','-Y Node','+Z Node','-Z Node');

% Plot the reflected solar radiation heat loads over one orbit for the hot case
figure(2);
set(2,'WindowStyle','Docked');
grid on
title('Hot Case Reflected Solar Radiation Heat Loads over a Single Orbit','fontsize', 15)
xlabel('Time from Deployment (s)', 'fontsize', 15);
ylabel('Heat Load (W)', 'fontsize',15);
hold on;
plot(t(1:orbital_period),Qdot_albedo_hot_plusX(1:orbital_period));
hold on;
plot(t(1:orbital_period),Qdot_albedo_hot_minusX(1:orbital_period));
hold on;
plot(t(1:orbital_period),Qdot_albedo_hot_plusY(1:orbital_period));
hold on;
plot(t(1:orbital_period),Qdot_albedo_hot_minusY(1:orbital_period));
hold on;
plot([0,orbital_period],[0,0]);
hold on;
plot(t(1:orbital_period),Qdot_albedo_hot_minusZ(1:orbital_period));
legend('+X Node','-X Node','+Y Node','-Y Node','+Z Node','-Z Node','Location','southeast');

% Plot the Earth IR Radiation heat loads over one orbit
figure(3);
set(3,'WindowStyle','Docked');
grid on
title('Earth IR Radiation Heat Loads over a Single Orbit','fontsize', 15)
xlabel('Time from Deployment (s)', 'fontsize', 15);
ylabel('Heat Load (W)', 'fontsize',15);
hold on;
plot([0,orbital_period],[Qdot_earth_IR_plusX,Qdot_earth_IR_plusX]);
hold on;
plot([0,orbital_period],[Qdot_earth_IR_minusZ,Qdot_earth_IR_minusZ]);
hold on;
plot([0,orbital_period],[0,0]);
legend('±X and ±Y Nodes','-Z Node','+Z Node','Location','southeast');

% Plot the internally generated heat loads over one orbit
figure(4);
set(4,'WindowStyle','Docked');
grid on
title('Internally Generated Power over a Single Orbit','fontsize', 15)
xlabel('Time from Deployment (s)', 'fontsize', 15);
ylabel('Heat Load (W)', 'fontsize',15);
hold on;
plot([0,orbital_period],[large_node_power,large_node_power]);
hold on;
plot([0,orbital_period],[small_node_power,small_node_power]);
legend('±X and ±Y Nodes','±Z Nodes','Location','southeast');

% Plot the combined heat loads over one orbit for the hot case
figure(5);
set(5,'WindowStyle','Docked');
grid on
title('Hot Case Combined Heat Loads over a Single Orbit','fontsize', 15)
xlabel('Time from Deployment (s)', 'fontsize', 15);
ylabel('Heat Load (W)', 'fontsize',15);
hold on;
plot(t(1:orbital_period),total_heat_load_plusX_hot(1:orbital_period));
hold on;
plot(t(1:orbital_period),total_heat_load_minusX_hot(1:orbital_period));
hold on;
plot(t(1:orbital_period),total_heat_load_plusY_hot(1:orbital_period));
hold on;
plot(t(1:orbital_period),total_heat_load_minusY_hot(1:orbital_period));
hold on;
plot(t(1:orbital_period),total_heat_load_plusZ_hot(1:orbital_period));
hold on;
plot(t(1:orbital_period),total_heat_load_minusZ_hot(1:orbital_period));
legend('+X Node','-X Node','+Y Node','-Y Node','+Z Node','-Z Node');

% Plot the diirect solar radiation heat loads over one orbit for the cold case
figure(6);
set(6,'WindowStyle','Docked');
grid on
title('Cold Case Direct Solar Radiation Heat Loads over a Single Orbit','fontsize', 15)
xlabel('Time from Deployment (s)', 'fontsize', 15);
ylabel('Heat Load (W)', 'fontsize',15);
hold on;
plot(t(1:orbital_period),Qdot_solar_cold_plusX(1:orbital_period));
hold on;
plot([0,orbital_period],[0,0]);
hold on;
plot(t(1:orbital_period),Qdot_solar_cold_plusY(1:orbital_period));
hold on;
plot(t(1:orbital_period),Qdot_solar_cold_minusY(1:orbital_period));
hold on;
plot(t(1:orbital_period),Qdot_solar_cold_plusZ(1:orbital_period));
hold on;
plot(t(1:orbital_period),Qdot_solar_cold_minusZ(1:orbital_period));
legend('+X Node','-X Node','+Y Node','-Y Node','+Z Node','-Z Node');

% Plot the reflected solar radiation heat loads over one orbit for the cold case
figure(7);
set(7,'WindowStyle','Docked');
grid on
title('Cold Case Reflected Solar Radiation Heat Loads over a Single Orbit','fontsize', 15)
xlabel('Time from Deployment (s)', 'fontsize', 15);
ylabel('Heat Load (W)', 'fontsize',15);
hold on;
plot(t(1:orbital_period),Qdot_albedo_cold_plusX(1:orbital_period));
hold on;
plot(t(1:orbital_period),Qdot_albedo_cold_minusX(1:orbital_period));
hold on;
plot(t(1:orbital_period),Qdot_albedo_cold_plusY(1:orbital_period));
hold on;
plot(t(1:orbital_period),Qdot_albedo_cold_minusY(1:orbital_period));
hold on;
plot([0,orbital_period],[0,0]);
hold on;
plot(t(1:orbital_period),Qdot_albedo_cold_minusZ(1:orbital_period));
legend('+X Node','-X Node','+Y Node','-Y Node','+Z Node','-Z Node','Location','southeast');

% Plot the combined heat loads over one orbit for the cold case
figure(8);
set(8,'WindowStyle','Docked');
grid on
title('Cold Case Combined Heat Loads over a Single Orbit','fontsize', 15)
xlabel('Time from Deployment (s)', 'fontsize', 15);
ylabel('Heat Load (W)', 'fontsize',15);
hold on;
plot(t(1:orbital_period),total_heat_load_plusX_cold(1:orbital_period));
hold on;
plot(t(1:orbital_period),total_heat_load_minusX_cold(1:orbital_period));
hold on;
plot(t(1:orbital_period),total_heat_load_plusY_cold(1:orbital_period));
hold on;
plot(t(1:orbital_period),total_heat_load_minusY_cold(1:orbital_period));
hold on;
plot(t(1:orbital_period),total_heat_load_plusZ_cold(1:orbital_period));
hold on;
plot(t(1:orbital_period),total_heat_load_minusZ_cold(1:orbital_period));
legend('+X Node','-X Node','+Y Node','-Y Node','+Z Node','-Z Node');

% Plot all nodal temperatures for the hot case
figure(9);
set(9,'WindowStyle','Docked');
grid on
title('Hot Case Temperatures over Ten Orbits','fontsize', 15);
xlabel('Time from Deployment (s)', 'fontsize', 15);
ylabel('Temperature (°C)', 'fontsize',15);
hold on;
plot(t,T_hot_plusX);
hold on
plot(t,T_hot_minusX);
hold on;
plot(t,T_hot_plusY);
hold on
plot(t,T_hot_minusY);
hold on;
plot(t,T_hot_plusZ);
hold on
plot(t,T_hot_minusZ);
legend('+X Node','-X Node','+Y Node','-Y Node','+Z Node','-Z Node','Location','northeast')

% Plot all nodal temperatures for the cold case
figure(10);
set(10,'WindowStyle','Docked');
grid on
title('Cold Case Temperatures over Ten Orbits','fontsize', 15);
xlabel('Time from Deployment (s)', 'fontsize', 15);
ylabel('Temperature (°C)', 'fontsize',15);
hold on;
plot(t,T_cold_plusX);
hold on
plot(t,T_cold_minusX);
hold on;
plot(t,T_cold_plusY);
hold on
plot(t,T_cold_minusY);
hold on;
plot(t,T_cold_plusZ);
hold on
plot(t,T_cold_minusZ);
legend('+X Node','-X Node','+Y Node','-Y Node','+Z Node','-Z Node','Location','northeast')