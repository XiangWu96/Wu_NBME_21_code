%% Monte Carlo code
% Please ensure that the script location.m is in the same directory as this code!

clear all
close all

% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
target_depths = 0.9; % In mm, depth from surface of cortex
targetRadius = 0.8; % In mm, e.g. 0.8 mm radius gives an injection volume of 2.5 µL
NumPackets = 1000; % Increase to average over a greater number of random trajectories
input_power = 8; % In mW/mm^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All depths and widths are in mm
skin_depth=0.6;
bone_depth=0.3;
tissue_depth=10;
tissue_width=25;
total_depth=skin_depth+bone_depth+tissue_depth; % Total depth of simulation domain

interface_position = [0 skin_depth skin_depth+bone_depth 0; skin_depth 0 skin_depth+bone_depth 0; skin_depth skin_depth+bone_depth 0 0; 0 0 0 0]; 
% Gives the y-position of the interface for refraction from row # (index_old) to
% column # (index_new). Index values: 1 = skin, 2 = bone, 3 = tissue

power = input_power/1E5; % in W/(0.1 mm^2), power incident on a pixel from a 1 W/cm^2 laser
sim_duration = 0.1; % In seconds, arbitary time (will be accounted for in heating code below)
energy = power*sim_duration; % Joules incident on pixel during entire sim duration

EnerPacket = energy/NumPackets; % Calculating the energy of each packet

mua = [0.0896 0.015 0.04052]; % Absorption coefficient, [skin bone brain], unit is mm^-1
mus = [3.79 16 5.90]; % Scattering coefficient, [skin bone brain], unit is mm^-1
g = [0.929 0.9 0.9]; % [skin bone brain], unitless
n = [1.38 1.56 1.3526 1.33]; % [skin bone brain air], refractive index

Cv = [3391 1794 3630]; % Specific heat, [skin bone brain], unit is J/(kg K)
density = [1109 1543 1046]; % Tissue density, [skin bone brain], unit is kg/m^3

Cv_minds = 4200; % unit is J/(kg K)
density_minds = 1000; % unit is kg/m^3
mua_MINDS = 10; % unit is mm^-1
conversion_efficiency = 0.71; % = 1 when the device uses 100% of the absorbed energy, between 0 and 1 otherwise

PDF_scalp = @(cos_theta) (1-g(1).^2)./(2.*(1+g(1).^2-2.*g(1).*cos_theta).^(3/2)); % Calculating the probability density function for scattering
PDF_skull = @(cos_theta) (1-g(2).^2)./(2.*(1+g(2).^2-2.*g(2).*cos_theta).^(3/2));
PDF_brain = @(cos_theta) (1-g(3).^2)./(2.*(1+g(3).^2-2.*g(3).*cos_theta).^(3/2));

int_resolution = 1E4; % 1E4 should work well, specifies the resolution of our CDF
count = 0;
CDF = zeros(int_resolution,3);
x_CDF = zeros(int_resolution,3);

for i = linspace(-1,1,int_resolution)
    count = count+1;
    
    CDF(count,1) = integral(PDF_scalp,-1,i);
    x_CDF(count,1) = i;
    
    CDF(count,2) = integral(PDF_skull,-1,i);
    x_CDF(count,2) = i;
    
    CDF(count,3) = integral(PDF_brain,-1,i);
    x_CDF(count,3) = i;
end

stepSize=0.1; % In mm, meshing the domain on a square grid

xPositions = linspace(-5,5,101); % For 5 mm laser radius (10 mm diameter) and 0.1 mm step size

width_tuning = 7;
width = round((width_tuning+tissue_width/2)/stepSize)-round((-width_tuning+tissue_width/2)/stepSize)+1; % Make the width_tuning values bigger to increase width of simulation

for q = 1:length(target_depths) % Iterating through different target depths (if desired)
    
targetY=skin_depth+bone_depth+target_depths(q); % In mm
sumTempMatrix = zeros(round(total_depth/stepSize),round(width)); % Generate a new matrix for each depth

% Start of MINDS finding code
MINDS_location = zeros(round(total_depth/stepSize),round(tissue_width/stepSize));
MINDS_X = 0; % Assuming that target is in the center of the simulation
MINDS_Y = targetY;

y_values = linspace(0,total_depth,1729); % Random number of values within, just trying to ensure we hit all points on the grid
x_values = linspace(-targetRadius-2,targetRadius+2,1729); % Arbitrary safe radius, same motive as above

for a = 1:length(y_values)
    for b = 1:length(x_values)
        if sqrt((y_values(a)-MINDS_Y)^2+(x_values(b)-MINDS_X)^2)<=targetRadius
            MINDS_location(ceil(y_values(a)/stepSize),ceil((x_values(b)+tissue_width/2)/stepSize)) = 1;
        end
    end
end

left = round((MINDS_X-width_tuning+tissue_width/2)/stepSize); % These are used
right = round((MINDS_X+width_tuning+tissue_width/2)/stepSize);

sumTempMINDS = MINDS_location(:,left:right);
% End of MINDS finding code

targetEnergy=0;

for m = 1:length(xPositions)

tempMatrix=zeros(round(total_depth/stepSize),round(tissue_width/stepSize));

targetX=xPositions(m); % mm

for i=1:NumPackets

    currentAngle=0; % Assuming normal incidence into the brain
    currentEnergy=EnerPacket;
    currentX=0; % mm
    currentY=0; % mm

    % Travel in the brain tissue
    flag = 0; % Flag goes to 1 when packet exits the simulation bounds, resets to 0 for each new packet
    while (currentY<total_depth) && (currentY>=0) && abs(currentX)<(tissue_width/2) && currentEnergy>(EnerPacket/100000000)
        index_start_scatter = location(currentY,skin_depth,bone_depth,stepSize); % Index where our packet begins next scattering event
        freeTravelDist=-log(rand)/mus(index_start_scatter); % mu_s is location dependent
        numSteps=floor(freeTravelDist/stepSize);
        remainder=freeTravelDist;
        for j=1:numSteps
            index_new = location((stepSize*cos(currentAngle)+currentY),skin_depth,bone_depth,stepSize); % Index of where our current step would take us, neglecting refraction
            index_old = location(currentY,skin_depth,bone_depth,stepSize); % Index of where we are now
            if index_new ~= index_old
                % Calculating the angle between ray and normal angle, the calculation is dependent on the quadrant of currentAngle
                if 0 <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < (pi/2)
                        diffFromNormalAngle = mod(currentAngle,2*pi);
                elseif (pi/2) <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < pi
                        diffFromNormalAngle = mod(pi-currentAngle,2*pi);
                elseif pi <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < (3*pi/2)
                        diffFromNormalAngle = mod(currentAngle-pi,2*pi);
                elseif (3*pi/2) <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < (2*pi)
                        diffFromNormalAngle = mod(2*pi-currentAngle,2*pi);
                end
                if n(index_new) < n(index_old) && diffFromNormalAngle > asin(n(index_new)/n(index_old)) % Condition for total internal reflection
                    if 0 <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < (pi/2)
                        newAngle = pi-diffFromNormalAngle;
                    elseif (pi/2) <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < pi
                        newAngle = diffFromNormalAngle;
                    elseif pi <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < (3*pi/2)
                        newAngle = 2*pi-diffFromNormalAngle;
                    elseif (3*pi/2) <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < (2*pi)
                        newAngle = pi+diffFromNormalAngle;
                    end
                else % Refraction occurs, not total internal reflection
                    if 0 <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < (pi/2)
                        newAngle = asin(n(index_old)/n(index_new)*sin(diffFromNormalAngle));
                    elseif (pi/2) <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < pi
                        newAngle = pi-asin(n(index_old)/n(index_new)*sin(diffFromNormalAngle));
                    elseif pi <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < (3*pi/2)
                        newAngle = pi+asin(n(index_old)/n(index_new)*sin(diffFromNormalAngle));
                    elseif (3*pi/2) <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < (2*pi)
                        newAngle = 2*pi-asin(n(index_old)/n(index_new)*sin(diffFromNormalAngle));
                    end
                end   
                s1 = abs(interface_position(index_old,index_new)-currentY)/cos(diffFromNormalAngle);
                s2 = stepSize-s1; % The two path lengths must add up to the total stepSize
                if s1>stepSize || s2>stepSize % Just a precautionary check
                    error('Step size error occurred')
                end
                x1 = s1*sin(currentAngle);
                x2 = s2*sin(newAngle);
                newX = currentX+x1+x2;
                y1 = interface_position(index_old,index_new)-currentY;
                y2 = s2*cos(newAngle);
                newY = currentY+y1+y2;
                currentAngle = newAngle;
            else % If we didn't actually cross an interface, it's a much simpler case
                newX = stepSize*sin(currentAngle)+currentX;
                newY = stepSize*cos(currentAngle)+currentY;
            end
            if newY<total_depth && newY>=1E-9 && abs(newX)<(tissue_width/2) % If our new coordinates are within the simulation boundaries, go there
                currentY=round(newY,10);
                currentX=round(newX,10);
                if sqrt((currentY-targetY)^2+(currentX-targetX)^2)<=targetRadius
                    tempMatrix(ceil(currentY/stepSize),ceil((currentX+tissue_width/2)/stepSize))=tempMatrix(ceil(currentY/stepSize),ceil((currentX+tissue_width/2)/stepSize))+conversion_efficiency.*currentEnergy*(1-exp(-mua_MINDS*stepSize))/(Cv_minds*density_minds*(stepSize*1E-3)^3);
                    targetEnergy=targetEnergy+conversion_efficiency.*currentEnergy*(1-exp(-mua_MINDS*stepSize));
                    currentEnergy=currentEnergy*exp(-mua_MINDS*stepSize);                    
                else
                    tempMatrix(ceil(currentY/stepSize),ceil((currentX+tissue_width/2)/stepSize))=tempMatrix(ceil(currentY/stepSize),ceil((currentX+tissue_width/2)/stepSize))+currentEnergy*(1-exp(-mua(location(currentY,skin_depth,bone_depth,stepSize))*stepSize))/(Cv(location(currentY,skin_depth,bone_depth,stepSize))*density(location(currentY,skin_depth,bone_depth,stepSize))*(stepSize*1E-3)^3);
                    currentEnergy=currentEnergy*exp(-mua(location(currentY,skin_depth,bone_depth,stepSize))*stepSize);
                end
                remainder=remainder-stepSize;
            else
                flag = 1; % If our new coordinates are not in the simulation boundaries, break out of the for loop, and turn on flag to break out of the while loop, restarting with the next packet
                break % Breaks the loop over # of steps if the packet has exited the boundaries
            end
        end 
        if flag == 1 % If we already exited the bounds without reaching the remainder, we don't want to also process the remainder step, break out of the while loop
            break
        end
        index_new = location((remainder*cos(currentAngle)+currentY),skin_depth,bone_depth,stepSize);
        index_old = location(currentY,skin_depth,bone_depth,stepSize);
            if index_new ~= index_old
                if 0 <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < (pi/2)
                        diffFromNormalAngle = mod(currentAngle,2*pi);
                elseif pi/2 <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < (pi)
                        diffFromNormalAngle = mod(pi-currentAngle,2*pi);
                elseif pi <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < (3*pi/2)
                        diffFromNormalAngle = mod(currentAngle-pi,2*pi);
                elseif 3*pi/2 <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < (2*pi)
                        diffFromNormalAngle = mod(2*pi-currentAngle,2*pi);
                end
                if n(index_new) < n(index_old) && diffFromNormalAngle > asin(n(index_new)/n(index_old)) % Condition for total internal reflection
                    if 0 <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < (pi/2)
                        newAngle = pi-diffFromNormalAngle;
                    elseif pi/2 <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < (pi)
                        newAngle = diffFromNormalAngle;
                    elseif pi <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < (3*pi/2)
                        newAngle = 2*pi-diffFromNormalAngle;
                    elseif 3*pi/2 <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < (2*pi)
                        newAngle = pi+diffFromNormalAngle;
                    end
                else % Normal refraction occurs
                    if 0 <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < (pi/2)
                        newAngle = asin(n(index_old)/n(index_new)*sin(diffFromNormalAngle));
                    elseif pi/2 <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < (pi)
                        newAngle = pi-asin(n(index_old)/n(index_new)*sin(diffFromNormalAngle));
                    elseif pi <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < (3*pi/2)
                        newAngle = pi+asin(n(index_old)/n(index_new)*sin(diffFromNormalAngle));
                    elseif 3*pi/2 <= mod(currentAngle,2*pi) && mod(currentAngle,2*pi) < (2*pi)
                        newAngle = 2*pi-asin(n(index_old)/n(index_new)*sin(diffFromNormalAngle));
                    end
                end   
                s1 = abs(interface_position(index_old,index_new)-currentY)/cos(diffFromNormalAngle);
                s2 = remainder-s1;
                if s1>remainder || s2>remainder || remainder>stepSize
                    error('Step size error occurred')
                end
                x1 = s1*sin(currentAngle);
                x2 = s2*sin(newAngle);
                newX = currentX+x1+x2;
                y1 = interface_position(index_old,index_new)-currentY;
                y2 = s2*cos(newAngle);
                newY = currentY+y1+y2;
                currentAngle = newAngle;
            else
                newX = remainder*sin(currentAngle)+currentX;
                newY = remainder*cos(currentAngle)+currentY;
            end
        if newY<total_depth && newY>=1E-9 && abs(newX)<(tissue_width/2)
            currentY=round(newY,10);
            currentX=round(newX,10);
            currentAngle=currentAngle+(2*randi([0 1])-1)*acos(interp1(CDF(:,location(currentY,skin_depth,bone_depth,stepSize)),x_CDF(:,location(currentY,skin_depth,bone_depth,stepSize)),rand)); % Calculating the new angle for the next scattering event
            if sqrt((currentY-targetY)^2+(currentX-targetX)^2)<=targetRadius
                tempMatrix(ceil(currentY/stepSize),ceil((currentX+tissue_width/2)/stepSize))=tempMatrix(ceil(currentY/stepSize),ceil((currentX+tissue_width/2)/stepSize))+conversion_efficiency.*currentEnergy*(1-exp(-mua_MINDS*remainder))/(Cv_minds*density_minds*(stepSize*1E-3)^3);
                targetEnergy=targetEnergy+conversion_efficiency.*currentEnergy*(1-exp(-mua_MINDS*remainder));
                currentEnergy=currentEnergy*exp(-mua_MINDS*remainder); 
            else
                tempMatrix(ceil(currentY/stepSize),ceil((currentX+tissue_width/2)/stepSize))=tempMatrix(ceil(currentY/stepSize),ceil((currentX+tissue_width/2)/stepSize))+currentEnergy*(1-exp(-mua(location(currentY,skin_depth,bone_depth,stepSize))*remainder))/(Cv(location(currentY,skin_depth,bone_depth,stepSize))*density(location(currentY,skin_depth,bone_depth,stepSize))*(stepSize*1E-3)^3);
                currentEnergy=currentEnergy*exp(-mua(location(currentY,skin_depth,bone_depth,stepSize))*remainder);
            end
        else
            break
        end
    end
end

left = round((targetX-width_tuning+tissue_width/2)/stepSize);
right = round((targetX+width_tuning+tissue_width/2)/stepSize);

sumTempMatrix = sumTempMatrix + tempMatrix(:,left:right);

percent_complete = num2str(round(100*m/length(xPositions),1));
text_output = strcat(percent_complete,'% complete');
disp(text_output)

end

end

% SumTempMatrix gives the abso

%% Heat diffusion code

time_step = 0.0001; % In s
video_timestep = 0.01; % In s, how much time resolution is needed for the final output matrix

times = [0 2 3]; % [Time when heating starts   Time when heating stops   Time when simulation stops]

waiting_Time = 0;

total_Time = ceil(max(times)); % In s, including waiting_Time

column = [36.4992163378730;36.5199045363406;36.5399632606677;36.5594180946987;36.5782938575621;36.5966146353056;36.6108152794713;36.6219976641139;36.6325634640147;36.6431148223471;36.6539939385512;36.6645522427838;36.6747991308846;36.6847437252471;36.6943948829289;36.7037612035252;36.7128510368112;36.7216724901611;36.7302334357497;36.7385415175441;36.7466041580897;36.7544285650989;36.7620217378466;36.7693904733794;36.7765413725431;36.7834808458348;36.7902151190843;36.7967502389692;36.8030920783708;36.8092463415728;36.8152185693099;36.8210141436690;36.8266382928487;36.8320960957806;36.8373924866169;36.8425322590877;36.8475200707334;36.8523604470138;36.8570577853000;36.8616163587504;36.8660403200761;36.8703337051983;36.8745004368006;36.8785443277801;36.8824690846003;36.8862783105487;36.8899755089012;36.8935640859983;36.8970473542332;36.9004285349566;36.9037107612999;36.9068970809192;36.9099904586634;36.9129937791675;36.9159098493748;36.9187414009894;36.9214910928616;36.9241615133083;36.9267551823701;36.9292745540078;36.9317220182398;36.9340999032228;36.9364104772768;36.9386559508569;36.9408384784743;36.9429601605664;36.9450230453196;36.9470291304454;36.9489803649113;36.9508786506290;36.9527258441002;36.9545237580225;36.9562741628565;36.9579787883549;36.9596393250560;36.9612574257424;36.9628347068662;36.9643727499422;36.9658731029097;36.9673372814654;36.9687667703670;36.9701630247104;36.9715274711796;36.9728615092728;36.9741665125031;36.9754438295779;36.9766947855552;36.9779206829795;36.9791228029980;36.9803024064580;36.9814607349861;36.9825990120510;36.9837184440100;36.9848202211407;36.9859055186582;36.9869754977199;36.9880313064168;36.9890740807547;36.9901049456234;36.9911250157575;36.9921353966871;36.9931371856815;36.9941314726851;36.9951193412465;36.9961018694432;36.9970801308004;36.9980551952060;36.9990281298233;37]; % Using the cooling parameter, this is the temperature distribution that results during steady state conditions, which we use as the starting point for our thermal simulations
cooling_parameter = 0.000021; % Empirically determined cooling parameter that gives a temperature distribution throughout the scalp / skull / brain that matches our experimental observations in awake mice

frames = floor(total_Time/video_timestep);

sumTempMatrix_normalized = sumTempMatrix.*time_step./sim_duration;

zero_matrix = zeros(size(sumTempMatrix_normalized,1),1);
sumTempMatrix_resized = [zero_matrix, sumTempMatrix_normalized, zero_matrix];

[depth_pixel, width_pixel] = size(sumTempMatrix_resized);

number_maps = round(total_Time./time_step); % Total number of maps we need to generate
Tem_old = zeros(depth_pixel,width_pixel); % Setting up matrices to be filled in below
Tem_new = zeros(depth_pixel,width_pixel);
Tem_final = zeros(depth_pixel,width_pixel,frames+1);
on_or_off = zeros(frames+1,1);
Tem_old(1:end) = repelem(column,1,width_pixel); % Initial condition for tissue

heating_step_heating = zeros(depth_pixel,width_pixel);
heating_step_heating(1:end,:) = sumTempMatrix_resized; % Initial condition

heating_step_cooling = zeros(depth_pixel, width_pixel); % No temperature increase is applied if the laser is off

zero_MINDS = zeros(size(sumTempMINDS,1),1);
sumTempMINDS = [zero_MINDS, sumTempMINDS, zero_MINDS];
binary_Tem_MINDS = sumTempMINDS;

Cv = [3391 1794 3630]; % [skin bone brain], unit is J/(kg K)
density = [1109 1543 1046]; % [skin bone brain], unit is kg/m^3
k = [0.37 0.32 0.51]; % Thermal conductivity, [skin bone brain], unit is W/m/K
alpha = k./Cv./density; % Thermal diffusivity, [skin bone brain], unit is m^2/s

f_blood = 1.12E-5; % In m^3/(kg s)
Cv_blood = 3617; % In J/(kg K)
density_blood = 1050; % In kg/m^3
T_blood = 37; % In C

% Assuming thermal diffusivity is constant in horizontal direction (same in
% the same tissue), and only changes along depth (different in different
% tissues). Thus only need to make it a function of depth

air_depth = 0;
depth_pixel_with_air = depth_pixel;

alpha_depth(round(air_depth./stepSize)+1:round((air_depth+skin_depth)./stepSize)) = alpha(1);
alpha_depth(round((air_depth+skin_depth)./stepSize)+1:round((air_depth+skin_depth+bone_depth)./stepSize)) = alpha(2);
alpha_depth(round((air_depth+skin_depth+bone_depth)./stepSize)+1:depth_pixel_with_air) = alpha(3);

Cv_depth(round(air_depth./stepSize)+1:round((air_depth+skin_depth)./stepSize)) = Cv(1);
Cv_depth(round((air_depth+skin_depth)./stepSize)+1:round((air_depth+skin_depth+bone_depth)./stepSize)) = Cv(2);
Cv_depth(round((air_depth+skin_depth+bone_depth)./stepSize)+1:depth_pixel_with_air) = Cv(3);

binary_depth(round(air_depth./stepSize)+1:round((air_depth+skin_depth)./stepSize)) = 1;
binary_depth(round((air_depth+skin_depth)./stepSize)+1:round((air_depth+skin_depth+bone_depth)./stepSize)) = 1;
binary_depth(round((air_depth+skin_depth+bone_depth)./stepSize)+1:depth_pixel_with_air) = 1;

count = 1;
old_heating = 0;
flag = 0;

for kk=2:1:number_maps
    
    kk/number_maps

    if kk*time_step < waiting_Time
        new_heating = 0;
        heating_step = heating_step_cooling;
        
    elseif kk*time_step >= waiting_Time && kk*time_step >= times(1) && kk*time_step <= times(2)
        new_heating = 1;
        heating_step = heating_step_heating;
        
    else
        new_heating = 0;
        heating_step = heating_step_cooling;
        
    end
    if kk ~= 2
    Tem_old = Tem_new;
    end
    for ii=1:1:depth_pixel_with_air-1
        for jj=1:1:width_pixel
            % Calcuate the temperature change due to vertical heat diffusion
            if ii==1 % At the first point
                delta_T_x = (time_step./(stepSize/1000).^2).*(1/2)*(alpha_depth(ii+1)+alpha_depth(ii))*(Tem_old(ii+1,jj)-Tem_old(ii,jj));
            elseif ii==depth_pixel_with_air % At the last point
                delta_T_x = (time_step./(stepSize/1000).^2).*(1/2)*(alpha_depth(ii-1)+alpha_depth(ii))*(Tem_old(ii-1,jj)-Tem_old(ii,jj));
            else % In the middle, basically the sum of the above two situations
                delta_T_x = (time_step./(stepSize/1000).^2).*(1/2)*(alpha_depth(ii+1)+alpha_depth(ii))*(Tem_old(ii+1,jj)-Tem_old(ii,jj)) + (time_step./(stepSize/1000).^2).*(1/2)*(alpha_depth(ii-1)+alpha_depth(ii))*(Tem_old(ii-1,jj)-Tem_old(ii,jj));
            end
            
            % Calculate the temperature change due to horizontal heat diffusion
            if jj==1
                delta_T_y = (time_step./(stepSize/1000).^2).*alpha_depth(ii).*(Tem_old(ii,jj+1)-Tem_old(ii,jj));
            elseif jj==width_pixel
                delta_T_y = (time_step./(stepSize/1000).^2).*alpha_depth(ii).*(Tem_old(ii,jj-1)-Tem_old(ii,jj));
            else
                delta_T_y = (time_step./(stepSize/1000).^2).*alpha_depth(ii).*(Tem_old(ii,jj+1)+Tem_old(ii,jj-1)-2.*Tem_old(ii,jj));
            end
            
            Tem_new(ii,jj) = Tem_old(ii,jj) + delta_T_x + delta_T_y + heating_step(ii,jj) + binary_depth(ii)*time_step*density_blood*Cv_blood*f_blood*(T_blood-Tem_old(ii,jj))/Cv_depth(ii); % Pennes bioheat equation
            
        end
    end
    
    Tem_new(end,:) = 37;
    Tem_new(:,1) = column;
    Tem_new(:,end) = column;
    
    Tem_new(1,:) = Tem_new(1,:)-cooling_parameter;
    
    if kk == 2
        Tem_final(:,:,1) = Tem_old; % Keeping the first frame
        on_or_off(1) = new_heating;
    elseif mod(kk,number_maps/frames) == 0
        count = count+1;
        Tem_final(:,:,count) = Tem_new; % Keeping a frame every 0.01 s
        on_or_off(count) = new_heating;
    end
    
    old_heating = new_heating;
    
end

% The Tem_final matrix can be used for plotting temperature vs. time as desired, and desired
% quantities can be extracted by changing the sample lines below:
cortex_max_heating = max(Tem_final(round((skin_depth+bone_depth)/stepSize+1),:,:));
MINDS_max_heating = max(max(binary_Tem_MINDS.*Tem_final));
