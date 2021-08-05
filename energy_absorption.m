%% Monte Carlo code
% Please ensure that the script location.m is in the same directory as this code!

clear all
close all

% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
target_depth = 5; % In mm, depth from surface of cortex
targetRadius = 0.1; % In mm
NumPackets = 1000;
input_power = 10; % In mW/mm^2
wavelength = 1064; % Options are 635 nm, 980 nm, and 1064 nm
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

power = input_power/1E5; % in W/(0.1 mm^2), power incident on our pixel from a 1 W/cm^2 laser
sim_duration = 0.1; % In seconds, arbitary time over which cooling is negligible
energy = power*sim_duration; % Joules incident on pixel during entire sim duration

EnerPacket = energy/NumPackets;

if wavelength == 635
mua = [0.0953 0.059 0.10783]; % Absorption coefficient, [skin bone brain], unit is mm^-1
mus = [4.725 33.5 8.9965]; % Scattering coefficient, [skin bone brain], unit is mm^-1
g = [0.8 0.9 0.898]; % [skin bone brain], unitless
n = [1.38 1.56 1.35 1.33]; % [skin bone brain air], refractive index

elseif wavelength == 980
mua = [0.101 0.022 0.07205]; % % Absorption coefficient, [skin bone brain], unit is mm^-1
mus = [4.08 18 6.331179282]; % Scattering coefficient, [skin bone brain], unit is mm^-1
g = [0.904 0.9 0.9]; % [skin bone brain], unitless
n = [1.38 1.56 1.35 1.33]; % [skin bone brain air], refractive index
   
elseif wavelength == 1064
mua = [0.0896 0.015 0.04052]; % Absorption coefficient, [skin bone brain], unit is mm^-1
mus = [3.79 16 5.90]; % Scattering coefficient, [skin bone brain], unit is mm^-1
g = [0.929 0.9 0.9]; % [skin bone brain], unitless
n = [1.38 1.56 1.3526 1.33]; % [skin bone brain air], refractive index
end

Cv = [3391 1794 3630]; % [skin bone brain], unit is J/(kg K)
density = [1109 1543 1046]; % [skin bone brain], unit is kg/m^3

Cv_minds = 4200; % unit is J/(kg K)
density_minds = 1000; % unit is kg/m^3
mua_MINDS = 10; % unit is mm^-1
conversion_efficiency = 0.71; % = 1 when the device uses 100% of the absorbed energy, between 0 and 1 otherwise

PDF_scalp = @(cos_theta) (1-g(1).^2)./(2.*(1+g(1).^2-2.*g(1).*cos_theta).^(3/2));
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

stepSize=0.01; % In mm

xPositions = linspace(-5,5,1001); % For 5 mm laser radius (10 mm diameter) and 0.1 mm step size

width_tuning = 7;
width = round((width_tuning+tissue_width/2)/stepSize)-round((-width_tuning+tissue_width/2)/stepSize)+1; % Kind of arbitrary, playing around with width of produced image. Make the width_tuning values bigger to increase width of simulation

for q = 1:length(target_depth) % Iterating through different target depths (if desired)
    
targetY=skin_depth+bone_depth+target_depth(q); % In mm
sumEnerMatrix = zeros(round(total_depth/stepSize),round(width)); % Generate a new matrix for each depth

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

enerMatrix=zeros(round(total_depth/stepSize),round(tissue_width/stepSize));

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
                    enerMatrix(ceil(currentY/stepSize),ceil((currentX+tissue_width/2)/stepSize))=enerMatrix(ceil(currentY/stepSize),ceil((currentX+tissue_width/2)/stepSize))+currentEnergy*(1-exp(-mua_MINDS*stepSize));
                    targetEnergy=targetEnergy+conversion_efficiency.*currentEnergy*(1-exp(-mua_MINDS*stepSize));
                    currentEnergy=currentEnergy*exp(-mua_MINDS*stepSize);                    
                else
                    enerMatrix(ceil(currentY/stepSize),ceil((currentX+tissue_width/2)/stepSize))=enerMatrix(ceil(currentY/stepSize),ceil((currentX+tissue_width/2)/stepSize))+currentEnergy*(1-exp(-mua(location(currentY,skin_depth,bone_depth,stepSize))*stepSize));
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
                enerMatrix(ceil(currentY/stepSize),ceil((currentX+tissue_width/2)/stepSize))=enerMatrix(ceil(currentY/stepSize),ceil((currentX+tissue_width/2)/stepSize))+currentEnergy*(1-exp(-mua_MINDS*remainder));
                targetEnergy=targetEnergy+conversion_efficiency.*currentEnergy*(1-exp(-mua_MINDS*remainder));
                currentEnergy=currentEnergy*exp(-mua_MINDS*remainder); 
            else
                enerMatrix(ceil(currentY/stepSize),ceil((currentX+tissue_width/2)/stepSize))=enerMatrix(ceil(currentY/stepSize),ceil((currentX+tissue_width/2)/stepSize))+currentEnergy*(1-exp(-mua(location(currentY,skin_depth,bone_depth,stepSize))*remainder));
                currentEnergy=currentEnergy*exp(-mua(location(currentY,skin_depth,bone_depth,stepSize))*remainder);
            end
        else
            break
        end
    end
end

left = round((targetX-width_tuning+tissue_width/2)/stepSize);
right = round((targetX+width_tuning+tissue_width/2)/stepSize);

sumEnerMatrix = sumEnerMatrix + enerMatrix(:,left:right);

percent_complete = num2str(round(100*m/length(xPositions),1));
text_output = strcat(percent_complete,'% complete');
disp(text_output)

end

end

% sumEnerMatrix matrix gives the absorbed energy distribution in tissue