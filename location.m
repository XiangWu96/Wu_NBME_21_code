function [index] = location(currentY,skin_depth,bone_depth,stepSize)

% Index key
% 1 = skin
% 2 = bone 
% 3 = tissue

if skin_depth == 0 && bone_depth == 0 % No scalp or skull
    if currentY < 0
    index = 4; % air
    else
    index = 3;
    end

 
elseif skin_depth == 0 % Transcranial-stimulation
    if currentY < 0
        index = 4; % air
    elseif ceil(currentY/stepSize) <= ceil(bone_depth/stepSize)
        index = 2; % bone
    else
        index = 3; % tissue
    end


else % Trans-scalp stimulation
    if currentY < 0
        index = 4; % air
    elseif ceil(currentY/stepSize) <= ceil(skin_depth/stepSize)
    index = 1; % skin
    elseif ceil(currentY/stepSize) <= ceil((skin_depth+bone_depth)/stepSize)
        index = 2; % bone
    else
        index = 3; % tissue
    end

end

end

