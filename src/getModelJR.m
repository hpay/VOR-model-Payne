function [Ehat, Phat] = getModelJR(K_or_F, I, head, target, light, sines)
% function [Ehat, Phat] = getModelJR(K_or_F, I, head, target, light, sines)
%
% Get model predictions 

% Make sure its ok to use impulse response version of simulation
if any(light)
    I.impulse_or_closedloop = 0; % Use closed loop version 
end

% Get impulse response if needed
if I.impulse_or_closedloop    
    if isfield(K_or_F,'Kb_eye')
        F = getImpulse(K_or_F, I, 0);
    else
        F = K_or_F;
    end
else
    K = K_or_F;
end

[Ehat, Phat] = deal(cell(size(head)));
for jj = 1:length(head)
    
    if I.impulse_or_closedloop
        
        % Run model using impulse response
        [Ehat{jj}, Phat{jj}] = modelImpulse(F, I, head{jj}, target{jj}, light(jj), sines(jj));
    else
        
        % Run model using closed loop
        [Ehat{jj}, Phat{jj}] = modelClosedloop(K, I, head{jj}, target{jj}, light(jj), sines(jj));
    end
    
end

if length(Ehat)==1
    Ehat = Ehat{1};
    Phat = Phat{1};
end