function [Ehat, Phat, varargout] = getFreq(KF, I, freqs)
% function [Ehat, Phat, E_gain, E_phase, P_gain, P_phase] = getFreq(K/F, I, RL_data)

% Get impulse response if needed
if I.impulse_or_closedloop
    if isfield(KF,'Kb_eye')
        F = getImpulse(KF, I, 1);
    else
        F = KF;
    end
else
    K = KF;
end


dt = KF.dt;

% Generate the gain and phase data for eye and pc
[E_gain, E_phase, P_gain, P_phase] = deal(NaN(size(freqs)));

[Ehat, Phat] = deal(cell(size(freqs)));
for jj = 1:length(freqs)
    
    T_temp = ceil(2*freqs(jj))/freqs(jj); % Max time - should be a multiple of sine period
    tt = dt*(1:round(T_temp/dt))'-dt;
    
    head_curr = sin(2*pi*freqs(jj)*tt);
    if I.impulse_or_closedloop
        [Ehat{jj}, Phat{jj}] = modelImpulse(F, I, head_curr, zeros(size(head_curr)), 0, freqs(jj));
    else
        [Ehat{jj}, Phat{jj}] = modelClosedloop(K, I, head_curr, zeros(size(head_curr)), 0, freqs(jj));
    end
    
    if nargout >2
        [E_gain(jj), E_phase(jj)] = fitsine(tt, Ehat{jj}, freqs(jj));
        [P_gain(jj), P_phase(jj)] = fitsine(tt, Phat{jj}, freqs(jj));
    end
    
    Ehat{jj} = Ehat{jj}(1:round(2/dt));
    Phat{jj} = Phat{jj}(1:round(2/dt));
    
end

if nargout>2
    varargout = {E_gain, E_phase, P_gain, P_phase};
end

end