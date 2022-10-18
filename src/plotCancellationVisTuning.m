function hf_canc = plotCancellationVisTuning(I,  K1, K2, K0, K1_vis, K2_vis, K0_vis)     
hf_canc = figure;

    sine_on = 0.5; % (Hz)
    light_on = 1;
    
    dt = K1.dt;
    tt = (1:2/dt)'*dt-dt;
    head_curr = sin(2*pi*.5*tt);
    target_curr = head_curr;
    
%         if I.impulse_or_closedloop
%             [Ehat1, Phat1] = modelImpulse(getImpulse(K1, I), head_curr, target_curr,  light_on, sine_on);
%             [Ehat0, Phat0] = modelImpulse(getImpulse(K0, I), head_curr, target_curr,  light_on, sine_on);
%             [Ehat2, Phat2] = modelImpulse(getImpulse(K2, I), head_curr, target_curr,  light_on, sine_on);
%         else
    [Ehat1, Phat1] = modelClosedloop(K1, I,  head_curr, target_curr, light_on, sine_on);
    [Ehat0, Phat0] = modelClosedloop(K0, I, head_curr, target_curr, light_on, sine_on);
    [Ehat2, Phat2] = modelClosedloop(K2, I, head_curr, target_curr, light_on, sine_on);
%         end
    subplot(221);
    h = plot(tt, Ehat0,'Color',I.learn_color(1,:)); hold on
    h(2) = plot(tt, Ehat1,'Color',I.learn_color(2,:)); hold on
    h(3) = plot(tt, Ehat2,'Color',I.learn_color(3,:)); hold on
   
    ylim([-1 1]); ylabel('Eye gain'); xlabel('Time (s)')
    legend(h, {'Low','Normal','High'})
    title('Cancellation: before tuning')
    axis square; fixticks
    
    subplot(223);
    h = plot(tt, Phat0,'Color',I.learn_color(1,:)); hold on
    h(2) = plot(tt, Phat1,'Color',I.learn_color(2,:)); hold on
    h(3) = plot(tt, Phat2,'Color',I.learn_color(3,:)); hold on
    
    ylim(1.5*[-1 1]); ylabel('Eye gain'); xlabel('Time (s)')
    legend(h, {'Low','Normal','High'})
    title('Cancellation: before tuning: Purkinje cell')
    axis square; fixticks
    
    % Run cancellation at 0.5 Hz after learning
%         if I.impulse_or_closedloop
%             [Ehat1, Phat1] = modelImpulse(getImpulse(K1_vis, I), head_curr, target_curr,  light_on, sine_on);
%             [Ehat0, Phat0] = modelImpulse(getImpulse(K0_vis, I), head_curr, target_curr,  light_on, sine_on);
%             [Ehat2, Phat2] = modelImpulse(getImpulse(K2_vis, I), head_curr, target_curr,  light_on, sine_on);
%         else
    [Ehat1, Phat1] = modelClosedloop(K1_vis, I, head_curr, target_curr, light_on, sine_on);
    [Ehat0, Phat0] = modelClosedloop(K0_vis, I, head_curr, target_curr, light_on, sine_on);
    [Ehat2, Phat2] = modelClosedloop(K2_vis, I, head_curr, target_curr, light_on, sine_on);
%         end
    subplot(222);
    
    h = plot(tt, Ehat0,'Color',I.learn_color(1,:)); hold on
    h(2) = plot(tt, Ehat1,'Color',I.learn_color(2,:));
    h(3) = plot(tt, Ehat2,'Color',I.learn_color(3,:));
    
    ylim([-1 1]); xlabel('Time (s)')
    title('Cancellation: after tuning')
    axis square; fixticks
    
    
    subplot(224);
    h = plot(tt, Phat0,'Color',I.learn_color(1,:)); hold on
    h(2) = plot(tt, Phat1,'Color',I.learn_color(2,:)); hold on
    h(3) = plot(tt, Phat2,'Color',I.learn_color(3,:)); hold on
    
    ylim(1.5*[-1 1]); xlabel('Time (s)')
    title('Cancellation: after tuning - Purkinje cell')
    axis square;
    fprintf('Eye canc: x0: %g   x1: %g   x2: %g\n', max(Ehat0), max(Ehat1), max(Ehat2));
    