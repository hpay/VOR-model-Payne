function [hf, hp] = plotFiltersLearn(K1, K2, color1, color2)
% Only plots the vestibular weights and their changes
    
% Plot resulting filter
hf = gcf;
dt = K1.dt;
ntf = length(K1.PH);
tf = (0:ntf-1)*dt;

% k_PH
hs = subplot(2, 2, 1);
hp = plot(tf, K1.PH,'Color',color1);     hold on
hp(2) = plot(tf, K2.PH,'Color',color2);     hold on
ylabel('weight')
title('k_{PH}')

% delta k_PH
hs(1,2) = subplot(2, 2, 2);
plot(tf, K2.PH - K1.PH,'Color',color2);   hold on;
title('\Deltak_{PH}')

% k_EH
hs(2,1) = subplot(2, 2, 3);
plot(tf, K1.EH,'Color',color1);  hold on;
plot(tf, K2.EH,'Color',color2);  hold on;
ylabel('weight')
title('k_{EH}');
xlabel('Time (s)')

% delta k_EH
hs(2,2) = subplot(2, 2, 4);
plot(tf, K2.EH - K1.EH,'Color',color2);  hold on;
title('\Delta k_{EH}');
xlabel('Time (s)')

drawnow

        
        
   