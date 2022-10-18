function  plotHeadBreak(I, K, hs)



% Generate smoothed steps
T_smooth = 0.025;   % (s)
tbefore = -.1;
tbreak = 0.5;       % (s) time at which head stops moving
tafter = 1;         % (s) Time to go out for step
amp = 1;            % (deg/s) Scale for step

dt = K.dt;
tt = tbefore:dt:tafter;

% Generate smoothed step
target_curr = double(tt>0)*amp;
target_curr = smooth(target_curr, round(T_smooth/dt));
target_curr = smooth(target_curr, round(T_smooth/dt));


head_curr = double(tt>0 & tt<tbreak)*amp;
head_curr = smooth(head_curr, round(T_smooth/dt));
head_curr = smooth(head_curr, round(T_smooth/dt));


% Run model
% K, I, head_vel, target_vel,  light, sine, pc_stim, N)
[E, P, Psource] = modelClosedloop(K, I, head_curr, target_curr, 1, 0);
G = E + head_curr;

% Plot results

drawnow
if exist('hs','var')
    axes(hs(1)); drawnow;
else
    hf = figure;
end
plot(tt, cumsum(head_curr)*dt,'k-','LineWidth',2); hold on;
plot(tt, cumsum(target_curr)*dt,'--','Color',.6*[1 1 1],'LineWidth',2); hold on;
plot(tt, -.5+cumsum(G)*dt,'k','LineWidth',1.5)

reterror = cumsum(target_curr)*dt - cumsum(G)*dt;
plot(tt, -1+ reterror,'Color',[.5 0 .7],'LineWidth',1.5)

title(sprintf('Feedback = %g',K.PF));
xlim(tt([1 end]))
fixticks
ylabel('Position (deg)'); drawnow

legend('Head','Target','Gaze','Retinal error','Location','eastoutside')
legend('boxoff');  drawnow



%% Plot the source of PC activity

if exist('hs','var')
    axes(hs(2)); drawnow;
else
    hf(2) = figure;
end
drawnow
plot(tt,head_curr,'k-','LineWidth',2); hold on;
plot(tt,target_curr,'--','Color',.6*[1 1 1],'LineWidth',2); hold on;

plot(tt, -1+G,'Color',0*[1 1 1],'LineWidth',1.5)
plot(tt, -2+P,'Color',.5*[1 1 1],'LineWidth',1.5)

plot(tt, -3 + Psource(:,1,1),'b');
plot(tt, -3 + sum(Psource(:,2,:),3),'g');
plot(tt, -3 + sum(Psource(:,3,:),3),'r');
xlim(tt([1 end]))
ylabel('Deg/s or sp/s')

fixticks
ylim([-4.5 1.5])
xlabel('Time (s)')
title(sprintf('Feedback = %g',K.PF)); drawnow

legend('Head','Target','Gaze','P total','P from head','P from vision','P from efference','Location','eastoutside')
legend('boxoff'); 
drawnow
