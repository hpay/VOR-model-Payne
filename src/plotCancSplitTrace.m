%% Plot PC activity during VORD and cancellation
function fH_traces =   plotCancSplitTrace(I, K1, K2, K0)
% Plot inputs to PC during a step of cancellation stim

% Parameters

% split_target_retina = 1; % plot target and retinal slip separately or not?

learnStr = {'Low', 'Normal','High'};

Psource = [];


T_smooth = 0.025;   % (s) % 0.01
tafter = 0.5;      % (s) Time to go out for step
tbefore = -.1;
amp = 1;           % (deg/s) Scale for step

dt = K1.dt;
tt = tbefore:dt:tafter;

% Generate smoothed step
head_curr = double(tt>0)*amp;
head_curr = smooth(head_curr, round(T_smooth/dt));
head_curr = smooth(head_curr, round(T_smooth/dt));
target_curr = head_curr;
% K, I, head_vel, target_vel,  light, sine, pc_stim, N)
[E(:,2), P(:,2), Psource(:,:,:,2)] = modelClosedloop(K1, I, head_curr, target_curr, 1, 0);
[E(:,3), P(:,3), Psource(:,:,:,3)] = modelClosedloop(K2, I, head_curr, target_curr, 1, 0);
[E(:,1), P(:,1), Psource(:,:,:,1)] = modelClosedloop(K0, I, head_curr, target_curr, 1, 0);

%% Plot the traces for each input separately
colors = [.8 .8 .8; .5 .5 .5; .2 .2 .2];
xt = -.1;

scale_small = [1 1 1 1];
% if K1.PF == 1
%     scale_small = [50 10 5 10 ];
% else
%     scale_small = [10  0.5 10 10];
% end

labels_source= {'Head','Visual','Efference'};
% offset = 10;
    fH_traces = figure;

for jj =1:3 % Loop of Low Normal High gain
    
    text(xt,0, 'Head','Color','r','Vert','bottom'); hold on
    plot(tt, Psource(:,1,1,jj), 'Color', colors(jj,:) ); hold on;
    
    text(xt,-2, sprintf('Retinal slip x %g', scale_small(1)),'Color','r','Vert','bottom')
    plot(tt, -2+scale_small(1)*Psource(:,strcmp(labels_source,'Visual'),1,jj),'Color',colors(jj,:)); hold on;
    
    text(xt,-4, sprintf('Target x %g', scale_small(2)),'Color','r','Vert','bottom')
    plot(tt, -4+scale_small(2)*Psource(:,strcmp(labels_source,'Visual'),2,jj),'Color',colors(jj,:)); hold on;
    
    text(xt,-6.32, sprintf('Efference copy x %g',scale_small(3)),'Color','r','Vert','bottom')
    plot(tt, -6+scale_small(3)*Psource(:,strcmp(labels_source,'Efference'),1,jj),'Color',colors(jj,:)); hold on;
    
    text(xt,-10, 'Purkinje total','Color','r','Vert','bottom')
    Phat_total = sum(sum(Psource(:,:,:,jj),2),3);
    %     Phat_total = P(:,jj);
    plot(tt, -10+Phat_total,'Color',colors(jj,:),'LineWidth',2); hold on;
%     plot(tt,-10+P(:,jj),'--','Color',colors(jj,:),'LineWidth',4); % These should be the same
    
    text(xt,-12, sprintf('Eye total x %i', scale_small(4)),'Color','r','Vert','bottom')
    plot(tt, -12+E(:,jj)*scale_small(4),'Color',colors(jj,:),'LineWidth',2); hold on;
    
    title(sprintf('Purkinje cell during cancellation: Positive feedback = %g, gain = %s', K1.PF, learnStr{jj}))
    fixticks
end


