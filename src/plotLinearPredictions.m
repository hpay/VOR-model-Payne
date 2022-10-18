function [err_eye, err_pc, Ehat_linear, Phat_linear]  = plotLinearPredictions(S_eye, S_pc, R_eye, R_pc, Kb_eye, Kb_pc)

Ehat_linear = S_eye*Kb_eye;
Phat_linear = S_pc*Kb_pc;

err_eye = sqrt(mean((Ehat_linear- R_eye).^2));
err_pc = sqrt(mean((Phat_linear- R_pc).^2));

tt_all = (0:length(Ehat_linear)-1);
figure;  subplot(211);
plot([tt_all(1) tt_all(end)],[0 0],'k--','Color',.7*[1 1 1],'Clipping','off'); hold on;
plot(tt_all, R_eye, 'Color',.2*[1 1 1 1],'Clipping','off','LineWidth',2); hold on;
plot(tt_all, Ehat_linear,'-k','Clipping','off');
ylabel('Eye (deg/s)');   set(gca,'XColor','w')
box off

subplot(212)
plot([tt_all(1) tt_all(end)],[0 0],'k--','Color',.7*[1 1 1]); hold on;
plot(tt_all, R_pc,'Color',[1 .7 .7 1],'Clipping','off','LineWidth',2); hold on;
plot(tt_all,Phat_linear,'r-','Clipping','off')
xlabel('Time (s)');    ylabel('PC (sp/s)'); box off

linkaxes;   ylim(30*[-1 1]);   xlim([0 tt_all(end)]);   
