
function post_processing_and_plots(step)

global DATA

% clean extra allocated memory
DATA.PCA(step+1:end)= [];
DATA.PCA_MJ(step+1:end)= [];
DATA.PCA_MJ_outliers(step+1:end)= [];
DATA.PCAt(step+1:end)= [];
DATA.realPCA(step+1:end)= [];
DATA.calcPCA(step+1:end)= [];
DATA.calcPCA_MJ(step+1:end)= [];
DATA.errorXX(step+1:end,:)= [];
DATA.stdXX(step+1:end,:)= [];


% plots - P(CA) VS time
figure; hold on; grid on;
title('Instant P(CA)');
xlabel('Time'); ylabel('P(CA)');
plot(DATA.PCA,'-b');
plot(DATA.PCA_MJ,'--g');
% plot(DATA.PCA_MJ_outliers,'-*k');
% plot(DATA.PCAt,'or');
% idx= find(DATA.PCAt == 0);
% for i= 1:length(idx)
%     line([idx(i),idx(i)],[0,1],'color','red');
% end
axis([0,step,0,1]);


% % plots - Error Vs Covariance - X & Y
% figure; hold on; grid on
% plot(errorXX(:,1),'b-');
% plot(errorXX(:,2),'r-');
% plot(stdXX(:,1),'b--','linewidth',5);
% plot(stdXX(:,2),'r--','linewidth',5);
% for i= 1:length(idx)
%     line([idx(i),idx(i)],[0,max(stdXX(:))],'color','black');
% end
% xlabel('time')
% ylabel('m')
% legend('error X','error Y','covariance X','covariance Y','location','southeast');

% % plots - Error Vs Covariance - phi
% figure; hold on; grid on
% plot(errorXX(:,3),'g-');
% plot(stdXX(:,3),'g--','linewidth',5);
% for i= 1:length(idx)
%     line([idx(i),idx(i)],[0,max(stdXX(:,3))],'color','black');
% end
% xlabel('time')
% ylabel('rad')
% legend('error phi','covariance phi','location','southeast');

% plots - realPCA Vs calcPCA
figure; hold on; grid on
title('Averaged P(CA)');
xlabel('Time'); ylabel('Averaged P(CA)');
% plot(DATA.realPCA,'-r','linewidth',2);
plot(DATA.calcPCA,'b-','linewidth',2);
plot(DATA.calcPCA_MJ,'g--','linewidth',2);
% legend('Real P(CA)','Average calculated P(CA)','Average calculated P(CA) - MJ','location','southeast');
axis([0,step,0,1]);



% find(errorXX(:,1) > stdXX(:,1))
% find(errorXX(:,2) > stdXX(:,2))
% find(errorXX(:,3) > stdXX(:,3))
































