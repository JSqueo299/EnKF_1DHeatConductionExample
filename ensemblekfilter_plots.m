function ensemblekfilter_plots(tPlot,q,w,v,y_meas,x_EnKF,x_tr,dt,solverRuns,T_FD)
   kPlot = tPlot/(dt*solverRuns);    % Index for specified time to access matrices
   
    blue = [0, 0.4470, 0.7410];
    red = [0.6350, 0.0780, 0.1840];
    purple = [0.4940, 0.1840, 0.5560];
    cyan = [0.3010, 0.7450, 0.9330];
    turq = [0, 0.75, 0.75];
    orange = [0.9100,0.4100,0.1700];
    green = [0.4660, 0.6740, 0.1880];
    yellow = [0.9290, 0.6940, 0.1250];
   
   [row,col] = find(y_meas(:,1));
    %% PLOTS
    figure(2);
    hold on;
    box on;
    plot(row,nonzeros(y_meas(:,kPlot)),'kx','linewidth',2);
    plot(T_FD(kPlot,:),':','Color',red,'linewidth',2);
    plot(x_EnKF(:,kPlot),'Color',turq,'linewidth',2);
    
    
%     str = sprintf('~/OpenFOAM/jns14008-5.x/1D_heatConductionEnKF/%.15g',tPlot);
%     cd(str);
%     
%     fid = fopen('T');
%     temp = textscan(fid,'%f', 'headerlines',22);
%     fclose(fid);
%     Tsolver(:,1) = temp{1};
%     plot(Tsolver(:,1)',':','Color',red,'linewidth',2);
    
%     plot(x_tr(:,kPlot),'r:','linewidth',2);
    legend('measurements','w/o filter', 'filter'); legend boxoff;

%% Frequency Tests
%     plot(y_meas(:,kPlot),'x','linewidth',2);
%     plot(x_EnKF(:,kPlot),'--','linewidth',2);
%     legend('Meas','Original k=80W/m*K','Meas F=1','EnKF F=1', 'Meas F=50','EnKF F=50','Meas F=100','EnKF F=100');
%     ylim([200 700]);  
%%
%     str1 = sprintf('q = %g, t = %g s, w = %g, v = %g', q,tPlot,w,v);
    str1 = sprintf('1D Heat Conduction with q = %g, w = %g, v = %g', q,w,v);
%     str1 = sprintf('1D Heat Transient Conduction with Varying Ensemble Size');
    title(str1);
    str2 = sprintf('T(t = %g s)', tPlot);
    ylabel(str2);
    xlabel('Cell ID');
    ylim([300 600]);
    hold off;

saveDirectory = '~/dataAssimilation/tests';
saveas(gcf, fullfile(saveDirectory, 'old_1DHeatConduction.tif'));
    
end
