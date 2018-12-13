% Create plots for different aspects of the 'gen_sim_data_ipas_v_info_parallel'
% sim_data. 

clear variables; clc;

FilterSpec = '*.mat';
[FileName, PathName, ~] = uigetfile(FilterSpec);
theFile = fullfile(PathName,FileName);
temp = load(theFile);
full_sims = temp.full_sims;
sim_meta = temp.sim_meta_data;

% Lists of chosen values in parameter study
grid_list = sim_meta.grid_list;
param_list = sim_meta.param_list;
sensor_list = sim_meta.sensor_list;

sim_fields = fieldnames(full_sims(1).data);

% Preallocate
ipas_iters        = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
ipas_pcost_var    = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
ipas_pcost_true   = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
ipas_pcost_exp    = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
ipas_pcost_inc    = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
ipas_time         = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
ipas_nmeas        = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
ipas_pcent_subopt = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));

FP_iters          = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
FP_pcost_var      = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
FP_pcost_true     = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
FP_pcost_exp      = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
FP_pcost_inc      = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
FP_time           = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
FP_nmeas          = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
FP_pcent_subopt   = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));

MSE_iters         = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
MSE_pcost_var     = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
MSE_pcost_true    = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
MSE_pcost_exp     = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
MSE_pcost_inc     = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
MSE_time          = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
MSE_nmeas         = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
MSE_pcent_subopt  = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));

RAND_iters        = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
RAND_pcost_var    = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
RAND_pcost_true   = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
RAND_pcost_exp    = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
RAND_pcost_inc    = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
RAND_time         = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
RAND_nmeas        = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));
RAND_pcent_subopt = zeros(length(grid_list), length(param_list), length(sensor_list),length(full_sims));


% Process data to use in plots
for n_sim = 1:length(full_sims)
    for n_grid = 1:length(grid_list)
        for n_p = 1:length(param_list)
            for n_s = 1:length(sensor_list)
                for fld = 1:numel(sim_fields)
                    if strcmp(sim_fields{fld}, 'ipas')
                        ipas_iters(n_grid, n_p, n_s, n_sim)      = full_sims(n_sim).data(n_grid, n_p, n_s).ipas.iterations;
                        ipas_pcost_var(n_grid, n_p, n_s, n_sim)  = full_sims(n_sim).data(n_grid, n_p, n_s).ipas.path_cost_var;
                        ipas_pcost_true(n_grid, n_p, n_s, n_sim) = full_sims(n_sim).data(n_grid, n_p, n_s).ipas.path_cost.true;
                        ipas_pcost_exp(n_grid, n_p, n_s, n_sim)  = full_sims(n_sim).data(n_grid, n_p, n_s).ipas.path_cost.expected;
                        ipas_pcost_inc(n_grid, n_p, n_s, n_sim)  = full_sims(n_sim).data(n_grid, n_p, n_s).ipas.path_cost.incurred;
                        ipas_time(n_grid, n_p, n_s, n_sim)       = full_sims(n_sim).data(n_grid, n_p, n_s).ipas.comp_time;
                        ipas_nmeas(n_grid, n_p, n_s, n_sim)      = full_sims(n_sim).data(n_grid, n_p, n_s).ipas.n_measurements;
                        
                        diff = (ipas_pcost_inc(n_grid, n_p, n_s, n_sim) - ipas_pcost_true(n_grid, n_p, n_s, n_sim))/ipas_pcost_true(n_grid, n_p, n_s, n_sim);
                        ipas_pcent_subopt(n_grid, n_p, n_s, n_sim) = diff * 100;
                    end
                    
                    if strcmp(sim_fields{fld}, 'FP')
                        FP_iters(n_grid, n_p, n_s, n_sim)      = full_sims(n_sim).data(n_grid, n_p, n_s).FP.iterations;
                        FP_pcost_var(n_grid, n_p, n_s, n_sim)  = full_sims(n_sim).data(n_grid, n_p, n_s).FP.path_cost_var;
                        FP_pcost_true(n_grid, n_p, n_s, n_sim) = full_sims(n_sim).data(n_grid, n_p, n_s).FP.path_cost.true;
                        FP_pcost_exp(n_grid, n_p, n_s, n_sim)  = full_sims(n_sim).data(n_grid, n_p, n_s).FP.path_cost.expected;
                        FP_pcost_inc(n_grid, n_p, n_s, n_sim)  = full_sims(n_sim).data(n_grid, n_p, n_s).FP.path_cost.incurred;
                        FP_time(n_grid, n_p, n_s, n_sim)       = full_sims(n_sim).data(n_grid, n_p, n_s).FP.comp_time;
                        FP_nmeas(n_grid, n_p, n_s, n_sim)      = full_sims(n_sim).data(n_grid, n_p, n_s).FP.n_measurements;
                        
                        diff = (FP_pcost_inc(n_grid, n_p, n_s, n_sim) - FP_pcost_true(n_grid, n_p, n_s, n_sim))/FP_pcost_true(n_grid, n_p, n_s, n_sim);
                        FP_pcent_subopt(n_grid, n_p, n_s, n_sim) = diff * 100;
                    end
                    
                    if strcmp(sim_fields{fld}, 'MSE')
                        MSE_iters(n_grid, n_p, n_s, n_sim)      = full_sims(n_sim).data(n_grid, n_p, n_s).MSE.iterations;
                        MSE_pcost_var(n_grid, n_p, n_s, n_sim)  = full_sims(n_sim).data(n_grid, n_p, n_s).MSE.path_cost_var;
                        MSE_pcost_true(n_grid, n_p, n_s, n_sim) = full_sims(n_sim).data(n_grid, n_p, n_s).MSE.path_cost.true;
                        MSE_pcost_exp(n_grid, n_p, n_s, n_sim)  = full_sims(n_sim).data(n_grid, n_p, n_s).MSE.path_cost.expected;
                        MSE_pcost_inc(n_grid, n_p, n_s, n_sim)  = full_sims(n_sim).data(n_grid, n_p, n_s).MSE.path_cost.incurred;
                        MSE_time(n_grid, n_p, n_s, n_sim)       = full_sims(n_sim).data(n_grid, n_p, n_s).MSE.comp_time;
                        MSE_nmeas(n_grid, n_p, n_s, n_sim)      = full_sims(n_sim).data(n_grid, n_p, n_s).MSE.n_measurements;
                        
                        diff = (MSE_pcost_inc(n_grid, n_p, n_s, n_sim) - MSE_pcost_true(n_grid, n_p, n_s, n_sim))/MSE_pcost_true(n_grid, n_p, n_s, n_sim);
                        MSE_pcent_subopt(n_grid, n_p, n_s, n_sim) = diff * 100;
                    end
                    
                    if strcmp(sim_fields{fld}, 'RAND')
                        RAND_iters(n_grid, n_p, n_s, n_sim)      = full_sims(n_sim).data(n_grid, n_p, n_s).RAND.iterations;
                        RAND_pcost_var(n_grid, n_p, n_s, n_sim)  = full_sims(n_sim).data(n_grid, n_p, n_s).RAND.path_cost_var;
                        RAND_pcost_true(n_grid, n_p, n_s, n_sim) = full_sims(n_sim).data(n_grid, n_p, n_s).RAND.path_cost.true;
                        RAND_pcost_exp(n_grid, n_p, n_s, n_sim)  = full_sims(n_sim).data(n_grid, n_p, n_s).RAND.path_cost.expected;
                        RAND_pcost_inc(n_grid, n_p, n_s, n_sim)  = full_sims(n_sim).data(n_grid, n_p, n_s).RAND.path_cost.incurred;
                        RAND_time(n_grid, n_p, n_s, n_sim)       = full_sims(n_sim).data(n_grid, n_p, n_s).RAND.comp_time;
                        RAND_nmeas(n_grid, n_p, n_s, n_sim)      = full_sims(n_sim).data(n_grid, n_p, n_s).RAND.n_measurements;
                        
                        diff = (RAND_pcost_inc(n_grid, n_p, n_s, n_sim) - RAND_pcost_true(n_grid, n_p, n_s, n_sim))/RAND_pcost_true(n_grid, n_p, n_s, n_sim);
                        RAND_pcent_subopt(n_grid, n_p, n_s, n_sim) = diff * 100;
                    end
                end
            end
        end
    end
end


% Plot Iterations
figure()
linestyles = cellstr(char('-','--'));
ncolors = length(param_list)/2;
% np_colors = distinguishable_colors(ncolors);
% np_colors = linspecer(ncolors,'qualitative');
% np_colors = linspecer(ncolors,'sequential');
np_colors = [
    0.00  0.00  1.00
    0.00  0.50  0.00
    1.00  0.00  0.00
    0.00  0.75  0.75];
axes('NextPlot','replacechildren', 'ColorOrder',np_colors);
for n_p = 1:length(param_list)
    plot(sensor_list, mean(squeeze(ipas_iters(1,n_p,:,:)),2),...
        'LineStyle',linestyles{ceil(n_p/(length(param_list)/2))},...
        'LineWidth',2)
    hold on;
    legendInfo{n_p} = ['Np = ' num2str(param_list(n_p)^2)];
end
set(gca, 'FontSize', 14)
title('IPAS Iterations to Convergence')
xlabel('No. of Sensors available, N_S')
xlim([min(sensor_list), max(sensor_list)])
ylabel('Iterations')
legend(legendInfo)

% IPAS Iterations log-curve
figure()
for n_p = 1:length(param_list)
    loglog(sensor_list, mean(squeeze(ipas_iters(1,n_p,:,:)),2))
    hold on;
    legendInfo{n_p} = ['Np = ' num2str(param_list(n_p)^2)];
end
set(gca, 'FontSize', 14)
title('IPAS Iterations to Convergence')
xlabel('No. of Sensors available, N_S')
xlim([min(sensor_list), max(sensor_list)])
ylabel('Iterations')
legend(legendInfo, 'Box', 'off')

% Plot Path cost variance
figure()
for n_p = 1:length(param_list)
    semilogy(sensor_list, mean(squeeze(ipas_pcost_var(1,n_p,:,:)),2))
    hold on;
    legendInfo{n_p} = ['Np = ' num2str(param_list(n_p)^2)];
end
legend(legendInfo)
title('IPAS Final Path Cost variance')
xlabel('No. of Sensors available')
ylabel('variance')

% Plot IPAS computation time
figure()
for n_p = 1:length(param_list)
    semilogy(sensor_list, mean(squeeze(ipas_time(1,n_p,:,:)),2))
    hold on;
    legendInfo{n_p} = ['Np = ' num2str(param_list(n_p)^2)];
end
legend(legendInfo)
title('IPAS Computation time')
xlabel('No. of Sensors available')
ylabel('Time, s')

% Plot IPAS number of measurements
figure()
ncolors = length(param_list);
np_colors = distinguishable_colors(ncolors);
% np_colors = linspecer(ncolors,'qualitative');
% np_colors = linspecer(ncolors,'sequential');
axes('NextPlot','replacechildren', 'ColorOrder',np_colors);
for n_p = 1:length(param_list)
    plot(sensor_list, mean(squeeze(ipas_nmeas(1,n_p,:,:)),2))
    hold on;
    legendInfo{n_p} = ['Np = ' num2str(param_list(n_p)^2)];
end
set(gca, 'FontSize', 14)
legend(legendInfo, 'Location','Best', 'Box', 'off')
title('IPAS - No. of measurements')
xlabel('No. of Sensors available')
xlim([min(sensor_list), max(sensor_list)])
ylabel('No. of measurements')

% Plot IPAS number of measurements X-axis: NP
figure()
nscolors = ceil(length(sensor_list)/2);
ns_colors = [
    0.00  0.00  1.00
    0.00  0.50  0.00
    1.00  0.00  0.00
    0.00  0.75  0.75
    0.75  0.00  0.75];
m00 = 1;
axes('NextPlot','replacechildren', 'ColorOrder',ns_colors);
for n_s = 10:10:sensor_list(end)
    plot(param_list.^2, mean(squeeze(ipas_nmeas(1,:,n_s,:)),2),...
        'LineStyle',linestyles{ceil(n_s/nscolors)},...
        'LineWidth',2)
    hold on;
    legendInfoNs{m00} = ['Ns = ' num2str(n_s)];
    m00 = m00 + 1;
end
set(gca, 'FontSize', 14)
legend(legendInfoNs, 'Location','Best', 'Box', 'off')
title('IPAS - No. of measurements')
xlabel('N_P')
xlim([min(param_list)^2, max(param_list)^2])
ylabel('No. of measurements')

% SUBOPTIMALITY OF DIFFERENT SENSOR PLACEMENT METHODS
% Plot IPAS percent suboptimality
figure()
axes('NextPlot','replacechildren', 'ColorOrder',np_colors);
for n_p = 1:length(param_list)
    plot(sensor_list, mean(squeeze(ipas_pcent_subopt(1,n_p,:,:)),2),...
        'LineStyle',linestyles{ceil(n_p/(length(param_list)/2))})
    hold on;
    legendInfo{n_p} = ['Np = ' num2str(param_list(n_p)^2)];
end
legend(legendInfo)
title('IPAS - Path Cost Suboptimality')
xlabel('No. of Sensors available, N_S')
xlim([min(sensor_list), max(sensor_list)])
ylabel('% (True - Inc)/True')

% Plot FP percent suboptimality
figure()
axes('NextPlot','replacechildren', 'ColorOrder',np_colors);
for n_p = 1:length(param_list)
    plot(sensor_list, mean(squeeze(FP_pcent_subopt(1,n_p,:,:)),2),...
        'LineStyle',linestyles{ceil(n_p/(length(param_list)/2))})
    hold on;
    legendInfo{n_p} = ['Np = ' num2str(param_list(n_p)^2)];
end
legend(legendInfo)
title('Frame Potential - Path Cost Suboptimality')
xlabel('No. of Sensors available, N_S')
xlim([min(sensor_list), max(sensor_list)])
ylabel('% (True - Inc)/True')

% Plot MSE percent suboptimality
if strcmp(sim_fields{fld}, 'MSE')
    figure()
    axes('NextPlot','replacechildren', 'ColorOrder',np_colors);
    for n_p = 1:length(param_list)
        plot(sensor_list, mean(squeeze(MSE_pcent_subopt(1,n_p,:,:)),2),...
            'LineStyle',linestyles{ceil(n_p/(length(param_list)/2))})
        hold on;
        legendInfo{n_p} = ['Np = ' num2str(param_list(n_p)^2)];
    end
    legend(legendInfo)
    title('Mean Squared Error - Path Cost Suboptimality')
    xlabel('No. of Sensors available, N_S')
    xlim([min(sensor_list), max(sensor_list)])
    ylabel('% (True - Inc)/True')
end

% Plot RAND percent suboptimality
figure()
axes('NextPlot','replacechildren', 'ColorOrder',np_colors);
for n_p = 1:length(param_list)
    plot(sensor_list, mean(squeeze(RAND_pcent_subopt(1,n_p,:,:)),2),...
        'LineStyle',linestyles{ceil(n_p/(length(param_list)/2))})
    hold on;
    legendInfo{n_p} = ['Np = ' num2str(param_list(n_p)^2)];
end
legend(legendInfo)
title('Random Placement - Path Cost Suboptimality')
xlabel('No. of Sensors available, N_S')
xlim([min(sensor_list), max(sensor_list)])
ylabel('% (True - Inc)/True')

% Compare variance of different techniques
figure()
num_p = length(param_list);
n_subplot = [1 1; 1 2; 1 3; 2 2; 2 3; 2 3; 3 3; 3 3; 3 3];
for n_p = 1:num_p
    subplot(n_subplot(num_p,1), n_subplot(num_p,2), n_p)
    semilogy(sensor_list, mean(squeeze(ipas_pcost_var(1,n_p,:,:)),2),'b','LineWidth',2)
    hold on
    semilogy(sensor_list, mean(squeeze(FP_pcost_var(1,n_p,:,:)),2),'r','LineWidth',2)
    if strcmp(sim_fields{fld}, 'MSE')
        semilogy(sensor_list, mean(squeeze(MSE_pcost_var(1,n_p,:,:)),2),'k','LineWidth',2)
    end
    semilogy(sensor_list, mean(squeeze(RAND_pcost_var(1,n_p,:,:)),2),'g','LineWidth',2)
    
    
    set(gca,'FontSize',16)
    title(['N_P = ',num2str(param_list(n_p)^2)])
    xlabel('N_S'); ylabel('Var(J(v))');
    xlim([min(sensor_list), max(sensor_list)])
%     ylim([0 (max(MI_time))]); xlim([0 ss]);
end

JOURNAL = true;
if JOURNAL
   % Np = 4
   n_p = 1;
   figure()
   semilogy(sensor_list, mean(squeeze(ipas_pcost_var(1,n_p,:,:)),2),'b','LineWidth',2)
   hold on
   semilogy(sensor_list, mean(squeeze(FP_pcost_var(1,n_p,:,:)),2),'r','LineWidth',2)
   if strcmp(sim_fields{fld}, 'MSE')
      semilogy(sensor_list, mean(squeeze(MSE_pcost_var(1,n_p,:,:)),2),'k','LineWidth',2)
   end
   semilogy(sensor_list, mean(squeeze(RAND_pcost_var(1,n_p,:,:)),2),'g','LineWidth',2)
   
   
   set(gca,'FontSize',16)
   title(['N_P = ',num2str(param_list(n_p)^2)])
   xlabel('$N_S$','Interpreter','latex'); 
   ylabel('$Var[\hat{J}(\mathbf{v}_{\bar{\ell}}^*)]$','Interpreter','latex');
   xlim([min(sensor_list), max(sensor_list)])
   legend('IPAS','FP')
   
   % Np = 25
   n_p = 4;
   figure()
   semilogy(sensor_list, mean(squeeze(ipas_pcost_var(1,n_p,:,:)),2),'b','LineWidth',2)
   hold on
   semilogy(sensor_list, mean(squeeze(FP_pcost_var(1,n_p,:,:)),2),'r','LineWidth',2)
   if strcmp(sim_fields{fld}, 'MSE')
      semilogy(sensor_list, mean(squeeze(MSE_pcost_var(1,n_p,:,:)),2),'k','LineWidth',2)
   end
   semilogy(sensor_list, mean(squeeze(RAND_pcost_var(1,n_p,:,:)),2),'g','LineWidth',2)
   
   
   set(gca,'FontSize',16)
   title(['N_P = ',num2str(param_list(n_p)^2)])
   xlabel('$N_S$','Interpreter','latex'); 
   ylabel('$Var[\hat{J}(\mathbf{v}_{\bar{\ell}}^*)]$','Interpreter','latex');
   xlim([min(sensor_list), max(sensor_list)])
   
   % Np = 49
   n_p = 6;
   figure()
   semilogy(sensor_list, mean(squeeze(ipas_pcost_var(1,n_p,:,:)),2),'b','LineWidth',2)
   hold on
   semilogy(sensor_list, mean(squeeze(FP_pcost_var(1,n_p,:,:)),2),'r','LineWidth',2)
   if strcmp(sim_fields{fld}, 'MSE')
      semilogy(sensor_list, mean(squeeze(MSE_pcost_var(1,n_p,:,:)),2),'k','LineWidth',2)
   end
   semilogy(sensor_list, mean(squeeze(RAND_pcost_var(1,n_p,:,:)),2),'g','LineWidth',2)
   
   
   set(gca,'FontSize',16)
   title(['N_P = ',num2str(param_list(n_p)^2)])
   xlabel('$N_S$','Interpreter','latex'); 
   ylabel('$Var[\hat{J}(\mathbf{v}_{\bar{\ell}}^*)]$','Interpreter','latex');
   xlim([min(sensor_list), max(sensor_list)])
   legend('IPAS','FP')
   
   % Np = 81
   n_p = 8;
   figure()
   semilogy(sensor_list, mean(squeeze(ipas_pcost_var(1,n_p,:,:)),2),'b','LineWidth',2)
   hold on
   semilogy(sensor_list, mean(squeeze(FP_pcost_var(1,n_p,:,:)),2),'r','LineWidth',2)
   if strcmp(sim_fields{fld}, 'MSE')
      semilogy(sensor_list, mean(squeeze(MSE_pcost_var(1,n_p,:,:)),2),'k','LineWidth',2)
   end
   semilogy(sensor_list, mean(squeeze(RAND_pcost_var(1,n_p,:,:)),2),'g','LineWidth',2)
   
   
   set(gca,'FontSize',16)
   title(['N_P = ',num2str(param_list(n_p)^2)])
   xlabel('$N_S$','Interpreter','latex'); 
   ylabel('$Var[\hat{J}(\mathbf{v}_{\bar{\ell}}^*)]$','Interpreter','latex');
   xlim([min(sensor_list), max(sensor_list)])
   legend('IPAS','FP')
end

% Compare computation time of different techniques
figure()
num_p = length(param_list);
n_subplot = [1 1; 1 2; 1 3; 2 2; 2 3; 2 3; 3 3; 3 3; 3 3];
for n_p = 1:num_p
    subplot(n_subplot(num_p,1), n_subplot(num_p,2), n_p)
    semilogy(sensor_list, mean(squeeze(ipas_time(1,n_p,:,:)),2),'b','LineWidth',2)
    hold on
    semilogy(sensor_list, mean(squeeze(FP_time(1,n_p,:,:)),2),'r','LineWidth',2)
    if strcmp(sim_fields{fld}, 'MSE')
        semilogy(sensor_list, mean(squeeze(MSE_time(1,n_p,:,:)),2),'k','LineWidth',2)
    end
    semilogy(sensor_list, mean(squeeze(RAND_time(1,n_p,:,:)),2),'g','LineWidth',2)
    
    
    set(gca,'FontSize',16)
    title(['N_P = ',num2str(param_list(n_p)^2)])
    xlabel('N_S'); ylabel('Time, s');
%     ylim([0 (max(MI_time))]); xlim([0 ss]);
end

%% GRID SIZE COMPARISONS
% 
% % Compare computation time of different techniques
% figure()
% num_p = length(param_list);
% n_subplot = [1 1; 1 2; 1 3; 2 2; 2 3; 2 3; 3 3; 3 3; 3 3];
% for n_p = 1:num_p
%     subplot(n_subplot(num_p,1), n_subplot(num_p,2), n_p)
%     semilogy(sensor_list, mean(squeeze(ipas_time(1,n_p,:,:)),2),'b','LineWidth',2)
%     hold on
%     semilogy(sensor_list, mean(squeeze(FP_time(1,n_p,:,:)),2),'r','LineWidth',2)
%     semilogy(sensor_list, mean(squeeze(MSE_time(1,n_p,:,:)),2),'k','LineWidth',2)
%     semilogy(sensor_list, mean(squeeze(RAND_time(1,n_p,:,:)),2),'g','LineWidth',2)
%     
%     
%     set(gca,'FontSize',16)
%     title(['N_P = ',num2str(param_list(n_p)^2)])
%     xlabel('N_S'); ylabel('Time, s');
% %     ylim([0 (max(MI_time))]); xlim([0 ss]);
% end
% suptitle(['Grid Size = ',num2str(grid_list(1)^2)]);
% % set(h,'FontSize',20,'FontWeight','normal')
% 
% figure()
% num_p = length(param_list);
% n_subplot = [1 1; 1 2; 1 3; 2 2; 2 3; 2 3; 3 3; 3 3; 3 3];
% for n_p = 1:num_p
%     subplot(n_subplot(num_p,1), n_subplot(num_p,2), n_p)
%     semilogy(sensor_list, mean(squeeze(ipas_time(2,n_p,:,:)),2),'b','LineWidth',2)
%     hold on
%     semilogy(sensor_list, mean(squeeze(FP_time(2,n_p,:,:)),2),'r','LineWidth',2)
%     semilogy(sensor_list, mean(squeeze(MSE_time(2,n_p,:,:)),2),'k','LineWidth',2)
%     semilogy(sensor_list, mean(squeeze(RAND_time(2,n_p,:,:)),2),'g','LineWidth',2)
%     
%     
%     set(gca,'FontSize',16)
%     title(['N_P = ',num2str(param_list(n_p)^2)])
%     xlabel('N_S'); ylabel('Time, s');
% %     ylim([0 (max(MI_time))]); xlim([0 ss]);
% end
% suptitle(['Grid Size = ',num2str(grid_list(2)^2)])


