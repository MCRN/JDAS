% DA_test.m
%
% test my EnKF against a free run without DA
%
% Andy Reagan
% 2013-10-15

clear all
close all
addpath(genpath('/Users/andyreagan/work/2013/2013-05data-assimilation/src'))

%% parameters
numRuns = 5;
obs_error_std = 0.05;
runTime = 50; %  model units
window = 1; %  window
modelname = 'lorenz63';

%%
if strcmp(modelname,'lorenz63')
    dim = 3;
else
    dim = 36000*6;
end


%% initialize storage
num_windows = ceil(runTime/window);
truth_vec = zeros(dim,num_windows+1);
ETKF_f_vec = zeros(dim,num_windows+1);
threedvar_f_vec = zeros(dim,num_windows+1);
direct_f_vec = zeros(dim,num_windows+1);
none_f_vec = zeros(dim,num_windows+1);


errors_DA = zeros(3,num_windows+1);
errors_3dvar = zeros(3,num_windows+1);
errors_direct = zeros(3,num_windows+1);
errors_none = zeros(3,num_windows+1);

%% loop over number of runs (smooth out errors)
for i=1:numRuns
    % make a truth run
    % if modelname ...
    truthmodel = lorenz63();
    % end
    truthmodel.init();
    truthmodel.window = window;
    % save the IC because I've give the DA,noDA runs good IC?
    truth_vec(:,1) = truthmodel.x;
    % generate the truth timeseries
    for j=1:num_windows
        truthmodel.run()
        truth_vec(:,j+1) = truthmodel.x;
    end
    
    % observing everything
    H = eye(3);
    R = obs_error_std*eye(3);
    
    % pass all of the observations, including time 0
    obs_pert = truth_vec + obs_error_std*randn(size(truth_vec));
    [ETKF_f_vec,~] = modelDAinterface(@lorenz63,'ETKF',obs_pert,0:window:runTime,H,R,window,0,0);
    [threedvar_f_vec,~] = modelDAinterface(@lorenz63,'3dvar',obs_pert,0:window:runTime,H,R,window,0,0,'climcov',sprintf('B%.2f.mat',window));
    [direct_f_vec,~] = modelDAinterface(@lorenz63,'direct',obs_pert,0:window:runTime,H,R,window,0,0);
    [none_f_vec,~] = modelDAinterface(@lorenz63,'none',obs_pert,0:window:runTime,H,R,window,0,0);
    % compute the running RMS error of each of those forecasts
    for j=1:num_windows+1
        % for each variable
        for k=1:dim
            errors_DA(k,j) = errors_DA(k,j)+rmse(truth_vec(k,1:j),ETKF_f_vec(k,1:j));
            errors_3dvar(k,j) = errors_3dvar(k,j)+rmse(truth_vec(k,1:j),threedvar_f_vec(k,1:j));
            errors_direct(k,j) = errors_direct(k,j)+rmse(truth_vec(k,1:j),direct_f_vec(k,1:j));
            errors_none(k,j) = errors_none(k,j)+rmse(truth_vec(k,1:j),none_f_vec(k,1:j));
        end
    end
    fprintf('run %d finished\n',i);
end

% correct for number of runs
errors_DA = errors_DA./numRuns;
errors_3dvar = errors_3dvar./numRuns;
% errors_direct = errors_direct./numRuns;
errors_none = errors_none./numRuns;

%% plot it

figure;

set(gcf,'DefaultAxesFontname','helvetica');
set(gcf,'DefaultLineColor','r');
set(gcf,'DefaultLineMarkerSize',5);
set(gcf,'DefaultLineMarkerEdgeColor','k');
set(gcf,'DefaultLineMarkerFaceColor','g');
set(gcf,'DefaultAxesLineWidth',0.5);
set(gcf,'PaperPositionMode','auto');

tmpsym = {'o','s','v','o','s','v'};
tmpcol = {'g','b','r','k','c','m'};
 
% plot grey lines
plot(0:window:runTime,errors_DA(1,:),'Color',0.7*[1 1 1]);
hold on;
plot(0:window:runTime,errors_direct(1,:),'Color',0.7*[1 1 1]);
plot(0:window:runTime,errors_3dvar(1,:),'Color',0.7*[1 1 1]);

% plot some marks
i=2;
tmph(i) = plot(0:window:runTime,errors_DA(1,:),'Marker',tmpsym{4},'MarkerFaceColor',tmpcol{4},'LineStyle','none');
legendcell{i} = 'ETKF';
i=1;
tmph(i) = plot(0:window:runTime,errors_direct(1,:),'Marker',tmpsym{5},'MarkerFaceColor',tmpcol{5},'LineStyle','none');
legendcell{i} = 'No Obs';
i=3;
tmph(i) = plot(0:window:runTime,errors_3dvar(1,:),'Marker',tmpsym{i},'MarkerFaceColor',tmpcol{i},'LineStyle','none');
legendcell{i} = '3D-Var';
% 
tmpxlab=xlabel('Time $t$', ...
    'fontsize',30,'verticalalignment','top','fontname','helvetica','interpreter','latex');
% set(tmpxlab,'position',get(tmpxlab,'position') - [0 .07 0]);

tmpylab=ylabel('RMS Error $\epsilon$','fontsize',30,'verticalalignment','bottom','fontname','helvetica','interpreter','latex');
% set(tmpylab,'position',get(tmpylab,'position') + [.05 4 0]);

tmplh = legend(tmph,legendcell,'location','northeast'); %,'No Obs'
set(tmplh,'position',get(tmplh,'position')+[-.1 -0.3 0 0])
% change font
tmplh = findobj(tmplh,'type','text');
set(tmplh,'FontSize',30);
% remove box:
legend boxoff

xlim([0 25]);

psprintcpdf_keeppostscript('DA_test_noname_w3dvar');

% xlim([0 70]);

% psprintcpdf_keeppostscript('DA_test_noname006');

