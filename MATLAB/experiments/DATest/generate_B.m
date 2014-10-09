% generate_B.m
%
% make some climatological B matrices for 3dvar
%
% Andy Reagan
% 2014-10-02

clear all
close all
addpath(genpath('/Users/andyreagan/work/2013/2013-05data-assimilation/src'))

%% parameters
numRuns = 100;

modelname = 'lorenz63';
dim = 3;

%% initialize storage

for window=0.05:0.05:1
    B = zeros(dim);
    %% loop over number of runs (smooth out errors)
    for i=1:numRuns
        model = lorenz63();
        model.init();
        model.window = window;
        initial_x = model.x;
        % generate the truth timeseries
        model.run()
        final_x = model.x;
        
        B = B + (final_x-initial_x)*(final_x-initial_x)';
        
        fprintf('run %d finished\n',i);
        % disp((final_x-initial_x)*(final_x-initial_x)');
    end
    fprintf('final B from window %f:\n',window)
    B = B./numRuns;
    disp(B);
    Binv = inv(B);
    save(sprintf('B%.2f.mat',window),'B','Binv')
end

