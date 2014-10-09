function X_a = EnKF(X_f,y_o,H,R,delta,varargin)

% EnKF.m
%
% generalized EnKF: compute the new analysis
% "perturbed observation" method
%
% INPUTS
%  X_f:     the whole forecast, columns are model state
%
% OUTPUTS
%  X_a:     the analysis matrix, columns are model state analysis
%
% written by Andy Reagan
% 2013-10-13

silent = 0;

i=1; 
while i<=length(varargin), 
  argok = 1; 
  if ischar(varargin{i}), 
    switch varargin{i},
        case 'silent',       silent  = 1;
        otherwise, argok=0; 
    end
  end
  if ~argok, 
    fprintf('invalid argument %s\n',varargin{1})
  end
  i = i+1; 
end


% set the error the randomly perturbing the analysis IC for ensemble members
% this is the std dev of a normal error
PertError = R(1,1); % or error_max

N = length(X_f(1,:));

% let x_f now be the average
x_f = mean(X_f,2);

% X_f_diff = X_f;
% for j=1:N
%     % this computes the difference from the mean without that forecast
%     % included
%     X_f_diff(:,j) = X_f(:,j)-mean([X_f(:,j+1:end) X_f(:,1:j-1)],2);
%     % this includes all of the forecasts in the mean
%     % X_f_diff(:,j) = X_f(:,j)-mean(X_f,2);
% end
% % for speed, could also try
X_f_diff = X_f - repmat(x_f,1,N);

% multiplicative inflation
X_f_diff = sqrt(1+delta).*X_f_diff;

p_f = 1/(N-1)*(X_f_diff*X_f_diff');

K = (p_f*H')/(R+H*p_f*H');

% now compute the analysis for each ensemble member
X_a=X_f;
for j=1:N
    % perturb y_o according to R
    pertY = y_o+randn(size(y_o))*PertError;

    % innovation
    d = pertY - H*(x_f+X_f_diff(:,j));

    % here we go
    X_a(:,j) = x_f+X_f_diff(:,j) + K*d;
end

if ~silent
  disp(K);
  disp(p_f);
end