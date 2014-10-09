function X_a = EnSRF(X_f,y_o,H,R,delta,varargin)

% EnSRF.m
%
% Ensemble Square Root Filter: compute the new analysis
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

% if there are any Nan's in there...fix it
% for i=1:length(p_f(:,1))
%   for j=1:length(p_f(1,:))
%     if isnan(p_f(i,j))
%       p_f(i,j) = 100;
%     end 
%   end
% end

% calming things down
% p_f(isnan(p_f)) = 100;
% p_f(isinf(p_f)) = 100;
% p_f(p_f>100) = 100;
% p_f(p_f<-100) = -100;

K = (p_f*H')/(R+H*p_f*H');
d = y_o - H*x_f;

if ~silent
    fprintf('======================\n')
end
% % calm the forecast error for initial stuff
% if norm(d) > 5 && p_f(1,1) > 5
%     fprintf('calming by %f\n',norm(d))
%     p_f = p_f/(norm(d)*100);
%     X_f_diff = X_f_diff/(norm(d)*100);
% end

x_a = x_f + K*d;

% now do the square root thing
if ~silent
    fprintf('the forecast X_f is\n')
    disp(X_f)
    fprintf('the obs err cov R is\n')
    disp(R)
    fprintf('the obs is\n')
    disp(y_o)
    fprintf('the forecast err cov p_f is\n')
    disp(p_f)
    disp(H*p_f*H'+R)
    fprintf('the analysis x_a is\n')
    disp(x_a)
end

% Kt = p_f*H'/(sqrtm(H*p_f*H'+R))')/(sqrtm(H*p_f*H'+R)+sqrtm(R));
Kt = p_f*H'*(inv(sqrtm(H*p_f*H'+R))')*inv(sqrtm(H*p_f*H'+R)+sqrtm(R));

X_a_diff = (eye(length(H(:,1)))-Kt*H)*X_f_diff;

% now compute the analysis for each ensemble member
X_a=repmat(x_a,1,N) + X_a_diff;

if ~silent
    fprintf('the analysis ensemble is\n')
    disp(X_a)
end

% if ~silent
%     disp(p_f);
% end