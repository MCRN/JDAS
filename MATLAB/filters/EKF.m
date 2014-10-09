function [x_a,p_a] = EKF(x_f,y_o,H,R,p_f,varargin)

% EKF
%
% generalized EKF: compute a new error covariance and forecast
%
% this one is different since it's only going to run on the lorenz63...
% so there is a lot of hard-coded function calling and I needed to pass in 
% t,window_len,x_a (prev anal to run the TLM)...
%
% INPUTS
%   x_f: the previous analysis, column vector
%   y_o: observations
%   H:   obs operator
%
% OUTPUTS
%  x_a -> the analysis
%
% written by Andy Reagan
% 13-06-13

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

% disp(size(p_f));
% disp(size(H));
% disp(size(R));

% now need to compute the new analysis
K = (p_f*H')/(R+H*p_f*H');

% innovation
d = y_o - H*x_f;

% analysis
x_a = x_f + K*d; 
p_a = (eye(length(x_f)) - K*H)*p_f;

if ~silent
  disp(p_f);
  disp(p_a);
end

% this is now passed
%save('P.mat','p_a');
