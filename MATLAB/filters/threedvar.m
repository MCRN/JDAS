function [x_a] = threedvar(x_f,y_o,H,R,B,varargin)

% 3D-Var
%
% Solve the 3D-Var cost function
%
% INPUTS
%   x_f: the previous analysis, column vector
%   y_o: observations
%   H:   obs operator
%   R:   obs cov
%   B:   model cov
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
% (B*H')/(R+H*B*H');
% note that R is easily invertible
% precompute intermediary matrix
V = (inv(B)+H'*inv(R)*H)\(H'*inv(R));

% innovation
d = y_o - H*x_f;

% here we go
x_a = x_f + V*d; 

if ~silent
  disp(B);
end

% this is now passed
%save('P.mat','p_a');
