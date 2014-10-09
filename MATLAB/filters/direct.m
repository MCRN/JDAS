function x_a = direct(x_f,y_o,H,varargin)

% direct.m
%
% directly insert the observations
%
% INPUTS
%  x_f:     the forecast
%
% OUTPUTS
%  x_a:     the analysis
%
% written by Andy Reagan
% 2013-10-14

d = y_o - H*x_f;

% disp(H)
% disp(d)
% disp(y_o)
% disp(x_f)

x_a = x_f + H'*d;


