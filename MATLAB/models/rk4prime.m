function [tout, yout, L] = rk4prime(FunFcn,J,params, tspan, y0, ssize)

% RK4	Integrates a system of ordinary differential equations using
%	the fourth order Runge-Kutta  method.  See also ODE45 and
%	ODEDEMO.M.
%	[t,y] = rk4('yprime', tspan, y0) integrates the system
%	of ordinary differential equations described by the M-file
%	yprime.m over the interval tspan=[t0,tfinal] and using initial
%	conditions Y0.
%	[t, y] = rk4(F, tspan, y0, ssize) uses step size ssize
%
% INPUT:
% F     - String containing name of user-supplied problem description.
%         Call: yprime = fun(t,y,params) where F = 'fun'.
%         t      - Time (scalar).
%         y      - Solution column-vector.
%         yprime - Returned derivative column-vector; yprime(i) = dy(i)/dt.
%         params - Parameters of the function, cell array
% tspan = [t0, tfinal], where t0 is the initial value of t, and tfinal is
%         the final value of t.
% y0    - Initial value column-vector.
% ssize - The step size to be used. (Default: ssize = (tfinal - t0)/100).
%
% OUTPUT:
% t  - Returned integration time points (column-vector).
% y  - Returned solution, one solution column-vector per tout-value.
%
% The result can be displayed by: plot(t,y).

% Initialization

t0=tspan(1);
tfinal=tspan(2);
pm = sign(tfinal - t0);  % Which way are we computing?
if nargin < 4, ssize = (tfinal - t0)/100; end
if ssize < 0, ssize = -ssize; end
h = pm*ssize;
t = t0;
y = y0(:);

% We need to compute the number of steps.

dt = abs(tfinal - t0);
N = floor(dt/ssize) + 1;
if (N-1)*ssize < dt
  N = N + 1;
end

% Initialize the output.

tout = zeros(N,1);
tout(1) = t;
yout = zeros(N,size(y,1));
yout(1,:) = y.';
k = 1;

% The main loop 

% initialize L
dim = length(y0);
L = eye(dim);
while k < N
  if pm*(t + h - tfinal) > 0 
    h = tfinal - t; 
    tout(k+1) = tfinal;
  else
    tout(k+1) = t0 +k*h;
  end
  k = k + 1;
  % Compute the slopes
  s1 = feval(FunFcn,t, y,params); s1 = s1(:);
  L1 = feval(J,t,y,params); 
  s2 = feval(FunFcn,t + h/2, y + h/2*s1,params); s2=s2(:);
  L2 = feval(J,t+h/2,y+h/2*s1,params)*(eye(dim)+h/2*L1);
  s3 = feval(FunFcn,t + h/2, y + h/2*s2,params); s3=s3(:);
  L3 = feval(J,t+h/2,y+h/2*s2,params)*(eye(dim)+h/2*L2);
  s4 = feval(FunFcn,t + h, y + h*s3,params); s4=s4(:);
  L4 = feval(J,t+h,y+h*s3,params)*(eye(dim)+h*L3);
  % Compute L here
  % L = L*(eye(dim) + h/6*(feval(J,t,y,params)+2*feval(J,t+h/2,y+h/2*s1,params)*(eye(dim)+h/2*feval(J,t,y,params))+2*feval(J,t+h/2,y+h/2*s2,params)*(eye(dim)+h/2*feval(J,t+h/2,y+h/2*s1,params))*(eye(dim)+h/2*feval(J,t,y,params))+feval(J,t+h,y+h*s3,params)*(eye(dim)+h*feval(J,t+h/2,y+h/2*s2,params))*(eye(dim)+h/2*feval(J,t+h/2,y+h/2*s1,params))*(eye(dim)+h/2*feval(J,t,y,params))));
  L = L*(eye(dim) + h/6*(L1 + 2*L2 + 2*L3 + L4));
  y = y + h*(s1 + 2*s2 + 2*s3 + s4)/6;
  t = tout(k);
  yout(k,:) = y.';
end;
