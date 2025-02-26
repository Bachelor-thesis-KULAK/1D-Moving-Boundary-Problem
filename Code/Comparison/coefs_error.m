function [error_wave, error_timederivative_wave] = coefs_error(f,g,R,L0,n_max,abstol,reltol,testsize)
%COEFS_ERROR Calculates the error bound on the wave amplitude due to how well the initial conditions are satisfied (see Appendix)
%Optionally returns the error bound on the timederivative as well
% INPUTS:
%   f: function_handle or interpolant object, initial position
%   g: function_handle or interpolant object, initial velocity
%   R: function_handle or interpolant object, the transformation R
%   L0: 1x1 double, initial length of the domain
%   optional n_max: double, number of eigenmodes (-n_max, n_max)
%   optional abstol: 1x1 double, absolute tolerance, used for integration
%   optional reltol: 1x1 double, relative tolerance, used for integration
%   optional testsize: 1x1 double, size of the array to look for min and max of w_error
% OUTPUTS:
%   error_wave: double, error bound on wave amplitude, same size as n_max 
%   error_timederivative_wave: double, error bound on time derivative of wave amplitude, same size as n_max
arguments
    f
    g
    R
    L0 (1,1) double
    n_max double = 1000
    abstol (1,1) double = 4e-12
    reltol (1,1) double = 4e-12
    testsize (1,1) double = 50000 % highest oscillating mode is exp(i*n*pi*R(xi)), so on average 50 points per half period
end

[n_max_U,~,iU] = unique(n_max); % sorted list of unique values, makes it easier to calculate
    
% we write the wave solution u(x,t) as u(x,t) = w(t+x) - w(t-x)
% timederivative is w'(t+x) - w'(t-x) = dw(t+x) - dw(t-x)

dR = derivative(R);
df = derivative(f);
num_coefs = max(n_max_U);
coefs = calc_coefs(df, g, R, L0, num_coefs, abstol, reltol);

x = linspace(0,L0,testsize).';

f_exact = f(x) + 0*x;
opts = odeset(RelTol=3e-14,AbsTol=1e-15);
[~,G_exact] = ode89(@(t,y) g(t), x, 0, opts);
if nargout>1
    df_exact = df(x) + 0*x;
    g_exact = g(x) + 0*x;
end

w_exact_left = (-f_exact+G_exact)/2;
w_exact_right = (f_exact+G_exact)/2;
w_exact = [w_exact_left(end:-1:2); w_exact_right]; % first part is for negative values, second for positive
if nargout>1, dw_exact = [df_exact-g_exact, df_exact+g_exact]/2; end

% Unvectorized (for loop) method is faster
X = [-x(end:-1:2); x];
w_approx = zeros(length(X),numel(n_max_U));
w_approx_now = zeros(size(X));
if nargout>1
    dw_approx = zeros(length(X),numel(n_max_U));
    dw_approx_now = zeros(size(X));
end
coef_index = 1;
for n = 1:num_coefs
    complex_exponent = exp(-1i*pi*n.*R(X));
    w_approx_now = w_approx_now - coefs(n) .* complex_exponent;
    if nargout>1, dw_approx_now = dw_approx_now + (1i*pi*n) .* coefs(n) .* dR(X) .* complex_exponent; end
    if n == n_max_U(coef_index)
        w_approx(:,coef_index) = w_approx_now;
        if nargout>1, dw_approx(:,coef_index) = dw_approx_now; end
        coef_index = coef_index + 1;
    end
end

% + Complex Conjugate (because we only looped over positive n):
w_approx = 2*real(w_approx);
if nargout>1, dw_approx = 2*real(dw_approx); end

error_wave = zeros(1,numel(n_max_U));
error_timederivative_wave = zeros(1,numel(n_max_U));

for i = 1:numel(n_max_U)
    w_approx_nm = w_approx(:,i);
    w = w_exact-w_approx_nm;
    error_wave(i) = max(w) - min(w);

    if nargout>1
        dw_approx_nm = dw_approx(:,i);
        dw = dw_exact-dw_approx_nm;
        error_timederivative_wave(i) = max(dw) - min(dw);
    end
end

error_wave = error_wave(iU); % select original order of n_max
if nargout>1, error_timederivative_wave = error_timederivative_wave(iU); end