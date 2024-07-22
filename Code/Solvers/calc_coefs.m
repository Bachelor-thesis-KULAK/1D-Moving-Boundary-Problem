function [coefs] = calc_coefs(df, g, R, L0, n_max, abstol, reltol)
%CALC_COEFS Calculates the coefficients in the expansion of eigenmodes exp(-i*pi*n*R(t-x)) - exp(-i*pi*n*R(t+x))
% INPUTS:
%   df: function_handle or interpolant object, the derivative of the inital position
%   g: function_handle or interpolant object, initial velocity
%   R: function_handle or interpolant object, the transformation
%   L0: double, initial length of domain at t = 0
%   optional n_max: double, number of coefficients to use in the series (from n = -nmax to n = nmax)
%   optional abstol: double, absolute tolerance for integral in calc_coefs
%   optional reltol: double, relative tolerance for integral in calc_coefs
% OUTPUTS:
%   coefs: 1 x n_max complex double, the coefficients in the series for n > 0
% EXAMPLE:
%   [L,R] = get_LR("linearL")
%   [f,g] = get_fg(L,"gauss")
%   coefs = calc_coefs(derivative(f),g,R,L(0),100,1e-11,1e-11)
arguments
    df
    g
    R
    L0 double
    n_max = 200
    abstol = 4e-12
    reltol = 4e-12
end
% f and g must be extended oddly on the domain [-L(0),L(0)] by letting f(-x) = -f(x) and g(-x) = -g(x)
% so f' must be extended to an even function on the domain [-L(0),L(0)]
print_indented_message("Calc coefs",true)
timer = tic;

df = @(x) df(abs(x)); % even extension of f'
g = @(x) sign(x).*g(abs(x)); % odd extension of g
n_array = 1:n_max;

% Calculate each coefficient
coefs = zeros(1,n_max);
for i = n_array
    integral_fun = @(x) (df(x) + g(x)) .* exp(1i*pi*n_array(i)*R(x));
    coefs(i) = integral(integral_fun, -L0, L0, RelTol=reltol, AbsTol=abstol);
end
coefs = -1i./(4*pi*n_array) .* coefs;

print_indented_message("Finished in " + string(toc(timer)))
end