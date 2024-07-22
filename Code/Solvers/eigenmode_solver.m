function [solObj] = eigenmode_solver(f,g,R,L,kwargs)
%EIGENMODE_SOLVER Returns the wave solution as a solutionObject
%Based on u(x,t) = Σ c_n ( exp(-i*pi*n*R(t-x)) - exp(-i*pi*n*R(t+x)) )
% INPUTS:
%   f: function_handle, inital position
%   g: function_handle, initial velocity
%   R: function_handle or interpolant object, the transformation
%   L: function_handle or interpolant object, length of the domain
%   optional nx: double, number of x-points
%   optional t: double, time array
%   optional n_max: double, number of coefficients to use in the series (from n = -nmax to n = nmax)
%   optional abstol: double, absolute tolerance for integral in calc_coefs
%   optional reltol: double, relative tolerance for integral in calc_coefs
%   optional coefs: double, if the coefs c_n are known: bypass calc_coefs
% OUTPUTS:
%   solObj: solutionObject (See solutionObject.m)
% EXAMPLE:
%   [L,R] = get_LR("linearL")
%   [f,g] = get_fg(L,"gauss")
%   sol = eigenmode_solver(f,g,R,L)
arguments
    f function_handle
    g function_handle
    R
    L
    kwargs.nx = 500
    kwargs.t = 0:(1/120):5
    kwargs.n_max = 100
    kwargs.abstol = 4e-12
    kwargs.reltol = 4e-12
    kwargs.coefs = []
end

print_indented_message("Eigenmode solver",true)
timer = tic;

[nx, t, n_max, abstol, reltol] = deal(kwargs.nx, kwargs.t, ...
    kwargs.n_max, kwargs.abstol, kwargs.reltol);

if isempty(kwargs.coefs) % if needed, calculate the coefficients
    df = derivative(f);
    coefs = calc_coefs(df,g,R,L(0),n_max,abstol,reltol);
else
    coefs = kwargs.coefs;
    n_max = numel(coefs); % overwrite n_max
end
s = linspace(0,1,nx);
[T,S] = meshgrid(t,s);
X = S .* L(t);

expR1 = exp(-1i*pi*R(T-X));
expR2 = exp(-1i*pi*R(T+X));

U1 = zeros(size(expR1));
U2 = zeros(size(expR1));

for n = 1:n_max 
    % Use Horner to calculate the wave solution
    coef = coefs(end+1-n);
    U1 = ( U1 + coef ) .* expR1;
    U2 = ( U2 + coef ) .* expR2;
end
U = U1-U2;

U = 2*real(U); % no longer looping over negative n-values

solObj = solutionObject(U,X,t);
print_indented_message("Finished in " + string(toc(timer)))

end