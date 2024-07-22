function [solObj,R] =  moore_solver(f,g,L,kwargs)
%MOORE_SOLVER Returns the wave solution of Moore's method, as a solutionObject
%Based on R(t+L(t)) = 2+R(t-L(t))
% INPUTS:
%   f: function_handle, inital position
%   g: function_handle, initial velocity
%   L: function_handle or interpolant object, length of the domain
%   optional nx: double, number of x-points
%   optional t: double, time array
%   optional n_max: double, number of coefficients to use in the series (from n = -nmax to n = nmax)
%   optional abstol: double, absolute tolerance for integral in calc_coefs
%   optional reltol: double, relative tolerance for integral in calc_coefs
%   optional Kmax: maximum number of gamma_l's to be calculated (more than 4 is not feasible except for simple examples)
%   optional R: if specified, R will not be calculated anymore and this R will be used instead
% OUTPUTS:
%   solObj: solutionObject (See solutionObject.m)
%   R: interpolant object, returns the transformation R
% EXAMPLE:
%   [L,R] = get_LR("linearL")
%   [f,g] = get_fg(L,"gauss")
%   sol = interpolation_solver(f,g,L,nx=1000,t=0:1/100:10,resolution=10000,method='spline')

arguments
    f function_handle
    g function_handle
    L
    kwargs.nx = 500
    kwargs.t = 0:(1/120):10
    kwargs.n_max = 100
    kwargs.abstol = 4e-12
    kwargs.reltol = 4e-12
    kwargs.Kmax = 20
    kwargs.R = []
end

print_indented_message("Moore solver",true)
timer = tic;

if isempty(kwargs.R)
    R = find_R_moore(L,t_max=max(kwargs.t),Kmax=kwargs.Kmax);
else
    R = kwargs.R;
end

% now, use eigenmode solver to get solution object
solObj = eigenmode_solver(f,g,R,L,nx=kwargs.nx,t=kwargs.t,n_max=kwargs.n_max,abstol=kwargs.abstol,reltol=kwargs.reltol);

print_indented_message("Finished in " + string(toc(timer)))
end