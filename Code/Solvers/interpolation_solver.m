function [solObj,R] = interpolation_solver(f,g,L,kwargs)
%INTERPOLATION_SOLVER Returns the wave solution of the interpolation method based on R, as a solutionObject
%By default the optimal number of gamma_l's is sought and selected
% INPUTS:
%   f: function_handle, inital position
%   g: function_handle, initial velocity
%   L: function_handle or interpolant object, length of the domain
%   optional resolution: double, the density of time points per unit of 1/L
%   optional nx: double, number of x-points
%   optional t: double, time array
%   optional n_max: double, number of coefficients to use in the series (from n = -nmax to n = nmax)
%   optional abstol: double, absolute tolerance for integral in calc_coefs
%   optional reltol: double, relative tolerance for integral in calc_coefs
%   optional method: 'linear', 'cubic' or 'spline' (= default), the method of interpolation
%   optional initialR: 'quadratic' or 'cubic' (= default) (or 'linear' or 'sigmoid-like'), the choice of R on [-L0,L0]
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
    kwargs.resolution = 1000
    kwargs.nx = 500
    kwargs.t = 0:(1/120):5
    kwargs.n_max = 100
    kwargs.abstol = 4e-12
    kwargs.reltol = 4e-12
    kwargs.method = "spline"
    kwargs.initialR = "cubic"
    kwargs.R = []
end

print_indented_message("Interpolation solver",true)
timer = tic;

if isempty(kwargs.R)
    R = find_R_interpolation(L,t_max=max(kwargs.t),resolution=kwargs.resolution,method=kwargs.method,initialR=kwargs.initialR);
else
    R = kwargs.R;
end

% now, use eigenmode solver to get solution object
solObj = eigenmode_solver(f,g,R,L,nx=kwargs.nx,t=kwargs.t,n_max=kwargs.n_max,abstol=kwargs.abstol,reltol=kwargs.reltol); % default 1000 nx values iguess

print_indented_message("Finished in " + string(toc(timer)))

end