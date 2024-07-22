function [solObj,w] = characteristics_solver(f,g,L,kwargs)
%CHARACTERISTICS_SOLVER Returns the wave solution of the interpolation-based method of characteristics, as a solutionObject
%Based on u(x,t) = w(t-x)-w(t+x) and w(t+L(t)) = w(t-L(t))
% INPUTS:
%   f: function_handle, inital position
%   g: function_handle, initial velocity
%   L: function_handle or interpolant object, length of the domain
%   optional nx: double, number of x-points
%   optional t: double, time array
%   optional resolution: double, the density of time points per unit of 1/L
%   optional method: 'linear', 'cubic' or 'spline' (= default), the method of interpolation
% OUTPUTS:
%   solObj: solutionObject (See solutionObject.m)
%   w: interpolant object, returns the function w in u(x,t) = w(t-x)-w(t+x)
% EXAMPLE:
%   [L,R] = get_LR("linearL")
%   [f,g] = get_fg(L,"gauss")
%   sol = characteristics_solver(f,g,L,nx=1000,t=0:1/100:10,resolution=10000,method='spline')

arguments
    f function_handle
    g function_handle
    L
    kwargs.nx = 500
    kwargs.t = 0:(1/120):5
    kwargs.resolution = 50000
    kwargs.method = "spline"
end

print_indented_message("Characteristics solver",true)
timer = tic;

t_max = max(kwargs.t);

w = find_w_interpolation(f,g,L,t_max=t_max,resolution=kwargs.resolution,method=kwargs.method);

t = kwargs.t; % t for solObj from here
s = linspace(0,1,kwargs.nx);
[~,S] = meshgrid(t,s);
X = S .* L(t);
solObj = solutionObject(w(t-X)-w(t+X), X, t);

print_indented_message("Finished in " + string(toc(timer)))

end