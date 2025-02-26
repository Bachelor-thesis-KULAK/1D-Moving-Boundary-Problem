function [R_PP] = find_R_interpolation_density(L,kwargs)
%FIND_R_INTERPOLATION_DENSITY Finds the transformation R on [-L(0), tmax+L(tmax)] from L using the interpolation method
%Based on R(t+L(t)) = 2+R(t-L(t))
% INPUTS:
%   L: function_handle or interpolant object, length of the domain
%   optional t_max: double, maximum time to calculate R for
%   optional resolution_density: double, the density of points per time to approximate the density for constructing the non-equidistant time array
%   optional resolution: double, the density of time points per unit of 1/L
%   optional method: 'linear', 'cubic' or 'spline' (= default), the method of interpolation
%   optional initialR: 'quadratic' or 'cubic' (= default) (or 'linear' or 'sigmoid-like'), the choice of R on [-L0,L0]
%   optional derivL0: specify the first (and second) derivative of L at t = 0 instead of calculating it
% OUTPUTS:
%   R_PP: interpolant object, returns the transformation R
% EXAMPLE:
%   [L,R] = get_LR("linearL")
%   R = find_R_interpolation(L,t_max=10,resolution=10000,method='spline')

arguments
    L
    kwargs.t_max {mustBeNonnegative} = 5
    kwargs.resolution_density {mustBePositive} = 1000
    kwargs.resolution {mustBePositive} = 10000 
    kwargs.method = "spline"
    kwargs.initialR = "quadratic"
    kwargs.derivL0 = []
end

print_indented_message("Find R Interpolation",true)
timer = tic;

t_max = kwargs.t_max;
kwargs.resolution = round(kwargs.resolution);

if ~any(kwargs.initialR == ["linear", "quadratic", "cubic", "sigmoid-like"])
    error("Unrecognized InitialR")
end

L0 = L(0);
if isempty(kwargs.derivL0)
    Ldot = derivative(L); Ldot0 = Ldot(0);
    if kwargs.initialR == "cubic"
        Lddot = derivative(Ldot); Lddot0 = Lddot(0);
    end
else
    Ldot0 = kwargs.derivL0(1);
    if numel(kwargs.derivL0) > 1
        Lddot0 = kwargs.derivL0(2);
    end
end

%% Initial domain
% First, find all values of R at t = 0
xi_0 = linspace(-L0, L0, 2*kwargs.resolution-1);
% the density of points in the time array is rho/L0 at t=0
% so we choose rho/L0 * (L0) = rho points between -L0 and 0 and rho points between 0 and L0
if kwargs.initialR == "linear"
    R_0_fun = @(xi_0) (xi_0+L0)/(L0);
elseif kwargs.initialR == "quadratic"
    deriv_left = (1+Ldot0)/L0; deriv_right = (1-Ldot0)/L0; % make sure derivative is continuous
    R_0_fun = @(xi_0) (xi_0 + L0) .* (deriv_left + 1/(4*L0) * (xi_0 + L0) * (deriv_right-deriv_left) ); % first block, used for interpolation for next blocks
elseif kwargs.initialR == "cubic"
    b = 3/L0 * (1+Ldot0+Ldot0^2+Ldot0^3) / (3+Ldot0^2+L0*Lddot0);
    c = -3/(2*L0^2) * (Ldot0+2*Ldot0^2+Ldot0^3-L0*Lddot0) / (3+Ldot0^2+L0*Lddot0);
    d = 1/(2*L0^3) * (2*Ldot0^2-L0*Lddot0) / (3+Ldot0^2+L0*Lddot0);
    R_0_fun = @(xi_0) b * (xi_0 + L0) + c * (xi_0 + L0).^2 + d * (xi_0 + L0).^3;
elseif kwargs.initialR == "sigmoid-like"
    R_0_fun = @(xi_0) 2/pi * atan(sinh(5*tan(pi/2*xi_0/L0))) + 1;
end

R_PP = interpolant(xi_0,R_0_fun(xi_0),kwargs.method);
if t_max == 0
    return
end

%% Find time breakpoints
% t0 = 0
% t1 such that t1-L(t1) = t0+L(t0)
% t2 such that t2-L(t2) = t1+L(t1) ...
t_breakpoints = 0;
while t_max - L(t_max) > t_breakpoints(end) + L(t_breakpoints(end))
    fun_zero = @(t) (t-L(t)) - (t_breakpoints(end) + L(t_breakpoints(end)));
    t_breakpoint = fzero(fun_zero,[0,t_max]);
    if t_max >= t_breakpoint % t_max should be bigger than all breakpoints
        t_breakpoints(end+1) = t_breakpoint;
    else
        break
    end
end
num_breakpoints = length(t_breakpoints);
if num_breakpoints > 1000
    warning("t_max is too high or length L gets too small, too many reflections")
end

%% Make t array with density 1/L*resolution
t = linspace(0,t_max,ceil(1+kwargs.resolution_density*t_max)); % initial equidistant array
dt = (t(end)-t(1))/(length(t)-1);
int_density = cumtrapz(dt,1./L(t)); % integral of density, L is efficiently evaluated in not many points
int_end = int_density(end);

dy = int_end / ceil(int_end*kwargs.resolution); % denominator is approximated number of wished points
y_breakpoints = interp1(t,int_density,t_breakpoints); % include these points
indices_breakpoints = ones(size(y_breakpoints));
y = cell(1,1+num_breakpoints); y{1} = 0; % start
for i = 1:length(y_breakpoints)
    start_point_y = y_breakpoints(i);
    if i ~= length(y_breakpoints)
        end_point_y = y_breakpoints(i+1);
    else
        end_point_y = int_end;
    end
    num_points = ceil( (end_point_y - start_point_y) / dy );
    if i ~= length(y_breakpoints)
        indices_breakpoints(i+1) = indices_breakpoints(i) + num_points;
    end
    y_new = linspace(start_point_y, end_point_y, 1+num_points); % choose equidistant array with wished number of points in transformed space
    y{i+1} = y_new(2:end); % index from 2 to avoid double points
end
y = horzcat(y{:});
t = interp1(int_density,t,y); % inverse of y is an array with given density 1/L
t(indices_breakpoints) = t_breakpoints; % these are the breakpoints

%% Make R
L_in_t = L(t);

R_PP_coefs = cell(1,1+num_breakpoints); R_PP_coefs{1} = R_PP.PP.coefs;
R_PP_breaks = cell(1,1+num_breakpoints); R_PP_breaks{1} = R_PP.PP.breaks;
complete = false;
% Copy the values of R at previously calculated points and add 2
for i = 1:num_breakpoints
    if i == num_breakpoints
        t_indices = t_breakpoints(i) <= t;
    else
        t_indices = t_breakpoints(i) <= t & t <= t_breakpoints(i+1);
    end
    t_arr = t(t_indices);
    if t_arr(1) >= t_max
        break % we are done
    elseif i == num_breakpoints || min(t_max - t_breakpoints(i+1), t_max - t_arr(end)) < 20*eps(t_max)
        t_arr(end) = t_max;
        complete = true;
    elseif numel(t_arr) <= 3
        t_arr = linspace(t_arr(1),t_max,4);
    end
    print_indented_message("Spline object " + i + ": with " + length(t_arr) + " points")
    L_arr = L_in_t(t_indices);
    new_R = 2 + R_PP(t_arr-L_arr);
    new_R_PP = interpolant(t_arr+L_arr, new_R, kwargs.method);
    common_break = R_PP.PP.breaks(end);
    if abs(new_R_PP(common_break) - R_PP(common_break)) > 1000*eps
        error('Interpolation objects are not consistent (this means that R is discontinuous)')
    end
    R_PP = new_R_PP;
    % Merge
    R_PP_coefs{i+1} = new_R_PP.PP.coefs;
    R_PP_breaks{i+1} = new_R_PP.PP.breaks(2:end);
    if complete, break, end
end

% Merge R's
R_PP = interpolant(mkpp(horzcat(R_PP_breaks{:}),vertcat(R_PP_coefs{:})));

print_indented_message("Finished in " + string(toc(timer)))

end