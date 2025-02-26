function [w_PP] = find_w_interpolation_density(f,g,L,kwargs)
%FIND_W_INTERPOLATION_DENSITY Finds the function w on [-L(0), tmax+L(tmax)] from L using the characteristics method
%Based on w(t+L(t)) = w(t-L(t))
% INPUTS:
%   f: function_handle, inital position
%   g: function_handle, initial velocity
%   L: function_handle or interpolant object, length of the domain
%   optional t_max: double, maximum time to calculate w for
%   optional resolution: double, the density of time points per unit of 1/L
%   optional method: 'linear', 'cubic' or 'spline' (= default), the method of interpolation
% OUTPUTS:
%   w_PP: interpolant object, returns the function w
% EXAMPLE:
%   [L,R] = get_LR("linearL")
%   [f,g] = get_fg(L,"gauss")
%   w = find_w_interpolation(f,g,L,t_max=10,resolution=10000,method='spline')

arguments
    f function_handle
    g function_handle
    L
    kwargs.t_max {mustBeNonnegative} = 5
    kwargs.resolution_density {mustBePositive} = 1000
    kwargs.resolution {mustBePositive} = 50000
    kwargs.method = "spline"
end

print_indented_message("Find w Interpolation",true)
timer = tic;

t_max = kwargs.t_max;
kwargs.resolution = round(kwargs.resolution);

%% Initial domain
% First, find all values of w at t = 0
xi_0_neg = linspace(-L(0), 0, kwargs.resolution); xi_0_neg = xi_0_neg(1:end-1); % we only want 0 once in our xi_0 array!
xi_0_pos = linspace(0, L(0), kwargs.resolution);
xi_0 = [xi_0_neg, xi_0_pos];

opts = odeset(RelTol=3e-14,AbsTol=1e-15);
[~,U0] = ode45(@(t,y) g(t), xi_0_pos, 0, opts);
U0 = reshape(U0, size(xi_0_pos));
w_0_neg = ( f(-xi_0_neg) - U0(end:-1:2) ) / 2;
w_0_pos = ( -f(xi_0_pos) - U0 ) / 2;
w_0 = [w_0_neg, w_0_pos];

w_PP = interpolant(xi_0,w_0,kwargs.method);
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

%% Make w
L_in_t = L(t);

w_PP_coefs = cell(1,1+num_breakpoints); w_PP_coefs{1} = w_PP.PP.coefs;
w_PP_breaks = cell(1,1+num_breakpoints); w_PP_breaks{1} = w_PP.PP.breaks;
complete = false;
% Copy the values of w at previously calculated points
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
    new_w = w_PP(t_arr-L_arr);
    new_w_PP = interpolant(t_arr+L_arr, new_w, kwargs.method);
    common_break = w_PP.PP.breaks(end);
    if abs(new_w_PP(common_break) - w_PP(common_break)) > 10^(-10)
        error('Interpolation objects are not consistent (this means that w is discontinuous)')
    end
    w_PP = new_w_PP;
    % Merge
    w_PP_coefs{i+1} = new_w_PP.PP.coefs;
    w_PP_breaks{i+1} = new_w_PP.PP.breaks(2:end);
    if complete, break, end
end

% Merge w's
w_PP = interpolant(mkpp(horzcat(w_PP_breaks{:}),vertcat(w_PP_coefs{:})));

print_indented_message("Finished in " + string(toc(timer)))

end