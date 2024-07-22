function [R_PP] = find_R_interpolation(L,kwargs)
%FIND_R_INTERPOLATION Finds the transformation R on [-L(0), tmax+L(tmax)] from L using the interpolation method
%Based on R(t+L(t)) = 2+R(t-L(t))
% INPUTS:
%   L: function_handle or interpolant object, length of the domain
%   optional t_max: double, maximum time to calculate R for
%   optional resolution: double, the density of time points per unit of 1/L
%   optional method: 'linear', 'cubic' or 'spline' (= default), the method of interpolation
%   optional initialR: 'quadratic' or 'cubic' (= default) (or 'linear' or 'sigmoid-like'), the choice of R on [-L0,L0]
%   optional margin: if specified, the function R will be calculated from -L(0)-margin to tmax+L(tmax)+margin (a little outside the domain specified)
%   optional derivL0: specify the first (and second) derivative of L at t = 0 instead of calculating it
% OUTPUTS:
%   R_PP: interpolant object, returns the transformation R
% EXAMPLE:
%   [L,R] = get_LR("linearL")
%   R = find_R_interpolation(L,t_max=10,resolution=10000,method='spline')

arguments
    L
    kwargs.t_max = 5 
    kwargs.resolution = 1000
    kwargs.method = "spline"
    kwargs.initialR = "cubic"
    kwargs.margin = 0
    kwargs.derivL0 = []
end

print_indented_message("Find R Interpolation",true)
timer = tic;

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

%% Find time breakpoints
% t0 = 0
% t1 such that t1-L(t1) = t0+L(t0)
% t2 such that t2-L(t2) = t1+L(t1) ...
t_breakpoints = 0;
t_maxmarg = kwargs.t_max + kwargs.margin;
while t_maxmarg - L(t_maxmarg) >= t_breakpoints(end) + L(t_breakpoints(end))
    fun_zero = @(t) (t-L(t)) - (t_breakpoints(end) + L(t_breakpoints(end)));
    t_breakpoints(end+1) = fzero(fun_zero,[0,t_maxmarg]);
end
t_breakpoints_before = 0;
t_min = -kwargs.margin;
while t_min + L(t_min) <= t_breakpoints_before(end) - L(t_breakpoints_before(end))
    fun_zero = @(t) (t+L(t)) - (t_breakpoints_before(end)-L(t_breakpoints_before(end)));
    t_breakpoints_before(end+1) = fzero(fun_zero,[t_min,0]);
end

%% Make t array with density 1/L*resolution
if ceil(1000*kwargs.t_max) < 2
    t = 0;
else
    t = linspace(0,kwargs.t_max,ceil(1000*kwargs.t_max)); % initial equidistant array
    dt = (t(end)-t(1))/(length(t)-1);
    int_density = cumtrapz(dt,1./L(t)); % integral of density, L is efficiently evaluated in not many points
    num_points = ceil( int_density(end)*kwargs.resolution );
    if num_points > 10^8
        error("length L gets too small (or resolution too high), running this program will crash")
    end
    y = linspace(0,int_density(end),num_points); % choose equidistant array with wished number of points in transformed space
    t = interp1(int_density,t,y,'spline'); % inverse of y is an array with given density 1/L
end

%% Make t_after array for xi > tmax+L(tmax)
if kwargs.margin > 0
    % Make t_after array
    t_after = linspace(kwargs.t_max,t_maxmarg,ceil(1000*kwargs.margin)); % initial equidistant array
    dt_after = (t_after(end)-t_after(1))/(length(t_after)-1);
    int_density = cumtrapz(dt_after,1./L(t_after)); % integral of density, L is efficiently evaluated in not many points
    num_points = ceil( int_density(end)*kwargs.resolution );
    if num_points > 10^8
        error("length L gets too small (or resolution too high), running this program will crash")
    end
    y = linspace(0,int_density(end),num_points); % choose equidistant array with wished number of points in transformed space
    t_after = interp1(int_density,t_after,y,'spline'); % inverse of y is an array with given density 1/L
    t = [t, t_after];
end

t = unique([t, t_breakpoints]);

%% Make t_before array for xi < -L(0)
if kwargs.margin > 0
    % Make t_before array
    t_before = linspace(-kwargs.margin,0,ceil(1000*kwargs.margin)); % initial equidistant array
    dt_before = (t_before(end)-t_before(1))/(length(t_before)-1);
    int_density = cumtrapz(dt_before,1./L(t_before)); % integral of density, L is efficiently evaluated in not many points
    num_points = ceil( int_density(end)*kwargs.resolution );
    if num_points > 10^8
        error("length L gets too small (or resolution too high), running this program will crash")
    end
    y = linspace(0,int_density(end),num_points); % choose equidistant array with wished number of points in transformed space
    t_before = interp1(int_density,t_before,y,'spline'); % inverse of y is an array with given density 1/L
    t_before = unique([t_before, t_breakpoints_before]);
end

%% Make R
L_in_t = L(t);

% First, find all values of R at t = 0
xi_0 = linspace(-L0, L0, 2*kwargs.resolution);
% the density of points in the time array is rho/L0 at t=0
% so we choose rho/L0 * (2L0) = 2 rho points
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
% Copy the values of R at previously calculated points and add 2
for i = 1:length(t_breakpoints)
    if i == length(t_breakpoints)
        t_indices = t_breakpoints(i) <= t;
    else
        t_indices = t_breakpoints(i) <= t & t <= t_breakpoints(i+1);
    end
    print_indented_message("Spline object " + i + ": from index " + find(t_indices,1,"first") + " to " + find(t_indices,1,"last"))
    t_arr = t(t_indices);
    if numel(t_arr) > 1
        L_arr = L_in_t(t_indices);
        new_R = 2 + R_PP(t_arr-L_arr);
        new_R_PP = interpolant(t_arr+L_arr, new_R, kwargs.method);
        common_break = R_PP.PP.breaks(end);
        if abs(new_R_PP(common_break) - R_PP(common_break)) > 10^(-10)
            error('Interpolation objects are not consistent (this means that R is discontinuous)')
        end
        % Merge pp objects between regions
        PP = mkpp([R_PP.PP.breaks, new_R_PP.PP.breaks(2:end)], [R_PP.PP.coefs; new_R_PP.PP.coefs]);
        R_PP = interpolant(PP);
    end
end

%% Make R_before
if kwargs.margin > 0
    L_in_t_before = L(t_before);

    % Copy the values of R at previously calculated points and add 2
    for i = 1:length(t_breakpoints_before)
        if i == length(t_breakpoints_before)
            t_indices = t_before <= t_breakpoints_before(i);
        else
            t_indices = t_breakpoints_before(i+1) <= t_before & t_before <= t_breakpoints_before(i);
        end
        print_indented_message(find(t_indices,1,"last") + " to " + find(t_indices,1,"first"))
        t_arr = t_before(t_indices);
        if numel(t_arr) > 1
            L_arr = L_in_t_before(t_indices);
            new_R = -2 + R_PP(t_arr+L_arr);
            new_R_PP = interpolant(t_arr-L_arr, new_R, kwargs.method);
            % Merge pp objects between regions
            PP = mkpp([new_R_PP.PP.breaks(1:end-1), R_PP.PP.breaks], [new_R_PP.PP.coefs; R_PP.PP.coefs]);
            R_PP = interpolant(PP);
        end
    end
end

print_indented_message("Finished in " + string(toc(timer)))

end