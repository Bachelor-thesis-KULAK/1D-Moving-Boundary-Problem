function [w_PP] = find_w_interpolation(f,g,L,kwargs)
%FIND_W_INTERPOLATION Finds the function w on [-L(0), tmax+L(tmax)] from L using the characteristics method
%Based on w(t+L(t)) = w(t-L(t))
% INPUTS:
%   f: function_handle, inital position
%   g: function_handle, initial velocity
%   L: function_handle or interpolant object, length of the domain
%   optional t_max: double, maximum time to calculate R for
%   optional resolution: double, the density of time points per unit of 1/L
%   optional method: 'linear', 'cubic' or 'spline' (= default), the method of interpolation
%   optional margin: if specified, the function R will be calculated from -L(0)-margin to tmax+L(tmax)+margin (a little outside the domain specified)
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
    kwargs.t_max = 5
    kwargs.resolution = 50000
    kwargs.method = "spline"
    kwargs.margin = 0
end

print_indented_message("Find w interpolation",true)
timer = tic;

t_max = kwargs.t_max;
kwargs.resolution = round(kwargs.resolution);

%% Find time breakpoints
% t0 = 0
% t1 such that t1-L(t1) = t0+L(t0)
% t2 such that t2-L(t2) = t1+L(t1) ...
t_breakpoints = 0;
t_maxmarg = t_max + kwargs.margin;
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
if ceil(1000*t_max) < 2
    t = 0;
else
    t = linspace(0,t_max,ceil(1000*t_max)); % initial equidistant array
    dt = (t(end)-t(1))/(length(t)-1);
    int_density = cumtrapz(dt,1./L(t)); % integral of density, L is efficiently evaluated in not many points
    num_points = ceil( int_density(end)*kwargs.resolution );
    if num_points > 10^8
        error("The resolution is too high (or the length L is too small), too many points")
    end
    y = linspace(0,int_density(end),num_points); % choose equidistant array with wished number of points in transformed space
    t = interp1(int_density,t,y,'spline'); % inverse of y is an array with given density 1/L
end

%% Make t_after array for xi > tmax+L(tmax)
if kwargs.margin > 0
    % Make t_after array
    t_after = linspace(t_max,t_maxmarg,ceil(1000*kwargs.margin)); % initial equidistant array
    dt_after = (t_after(end)-t_after(1))/(length(t_after)-1);
    int_density = cumtrapz(dt_after,1./L(t_after)); % integral of density, L is efficiently evaluated in not many points
    num_points = ceil( int_density(end)*kwargs.resolution );
    if num_points > 10^8
        error("The resolution is too high (or the length L is too small), too many points")
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

%% Make w
L_in_t = L(t);

% First, find all values of R at t = 0
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
% Copy the values of w at previously calculated points
for i = 1:length(t_breakpoints)
    if i == length(t_breakpoints)
        t_indices = t_breakpoints(i) <= t;
    else
        t_indices = t_breakpoints(i) <= t & t <= t_breakpoints(i+1);
    end
    t_arr = t(t_indices);
    if numel(t_arr) > 1
        print_indented_message("Spline object " + i + ": from index " + find(t_indices,1,"first") + " to " + find(t_indices,1,"last"))
        L_arr = L_in_t(t_indices);
        new_w = w_PP(t_arr-L_arr);
        new_w_PP = interpolant(t_arr+L_arr, new_w, kwargs.method);
        if size(w_PP.PP.coefs,2) ~= size(new_w_PP.PP.coefs,2)
            error('Resolution is too small (too little points in this region)')
        end
        % Merge pp objects between regions
        PP = mkpp([w_PP.PP.breaks, new_w_PP.PP.breaks(2:end)], [w_PP.PP.coefs; new_w_PP.PP.coefs]);
        w_PP = interpolant(PP);
    end
end

%% Make w_before
if kwargs.margin > 0
    L_in_t_before = L(t_before);
    
    % Copy the values of w at previously calculated points
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
            new_w = w_PP(t_arr+L_arr);
            new_w_PP = interpolant(t_arr-L_arr, new_w, kwargs.method);
            if size(w_PP.PP.coefs,2) ~= size(new_w_PP.PP.coefs,2)
                error('Resolution is too small (too little points in this region)')
            end
            % Merge pp objects between regions
            PP = mkpp([new_w_PP.PP.breaks(1:end-1), w_PP.PP.breaks], [new_w_PP.PP.coefs; w_PP.PP.coefs]);
            w_PP = interpolant(PP);
        end
    end
end

print_indented_message("Finished in " + string(toc(timer)))

end