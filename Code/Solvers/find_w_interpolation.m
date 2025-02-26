function [w_PP] = find_w_interpolation(f,g,L,kwargs)
%FIND_W_INTERPOLATION Finds the function w on [-L(0), tmax+L(tmax)] from L using the characteristics method
%Based on w(t+L(t)) = w(t-L(t))
% INPUTS:
%   f: function_handle, inital position
%   g: function_handle, initial velocity
%   L: function_handle or interpolant object, length of the domain
%   optional t_max: double, maximum time to calculate w for
%   optional resolution_t: double, number of points per time in the first interpolation t-L(t) -> t
%   optional resolution: double, number of points in region xi = -L0 to xi = L0 in the final interpolation object
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
    kwargs.resolution_t {mustBePositive} = 10000
    kwargs.resolution {mustBePositive} = 200000
    kwargs.method = "spline"
end

print_indented_message("Find w Interpolation",true)
timer = tic;

t_max = kwargs.t_max;
resolution_t = round(kwargs.resolution_t);
points_xi = round(kwargs.resolution); % between 0 and L(0)
L0 = L(0);

%% Initial domain
% First, find all values of w at t = 0
index0 = floor(points_xi/2);
isodd = points_xi - 2*index0;
xi_0 = linspace(-L0,L0,points_xi);
xi_0_pos = xi_0(index0+1:end);

opts = odeset(RelTol=3e-14,AbsTol=1e-15);
[~,U0] = ode45(@(t,y) g(t), xi_0_pos, 0, opts);
U0 = reshape(U0, size(xi_0_pos));
f0 = f(xi_0_pos);
w_0 = [f0(end:-1:1+isodd) - U0(end:-1:1+isodd), -f0 - U0]/2;

w_PP = interpolant(xi_0,w_0,kwargs.method);
if t_max == 0
    return
end
sz = size(w_PP.PP.coefs,2);

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

%% Make w
t = 0:1/resolution_t:(t_max+1/resolution_t);
I = interpolant(t-L(t),t,kwargs.method);

w_PP_coefs = cell(1,1+num_breakpoints); w_PP_coefs{1} = w_PP.PP.coefs;
w_PP_breaks = cell(1,1+num_breakpoints); w_PP_breaks{1} = w_PP.PP.breaks;
% Copy the values of w at previously calculated points
xi_I_max = t_max-L(t_max);
xi_new = xi_0;
complete = false;
for i = 1:num_breakpoints
    % find t_new
    if xi_new(end) > xi_I_max
        xi_new = xi_new(xi_new <= xi_I_max);
        if xi_new(1) >= xi_I_max
            break % we are done
        elseif numel(xi_new) <= 3 % too little points spaced far enough
            xi_new = linspace(xi_new(1),xi_I_max,4);
        elseif xi_new(end)-xi_new(end-1) < 2*(xi_I_max-xi_new(end))
            xi_new = [xi_new, xi_I_max]; % add new point
        else
            xi_new(end) = xi_I_max; % move point
        end
    end
    t_new = I(xi_new);
    print_indented_message("Spline object " + i + ": with " + length(t_new) + " points")
    t_new(1) = t_breakpoints(i); 
    % Force first time point on the breakpoint, since
    % by interpolation errors, it may not be exactly t_breakpoints(i)
    if i == num_breakpoints || min(t_max - t_breakpoints(i+1), t_max - t_new(end)) < 20*eps(t_max)
        t_new(end) = t_max;
        complete = true;
    else
        t_new(end) = t_breakpoints(i+1);
        % Force last time point on the breakpoint
    end

    L_new = L(t_new);
    new_w = w_PP(t_new-L_new);
    xi_new = t_new+L_new;
    new_w_PP = interpolant(xi_new, new_w, kwargs.method);
    common_break = w_PP.PP.breaks(end);
    if abs(new_w_PP(common_break) - w_PP(common_break)) > 10^(-10)
        error('Interpolation objects are not consistent (this means that w is discontinuous)')
    end
    w_PP = new_w_PP;
    % Merge
    % w_PP_coefs{i+1} = [zeros(new_sz(1), sz(2)-new_sz(2)), new_w_PP.PP.coefs];
    w_PP_coefs{i+1} = new_w_PP.PP.coefs;
    w_PP_breaks{i+1} = new_w_PP.PP.breaks(2:end);

    if complete, break, end
end

% Merge w's
w_PP = interpolant(mkpp(horzcat(w_PP_breaks{:}),vertcat(w_PP_coefs{:})));

print_indented_message("Finished in " + string(toc(timer)))

end