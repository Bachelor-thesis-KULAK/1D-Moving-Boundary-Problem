function [R_PP] = find_R_interpolation(L,kwargs)
%FIND_R_INTERPOLATION Finds the transformation R on [-L(0), tmax+L(tmax)] from L using the interpolation method
%Based on R(t+L(t)) = 2+R(t-L(t))
%
% Ensures every domain (copy of R) recieves equal amount of t-points
%
% INPUTS:
%   L: function_handle or interpolant object, length of the domain
%   optional t_max: double, maximum time to calculate R for
%   optional resolution_t: double, number of points per time in the first interpolation t-L(t) -> t
%   optional resolution: double, number of points in region xi = -L0 to xi = L0 in the final interpolation object
%   optional method: 'linear', 'cubic' or 'spline' (= default), the method of interpolation
%   optional initialR: 'quadratic' or 'cubic' (= default) (or 'linear' or 'sigmoid-like'), the choice of R on [-L0,L0]
%   optional derivL0: specify the first (and second) derivative of L at t = 0 instead of calculating it
% OUTPUTS:
%   R_PP: interpolant object, returns the transformation R
% EXAMPLE:
%   [L,R] = get_LR("linearL")
%   R = find_R_interpolation(L,t_max=10,resolution_xi=10000,method='spline')

arguments
    L
    kwargs.t_max {mustBeNonnegative} = 5 
    kwargs.resolution_t {mustBePositive} = 10000
    kwargs.resolution {mustBePositive} = 20000
    kwargs.method = "spline"
    kwargs.initialR = "quadratic"
    kwargs.derivL0 = []
end

print_indented_message("Find R Interpolation",true)
timer = tic;

t_max = kwargs.t_max;
resolution_t = round(kwargs.resolution_t);
points_xi = round(kwargs.resolution); % between 0 and L(0)

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
xi_0 = linspace(-L0, L0, points_xi);
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

%% Make R
t = 0:1/resolution_t:(t_max+1/resolution_t);
I = interpolant(t-L(t),t,kwargs.method);

R_PP_coefs = cell(1,1+num_breakpoints); R_PP_coefs{1} = R_PP.PP.coefs;
R_PP_breaks = cell(1,1+num_breakpoints); R_PP_breaks{1} = R_PP.PP.breaks;
% Copy the values of R at previously calculated points and add 2
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
    new_R = 2 + R_PP(t_new-L_new);
    xi_new = t_new+L_new;
    if numel(unique(xi_new)) ~= numel(xi_new)
        error('Help')
    end
    new_R_PP = interpolant(xi_new, new_R, kwargs.method);
    common_break = R_PP.PP.breaks(end);
    if abs(new_R_PP(common_break) - R_PP(common_break)) > 10^(-10)
        error('Interpolation objects are not consistent (this means that R is discontinuous)')
    end
    R_PP = new_R_PP;
    % Merge
    % R_PP_coefs{i+1} = [zeros(new_sz(1), sz(2)-new_sz(2)), new_R_PP.PP.coefs];
    R_PP_coefs{i+1} = new_R_PP.PP.coefs;
    R_PP_breaks{i+1} = new_R_PP.PP.breaks(2:end);

    if complete, break, end
end

% Merge R's
R_PP = interpolant(mkpp(horzcat(R_PP_breaks{:}),vertcat(R_PP_coefs{:})));

print_indented_message("Finished in " + string(toc(timer)))

end