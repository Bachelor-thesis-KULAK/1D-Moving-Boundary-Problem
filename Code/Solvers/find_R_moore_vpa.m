function [R] = find_R_moore_vpa(L,kwargs)
%FIND_R_MOORE Finds the function R on [-L(0), tmax+L(tmax)] using Moore's method and numerical derivatives with high precision vpa
%By default the optimal number of gamma_l's is sought and selected
% INPUTS:
%   L: function_handle or interpolant object, length of the domain
%   optional t_max: double, maximum time to calculate R for
%   optional density: double, density of xi-points in the interpolant object
%   optional Kmax: maximum number of gamma_l's to be calculated (more than 4 is not feasible except for simple examples)
%   optional margin: if specified, the function R will be calculated from -L(0)-margin to tmax+L(tmax)+margin (a little outside the domain specified)
%   optional K_array: if specified, ignore the minimization process of the error and instead
%           return a cell array with an interpolant R for each value Kmax in the K_array
% OUTPUTS:
%   R: interpolant object, returns the function R
% EXAMPLES:
%   [L,R] = get_LR("linearL")
%   R = find_R_moore(L,Kmax=4) 
%       finds the optimal number of terms (between 1 and 5) and returns an R with the optimal number
%   R = find_R_moore(L,K_array=[1,2,3,4])
%       returns a 1x4 cell array with an R for Kmax = 1, 2, 3 and 4
arguments
    L
    kwargs.t_max = 5
    kwargs.density = 1000
    kwargs.Kmax = 6
    kwargs.margin = 0
    kwargs.K_array = []
end

print_indented_message("Find R Moore",true)
timer = tic;
reset(symengine)
digits(2000)
min_xi = -vpa(kwargs.margin) - L(vpa(0)-vpa(kwargs.margin));
max_xi = vpa(kwargs.t_max)+vpa(kwargs.margin)+L(vpa(kwargs.t_max)+vpa(kwargs.margin));
xi0 = linspace(min_xi,max_xi,kwargs.density*(max_xi-min_xi));

if ~isempty(kwargs.K_array)
    Kmax = max(kwargs.K_array,[],"all");
else
    Kmax = kwargs.Kmax; % make variable
end

d = vpa(10)^(-20);
xi = cell(1,2*Kmax+1);
Lxi = cell(size(xi));
for i = 1:numel(xi)
    xi{i} = xi0 + d*((-Kmax)+(i-1));
    Lxi{i} = L(xi{i});
end

tt = linspace(0,kwargs.t_max,1e5); % test array for minimization

gamma = cell(1,Kmax+1);
derivs = cell(1,Kmax+1);
gamma{1} = cellfun(@(x) 1./x, Lxi, 'UniformOutput', false);
dR = 0;

if ~isempty(kwargs.K_array)  % store the dR
    dR_arr = cell(1,Kmax+1);
end

k = 0;
while k <= Kmax
    print_indented_message("Iteration " + string(k)) % print which iteration
    % Calculate next gamma_k
    for i = 1:k
        % Avoid matrix slicing, which takes ages
        deriv = derivs{k+1-i};
        gamm = gamma{k+1};
        new_deriv = cell(size(xi));
        new_gamm = cell(size(xi));
        for j = (1+k):(2*Kmax+1-k)
            new_deriv{j} = ( deriv{j-1}+deriv{j+1}-2*deriv{j} ) / d^2; % 2nd derivative
            new_gamm{j} = -Lxi{j}.^(2*i) / vpa(factorial(sym(2*i+1))) .* new_deriv{j};
            if ~isempty(gamm)
                new_gamm{j} = new_gamm{j} + gamm{j};
            end
        end
        derivs{k+1-i} = new_deriv;
        gamma{k+1} = new_gamm;
    end
    derivs{k+1} = gamma{k+1};

    dR = dR + gamma{k+1}{Kmax+1};
    dR_interp = interpolant(double(xi0),double(dR));

    if isempty(kwargs.K_array)
        R = dR_interp.integral;
        % minimize if no K_array is given
        new_error = rmse(R(tt+L(tt))-R(tt-L(tt)),2);
    
        % check for convergence
        if (k > 0 && new_error >= prev_error) || k == Kmax
            if k > 0 && new_error >= prev_error
                R = prev_R;
                best_k = k-1;
            else
                best_k = k;
            end
            print_indented_message("Best k: " + best_k)
            print_indented_message("Finished in " + string(toc(timer)))
            return
        end
        prev_error = new_error;
        prev_R = R;
    else
        dR_arr{k+1} = dR_interp;
    end

    k = k+1;
end


R = cell(size(kwargs.K_array));
for j = 1:numel(kwargs.K_array)
    dR_interp = dR_arr{kwargs.K_array(j)+1};
    R{j} = dR_interp.integral;
end

print_indented_message("Finished in " + string(toc(timer)))
end