%% This is the main script for the simulations and the implementation of the models.
%------------------------------------------------------------------------------------%

addpath('./')
set(0,'defaultfigureposition',[488,242,560*3/4,420*3/4])
set(0, 'DefaultLineLineWidth', 1);
set(0, 'DefaultLineMarkerSize', 5);
warning('on','backtrace')

% Set or get the silent flag (a persistent variable). 
% The silent flag is used in print_indented_message to log messages to the console.
% The silent flag can be:
%   0 (default): print_indented_message prints all messages that it receives to the console
%   1: print_indented_message only prints the bold hyperlinked function names to the console
%   2: print_indented_message is completely silent and prints nothing to the console
silent_flag(0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Moore series: convergence for linear L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the coefficients c_l
l_arr = 0:100;
c_l = ones(size(l_arr)); % c_0 = 1
for l = l_arr(2:end)
    c_l(l + 1) = sum(-c_l( l-(1:l) + 1) ./ (2*(1:l)+1));
end

% Plot the coefficients c_l
figure
plot(l_arr,c_l)
xlabel('$l$', interpreter='latex'), ylabel('$c_l$', interpreter='latex')
xlim([0,15])

% Plot the absolute value of the coefficients c_l
figure
semilogy(l_arr, abs(c_l))
xlabel('$l$', interpreter='latex'), ylabel('$|c_l|$', interpreter='latex')

% Plot the compliance of the boundary condition for 3 values of the velocity v
v_arr = [0.5, 0.75, 0.9];
figure
linespec = ["-.","--","-"];
for i = 1:numel(v_arr)
    v = v_arr(i);
    diffalpha = cumsum(c_l .* v.^(2*(l_arr)-1));
    semilogy(l_arr, abs(diffalpha * log((1+v)/(1-v)) - 2), linespec(i) ), hold on
    ylim([10^(-15),1])
end
xlabel('$n$', interpreter='latex')
ylabel('$|\sum_{l=0}^n \left( \alpha_{l}(t+L(t)) - \alpha_{l}(t-L(t)) \right)-2|$', interpreter='latex')
legend("$v = " + string(v_arr) + "$", interpreter='latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Moore series: convergence for exponential L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the coefficients c_l
l_arr = 0:80;
c_l = ones(size(l_arr)); % c_0 = 1
for l = l_arr(2:end)
    terms = - (2*l-2*(1:l)-1).^(2*(1:l)) ./ factorial(2*(1:l)+1) .* c_l(l-(1:l) + 1);
    c_l(l + 1) = sum(terms);
end

% Calculate the coefficients ctilde_l
ctilde_l = 10.^( log10(abs(c_l)) - (2*l_arr-1) .* log10(abs(2*l_arr-1)) );

% Plot the coefficients c_l and add a tangent line
figure
semilogy(l_arr, abs(c_l))
xlabel('$l$', interpreter='latex'), ylabel('$|c_l|$', interpreter='latex')
slope = log10(abs(c_l(42))) - log10(abs(c_l(41)));
interc = log10(abs(c_l(41))) - 40*slope;
hold on, semilogy(l_arr, 10.^(l_arr*slope+interc),'--')
legend("$|c_l|$", "tangent line at $l = 40$", interpreter='latex')

% Plot the coefficients ctilde_l
figure
semilogy(l_arr, ctilde_l)
xlabel('$l$', interpreter='latex'), ylabel('$|\tilde{c}_l|$', interpreter='latex')

% Plot the terms alpha_l for xi = -1 and 1
k = 0.25;
l_arr2 = 0:20;
c_l2 = c_l(l_arr2+1);
figure(Position=[488,242,560*3/4,420*3/4])
alph1 = c_l2 .* k.^(2*l_arr2-1) .* exp((1-2*l_arr2)*k) ./ (1-2*l_arr2);
alph2 = c_l2 .* k.^(2*l_arr2-1) .* exp(-(1-2*l_arr2)*k) ./ (1-2*l_arr2);
semilogy(l_arr2, abs(alph1)), hold on, semilogy(l_arr2, abs(alph2))
xlabel('$l$', interpreter='latex'), ylabel('$|\alpha_{l}|$', interpreter='latex')
legend("$|\alpha_{l}(" + [1, -1] + ")|$", interpreter='latex')

% Plot the compliance of the boundary condition for 4 values of k and 2 values of time t
figure
k_arr = flip([0.1,0.25,0.5,1]);
time = [0,0.5];
col = colororder;
col = flip([col(4,:) * 0.7; col(1:3,:)]); % Custom colors
linespec = ["-", "--"];
for j = 1:numel(time)
    t = time(j);
    for i = 1:numel(k_arr)
        k = k_arr(i);
        diffalph = -c_l2 .* k.^(2*l_arr2-1) .* exp((1-2*l_arr2)*k*t) .* 2.*sinh((1-2*l_arr2)*k*exp(-k*t)) ./ (2*l_arr2-1);
        semilogy(l_arr2,abs(cumsum(diffalph)-2),linespec(j),Color=col(i,:)), hold on
    end
end
ylim([10^(-15), 10^15])
xlabel('$n$', interpreter='latex')
ylabel('$|\sum_{l=0}^n \left( \alpha_{l}(t+L(t)) - \alpha_{l}(t-L(t)) \right)-2|$', interpreter='latex')
legend("$k=" + k_arr + "$", interpreter='latex')
add_legend("$t="+time+"$", ["k-","k--"]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Error on initial conditions, using coefs error (error bound on amplitude)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_from_file = true;
switch load_from_file
    case true
        load('err_initial_conditions.mat')
    case false
        n_coefs = floor(10.^(linspace(log10(4),log10(1000),40)));
        boundary_conditions = "sinhR";
        initial_conditions = ["gauss", "sin"];

        % initialisation and figures
        err = zeros([length(n_coefs), 3, ...
                    length(initial_conditions), ...
                    length(boundary_conditions)]);

        
        for j = 1:length(boundary_conditions)    
            [L,R] = get_LR(boundary_conditions(j));
            Rint = find_R_interpolation(L,initialR='cubic',t_max=0);
            Rmoore = find_R_moore(L,Kmax=4,t_max=5);

            for i = 1:length(initial_conditions)
                % get parameters for boundary+startcondition 
                [f,g] = get_fg(L,initial_conditions(i));
        
                [err(:,1,i,j),~] = coefs_error(f,g,R,L(0),n_coefs);
                [err(:,2,i,j),~] = coefs_error(f,g,Rint,L(0),n_coefs);
                [err(:,3,i,j),~] = coefs_error(f,g,Rmoore,L(0),n_coefs);
            end
        end
end
linespec = ["r-s","g-o","b-v"; 
            "r--s","g--o","b--v"];
figure
length_n = length(n_coefs);
for i = 1:length(initial_conditions)
    loglog(n_coefs,err(:,1,i,1,1), linespec(i,1), MarkerIndices=i:6:length_n),  hold on
    loglog(n_coefs,err(:,2,i,1,1), linespec(i,2), MarkerIndices=i+2:6:length_n),  
    loglog(n_coefs,err(:,3,i,1,1), linespec(i,3), MarkerIndices=i+4:6:length_n)
end
%title("Error bound on amplitude", interpreter='latex');
xlabel("Number of coefficients", interpreter='latex');ylabel("Error", interpreter='latex')
grid on

leg1 = add_legend(["Exact $R$","Interpolation","Moore"],linespec(1,:));
leg2 = add_legend(["Gaussian", "Sine"],["k","k--"]);
title(leg1, "Method", Interpreter='latex');
title(leg2, "Initial condition", Interpreter='latex')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Error on boundary condition as function of resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_from_file = true;
switch load_from_file
    case true
        load('err_boundary_conditions')
    case false
        As = [2,0.1];
        reso = 10.^linspace(1.5,4.5,50);
        t_max = 4;
        tt = linspace(0,t_max,179416);

        err_interpol = zeros(length(reso), numel(As));
        err_char = zeros(length(reso), numel(As));
        for j = 1:numel(As)
            L = get_LR('sinhR',A=As(j));
            [f,g] = get_fg(L, 'gauss');
            for i = 1:length(reso)
                R = find_R_interpolation(L, t_max=t_max, resolution=reso(i));
                w = find_w_interpolation(f,g,L, t_max=t_max, resolution=reso(i));
                err_interpol(i,j) = rmse(R(tt+L(tt))-R(tt-L(tt)), 2); % Root Mean Squared Error
                err_char(i,j) = rmse(w(tt+L(tt))-w(tt-L(tt)), 0);
            end
        end
end

figure
col = colororder;
col = [col(1,:)*0.7; col(2,:)*1.15];
xlabel('Resolution $\rho$', Interpreter='latex')

colororder(col)
yyaxis left
loglog(reso, err_interpol(:,1), Color=col(1,:)), hold on
loglog(reso, err_interpol(:,2), '--', Color=col(1,:)), hold off
ylabel('rms $|R(t+L(t))-R(t-L(t))-2|$', Interpreter='latex')
ylim([1e-15,5e-3])

yyaxis right
loglog(reso, err_char(:,1), Color=col(2,:)), hold on
loglog(reso, err_char(:,2), '--', Color=col(2,:))
ylabel('rms $|w(t+L(t))-w(t-L(t))|$', Interpreter='latex')
ylim([1e-15,5e-3])
xlim('tight')
grid on

add_legend("$A = " + As + "$",["k-","k--"]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Error on boundary condition as a function of velocity (for Moore)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_from_file = true; % Load from file, this takes really long to run
switch load_from_file
    case true
        load('err_boundary_conditions_moore')
    case false
        t_max = 4;
        tt = linspace(0,t_max,17416);
        K_cell = {[2,3,4,5], [5,10,25,100], [5,10,15]};
        v_cell = {10.^linspace(log10(0.05),log10(0.95),20).', 10.^linspace(log10(0.05),log10(0.95),200).', ...
            10.^linspace(log10(0.05),log10(0.95),200).'}; % Choose array of max velocities
        err_cell = arrayfun(@(i) zeros(numel(v_cell{i}),numel(K_cell{i})), 1:numel(K_cell), UniformOutput=false); % sinhR, linearL and expL
        
        %%% Calculations for sinhR
        % Find relation between parameter A en max velocity
        L = get_LR('sinhR',A=sym('A'));
        dL = derivative(L);
        A_arr = 10.^linspace(-3,2,10000);
        I = interpolant( arrayfun(@(a) abs(dL(fminsearch(@(t) -abs(dL(t,a)),0),a)), A_arr) , A_arr );
        
        K_arr = K_cell{1}; v = v_cell{1};
        A_values = I(v); % Corresponding parameter values
        % Find R for each parameter K and A
        R_moo_sinhR = cell(length(K_arr),length(A_values));
        for i = 1:length(A_values) 
            L = get_LR('sinhR',A=A_values(i));
            R_moo_sinhR(i,:) = find_R_moore_vpa(L,K_array=K_arr,t_max=t_max);
            for j = 1:length(K_arr)
                R = R_moo_sinhR{i,j};
                err_cell{1}(i,j) = rmse(R(tt+L(tt))-R(tt-L(tt)), 2); % Root Mean Squared Error
            end
        end
        
        %%% Calculations for linearL
        L = get_LR('linearL',v=sym('v'));
        K_arr = K_cell{2}; v = v_cell{2};
        R_moo_linearL = find_R_moore(L,K_array=K_arr,param_array=v,t_max=t_max); 
        for i = 1:length(v)
            L = get_LR('linearL',v=v(i));
            for j = 1:length(K_arr)
                R = R_moo_linearL{i,j};
                err_cell{2}(i,j) = rmse(R(tt+L(tt))-R(tt-L(tt)), 2); % Root Mean Squared Error
            end
        end
        
        %%% Calculations for expL
        L = get_LR('expL',k=sym('k')); % max speed at t=0
        K_arr = K_cell{3}; v = v_cell{3};
        R_moo_expL = find_R_moore(L,K_array=K_arr,param_array=v,t_max=t_max);
        for i = 1:length(v)
            L = get_LR('expL',k=v(i));
            for j = 1:length(K_arr)
                R = R_moo_expL{i,j};
                err_cell{3}(i,j) = rmse(R(tt+L(tt))-R(tt-L(tt)), 2); % Root Mean Squared Error
            end
        end
end
figure
% Here follows our finest selection of colors
bluepurple = [0.494,0.184,0.556; 0.062,0.258,0.501; 0,0.447,0.741; 0.329,0.713,1];
yellowred = [0.929,0.694,0.125; 0.96,0.466,0.16; 0.717,0.192,0.172; 0.372,0.105,0.031];
green = [0.286,0.858,0.25; 0.231,0.666,0.196; 0.007,0.345,0.054];
colors = {bluepurple, yellowred, green};
textlocations = {[0.055 1.34e-07; 0.055 4.09e-09; 0.053 1.27e-10; 0.061 3.65e-12], ...
                [0.337 2.07e-08; 0.394 4.33e-12; 0.595 6.50e-11; 0.749 6.50e-13], ...
                [0.366 5.71e-06; 0.378 0.040; 0.279 0.409]};
linespec = ["--","-","-."];

% Plot
for type = 1:3 % sinhR, linearL, expL
    K_arr = K_cell{type};
    for j = length(K_arr):-1:1
        col = colors{type}(j,:);
        loglog(v_cell{type},err_cell{type}(:,j),linespec(type),Color=col), hold on
        text(textlocations{type}(j,1), textlocations{type}(j,2), num2str(K_arr(j)),Interpreter="latex",Color=col)
    end
end

xlabel('Maximum $|\dot{L}|$', interpreter='latex'), ylabel('rms $|R(t+L(t))-R(t-L(t))-2|$', interpreter='latex')
ylim([10^-15,10^3]), grid on

leg = add_legend(["Hyperbolic sine $R$","Linear $L$","Exponential $L$"], linespec+"k");
set(leg, Position=[0.1586,0.7489,0.3468,0.1529])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RMSE error as function of time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_from_file = true;
switch load_from_file
    case true
        load('rms_error.mat')
    case false
        nx = 500;
        t_max = 5;
        time_arr = 0:1/120:t_max;
        boundary_conditions = ["linearL","sinhR"];
        initial_conditions = ["gauss","sin"];
        err = cell(length(boundary_conditions), length(initial_conditions));

        % make all permutations of boundary and initial conditions
        for j = 1:length(boundary_conditions)
            [L,R] = get_LR(boundary_conditions(j));
            dR = derivative(R);

            R_int = find_R_interpolation(L,t_max=t_max);
            R_moo = find_R_moore(L,t_max=t_max);
            for i = 1:length(initial_conditions)
                [f,g] = get_fg(L,initial_conditions(i));
                df = derivative(f);
                coefs = calc_coefs(df,g,R,L(0),300);
                coefs = reshape(coefs, [1,1,numel(coefs)]);
                n = reshape(1:numel(coefs), [1,1,numel(coefs)]);
                % Change f and g to finite combination of eigenmodes:
                f = @(x) 2*real( sum(coefs .* (exp(-1i*pi*n .* R(-x)) - exp(-1i*pi*n .* R(x))), 3) );
                g = @(x) 2*real( sum(coefs .* ( -1i*pi*n.*dR(-x).*exp(-1i*pi*n.*R(-x)) - (-1i)*pi*n.*dR(x).*exp(-1i*pi*n.*R(x)) ), 3) );
                
                u_ana = eigenmode_solver(f,g,R,L,nx=nx,t=time_arr, coefs=coefs);
                u_int = interpolation_solver(f,g,L,nx=nx,t=time_arr,n_max=200,R=R_int);
                u_cha = characteristics_solver(f,g,L,nx=nx,t=time_arr);
                u_moo = moore_solver(f,g,L,nx=nx,t=time_arr,n_max=200,R=R_moo);

                % RMSE error
                error_int = rms_error(u_int,u_ana);
                error_cha = rms_error(u_cha,u_ana);
                error_moo = rms_error(u_moo,u_ana);
        
                err{j,i} = {error_int, error_cha, error_moo}; %, error_bac};
            end
        end
end

linespec = ["-s","-o","-v"; 
            "--s","--o","--v"];
length_t = length(time_arr);

figs = cell(1,size(err,1));
for j = 1:size(err,1)
    figs{j} = figure;
    for i = 1:size(err,2)
        error_int = err{j,i}{1};
        error_cha = err{j,i}{2};
        error_moo = err{j,i}{3};
        
        col = colororder;
        col = col(1:4,:);
        semilogy(time_arr,error_int,linespec(i,1),Color=col(1,:),MarkerIndices=1+(i-1)*30:60:length_t), hold on,
        semilogy(time_arr,error_cha,linespec(i,2),Color=col(2,:),MarkerIndices=11+(i-1)*30:60:length_t),
        semilogy(time_arr,error_moo,linespec(i,3),Color=col(3,:),MarkerIndices=21+(i-1)*30:60:length_t),
    end
    ylim([10^-15,10^-1])
    xlabel('Time',Interpreter='latex');
    ylabel("Root mean squared error",Interpreter='latex');
end
figure(figs{1})

leg1 = add_legend(["Interpolation","Characteristics","Moore"], linespec(1,:), col);
leg2 = add_legend(["Gaussian","Sine"], ["k-","k--"]);
title(leg1, "Method", Interpreter='latex')
title(leg2, "Initial condition", Interpreter='latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run time: R_interpolation vs R_backtrace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_from_file = true;
switch load_from_file
    case true
        load('R_compute_time.mat')
    case false
        old_flag = silent_flag(2); % log nothing to the console
        
        t_max_arr = [0.5,5];
        L = get_LR('sinhR');
        dL = derivative(L);
        ddL = derivative(dL);
        derivL0 = [dL(0), ddL(0)];
        
        evaluate = @(f,x) f(x);
        find_R = {@(t_max,xi) evaluate(find_R_interpolation(L,t_max=t_max,derivL0=derivL0), xi), ...
                    @(t_max,xi) evaluate(find_R_backtrace(L,derivL0=derivL0), xi)};
        compute_times = cell(numel(t_max_arr),2);
        
        for i = 1:numel(t_max_arr)
            inputsize = 10.^linspace(1,8,50);
            compute_time = zeros(size(inputsize));
            t_max = t_max_arr(i);
            figure
            for j = 1:2
                for k = 1:numel(inputsize)
                    xi = linspace(t_max-L(t_max),t_max+L(t_max),inputsize(k));
                    f = find_R{j};
                    f(0,0); tic, [~] = toc; tic, [~] = toc; % Warm up
                    tic, f(t_max,xi); t_one = toc;
                    num = 1; t_tot = t_one;
                    while t_tot < 0.5
                        num = ceil(1/t_one); % run for roughly 1 second
                        tic
                        for counter = 1:num
                            f(t_max,xi);
                        end
                        t_tot = toc;
                        t_one = t_tot/num;
                    end
                    compute_time(k) = t_one;
                    disp(k)
                    if t_one > 10, break, end
                end
                compute_times{i,j} = {inputsize(1:k),compute_time(1:k)};
            end
        end
        silent_flag(old_flag); % set silent flag back
end

figure
linespec = ["-","--"];
col = colororder;
col(1,:) = col(1,:)*0.7; % DARK
col(2,:) = col(2,:)*1.15; % LIGHT
for i = 1:numel(t_max_arr)
    for j = 1:2
        loglog(compute_times{i,j}{1},compute_times{i,j}{2},linespec(i),Color=col(j,:)), hold on, grid on
        xlim([10,10^8])
        ylim([10^(-4),10])
    end
end
xlabel('Input size of $\xi$', Interpreter='latex')
ylabel('Computation time (s)', Interpreter='latex')

add_legend(["Interpolation", "Backtrace"], col);
add_legend("$t_0 = " + t_max_arr + "$", linespec+"k")