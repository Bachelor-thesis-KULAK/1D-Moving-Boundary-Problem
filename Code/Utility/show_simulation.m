function show_simulation(waves, kwargs)
%SHOW_SIMULATION Plots wave solutions at a specified fps or at specified time points
% INPUTS:
%   waves: solutionObject or cell array of solutionObject, waves to be plotted
%   optional fps: double, frames per second to be plotted (and fps saved to .mp4 if save is true)
%   optional play_ratio: double, specifies how much slower the video plays, 
%               1 time unit in the simulation corresponds to play_ratio seconds
%   optional size: 1x2 or 2x1 double, width (1st index) and height (2nd index) of the figure on which the animation plays
%   optional times: double, specify the times that needs to be plotted (instead of fps)
%   optional title: str, title of the plot
%   optional limits: double or str, specify the limits, can be 
%               "auto": fixed limits based on the max amplitude,
%               1x1 or 1x2 double: specifying the max absolute value (1x1) or limits (1x2)
%               Nx1 or Nx2 double: specifying the max absolute value or limits for each time point 
%   optional legend: str array, specify the legend titles
%   optional save: logical, saves the animation as .mp4 if true
%   optional save_title: str, name of the .mp4 if save is true
%   optional linewidth: double, specify linewidth
%   optional ordering: double, specify ordening of plotting the waves
% OUTPUTS:
%   none (optionally saves a .mp4)
%EXAMPLES:
%    1 wave:
%       plt(waves={solObj})
%       plt(waves=solObj)
%   multiple waves:
%       plt(waves={solObj1, solObj2, solObj3})
%   time array:
%       plt(waves=solObj, times=[2,3])
%   fps:
%       plt(waves=solObj, fps=1500)

arguments
    waves
    kwargs.fps double = 60
    kwargs.play_ratio double = 2 % play 2 times slower
    kwargs.size
    kwargs.times
    kwargs.title string
    kwargs.limits = "auto" % "auto", "dynamic" or specify
    kwargs.legend
    kwargs.linewidth = 1
    kwargs.ordering
    kwargs.save = false
    kwargs.save_title = "_"
end

if isempty(waves) % nothing to do
    return
end

if isa(waves, 'solutionObject') % waves can be specified as a solutionObject or as a cell array
    waves = {waves};
end

% find times to plot
if ~isfield(kwargs,'times')
    mintimes = cellfun(@(wave) min(wave.t), waves, UniformOutput=false);
    maxtimes = cellfun(@(wave) max(wave.t), waves, UniformOutput=false);
    kwargs.times = max([mintimes{:}]):(1/(kwargs.play_ratio*kwargs.fps)):min([maxtimes{:}]);
end
tmax = max(kwargs.times);
tmin = min(kwargs.times);
num_frames = length(kwargs.times);
if num_frames == 0 % nothing to do
    return
end

% trim waves to plotted times
new_waves = cellfun(@(wave) wave.reduceToTimes(kwargs.times), waves, UniformOutput=false);

% calculate x_max after trimming t values
x_max = max(new_waves{1}.x,[],"all"); 

%% Initialise figure
% Initialise figure with no menubar and no numbertitle + invisible
fig = figure(menubar="none", NumberTitle="off", Visible="off");
if isfield(kwargs,'size'), fig.Position(3:4) = kwargs.size; end
movegui(fig,'center')
set(fig, color='w') % make the background color fully white

% ax is axis for plotting
ax = axes(fig, Unit="normalized");
box on % make sure there is a full border around the figure
hold(ax, "on")

if num_frames == 1  % snapshot
    ax2 = gca;
else
    % ax2 is invisible axis for the colorbar (indicating the time)
    ax2 = copyobj(ax, fig); % copy the whole axis, to disconnect the colormap and the other plots (faster)
    set(ax2, Visible="off")
    % Set a colorbar
    cb = colorbar(ax2);
    colormap(ax2,[1,1,1]) % Set color of colormap to black
    ylabel(cb,"Time in simulation",interpreter='latex')
    num_ticks = length(cb.Ticks);
    tickLabels = string(round(linspace(tmin,tmax,num_ticks),1));
    cb.TickLabels = tickLabels;
    % Set positions
    ax_position_after_colorbar = get(ax2, 'Position');
    cb_position = get(cb, 'Position');
    new_width = 0.025;
    cb_position(3) = new_width;
    set(cb, 'Position', cb_position);
    set(ax, 'Position', ax_position_after_colorbar);
    set(ax2, 'Position', ax_position_after_colorbar);
end

if isfield(kwargs,"legend")
    % ax3 is axis for the legend
    ax3 = copyobj(ax, fig); % copy the whole axis, to disconnect the legend and the other plots
    set(ax3, Visible="off")
    hold(ax3, "on")
end

% Set xlabel, ylabel and title
xlabel(ax,"Position",interpreter='latex')
ylabel(ax,"Amplitude",interpreter='latex')
if isfield(kwargs,"title")
    title(ax,kwargs.title, interpreter = 'latex')
end

% Set xlimit and set ylimit for first frame
xlim(ax,[0,x_max])

yLimChangeIndices = zeros(1,num_frames); % per frame: 1 if new limits should be set, 0 otherwise
if class(kwargs.limits) == "char" || class(kwargs.limits) == "string"
    if kwargs.limits == "auto"
        maxValues = cellfun(@(wave) max(abs(wave.u(:))), new_waves);
        MAX = max(maxValues);
        ylim(ax,[-MAX, MAX])
    else
        error("Unknown limit mode")
    end
else
    if size(kwargs.limits,2) == 1  % only absolute value is specified
        kwargs.limits = [-kwargs.limits, kwargs.limits];
    end
    if size(kwargs.limits,1) == 1
        ylim(ax,kwargs.limits)
    elseif size(kwargs.limits,1) == numel(kwargs.times)
        maxValues = kwargs.limits;
        yLimChangeIndices = [1; any(diff(maxValues,1,1) ~= 0,2)];
        ylim(ax,maxValues(1,:))
    else
        error("The parameter 'limits' does not have the correct dimension" )
    end
end

%% Plot first frame
% First frame: initialize line objects
lineArray = gobjects(length(new_waves), 1); % initialize array without type
for j = 1:length(new_waves)
    wave = new_waves{j};
    lineArray(j) = plot(ax, wave.x(:,1), wave.u(:,1), lineWidth=kwargs.linewidth);
    if isfield(kwargs,"legend") && num_frames ~= 1
        plot(ax3, NaN, lineWidth=kwargs.linewidth) % also plot on ax3 for the legend
    end
end

% Mofify ordening of line objects
if isfield(kwargs,'ordering')
    lines = get(ax,"Children");
    set(ax,"Children", lines(kwargs.ordering))
end
if isfield(kwargs,"legend")
    legend(ax3,kwargs.legend,'AutoUpdate','off')
end

if num_frames == 1
    return % No other things to plot
end

% Plot the colorbar (indicating the time)
resolution_bar = 100;
cm = ones(resolution_bar+1,3);
current_t_fun = @(T) round( (T-tmin) / ((tmax-tmin)/resolution_bar) )+1;
current_t = current_t_fun(kwargs.times(1));
cm(1:current_t,:) = 0; % black
colormap(ax2,cm)

set(fig,Visible="on")

%% Videowriter
% Initiate video file
if kwargs.save 
    if isfield(kwargs,'save_title')
        file_name = kwargs.save_title;
    else
        file_name = "_";
    end
    if ~exist("./Vid","dir"); mkdir("./Vid_temp"); end
    v = VideoWriter("./Vid_temp/" + file_name);
    v.FrameRate = kwargs.fps;
    open(v)

    % Write first frame
    cdata = getframe(fig);
    writeVideo(v,cdata)
end

%% Animation loop
try
    for frame = 2:num_frames       
        % Update XData and YData from line objects
        for j = 1:length(new_waves)
            wave = new_waves{j};
            L = lineArray(j);
            L.YData = wave.u(:,frame); % change YData each loop
            if size(wave.x,2) ~= 1
                L.XData = wave.x(:,frame); % change XData as well (if variable)
            end
        end
        
        % Update y limits
        if yLimChangeIndices(frame)
            M = maxValues(:,frame);
            ylim(ax,M)
        end

        % Update colorbar
        current_t = current_t_fun(kwargs.times(frame));
        cm(1:current_t,:) = 0; % black
        colormap(ax2,cm)
        if num_ticks ~= length(cb.Ticks)
            num_ticks = length(cb.Ticks);
            tickLabels = string(round(linspace(tmin,tmax,num_ticks),1));
            cb.TickLabels = tickLabels;
        end
    
        drawnow
        if kwargs.save % add frame to vid
            cdata = getframe(fig);
            writeVideo(v,cdata)
        end
    
    end
catch % catch error when closing figure
end

if ishandle(fig)
    close(fig)
end
if kwargs.save
    close(v) % close if initiated
end

end