function pptool(operation, arg1, arg2),
% A Macro for efficient phase plane analysis

% Mikael Johansson 1995

global pp_figure pp_axes;
global f1_handle f2_handle res_handle hand_mat;
global eq_menu;
global c_x1 c_x2 half_width;

if nargin == 0,
	operation = 'init';
end;

if strcmp(operation, 'draw_phase_plane'),
	% Obtain function values and grid size
	hl = size(hand_mat, 1);
	set(pp_figure, 'Color', get(pp_figure, 'Color'));
	fun1 = get(f1_handle, 'String');
	fun2 = get(f2_handle, 'String');
	if (~isempty(fun1)) & (~isempty(fun2)),
	if nargin == 2, 
		res = 0;
		currpoint = get(pp_axes, 'CurrentPoint');
		x1_vec = currpoint(1,1);
		x2_vec = currpoint(1,2);
	else
		if hl,
			clear hand_mat;
			global hand_mat;
		end;
		cla;
		res  = get(res_handle, 'Value');
		ax1min = 0.9*(c_x1-half_width);
		ax1max = 0.9*(c_x1+half_width);
		ax2min = 0.9*(c_x2-half_width);
		ax2max = 0.9*(c_x2+half_width);
		x1_vec = linspace(ax1min, ax1max, 2*res+1);
		x2_vec = linspace(ax2min, ax2max, 2*res+1);
	end;
	hold on;
	for ii=1:2*res+1,
		for jj=1:2*res+1,
			hand_mat(jj,ii) = plot(NaN, NaN, ...
						'EraseMode', 'Back', ...
						'LineWidth', 1);
		end;
	end;
	t_0 = 0;
	t_f = 2; % 20 ->2
	pow = 1/3;
	tol = 1E-3;
	t = t_0;
	h_max = (t_f-t)/5;
	h_min = (t_f-t)*1e-10;
	h_initial = (t_f-t)/20;
	ax = 0.8*[0 1 0 0 1];
	ay = 0.15*[0 0 -1 1 0];
	ad = ax+sqrt(-1)*ay;			
	x_d = get(gca, 'XLim');
	scale = 0.75*(x_d(2)-x_d(1))/2;
	set(gcf, 'Pointer', 'watch');
	for ii=1:2*res+1,
		for jj=1:2*res+1,
			x_data = [x1_vec(jj)]; 
			y_data = [x2_vec(ii)];
			x_0 = [x1_vec(jj) x2_vec(ii)];
			%% Simulate ...
			t = t_0;
			h = h_initial;
			while h >= h_min,
				if t+h > t_f,
					h = t_f -t;
				end;
				x1 = x_0(1);
				x2 = x_0(2);
				x_1 = [eval(fun1) eval(fun2)];
				dx = h*x_1;
				x1 = x_0(1)+dx(1);
				x2 = x_0(2)+dx(2);
				x_2 = [eval(fun1) eval(fun2)];
				dx = h*(x_1+x_2)/4;
				x1 = x_0(1)+dx(1);
				x2 = x_0(2)+dx(2);
				x_3 = [eval(fun1) eval(fun2)];
				delta = norm(h*(x_1-2*x_3+x_2)/3, 'inf');
				tau = tol*max(norm(x_0, 'inf'), 1.0);
				ts = t;
				xs = x_0;
				if delta <= tau,
					t = t+h;
					x_0 = x_0 + h*(x_1+4*x_3+x_2)/6;
					%% Update handles
					x_data = [x_data x_0(1)];
					y_data = [y_data x_0(2)];
					set(hand_mat(jj, ii), ...
					'XData', x_data, ...
					'YData', y_data);
					drawnow;
				end;
				if delta ~= 0.0,
					h = min(h_max, 0.9*h*(tau/delta)^pow);
				end;
			end;
			index = round(0.075*length(x_data));
			if length(x_data) > 1,
				px = x_data(index:index+1);
				py = y_data(index:index+1);
				dd = (px(2)-px(1))+sqrt(-1)*(py(2)-py(1));
				dd = 0.1*scale*dd/max(abs(dd),0.001);
				aa = ad.*dd+ones(1,5)*(px(1)+sqrt(-1)*py(1));
				plot(real(aa), imag(aa), 'LineWidth', 2, ...
				'EraseMode', 'XOR');
			end;
		end;	
	end;
	set(gcf, 'Pointer', 'Arrow');
	end;
elseif strcmp(operation, 'init'),
	c_x1 = 0;
	c_x2 = 0;
	half_width = 1;
	pp_figure = figure('Name', 'Phase Plane Analysis Tool', ...
			'NumberTitle', 'Off', ...
			'BackingStore', 'Off', ...
	'WindowButtonDownFcn', 'pptool(''draw_phase_plane'', ''mouse'');');
%XX
	pp_axes  = axes('Position', [0.25 0.3 0.5 0.6], ...
			'Units', 'Normalized', ...
			'Box', 'On', ...
			'DrawMode', 'Fast', ...
			'XLim', [-1 1], ...
			'YLim', [-1 1]);
	title('Phase Plane');
	xlabel('x1');
	ylabel('x2');
	grid;

	uicontrol('Style', 'Frame', ...
		'Position', [0.025 0.025 0.95 0.2], ...
		'Units', 'Normalized');
	uicontrol('Style', 'Text', ...
		'String', 'State Equations', ...
		'Position', [0.2 0.15 0.275 0.05], ...
		'Units', 'Normalized');
	uicontrol('Style', 'Text', ...
		'Position', [0.05 0.10 0.15 0.05], ...
		'Units', 'Normalized', ...
		'String', 'f1(x1, x2) = ');
	uicontrol('Style', 'Text', ...
		'Position', [0.05 0.04 0.15 0.05], ...
		'Units', 'Normalized', ...
		'String', 'f2(x1, x2) = ');
	f1_handle = uicontrol('Style', 'Edit', ...
		'Position', [0.2 0.10 0.275 0.05], ...
		'Units', 'Normalized');
	f2_handle = uicontrol('Style', 'Edit', ...
		'Position', [0.2 0.04 0.275 0.05], ...
		'Units', 'Normalized');
	uicontrol('Style', 'Text', ...
		'String', 'Resolution', ...
		'Position', [0.5 0.10 0.15 0.05], ...
		'Units', 'Normalized');
	res_handle = uicontrol('Style', 'Popup', ...
		'String', 'Low|Medium|High', ...
		'Position', [0.65 0.10 0.15 0.05], ...
		'Units', 'Normalized');
	uicontrol('Style', 'Text', ...
		'String', 'Zoom', ...
		'Position', [0.5 0.04 0.2 0.05], ...
		'Units', 'Normalized');
	uicontrol('Style', 'Push', ...
		'String', '-', ...
		'Position', [0.65 0.04 0.05 0.05], ...
		'Units', 'Normalized', ...
		'Callback', 'pptool(''zoom_up'');');
	uicontrol('Style', 'Push', ...
		'String', '+', ...
		'Position', [0.7 0.04 0.05 0.05], ...
		'Units', 'Normalized', ...
		'Callback', 'pptool(''zoom_down'');');
	uicontrol('Style', 'Push', ...
		'String', 'Draw', ...
		'Position', [0.85 0.1 0.1 0.05], ...
		'Units', 'Normalized', ...
		'Callback', 'pptool(''draw_phase_plane'');');
	file_menu = uimenu('Label', 'File');
	uimenu(file_menu, 'Label', 'Print', 'Accelerator', 'P', ...
		'Callback', 'pptool(''print'');');
	uimenu(file_menu, 'Label', 'Close Window', ...
		'Separator', 'On', ...
		'Accelerator', 'W', ...
		'Callback', 'pptool(''close'');');

elseif strcmp(operation, 'zoom_up'),
	half_width = half_width*2;
	set(pp_axes, 'XLim', [c_x1-half_width, c_x1+half_width], ...
		'YLim', [c_x2-half_width, c_x2+half_width]);
elseif strcmp(operation, 'zoom_down'),
	half_width = half_width/2;
	set(pp_axes, 'XLim', [c_x1-half_width, c_x1+half_width], ...
		'YLim', [c_x2-half_width, c_x2+half_width]);	
elseif strcmp(operation, 'print');
	hl = size(hand_mat, 1);
	if hl,
		print;
	else
		fprintf('Figure is empty, there is nothing to print.\n');
	end;
elseif strcmp(operation, 'close'),
	close(pp_figure);
	clear pp_figure;
else
	operation
end;