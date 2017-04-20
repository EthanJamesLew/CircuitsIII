function varargout = noise_reduction(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @noise_reduction_OpeningFcn, ...
                   'gui_OutputFcn',  @noise_reduction_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end


function noise_reduction_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

syms s t Vcc;
vcc_lh_equ = 'Vcc/1031.25+(Vcc-(4.8529/s))/(1/(s*10*10^-9))+(Vcc-(5/s+1.471*10^(-3)*.01))/(100+0.01*s)=0';
vcc_hl_equ = 'Vcc/3300+(Vcc-(4.558/s))/(1/(s*10*10^(-9)))+(Vcc-(5/s+4.420*10^(-3)*.01))/(100+0.01*s)=0';

vcc_lh_equ_c = 'Vcc/1031.25+(Vcc-(4.7688/s))/(1/(s*10*10^-9))+(Vcc-(5/s+1.471*10^(-3)*.005))/(50+0.005*s)=0';
vcc_hl_equ_c = 'Vcc/3300.00+(Vcc-(4.92/s))/(1/(s*10*10^(-9)))+(Vcc-(5/s+4.624*10^(-3)*.005))/(50+0.005*s)=0';

Vcc_lh_s(s) = solve(vcc_lh_equ, 'Vcc');
Vcc_hl_s(s) = solve(vcc_hl_equ, 'Vcc');

Vcc_lh_s_c(s) = solve(vcc_lh_equ_c, 'Vcc');
Vcc_hl_s_c(s) = solve(vcc_hl_equ_c, 'Vcc');

Vcc_lh_t(t)=ilaplace(Vcc_lh_s);
Vcc_hl_t(t)=ilaplace(Vcc_hl_s);

Vcc_lh_t_c(t)=ilaplace(Vcc_lh_s_c);
Vcc_hl_t_c(t)=ilaplace(Vcc_hl_s_c);

handles.vcc_lh_t=simplify(Vcc_lh_t);
handles.vcc_hl_t=simplify(Vcc_hl_t);

handles.vcc_lh_c=simplify(Vcc_lh_t_c);
handles.vcc_hl_c=simplify(Vcc_hl_t_c);

%digits(4);
%fprintf('Low to High Solution: ');pretty(vpa(Vcc_lh_t_c));
%fprintf('High to Low Solution: ');pretty(vpa(Vcc_hl_t_c));
%digits(32);

handles.tmax = 1/4000.;
handles.x = 0:handles.tmax/1000:handles.tmax;

Array=csvread('LtSpice_data.txt');
col1 = Array(:, 1);
col2 = Array(:, 2);

f_osc_corr = 'osc_correction.csv';
f_osc_no_corr = 'osc_no_correction.csv';
raw_corr = csvread(f_osc_corr);
corr_t = raw_corr(:,1)-(.4*10^-5);
corr_v = raw_corr(:,2);

corr_t_1 = corr_t(1:2500/2);
corr_t_2 = corr_t(2500/2+1: 2500);

corr_v_1 = corr_v(1:2500/2);
corr_v_2 = corr_v(2500/2+1 : 2500);

raw_corr = csvread(f_osc_no_corr);
corr_t = raw_corr(:,1)-(.4*10^-5);
corr_v = raw_corr(:,2);
corr_e_1 = corr_t(1:2500/2);
corr_e_2 = corr_t(2500/2+1: 2500);

corr_ev_1 = corr_v(1:2500/2);
corr_ev_2 = corr_v(2500/2+1 : 2500);

corr_v_interp1 = interp1(transpose(corr_t_2),transpose(corr_v_2), handles.x);
corr_v_interp2 = interp1(transpose(corr_t_1-corr_t_1(13)),transpose(corr_v_1), handles.x);

hold(handles.axes_lh,'on');
hold(handles.axes_hl,'on');
title(handles.axes_hl, 'High to Low', 'Interpreter', 'latex');
title(handles.axes_lh, 'Low to High', 'Interpreter', 'latex');

xlabel(handles.axes_hl,'Time (s)', 'Interpreter', 'latex');
ylabel(handles.axes_hl,'Voltage (v)','Interpreter', 'latex');

xlabel(handles.axes_lh,'Time (s)', 'Interpreter', 'latex');
ylabel(handles.axes_lh,'Voltage (v)','Interpreter', 'latex');

xlim(handles.axes_lh, [0 handles.tmax]);
xlim(handles.axes_hl, [0 handles.tmax]);

handles.plot_vcc_lh_t  = plot(handles.axes_lh,handles.x, double(handles.vcc_lh_t(handles.x)), 'linewidth', 2);
handles.plot_vcc_hl_t  = plot(handles.axes_hl,handles.x, double(handles.vcc_hl_t(handles.x)),'linewidth', 2);

handles.plot_vcc_lh_c  = plot(handles.axes_lh,handles.x, double(handles.vcc_lh_c(handles.x)), 'linewidth', 2);
handles.plot_vcc_hl_c  = plot(handles.axes_hl,handles.x, double(handles.vcc_hl_c(handles.x)),'linewidth', 2);

handles.plot_lt_lh =plot(handles.axes_lh,col1(1:320), col2(1:320),'r','linewidth', 1);
handles.plot_lt_hl =plot(handles.axes_hl, col1(321:642)-col1(321),  col2(321:642),'r','linewidth', 1);

handles.plot_exp_hl = plot(handles.axes_lh, corr_t_2, (corr_v_2), 'linewidth', 2);
handles.plot_exp_lh = plot(handles.axes_hl, corr_t_1-corr_t_1(13), (corr_v_1), 'linewidth', 2);

handles.plot_e_hl = plot(handles.axes_lh, corr_e_2, corr_ev_2, 'linewidth', 2);
handles.plot_e_lh = plot(handles.axes_hl, corr_e_1-corr_e_1(13), (corr_ev_1), 'linewidth', 2);

handles.err_lh = corr_v_interp1 -double(handles.vcc_lh_c(handles.x)) ;
handles.err_hl =corr_v_interp2-double(handles.vcc_hl_c(handles.x)) ;
set(handles.plot_vcc_lh_t,'visible','off')
set(handles.plot_vcc_hl_t,'visible','off')
set(handles.plot_lt_lh,'visible','off')
set(handles.plot_lt_hl,'visible','off')
set(handles.plot_vcc_lh_c,'visible','off')
set(handles.plot_vcc_hl_c,'visible','off')
set(handles.plot_exp_hl,'visible','off')
set(handles.plot_exp_lh,'visible','off')
set(handles.plot_e_hl,'visible','off')
set(handles.plot_e_lh,'visible','off')

guidata(hObject, handles);
handles = guidata(hObject);
grid(handles.axes_lh,'on');
grid(handles.axes_hl,'on');

grid(handles.axes_lh,'minor');
grid(handles.axes_hl,'minor');

hold(handles.axes_lh,'off');
hold(handles.axes_hl,'off');




function varargout = noise_reduction_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;



% --- Executes on button press in checkbox_hand.
function checkbox_hand_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
hold(handles.axes_lh,'on');
hold(handles.axes_hl,'on');
check = get(handles.checkbox_hand, 'Value');
if check
    axes(handles.axes_lh)
    set(handles.plot_vcc_lh_t,'visible','on')
    axes(handles.axes_hl)
    set(handles.plot_vcc_hl_t,'visible','on')
else
    axes(handles.axes_lh)
    set(handles.plot_vcc_lh_t,'visible','off')
    axes(handles.axes_hl)
    set(handles.plot_vcc_hl_t,'visible','off')
end
hold(handles.axes_lh,'off');
hold(handles.axes_hl,'off');

% --- Executes on button press in checkbox_lt_sim.
function checkbox_lt_sim_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
hold(handles.axes_lh,'on');
hold(handles.axes_hl,'on');
check = get(handles.checkbox_lt_sim, 'Value');
if check
    axes(handles.axes_lh)
    set(handles.plot_lt_lh,'visible','on')
    axes(handles.axes_hl)
    set(handles.plot_lt_hl,'visible','on')
else
    axes(handles.axes_lh)
    set(handles.plot_lt_lh,'visible','off')
    axes(handles.axes_hl)
    set(handles.plot_lt_hl,'visible','off')
end
hold(handles.axes_lh,'off');
hold(handles.axes_hl,'off');


% --- Executes on button press in checkbox_lt_adj.
function checkbox_lt_adj_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
hold(handles.axes_lh,'on');
hold(handles.axes_hl,'on');
check = get(handles.checkbox_lt_adj, 'Value');
if check
    axes(handles.axes_lh)
    set(handles.plot_vcc_lh_c,'visible','on')
    axes(handles.axes_hl)
    set(handles.plot_vcc_hl_c,'visible','on')
else
    axes(handles.axes_lh)
    set(handles.plot_vcc_lh_c,'visible','off')
    axes(handles.axes_hl)
    set(handles.plot_vcc_hl_c,'visible','off')
end
hold(handles.axes_lh,'off');
hold(handles.axes_hl,'off');

% --- Executes on button press in checkbox_exp_c.
function checkbox_exp_c_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
hold(handles.axes_lh,'on');
hold(handles.axes_hl,'on');
check = get(handles.checkbox_exp_c, 'Value');
if check
    axes(handles.axes_lh)
    set(handles.plot_exp_lh,'visible','on')
    axes(handles.axes_hl)
    set(handles.plot_exp_hl,'visible','on')
else
    axes(handles.axes_lh)
    set(handles.plot_exp_lh,'visible','off')
    axes(handles.axes_hl)
    set(handles.plot_exp_hl,'visible','off')
end
hold(handles.axes_lh,'off');
hold(handles.axes_hl,'off');

function checkbox_exp_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
hold(handles.axes_lh,'on');
hold(handles.axes_hl,'on');
check = get(handles.checkbox_exp, 'Value');
if check
    axes(handles.axes_lh)
    set(handles.plot_e_lh,'visible','on')
    axes(handles.axes_hl)
    set(handles.plot_e_hl,'visible','on')
else
    axes(handles.axes_lh)
    set(handles.plot_e_lh,'visible','off')
    axes(handles.axes_hl)
    set(handles.plot_e_hl,'visible','off')
end
hold(handles.axes_lh,'off');
hold(handles.axes_hl,'off');


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
fig = figure('Position', [100, 100, 800, 600]);
subplot(2,1,1) 
hold on
xlim( [0 handles.tmax]);
plot( handles.x, sgolayfilt(handles.err_hl,1,31), 'linewidth', 2);
title( 'High to Low Difference', 'Interpreter', 'latex');

xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Voltage (v)','Interpreter', 'latex');
hold off

subplot(2,1,2) 
hold on
xlim( [0 handles.tmax]);
title( 'Low to High Difference', 'Interpreter', 'latex');
plot( handles.x, sgolayfilt(handles.err_lh,1,31), 'linewidth', 2);
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Voltage (v)','Interpreter', 'latex');
hold off