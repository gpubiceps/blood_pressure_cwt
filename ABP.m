function varargout = ABP(varargin)
% ABP MATLAB code for ABP.fig
%      ABP, by itself, creates a new ABP or raises the existing
%      singleton*.
%
%      H = ABP returns the handle to a new ABP or the handle to
%      the existing singleton*.
%
%      ABP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ABP.M with the given input arguments.
%
%      ABP('Property','Value',...) creates a new ABP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ABP_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ABP_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ABP

% Last Modified by GUIDE v2.5 17-Apr-2020 00:51:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ABP_OpeningFcn, ...
                   'gui_OutputFcn',  @ABP_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before ABP is made visible.
function ABP_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ABP (see VARARGIN)

% Choose default command line output for ABP
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ABP wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ABP_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
global Fs pressure avg_pr_x avg_pr_y A B C D file;

checkbox = get(handles.checkbox1,'Value'); % Значение чекбокса
slider = get(handles.slider1,'Value'); % Значение ползунка
window = get(handles.edit2,'string'); % Ширина окна
window = str2double(window);
len_signal = length(pressure);

left_bound = round((len_signal-window*Fs)*slider)+1;
right_bound = left_bound+window*Fs - 1;

T = 1/Fs;
tmax = right_bound/Fs;
t = left_bound/Fs:T:tmax; % Переменная времени

axes(handles.axes1);
cla('reset');
plot(t, pressure(left_bound:right_bound), 'b');
xlabel('t, c');
ylabel('Давление, мм рт. ст.');
title(file);
y_max = max(pressure(left_bound:right_bound))*1.1;
y_min = min(pressure(left_bound:right_bound))*0.8;
xlim([t(1) t(end)])
ylim([y_min y_max])
grid on
hold on

% Построение диастолы
index_start = find(A >= left_bound, 1, 'first');
index_stop = find(A <= right_bound, 1, 'last');
peaks_index = A(index_start: index_stop);
peaks_t = peaks_index./Fs;
for i=1:length(peaks_t)
    plot(peaks_t(i), pressure(peaks_index(i)), 'o', 'MarkerSize', 5, 'Color', 'k');
end
% Среднее диастолическое давление
mean_dia = mean(pressure(A(index_start: index_stop)));

% Построение систолы
index_start = find(B >= left_bound, 1, 'first');
index_stop = find(B <= right_bound, 1, 'last');
peaks_index = B(index_start: index_stop);
peaks_t = peaks_index./Fs;
for i=1:length(peaks_t)
    plot(peaks_t(i), pressure(peaks_index(i)), 'o', 'MarkerSize', 5, 'Color', 'k');
end
% Среднее систолическое давление
mean_sist = mean(pressure(B(index_start: index_stop)));

% Построение дикротического зубца
index_start = find(C >= left_bound, 1, 'first');
index_stop = find(C <= right_bound, 1, 'last');
peaks_index = C(index_start: index_stop);
peaks_t = peaks_index./Fs;
for i=1:length(peaks_t)
    plot(peaks_t(i), pressure(peaks_index(i)), 'o', 'MarkerSize', 5, 'Color', 'k');
end

% Построение дикротического пика
index_start = find(D >= left_bound, 1, 'first');
index_stop = find(D <= right_bound, 1, 'last');
peaks_index = D(index_start: index_stop);
peaks_t = peaks_index./Fs;
for i=1:length(peaks_t)
    plot(peaks_t(i), pressure(peaks_index(i)), 'o', 'MarkerSize', 5, 'Color', 'k');
end

index_start = find(avg_pr_x >= left_bound, 1, 'first');
index_stop = find(avg_pr_x <= right_bound, 1, 'last');
slice_pulse = avg_pr_x(index_start: index_stop);
slice_pulse_t = slice_pulse./Fs;
slice_pulse_ampl = avg_pr_y(index_start: index_stop);
% Среднее пульсовое давление
mean_pul = mean(slice_pulse_ampl);
if checkbox == 1
    plot(slice_pulse_t, slice_pulse_ampl, 'r--');
end

% Частота пульса
BPM = length(B)/(length(pressure)/Fs)*60;

%Заполнение таблицы
param = [{round(mean_sist)} {'115 - 140'};
    {round(mean_dia)} {'75 - 90'};
    {round(mean_pul)} {'40 - 50'};
    {round(BPM)} {'60 - 90'};];

set(handles.uitable1, 'Data', param);
set(handles.uitable1, 'FontSize', 9);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
global Fs pressure avg_pr_x avg_pr_y A B C D file;

file = uigetfile('*.mat');
pressure = load(file);
pressure = pressure.val;
Gain = 11.91; % Усиление
Baseline = -1227; % Нулевая линия
pressure = (pressure(70:350000) - Baseline)/Gain; % Нормирование по усилению и нулевой линии

Fs = get(handles.edit1,'string'); % Частота дискретизации
Fs = str2double(Fs);

slider = get(handles.slider1,'Value'); % Значение ползунка
window = get(handles.edit2,'string'); % Ширина окна
window = str2double(window);
len_signal = length(pressure);

left_bound = round((len_signal-window*Fs)*slider)+1;
right_bound = left_bound+window*Fs;

T = 1/Fs;
tmax = right_bound/Fs;
t = left_bound/Fs:T:tmax; % Переменная времени

n = 2; % Порядок фильтра
fc = 20; % Частота среза
Wn = fc*2/Fs; % Нормирование частоты среза
[b, a] = butter(n, Wn); % Коэффициенты фильтра Баттерворта
pressure = filtfilt(b, a, pressure); % Фильтрация в две стороны для устранения временной задержки

% График
axes(handles.axes1);
plot(t, pressure(left_bound:right_bound), 'b');
xlabel('t, c');
ylabel('Давление, мм рт. ст.');
title(file);
y_max = max(pressure(left_bound:right_bound))*1.1;
y_min = min(pressure(left_bound:right_bound))*0.8;
xlim([t(1) t(end)])
ylim([y_min y_max])
grid on
hold on

% Непрерывное вейвлет преобразование для выделения характерных точек
name = 'sym6'; % Название вейвлета
scale = 120; % Масштаб вейвлета, который выбирался исходя из скейлограммы
wt = cwt(pressure, scale, name); % Непрерывное вейвлет преобразование

% Поиск пересечений с нулем коэффициентов НВП
A = [];
C_bad = [];

for i=1:length(wt)-1
    if wt(i) <= 0 && wt(i+1) > 0
        A(end+1) = i;
    elseif wt(i) >= 0 && wt(i+1) < 0
        C_bad(end+1) = i;
    end
end

% Нахождение локальных максимумов сигнала
B = [];
D = [];

for i=1:length(A)
    for j=A(i):A(i)+Fs-2
        if pressure(j) < pressure(j+1) && pressure(j+1) > pressure(j+2)
            B(end+1) = j+1;
            break
        end
    end
end

for i=1:length(C_bad)
    for j=C_bad(i):C_bad(i)+Fs-2
        if length(pressure) < j+2
            break
        end
        if pressure(j) < pressure(j+1) && pressure(j+1) > pressure(j+2)
            D(end+1) = j+1;
            break
        end
    end
end

% Уточнение дикротического зубца
C = [];
for i=1:length(B)
    for j=B(i):B(i)+Fs-2
        if length(pressure) < j+2
            break
        end
        if pressure(j) > pressure(j+1) && pressure(j+1) < pressure(j+2)
            C(end+1) = j+1;
            break
        end
    end
end

% Рассчет среднего давления
avg_pr_y = [];
avg_pr_x = [];
for i=1:length(B)
    avg_pr_y(end+1) = (pressure(B(i)) - pressure(A(i)))/2 + pressure(A(i));
    avg_pr_x(end+1) = B(i);
end

% Построение пульсового давления
checkbox = get(handles.checkbox1,'Value'); % Значение чекбокса


index_start = find(avg_pr_x >= left_bound, 1, 'first');
index_stop = find(avg_pr_x <= right_bound, 1, 'last');
slice_pulse = avg_pr_x(index_start: index_stop);
slice_pulse_t = slice_pulse./Fs;
slice_pulse_ampl = avg_pr_y(index_start: index_stop);
if checkbox == 1
    plot(slice_pulse_t, slice_pulse_ampl, 'r--');
end
% Среднее пульсовое давление
mean_pul = mean(slice_pulse_ampl);

% Построение диастолы
index_start = find(A >= left_bound, 1, 'first');
index_stop = find(A <= right_bound, 1, 'last');
peaks_index = A(index_start: index_stop);
peaks_t = peaks_index./Fs;
for i=1:length(peaks_t)
    plot(peaks_t(i), pressure(peaks_index(i)), 'o', 'MarkerSize', 5, 'Color', 'k');
end
% Среднее диастолическое давление
mean_dia = mean(pressure(A(index_start: index_stop)));

% Построение систолы
index_start = find(B >= left_bound, 1, 'first');
index_stop = find(B <= right_bound, 1, 'last');
peaks_index = B(index_start: index_stop);
peaks_t = peaks_index./Fs;
for i=1:length(peaks_t)
    plot(peaks_t(i), pressure(peaks_index(i)), 'o', 'MarkerSize', 5, 'Color', 'k');
end
% Среднее систолическое давление
mean_sist = mean(pressure(B(index_start: index_stop)));

% Построение дикротического зубца
index_start = find(C >= left_bound, 1, 'first');
index_stop = find(C <= right_bound, 1, 'last');
peaks_index = C(index_start: index_stop);
peaks_t = peaks_index./Fs;
for i=1:length(peaks_t)
    plot(peaks_t(i), pressure(peaks_index(i)), 'o', 'MarkerSize', 5, 'Color', 'k');
end

% Построение дикротического пика
index_start = find(D >= left_bound, 1, 'first');
index_stop = find(D <= right_bound, 1, 'last');
peaks_index = D(index_start: index_stop);
peaks_t = peaks_index./Fs;
for i=1:length(peaks_t)
    plot(peaks_t(i), pressure(peaks_index(i)), 'o', 'MarkerSize', 5, 'Color', 'k');
end

% Частота пульса
BPM = length(B)/(length(pressure)/Fs)*60;

%Заполнение таблицы
param = [{round(mean_sist)} {'115 - 140'};
    {round(mean_dia)} {'75 - 90'};
    {round(mean_pul)} {'40 - 50'};
    {round(BPM)} {'60 - 90'};];

set(handles.uitable1, 'Data', param);
set(handles.uitable1, 'FontSize', 9);


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% Построение пульсового давления
global Fs pressure avg_pr_x avg_pr_y A B C D file;

checkbox = get(handles.checkbox1,'Value'); % Значение чекбокса
slider = get(handles.slider1,'Value'); % Значение ползунка
window = get(handles.edit2,'string'); % Ширина окна
window = str2double(window);
len_signal = length(pressure);

left_bound = round((len_signal-window*Fs)*slider)+1;
right_bound = left_bound+window*Fs;

if checkbox == 1
    index_start = find(avg_pr_x >= left_bound, 1, 'first');
    index_stop = find(avg_pr_x <= right_bound, 1, 'last');
    slice_pulse = avg_pr_x(index_start: index_stop);
    slice_pulse_t = slice_pulse./Fs;
    slice_pulse_ampl = avg_pr_y(index_start: index_stop);

    plot(slice_pulse_t, slice_pulse_ampl, 'r--');
else
    axes(handles.axes1);
    cla('reset');
    
    T = 1/Fs;
    tmax = right_bound/Fs;
    t = left_bound/Fs:T:tmax; % Переменная времени
    axes(handles.axes1);
    plot(t, pressure(left_bound:right_bound), 'b');
    xlabel('t, c');
    ylabel('Давление, мм рт. ст.');
    title(file);
    y_max = max(pressure(left_bound:right_bound))*1.1;
    y_min = min(pressure(left_bound:right_bound))*0.8;
    xlim([t(1) t(end)])
    ylim([y_min y_max])
    grid on
    hold on
    % Построение диастолы
    index_start = find(A >= left_bound, 1, 'first');
    index_stop = find(A <= right_bound, 1, 'last');
    peaks_index = A(index_start: index_stop);
    peaks_t = peaks_index./Fs;
    for i=1:length(peaks_t)
        plot(peaks_t(i), pressure(peaks_index(i)), 'o', 'MarkerSize', 5, 'Color', 'k');
    end

    % Построение диастолы
    index_start = find(B >= left_bound, 1, 'first');
    index_stop = find(B <= right_bound, 1, 'last');
    peaks_index = B(index_start: index_stop);
    peaks_t = peaks_index./Fs;
    for i=1:length(peaks_t)
        plot(peaks_t(i), pressure(peaks_index(i)), 'o', 'MarkerSize', 5, 'Color', 'k');
    end

    % Построение дикротического зубца
    index_start = find(C >= left_bound, 1, 'first');
    index_stop = find(C <= right_bound, 1, 'last');
    peaks_index = C(index_start: index_stop);
    peaks_t = peaks_index./Fs;
    for i=1:length(peaks_t)
        plot(peaks_t(i), pressure(peaks_index(i)), 'o', 'MarkerSize', 5, 'Color', 'k');
    end

    % Построение дикротического пика
    index_start = find(D >= left_bound, 1, 'first');
    index_stop = find(D <= right_bound, 1, 'last');
    peaks_index = D(index_start: index_stop);
    peaks_t = peaks_index./Fs;
    for i=1:length(peaks_t)
        plot(peaks_t(i), pressure(peaks_index(i)), 'o', 'MarkerSize', 5, 'Color', 'k');
    end
end
