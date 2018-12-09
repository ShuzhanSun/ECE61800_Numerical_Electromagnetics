function varargout = EE618Proj2_GUI(varargin)
% EE618PROJ2_GUI MATLAB code for EE618Proj2_GUI.fig
%      EE618PROJ2_GUI, by itself, creates a new EE618PROJ2_GUI or raises the existing
%      singleton*.
%
%      H = EE618PROJ2_GUI returns the handle to a new EE618PROJ2_GUI or the handle to
%      the existing singleton*.
%
%      EE618PROJ2_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EE618PROJ2_GUI.M with the given input arguments.
%
%      EE618PROJ2_GUI('Property','Value',...) creates a new EE618PROJ2_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EE618Proj2_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EE618Proj2_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EE618Proj2_GUI

% Last Modified by GUIDE v2.5 03-Apr-2018 20:07:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EE618Proj2_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @EE618Proj2_GUI_OutputFcn, ...
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

% Clear all 
% close all
clc

% --- Executes just before EE618Proj2_GUI is made visible.
function EE618Proj2_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EE618Proj2_GUI (see VARARGIN)

% Choose default command line output for EE618Proj2_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EE618Proj2_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%Initial code with static Picture%%%%%%%%%%%%%%%%%%
axes(handles.axes1);
imshow('RW_Uniform.PNG');
a = 0.06;
b = 0.03;
c = 1;
d = 25;
e = 25;
set(handles.edit1,'string',a);
set(handles.edit2,'string',b);
set(handles.edit3,'string',c);
set(handles.edit4,'string',d);
set(handles.edit5,'string',e);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Outputs from this function are returned to the command line.
function varargout = EE618Proj2_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%TE10%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dflag = 0;
iflag = 1;
cla reset;
axes(handles.axes1);
imshow('RW_Uniform.PNG');
wd = str2double(get(handles.edit1,'String'));
ht = str2double(get(handles.edit2,'String'));
Nx = str2double(get(handles.edit4,'String'));
Ny = str2double(get(handles.edit5,'String'));
epsi_r=str2double(get(handles.edit3,'String'));
mode=2;
%[DT]=MeshGen(wd,ht,Nx,Ny)
% [wavenumber_z,any_K_z]=PropConsTE20(wd,ht,Nx,Ny);
[wavenumber_z]=TE_PropConsFieldDis(wd,ht,Nx,Ny,epsi_r,mode,handles.axes2);
freq = 0.1e9:0.01e9:10e9;
any_K_z= sqrt((2*pi*freq/(1/sqrt((4*pi*1.0e-7)*(8.854e-12)*epsi_r))).^2-(((1*pi)./wd).^2)-(((0*pi)./ht).^2));
%% Plot Propagation Constant for TE10
axes(handles.axes3);
plot(freq/1e9,wavenumber_z, freq/1e9,any_K_z,':','linewidth',5);
%axis([0 7 0 50]);
hold on; 
% plot(freq/1e9,any_K_z,'--mo');
% hold on; 
set(gca,'fontsize',16);
xlabel('frequency ( GHz )','FontSize',16);
ylabel('Progagation constant ( m^{-1} )','FontSize',16);
legend('Numerical Sol', 'Analytical Sol','Location','northwest');
title('Propagation constant for TE10 mode');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%TE20%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dflag = 0;
iflag = 1;
cla reset;
wd = str2double(get(handles.edit1,'String'));
ht = str2double(get(handles.edit2,'String'));
Nx = str2double(get(handles.edit4,'String'));
Ny = str2double(get(handles.edit5,'String'));
epsi_r=str2double(get(handles.edit3,'String'));
mode=3;
%[DT]=MeshGen(wd,ht,Nx,Ny)
% [wavenumber_z,any_K_z]=PropConsTE20(wd,ht,Nx,Ny);
[wavenumber_z]=TE_PropConsFieldDis(wd,ht,Nx,Ny,epsi_r,mode,handles.axes2);
freq = 0.1e9:0.01e9:10e9;
any_K_z= sqrt((2*pi*freq/(1/sqrt((4*pi*1.0e-7)*(8.854e-12)*epsi_r))).^2-(((2*pi)./wd).^2)-(((0*pi)./ht).^2));
%% Plot Propagation Constant for TE20
axes(handles.axes3);
plot(freq/1e9,wavenumber_z, freq/1e9,any_K_z,':','linewidth',5);
hold on; 
% plot(freq/1e9,any_K_z,'--mo');
% hold on; 
set(gca,'fontsize',16);
xlabel('frequency ( GHz )','FontSize',16);
ylabel('Progagation constant ( m^{-1} )','FontSize',16);
legend('Numerical Sol', 'Analytical Sol','Location','northwest');
title('Propagation constant for TE20 mode');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%TE01%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dflag = 0;
iflag = 1;
cla reset;
wd = str2double(get(handles.edit1,'String'));
ht = str2double(get(handles.edit2,'String'));
Nx = str2double(get(handles.edit4,'String'));
Ny = str2double(get(handles.edit5,'String'));
epsi_r=str2double(get(handles.edit3,'String'));
mode=4;
%[DT]=MeshGen(wd,ht,Nx,Ny)
% [wavenumber_z,any_K_z]=PropConsTE20(wd,ht,Nx,Ny);
[wavenumber_z]=TE_PropConsFieldDis(wd,ht,Nx,Ny,epsi_r,mode,handles.axes2);
freq = 0.1e9:0.01e9:10e9;
any_K_z= sqrt((2*pi*freq/(1/sqrt((4*pi*1.0e-7)*(8.854e-12)*epsi_r))).^2-(((0*pi)./wd).^2)-(((1*pi)./ht).^2));
%% Plot Propagation Constant for TE01
axes(handles.axes3);
plot(freq/1e9,wavenumber_z, freq/1e9,any_K_z,':','linewidth',5);
hold on; 
% plot(freq/1e9,any_K_z,'--mo');
% hold on; 
set(gca,'fontsize',16);
xlabel('frequency ( GHz )','FontSize',16);
ylabel('Progagation constant ( m^{-1} )','FontSize',16);
legend('Numerical Sol', 'Analytical Sol','Location','northwest');
title('Propagation constant for TE01 mode');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%TE11%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dflag = 0;
iflag = 1;
cla reset;
wd = str2double(get(handles.edit1,'String'));
ht = str2double(get(handles.edit2,'String'));
Nx = str2double(get(handles.edit4,'String'));
Ny = str2double(get(handles.edit5,'String'));
epsi_r=str2double(get(handles.edit3,'String'));
mode=5;
%[DT]=MeshGen(wd,ht,Nx,Ny)
% [wavenumber_z,any_K_z]=PropConsTE20(wd,ht,Nx,Ny);
[wavenumber_z]=TE_PropConsFieldDis(wd,ht,Nx,Ny,epsi_r,mode,handles.axes2);
freq = 0.1e9:0.01e9:10e9;
any_K_z= sqrt((2*pi*freq/(1/sqrt((4*pi*1.0e-7)*(8.854e-12)*epsi_r))).^2-(((1*pi)./wd).^2)-(((1*pi)./ht).^2));
%% Plot Propagation Constant for TE11
axes(handles.axes3);
plot(freq/1e9,wavenumber_z, freq/1e9,any_K_z,':','linewidth',5);
hold on; 
% plot(freq/1e9,any_K_z,'--mo');
% hold on; 
set(gca,'fontsize',16);
xlabel('frequency ( GHz )','FontSize',16);
ylabel('Progagation constant ( m^{-1} )','FontSize',16);
legend('Numerical Sol', 'Analytical Sol','Location','northwest');
title('Propagation constant for TE11 mode');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%TM11%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dflag = 0;
iflag = 1;
cla reset;
wd = str2double(get(handles.edit1,'String'));
ht = str2double(get(handles.edit2,'String'));
Nx = str2double(get(handles.edit4,'String'));
Ny = str2double(get(handles.edit5,'String'));
epsi_r=str2double(get(handles.edit3,'String'));
%[DT]=MeshGen(wd,ht,Nx,Ny)
[wavenumber_z,any_K_z]=TM_PropCons(wd,ht,Nx,Ny);
[kc_TM]=TM_FieldDis(wd,ht,Nx,Ny,epsi_r,handles.axes2);
% [wavenumber_z,any_K_z]=PropConsFieldDisTM(wd,ht,Nx,Ny,epsi_r,6,handles.axes2);
freq = 0.1e9:0.01e9:10e9;
%% Plot Propagation Constant for TE21
axes(handles.axes3);
plot(freq/1e9,wavenumber_z, freq/1e9,any_K_z,':','linewidth',5);
hold on; 
% plot(freq/1e9,any_K_z,'--mo');
% hold on; 
set(gca,'fontsize',16);
xlabel('frequency ( GHz )','FontSize',16);
ylabel('Progagation constant ( m^{-1} )','FontSize',16);
legend('Numerical Sol', 'Analytical Sol','Location','northwest');
title('Propagation constant for TM11 mode');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%TE21%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dflag = 0;
iflag = 1;
cla reset;
wd = str2double(get(handles.edit1,'String'));
ht = str2double(get(handles.edit2,'String'));
Nx = str2double(get(handles.edit4,'String'));
Ny = str2double(get(handles.edit5,'String'));
epsi_r=str2double(get(handles.edit3,'String'));
mode=6;
%[DT]=MeshGen(wd,ht,Nx,Ny)
% [wavenumber_z,any_K_z]=PropConsTE20(wd,ht,Nx,Ny);
[wavenumber_z]=TE_PropConsFieldDis(wd,ht,Nx,Ny,epsi_r,mode,handles.axes2);
freq = 0.1e9:0.01e9:10e9;
any_K_z= sqrt((2*pi*freq/(1/sqrt((4*pi*1.0e-7)*(8.854e-12)*epsi_r))).^2-(((2*pi)./wd).^2)-(((1*pi)./ht).^2));
%% Plot Propagation Constant for TE21
axes(handles.axes3);
plot(freq/1e9,wavenumber_z, freq/1e9,any_K_z,':','linewidth',5);
hold on; 
% plot(freq/1e9,any_K_z,'--mo');
% hold on; 
set(gca,'fontsize',16);
xlabel('frequency ( GHz )','FontSize',16);
ylabel('Progagation constant ( m^{-1} )','FontSize',16);
legend('Numerical Sol', 'Analytical Sol','Location','northwest');
title('Propagation constant for TE21 mode');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%Mesh Generator%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dflag = 0;
iflag = 1;
axes(handles.axes2)
cla reset;
wd = str2double(get(handles.edit1,'String'));
ht = str2double(get(handles.edit2,'String'));
Nx = str2double(get(handles.edit4,'String'));
Ny = str2double(get(handles.edit5,'String'));
%[ele2node]=FEM2DMesh(a,b)
[DT]=MeshGen(wd,ht,Nx,Ny);
%[N_ele,x_node, y_node, ele2node]=FEM2DMesh(wd,ht);
%% Plot 2D patches of the computation domain
% %figure;
triplot(DT);
axis tight;
set(gca,'fontsize',16);
xlabel('Width ( m )','fontsize',16);
ylabel('Height ( m )','fontsize',16);
title('Generated Mesh');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Refresh%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(handles.axes2);
cla reset;
axes(handles.axes3);
cla reset;
a = 0.06;
b = 0.03;
c = 1;
d = 25;
e = 25;
set(handles.edit1,'string',a);
set(handles.edit2,'string',b);
set(handles.edit3,'string',c);
set(handles.edit4,'string',d);
set(handles.edit5,'string',e);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
