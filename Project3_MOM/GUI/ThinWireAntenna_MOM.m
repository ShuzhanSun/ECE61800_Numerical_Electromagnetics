function varargout = ThinWireAntenna_MOM(varargin)
% THINWIREANTENNA_MOM MATLAB code for ThinWireAntenna_MOM.fig
%      THINWIREANTENNA_MOM, by itself, creates a new THINWIREANTENNA_MOM or raises the existing
%      singleton*.
%
%      H = THINWIREANTENNA_MOM returns the handle to a new THINWIREANTENNA_MOM or the handle to
%      the existing singleton*.
%
%      THINWIREANTENNA_MOM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in THINWIREANTENNA_MOM.M with the given input arguments.
%
%      THINWIREANTENNA_MOM('Property','Value',...) creates a new THINWIREANTENNA_MOM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ThinWireAntenna_MOM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ThinWireAntenna_MOM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ThinWireAntenna_MOM

% Last Modified by GUIDE v2.5 24-Apr-2018 00:04:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ThinWireAntenna_MOM_OpeningFcn, ...
                   'gui_OutputFcn',  @ThinWireAntenna_MOM_OutputFcn, ...
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


% --- Executes just before ThinWireAntenna_MOM is made visible.
function ThinWireAntenna_MOM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ThinWireAntenna_MOM (see VARARGIN)

% Choose default command line output for ThinWireAntenna_MOM
handles.output = hObject;

axes(handles.axes_DemoFig);
imshow('DemoFig.PNG');

pushbutton_RefreshG_Callback(hObject, eventdata, handles);
pushbutton_RefreshI_Callback(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ThinWireAntenna_MOM wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ThinWireAntenna_MOM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_2a_Callback(hObject, eventdata, handles)
% hObject    handle to edit_2a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_2a as text
%        str2double(get(hObject,'String')) returns contents of edit_2a as a double
handles=guidata(hObject);
contents = (get(handles.popup_LBy2a,'String'));
LBy2a = str2double(contents{get(handles.popup_LBy2a,'Value')});
L = LBy2a*str2double(get(handles.edit_2a,'String'));
if L < 0
    msg = 'Error. 2a should be positve!';
    errordlg(msg);
    error(msg);
end
set(handles.edit_L,'String',num2str(L));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_2a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_2a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_L_Callback(hObject, eventdata, handles)
% hObject    handle to edit_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_L as text
%        str2double(get(hObject,'String')) returns contents of edit_L as a double


% --- Executes during object creation, after setting all properties.
function edit_L_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_NSeg_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NSeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NSeg as text
%        str2double(get(hObject,'String')) returns contents of edit_NSeg as a double
handles=guidata(hObject);

NSeg = str2double(get(handles.edit_NSeg,'String'));
if NSeg-floor(NSeg)~=0.0  || NSeg<=0
    msg = 'Error. Number of segments should be positve integer!';
    errordlg(msg);
    error(msg);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_NSeg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NSeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_LBy2a.
function popup_LBy2a_Callback(hObject, eventdata, handles)
% hObject    handle to popup_LBy2a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_LBy2a contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_LBy2a
handles=guidata(hObject);

contents = (get(hObject,'String'));
LBy2a = str2double(contents{get(hObject,'Value')});
L = LBy2a*str2double(get(handles.edit_2a,'String'));
set(handles.edit_L,'String',num2str(L));

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popup_LBy2a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_LBy2a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_Pattern.
function popup_Pattern_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Pattern (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_Pattern contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Pattern


% --- Executes during object creation, after setting all properties.
function popup_Pattern_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Pattern (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_maxLByLamb_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maxLByLamb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maxLByLamb as text
%        str2double(get(hObject,'String')) returns contents of edit_maxLByLamb as a double
handles=guidata(hObject);

maxLByLamb = str2double(get(handles.edit_maxLByLamb,'String'));
if maxLByLamb <= 0
    msg = 'Error. Max L/\lambda should be larger than 0.1!';
    errordlg(msg);
    error(msg);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_maxLByLamb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maxLByLamb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_RefreshG.
function pushbutton_RefreshG_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_RefreshG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);

axes(handles.axes_Admittance);
cla reset;

guidata(hObject, handles);


% --- Executes on button press in pushbutton_StopG.
function pushbutton_StopG_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_StopG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);

setappdata(handles.axes_Admittance,'StopPlot',1);

% updata handles structure
guidata(hObject,handles);


% --- Executes on button press in pushbutton_PlotG.
function pushbutton_PlotG_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_PlotG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);

Antenna(handles);
guidata(hObject,handles);

guidata(hObject, handles);




function edit_LByLamb_Callback(hObject, eventdata, handles)
% hObject    handle to edit_LByLamb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_LByLamb as text
%        str2double(get(hObject,'String')) returns contents of edit_LByLamb as a double
handles=guidata(hObject);

LByLamb = str2double(get(handles.edit_LByLamb,'String'));
if LByLamb <= 0
    msg = 'Error. L/\lambda should be positve!';
    errordlg(msg);
    error(msg);
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_LByLamb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_LByLamb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_RefreshI.
function pushbutton_RefreshI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_RefreshI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);

axes(handles.axes_Current);
cla reset;

guidata(hObject, handles);


% --- Executes on button press in pushbutton__plusQuaterLamb.
function pushbutton__plusQuaterLamb_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton__plusQuaterLamb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);

LByLamb = str2double(get(handles.edit_LByLamb,'String'));
LByLamb = LByLamb+0.25;
set(handles.edit_LByLamb,'String',num2str(LByLamb));

guidata(hObject, handles);


% --- Executes on button press in pushbutton_PlotI.
function pushbutton_PlotI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_PlotI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);

Antenna_PlotCurrent(handles);

guidata(hObject, handles);


% --- Executes on button press in pushbutton_minusQuaterLamb.
function pushbutton_minusQuaterLamb_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_minusQuaterLamb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);

LByLamb = str2double(get(handles.edit_LByLamb,'String'));
LByLamb = LByLamb-0.25;
set(handles.edit_LByLamb,'String',num2str(LByLamb));

guidata(hObject, handles);


% --- Executes on button press in pushbutton_resetPara.
function pushbutton_resetPara_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_resetPara (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);

set(handles.edit_2a,'String','0.01');
set(handles.edit_L,'String','0.742');
set(handles.edit_NSeg,'String','49');
set(handles.popup_LBy2a,'Value',1);
set(handles.popup_Pattern,'Value',1);

guidata(hObject, handles);
