function varargout = worms(varargin)
% WORMS M-file for worms.fig
%      WORMS, by itself, creates a new WORMS or raises the existing
%      singleton*.
%
%      H = WORMS returns the handle to a new WORMS or the handle to
%      the existing singleton*.
%
%      WORMS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WORMS.M with the given input arguments.
%
%      WORMS('Property','Value',...) creates a new WORMS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before worms_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to worms_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help worms

% Last Modified by GUIDE v2.5 01-Feb-2010 00:20:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @worms_OpeningFcn, ...
                   'gui_OutputFcn',  @worms_OutputFcn, ...
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


% --- Executes just before worms is made visible.
function worms_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to worms (see VARARGIN)

% Choose default command line output for worms
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes worms wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = worms_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function threshSlider_Callback(hObject, eventdata, handles)
% hObject    handle to threshSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
thresh = get(hObject,'Value');
thresh = floor(thresh);
set(handles.threshValue,'String',num2str(thresh));
im = getappdata(handles.imageAxes, 'ImageData');
axes(handles.imageAxes);
imagesc(double(im > thresh));
colormap gray;


% --- Executes during object creation, after setting all properties.
function threshSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in selectButton.
function selectButton_Callback(hObject, eventdata, handles)
% hObject    handle to selectButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[Filename,Pathname,Filterindex] = ...
    uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files';...
             '*.*','All Files' },'Select file',...
             pwd');
if (Filename == 0)
    % nothing...
else
    set(handles.filename, 'String', Filename);
    im=imread([Pathname,'/',Filename]);
    setappdata(handles.imageAxes, 'ImageData', double(rgb2gray(im)));
    axes(handles.imageAxes);
    imagesc(double(rgb2gray(im)));
    colormap(gray);
end


function filename_Callback(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filename as text
%        str2double(get(hObject,'String')) returns contents of filename as a double


% --- Executes during object creation, after setting all properties.
function filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function threshValue_Callback(hObject, eventdata, handles)
% hObject    handle to threshValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshValue as text
%        str2double(get(hObject,'String')) returns contents of threshValue as a double
thresh = str2double(get(hObject,'String'));
set(handles.threshSlider,'Value',thresh);
thresh = floor(thresh);
set(handles.threshValue,'String',num2str(thresh));
im = getappdata(handles.imageAxes, 'ImageData');
axes(handles.imageAxes);
imagesc(double(im > thresh));
colormap gray;


% --- Executes during object creation, after setting all properties.
function threshValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over axes background.
function imageAxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to imageAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in resetViewButton.
function resetViewButton_Callback(hObject, eventdata, handles)
% hObject    handle to resetViewButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
im = getappdata(handles.imageAxes, 'ImageData');
axes(handles.imageAxes);
imagesc(im);
colormap gray;
