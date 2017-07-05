function varargout = gui_spm_iTT(varargin)
% GUI_SPM_iTT M-file for gui_spm_iTT.fig
%      GUI_SPM_iTT, by itself, launches the iTT Configuration GUI
%
%      TT = GUI_SPM_iTT returns the configuration object
% 
% After accepting changes, GUI writes changes into the "config.ini"
% configuration file.
%_______________________________________________________________________
% Copyright (C) 2011 Biomedizinische NMR Forschungs GmbH

% Tibor Auer
% $Id: gui_spm_iTT.m 2011-07-06 $ 

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_spm_iTT_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_spm_iTT_OutputFcn, ...
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


% --- Executes just before gui_spm_iTT is made visible.
function gui_spm_iTT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_spm_iTT (see VARARGIN)

% loading config.ini and passing configuration object to the GUI
%--------------------------------------------------------------------------
handles.iTT = IniFile(fullfile(spm('dir'),'toolbox','iTT','config.ini'));
handles.output = handles.iTT;

% setting starting values based on the current configuration
%--------------------------------------------------------------------------
set(handles.e_ut,'String',num2str(handles.iTT.thresholds.ut));
set(handles.e_lt,'String',num2str(handles.iTT.thresholds.lt));
set(handles.ch_sh,'Value',handles.iTT.config.show);
set(handles.ch_getSPM,'Value',handles.iTT.config.check_getSPM);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_spm_iTT wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_spm_iTT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btn_OK.
function btn_OK_Callback(hObject, eventdata, handles)
% hObject    handle to btn_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% writing changes into the configuration file
%--------------------------------------------------------------------------
handles.iTT.Close;

close(handles.figure1);

% --- Executes on button press in btn_Cancel.
function btn_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% exiting without finilizing changes
%--------------------------------------------------------------------------
close(handles.figure1);

function e_ut_Callback(hObject, eventdata, handles)
% hObject    handle to e_ut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.iTT.thresholds.ut = str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function e_ut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e_ut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function e_lt_Callback(hObject, eventdata, handles)
% hObject    handle to e_lt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.iTT.thresholds.lt = str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function e_lt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e_lt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ch_sh.
function ch_sh_Callback(hObject, eventdata, handles)
% hObject    handle to ch_sh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.iTT.config.show = get(hObject,'Value');
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in btn_restore.
function btn_restore_Callback(hObject, eventdata, handles)
% hObject    handle to btn_restore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if cmd_restore
    msgbox('Original spm_getSPM.m files have been restored.','Restore','help')
else
    msgbox('No restore have been found. Possibly no need?','Restore','help')
end


% --- Executes on button press in ch_getSPM.
function ch_getSPM_Callback(hObject, eventdata, handles)
% hObject    handle to ch_getSPM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.iTT.config.check_getSPM = get(hObject,'Value');
% Update handles structure
guidata(hObject, handles);
