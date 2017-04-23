function varargout = Main(varargin)
% MAIN M-file for Main.fig
%      MAIN, by itself, creates a new MAIN or raises the existing
%      singleton*.
%
%      H = MAIN returns the handle to a new MAIN or the handle to
%      the existing singleton*.
%
%      MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN.M with the given input arguments.
%
%      MAIN('Property','Value',...) creates a new MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Main

% Last Modified by GUIDE v2.5 11-Mar-2010 02:50:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Main_OpeningFcn, ...
                   'gui_OutputFcn',  @Main_OutputFcn, ...
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


% --- Executes just before Main is made visible.
function Main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Main (see VARARGIN)

% Choose default command line output for Main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes Main wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in Pop.
function Pop_Callback(hObject, eventdata, handles)
% hObject    handle to Pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Pop


% --- Executes during object creation, after setting all properties.
function Pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Display.
function Display_Callback(hObject, eventdata, handles)
% hObject    handle to Display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

load TS300_M20_P5
L=20;
mesh=20;
dz=L/mesh;
D1_inner=8*10^-3;%m
D2_inner=10*10^-3;
D1_outer=14*10^-3;%m
D2_outer=16*10^-3;
dt=0.3;
Pop_value=get(handles.Pop,'value');
Pop_value_2=get(handles.Pop_2,'value');
axes(handles.axes1)
clc
hold off;
axes(handles.axes2)
hold off;
TimeSteps=300;
TimeStep=300;
%Par1=zeros(TimeStep,mesh+1);
%Par2=zeros(TimeStep,mesh+1);
switch Pop_value
    case 1
          for i=1:TimeSteps
            Par1(i,:)=Temp_Time{i}.inner.T;
          end
    case 2
        for i=1:TimeSteps
            Par1(i,:)=Temp_Time{i}.inner.P;
        end
    case 3
        for i=1:TimeSteps
            Par1(i,:)=Temp_Time{i}.outer.m;
        end
    case 4 
        for i=1:TimeSteps
            Par1(i,:)=Temp_Time{i}.outer.T;
        end
    case 5
        for i=1:TimeSteps
            Par1(i,:)=Temp_Time{i}.outer.P;
        end
    case 6
        for i=1:TimeSteps
            Par1(i,:)=Temp_Time{i}.outer.x;
        end
    case 7
        for i=1:TimeSteps
            Par1(i,:)=1-Temp_Time{i}.outer.void;
        end
    case 8 
        for i=1:TimeSteps
            Par1(i,:)=Temp_Time{i}.outer.vg;
        end
    case 9 
        for i=1:TimeSteps
            Par1(i,:)=Temp_Time{i}.outer.vl;
        end
end

switch Pop_value_2
    case 1
          for i=1:TimeSteps
            Par2(i,:)=Temp_Time{i}.inner.T;
          end
    case 2
        for i=1:TimeSteps
            Par2(i,:)=Temp_Time{i}.inner.P;
        end
    case 3
        for i=1:TimeSteps
            Par2(i,:)=Temp_Time{i}.outer.m;
        end
    case 4 
        for i=1:TimeSteps
            Par2(i,:)=Temp_Time{i}.outer.T;
        end
    case 5
        for i=1:TimeSteps
            Par2(i,:)=Temp_Time{i}.outer.P;
        end
    case 6
        for i=1:TimeSteps
            Par2(i,:)=Temp_Time{i}.outer.x;
        end
    case 7
        for i=1:TimeSteps
            Par2(i,:)=1-Temp_Time{i}.outer.void;
        end
    case 8 
        for i=1:TimeSteps
            Par2(i,:)=Temp_Time{i}.outer.vg;
        end
    case 9 
        for i=1:TimeSteps
            Par2(i,:)=Temp_Time{i}.outer.vl;
        end
end


%figure('name','Water countou map')

z=zeros(1,mesh);
for ii=1:TimeStep
    axes(handles.axes1)
    if get(handles.Contour,'value')
    y_mesh=200;
    x=[0:mesh-1]*dz;
    y=[1:y_mesh];
    [X,Y]=meshgrid(x,y);
    inner_fluid_mesh=1:floor(D1_inner*y_mesh/(D2_outer));
    inner_wall_mesh=ceil(D1_inner*y_mesh/(D2_outer)):floor(D2_inner*y_mesh/(D2_outer));
    outer_fluid_mesh=ceil(D2_inner*y_mesh/(D2_outer)):floor(D1_outer*y_mesh/(D2_outer));
    outer_wall_mesh=ceil(D1_outer*y_mesh/(D2_outer)):floor(D2_outer*y_mesh/(D2_outer));
    switch Pop_value
        case 7
            %[X,Y]=meshgrid(x(1:end-1),y(1:end-1));
            for i=inner_fluid_mesh
                Z(i,1:mesh)=mean(Par1(ii,1:end-1));
            end
            for i=inner_wall_mesh
                Z(i,1:mesh)=mean(Par1(ii,1:end-1));
            end
            for i=outer_fluid_mesh
                Z(i,1:mesh)=Par1(ii,1:end-1);
            end
            for i=outer_wall_mesh
                Z(i,1:mesh)=mean(Par1(ii,1:end-1));
            end
        otherwise
                for i=inner_fluid_mesh
                    Z(i,1:mesh)=mean(Par1(ii,1:end-1));
                end
                for i=inner_wall_mesh
                    Z(i,1:mesh)=mean(Par1(ii,1:end-1));
                end
                for i=outer_fluid_mesh
                    Z(i,1:mesh)=Par1(ii,1:end-1);
                end
                for i=outer_wall_mesh
                    Z(i,1:mesh)=mean(Par1(ii,1:end-1));
                end
      end
    surf(X,Y,Z);%title('subcooled water temperature vs. length');
    xlim([0,(mesh-1)*dz]);ylim([ceil(D2_inner*y_mesh/(D2_outer)),floor(D1_outer*y_mesh/(D2_outer))]);
    view(0,90);xlabel('length of pipe (m)');ylabel('temperature(c)')
    colorbar 
    %disp(ii)
    %legend(['Elapsed Time = ' num2str(0.3*ii) ' seconds'],'location','northoutside');
    %pause(dt);
    else
    A=rand(1,3);
    switch Pop_value
        case 7
       plot(0:mesh-1,Par1(ii,:),'color',A);
        otherwise
       plot(0:mesh,Par1(ii,:),'color',A);   
    end
    xlabel('length (m)');ylabel('temperatur of water(c)');
    %legend(['Elapsed Time = ' num2str(0.3*ii) ' seconds'],'location','northoutside');
    %pause(dt);
    hold on;
    %disp(ii)
    end
    
    
    axes(handles.axes2)
    if get(handles.Contour_2,'value')
   
    y_mesh=200;
    x=[0:mesh-1]*dz;
    y=[1:y_mesh];
    [X,Y]=meshgrid(x,y);
    inner_fluid_mesh=1:floor(D1_inner*y_mesh/(D2_outer));
    inner_wall_mesh=ceil(D1_inner*y_mesh/(D2_outer)):floor(D2_inner*y_mesh/(D2_outer));
    outer_fluid_mesh=ceil(D2_inner*y_mesh/(D2_outer)):floor(D1_outer*y_mesh/(D2_outer));
    outer_wall_mesh=ceil(D1_outer*y_mesh/(D2_outer)):floor(D2_outer*y_mesh/(D2_outer));
    switch Pop_value_2
  
        case 7
                          %    [X,Y]=meshgrid(x(1:end-1),y(1:end-1));
            for i=inner_fluid_mesh
                Z(i,1:mesh)=mean(Par2(ii,1:end-1));
            end
            for i=inner_wall_mesh
                Z(i,1:mesh)=mean(Par2(ii,1:end-1));
            end
            for i=outer_fluid_mesh
                Z(i,1:mesh)=Par2(ii,1:end-1);
            end
            for i=outer_wall_mesh
                Z(i,1:mesh)=mean(Par2(ii,1:end-1));
            end
        otherwise
                for i=inner_fluid_mesh
                    Z(i,1:mesh)=mean(Par2(ii,1:end-1));
                end
                for i=inner_wall_mesh
                    Z(i,1:mesh)=mean(Par2(ii,1:end-1));
                end
                for i=outer_fluid_mesh
                    Z(i,1:mesh)=Par2(ii,1:end-1);
                end
                for i=outer_wall_mesh
                    Z(i,1:mesh)=mean(Par2(ii,1:end-1));
                end
    end
    surf(X,Y,Z);%title('subcooled water temperature vs. length');
    xlim([0,(mesh-1)*dz]);ylim([ceil(D2_inner*y_mesh/(D2_outer)),floor(D1_outer*y_mesh/(D2_outer))]);
    view(0,90);%xlabel('length of pipe (m)');ylabel('temperature(c)')
    colorbar 
    %disp(ii)
    legend(['Elapsed Time = ' num2str(0.3*ii) ' seconds'],'location','northoutside');
    pause(dt);
    else
    A=rand(1,3);
    switch Pop_value_2
        case 7
       plot(0:mesh-1,Par2(ii,:),'color',A);
        otherwise
       plot(0:mesh,Par2(ii,:),'color',A);   
    end
    %xlabel('length (m)');
    legend(['Elapsed Time = ' num2str(0.3*ii) ' seconds'],'location','northoutside');
    pause(dt);hold on;
    %disp(ii)
    end        
end







% --- Executes on button press in Uniform.

% --- Executes on button press in D2.
function D2_Callback(hObject, eventdata, handles)
% hObject    handle to D2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.D2,'value',1)
set(handles.Contour,'value',0)
% Hint: get(hObject,'Value') returns toggle state of D2


% --- Executes on button press in Contour.
function Contour_Callback(hObject, eventdata, handles)
% hObject    handle to Contour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.D2,'value',0)
set(handles.Contour,'value',1)
% Hint: get(hObject,'Value') returns toggle state of Contour


% --- Executes on selection change in Pop_2.
function Pop_2_Callback(hObject, eventdata, handles)
% hObject    handle to Pop_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Pop_2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Pop_2


% --- Executes during object creation, after setting all properties.
function Pop_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pop_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in D2_2.
function D2_2_Callback(hObject, eventdata, handles)
% hObject    handle to D2_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.D2_2,'value',1)
set(handles.Contour_2,'value',0)
% Hint: get(hObject,'Value') returns toggle state of D2_2


% --- Executes on button press in Contour_2.
function Contour_2_Callback(hObject, eventdata, handles)
% hObject    handle to Contour_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.D2_2,'value',0)
set(handles.Contour_2,'value',1)
% Hint: get(hObject,'Value') returns toggle state of Contour_2


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)