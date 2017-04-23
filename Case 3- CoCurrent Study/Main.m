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

% Last Modified by GUIDE v2.5 10-Mar-2010 22:00:49

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

if get(handles.Uniform,'value')
    load Uniform_mesh100 
else
    load Pulse10_mesh100
end
L=7;
mesh=100;
dz=L/mesh;
D1_inner=8*10^-3;%m
D2_inner=10*10^-3;
D1_outer=14*10^-3;%m
D2_outer=16*10^-3;
dt=0.3;
Pop_value=get(handles.Pop,'value');
axes(handles.axes1)
hold off;
clear axes1
TimeStep=100;
if get(handles.D2,'value')
    switch Pop_value
        case 1
            for i=1:TimeStep
                A=rand(1,3);
                plot(1:mesh+1,Temp_Time{i}.inner.T,'color',A);
                xlabel('length (m)');ylabel('temperatur of water(c)');
                legend(['Elapsed Time = ' num2str(0.3*i) ' seconds'],'location','northoutside');
                pause(dt);hold on;
                disp(i)
            end
        case 2
            for i=1:TimeStep
                A=rand(1,3);
                plot(1:mesh+1,Temp_Time{i}.inner.P,'color',A);
                xlabel('length (m)');ylabel('Pressure of water(c)');
                legend(['Elapsed Time = ' num2str(0.3*i) ' seconds'],'location','northoutside');
                pause(dt);hold on;
                disp(i)
            end
        case 3
            for i=1:TimeStep
                A=rand(1,3);
                plot(1:mesh+1,Temp_Time{i}.outer.m,'color',A);
                xlabel('length (m)');ylabel('Mass Rate(m/s)');
                legend(['Elapsed Time = ' num2str(0.3*i) ' seconds'],'location','northoutside');
                pause(dt);hold on;
                disp(i)
            end
        case 4
            for i=1:TimeStep
                A=rand(1,3);
                plot(1:mesh+1,Temp_Time{i}.outer.T,'color',A);
                xlabel('length (m)');ylabel('Temperature of Steam(c)');
                legend(['Elapsed Time = ' num2str(0.3*i) ' seconds'],'location','northoutside');
                pause(dt);hold on;
                disp(i)
            end
        case 5
            for i=1:TimeStep
                A=rand(1,3);
                plot(1:mesh+1,Temp_Time{i}.outer.P,'color',A);
                xlabel('length (m)');ylabel('Pressure of Steam(kpa)');
                legend(['Elapsed Time = ' num2str(0.3*i) ' seconds'],'location','northoutside');
                pause(dt);hold on;ylime([0.8 1]);
                disp(i)
            end
        case 6
            for i=1:TimeStep
                A=rand(1,3);
                plot(1:mesh+1,Temp_Time{i}.outer.x,'color',A);
                xlabel('length (m)');ylabel('x');
                legend(['Elapsed Time = ' num2str(0.3*i) ' seconds'],'location','northoutside');
                pause(dt);hold on;
                disp(i)
            end
            
        case 7
            for i=1:TimeStep
                A=rand(1,3);
                plot(1:mesh,1-Temp_Time{i}.outer.void,'color',A);
                xlabel('length (m)');ylabel('void');
                legend(['Elapsed Time = ' num2str(0.3*i) ' seconds'],'location','northoutside');
                pause(dt);hold on;
                disp(i)
                ylim([0 0.1])
            end
        case 8
            for i=1:TimeStep
                A=rand(1,3);
                plot(1:mesh+1,Temp_Time{i}.outer.vg,'color',A);
                xlabel('length (m)');ylabel('velocity of water(m/s)');
                legend(['Elapsed Time = ' num2str(0.3*i) ' seconds'],'location','northoutside');
                pause(dt);hold on;
                disp(i)
            end
        case 9
            for i=1:TimeStep
                A=rand(1,3);
                plot(1:mesh+1,Temp_Time{i}.outer.vl,'color',A);
                xlabel('length (m)');ylabel('velocity of condensate(m/s)');
                legend(['Elapsed Time = ' num2str(0.3*i) ' seconds'],'location','northoutside');
                pause(dt);hold on;
                disp(i)
            end
    end
else
    switch Pop_value
        case 1
            for ii=1:TimeStep
                y_mesh=200;
                x=[0:mesh-1]*dz;
                y=[1:y_mesh];
                [X,Y]=meshgrid(x,y);
                inner_fluid_mesh=1:floor(D1_inner*y_mesh/(D2_outer));
                inner_wall_mesh=ceil(D1_inner*y_mesh/(D2_outer)):floor(D2_inner*y_mesh/(D2_outer));
                outer_fluid_mesh=ceil(D2_inner*y_mesh/(D2_outer)):floor(D1_outer*y_mesh/(D2_outer));
                outer_wall_mesh=ceil(D1_outer*y_mesh/(D2_outer)):floor(D2_outer*y_mesh/(D2_outer));
                for i=inner_fluid_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.inner.T(1:end-1));
                end
                for i=inner_wall_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.inner.T(1:end-1));
                end
                for i=outer_fluid_mesh
                    Z(i,1:mesh)=Temp_Time{ii}.inner.T(1:end-1);
                end
                for i=outer_wall_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.inner.T(1:end-1));
                end
                surf(X,Y,Z);title('subcooled water temperature vs. length');
                xlim([0,(mesh-1)*dz]);ylim([ceil(D2_inner*y_mesh/(D2_outer)),floor(D1_outer*y_mesh/(D2_outer))]);
                view(0,90);xlabel('length of pipe (m)');
                colorbar 
                disp(ii)
                legend(['Elapsed Time = ' num2str(0.3*ii) ' seconds'],'location','northoutside');
                pause(dt);
            end
        case 2
            for ii=1:TimeStep
                y_mesh=200;
                x=[0:mesh-1]*dz;
                y=[1:y_mesh];
                [X,Y]=meshgrid(x,y);
                inner_fluid_mesh=1:floor(D1_inner*y_mesh/(D2_outer));
                inner_wall_mesh=ceil(D1_inner*y_mesh/(D2_outer)):floor(D2_inner*y_mesh/(D2_outer));
                outer_fluid_mesh=ceil(D2_inner*y_mesh/(D2_outer)):floor(D1_outer*y_mesh/(D2_outer));
                outer_wall_mesh=ceil(D1_outer*y_mesh/(D2_outer)):floor(D2_outer*y_mesh/(D2_outer));
                for i=inner_fluid_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.inner.P(1:end-1));
                end
                for i=inner_wall_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.inner.P(1:end-1));
                end
                for i=outer_fluid_mesh
                    Z(i,1:mesh)=Temp_Time{ii}.inner.P(1:end-1);
                end
                for i=outer_wall_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.inner.P(1:end-1));
                end
                surf(X,Y,Z);title('subcooled water Pressure vs. length');
                xlim([0,(mesh-1)*dz]);ylim([ceil(D2_inner*y_mesh/(D2_outer)),floor(D1_outer*y_mesh/(D2_outer))]);
                view(0,90);xlabel('length of pipe (m)');
                colorbar 
                disp(ii)
                legend(['Elapsed Time = ' num2str(0.3*ii) ' seconds'],'location','northoutside');
                pause(dt);
            end
        case 3
            for ii=1:TimeStep
                y_mesh=200;
                x=[0:mesh-1]*dz;
                y=[1:y_mesh];
                [X,Y]=meshgrid(x,y);
                inner_fluid_mesh=1:floor(D1_inner*y_mesh/(D2_outer));
                inner_wall_mesh=ceil(D1_inner*y_mesh/(D2_outer)):floor(D2_inner*y_mesh/(D2_outer));
                outer_fluid_mesh=ceil(D2_inner*y_mesh/(D2_outer)):floor(D1_outer*y_mesh/(D2_outer));
                outer_wall_mesh=ceil(D1_outer*y_mesh/(D2_outer)):floor(D2_outer*y_mesh/(D2_outer));
                for i=inner_fluid_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.outer.m(1:end-1));
                end
                for i=inner_wall_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.outer.m(1:end-1));
                end
                for i=outer_fluid_mesh
                    Z(i,1:mesh)=Temp_Time{ii}.outer.m(1:end-1);
                end
                for i=outer_wall_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.outer.m(1:end-1));
                end
                surf(X,Y,Z);title('Mass flow rate of steam vs. length');
                xlim([0,(mesh-1)*dz]);ylim([ceil(D2_inner*y_mesh/(D2_outer)),floor(D1_outer*y_mesh/(D2_outer))]);
                view(0,90);xlabel('length of pipe (m)');
                colorbar 
                disp(ii)
                legend(['Elapsed Time = ' num2str(0.3*ii) ' seconds'],'location','northoutside');
                pause(dt);
            end
        case 4
            for ii=1:TimeStep
                y_mesh=200;
                x=[0:mesh-1]*dz;
                y=[1:y_mesh];
                [X,Y]=meshgrid(x,y);
                inner_fluid_mesh=1:floor(D1_inner*y_mesh/(D2_outer));
                inner_wall_mesh=ceil(D1_inner*y_mesh/(D2_outer)):floor(D2_inner*y_mesh/(D2_outer));
                outer_fluid_mesh=ceil(D2_inner*y_mesh/(D2_outer)):floor(D1_outer*y_mesh/(D2_outer));
                outer_wall_mesh=ceil(D1_outer*y_mesh/(D2_outer)):floor(D2_outer*y_mesh/(D2_outer));
                for i=inner_fluid_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.outer.T(1:end-1));
                end
                for i=inner_wall_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.outer.T(1:end-1));
                end
                for i=outer_fluid_mesh
                    Z(i,1:mesh)=Temp_Time{ii}.outer.T(1:end-1);
                end
                for i=outer_wall_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.outer.T(1:end-1));
                end
                surf(X,Y,Z);title('saturated steam temperature vs. length');
                xlim([0,(mesh-1)*dz]);ylim([ceil(D2_inner*y_mesh/(D2_outer)),floor(D1_outer*y_mesh/(D2_outer))]);
                view(0,90);xlabel('length of pipe (m)');
                colorbar 
                disp(ii)
                legend(['Elapsed Time = ' num2str(0.3*ii) ' seconds'],'location','northoutside');
                pause(dt);
            end
        case 5
            for ii=1:TimeStep
                y_mesh=200;
                x=[0:mesh-1]*dz;
                y=[1:y_mesh];
                [X,Y]=meshgrid(x,y);
                inner_fluid_mesh=1:floor(D1_inner*y_mesh/(D2_outer));
                inner_wall_mesh=ceil(D1_inner*y_mesh/(D2_outer)):floor(D2_inner*y_mesh/(D2_outer));
                outer_fluid_mesh=ceil(D2_inner*y_mesh/(D2_outer)):floor(D1_outer*y_mesh/(D2_outer));
                outer_wall_mesh=ceil(D1_outer*y_mesh/(D2_outer)):floor(D2_outer*y_mesh/(D2_outer));
                for i=inner_fluid_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.outer.P(1:end-1));
                end
                for i=inner_wall_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.outer.P(1:end-1));
                end
                for i=outer_fluid_mesh
                    Z(i,1:mesh)=Temp_Time{ii}.outer.P(1:end-1);
                end
                for i=outer_wall_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.outer.P(1:end-1));
                end
                surf(X,Y,Z);title('saturated steam pressure vs. length');
                xlim([0,(mesh-1)*dz]);ylim([ceil(D2_inner*y_mesh/(D2_outer)),floor(D1_outer*y_mesh/(D2_outer))]);
                view(0,90);xlabel('length of pipe (m)');
                colorbar 
                disp(ii)
                legend(['Elapsed Time = ' num2str(0.3*ii) ' seconds'],'location','northoutside');
                pause(dt);
            end
        case 6
            for ii=1:TimeStep
                y_mesh=200;
                x=[0:mesh-1]*dz;
                y=[1:y_mesh];
                [X,Y]=meshgrid(x,y);
                inner_fluid_mesh=1:floor(D1_inner*y_mesh/(D2_outer));
                inner_wall_mesh=ceil(D1_inner*y_mesh/(D2_outer)):floor(D2_inner*y_mesh/(D2_outer));
                outer_fluid_mesh=ceil(D2_inner*y_mesh/(D2_outer)):floor(D1_outer*y_mesh/(D2_outer));
                outer_wall_mesh=ceil(D1_outer*y_mesh/(D2_outer)):floor(D2_outer*y_mesh/(D2_outer));
                for i=inner_fluid_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.outer.x(1:end-1));
                end
                for i=inner_wall_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.outer.x(1:end-1));
                end
                for i=outer_fluid_mesh
                    Z(i,1:mesh)=Temp_Time{ii}.outer.x(1:end-1);
                end
                for i=outer_wall_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.outer.x(1:end-1));
                end
                surf(X,Y,Z);title('Quality vs. length');
                xlim([0,(mesh-1)*dz]);ylim([ceil(D2_inner*y_mesh/(D2_outer)),floor(D1_outer*y_mesh/(D2_outer))]);
                view(0,90);xlabel('length of pipe (m)');
                colorbar 
                disp(ii)
                legend(['Elapsed Time = ' num2str(0.3*ii) ' seconds'],'location','northoutside');
                pause(dt);
            end
        case 7
            for ii=1:TimeStep
                y_mesh=200;
                x=[0:mesh-1]*dz;
                y=[1:y_mesh];
                [X,Y]=meshgrid(x,y);
                inner_fluid_mesh=1:floor(D1_inner*y_mesh/(D2_outer));
                inner_wall_mesh=ceil(D1_inner*y_mesh/(D2_outer)):floor(D2_inner*y_mesh/(D2_outer));
                outer_fluid_mesh=ceil(D2_inner*y_mesh/(D2_outer)):floor(D1_outer*y_mesh/(D2_outer));
                outer_wall_mesh=ceil(D1_outer*y_mesh/(D2_outer)):floor(D2_outer*y_mesh/(D2_outer));
                for i=inner_fluid_mesh
                    Z(i,1:mesh)=1-mean(Temp_Time{2}.outer.void(1:end-1));
                end
                for i=inner_wall_mesh
                    Z(i,1:mesh)=1-mean(Temp_Time{2}.outer.void(1:end-1));
                end
                for i=outer_fluid_mesh
                    Z(i,1:mesh-1)=1-Temp_Time{ii}.outer.void(1:end-1);
                end
                for i=outer_wall_mesh
                    Z(i,1:mesh)=1-mean(Temp_Time{2}.outer.void(1:end-1));
                end
                surf(X,Y,Z);title('Liquid Film vs. length');
                xlim([0,(mesh-1)*dz]);ylim([ceil(D2_inner*y_mesh/(D2_outer)),floor(D1_outer*y_mesh/(D2_outer))]);
                view(0,90);xlabel('length of pipe (m)');ylabel('temperature(c)')
                colorbar 
                disp(ii)
                legend(['Elapsed Time = ' num2str(0.3*ii) ' seconds'],'location','northoutside');
                pause(dt);
            end
        case 8
            for ii=1:TimeStep
                y_mesh=200;
                x=[0:mesh-1]*dz;
                y=[1:y_mesh];
                [X,Y]=meshgrid(x,y);
                inner_fluid_mesh=1:floor(D1_inner*y_mesh/(D2_outer));
                inner_wall_mesh=ceil(D1_inner*y_mesh/(D2_outer)):floor(D2_inner*y_mesh/(D2_outer));
                outer_fluid_mesh=ceil(D2_inner*y_mesh/(D2_outer)):floor(D1_outer*y_mesh/(D2_outer));
                outer_wall_mesh=ceil(D1_outer*y_mesh/(D2_outer)):floor(D2_outer*y_mesh/(D2_outer));
                for i=inner_fluid_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.outer.vg(1:end-1));
                end
                for i=inner_wall_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.outer.vg(1:end-1));
                end
                for i=outer_fluid_mesh
                    Z(i,1:mesh)=Temp_Time{ii}.outer.vg(1:end-1);
                end
                for i=outer_wall_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.outer.vg(1:end-1));
                end
                surf(X,Y,Z);title('Steam velocity (m/s) vs. length');
                xlim([0,(mesh-1)*dz]);ylim([ceil(D2_inner*y_mesh/(D2_outer)),floor(D1_outer*y_mesh/(D2_outer))]);
                view(0,90);xlabel('length of pipe (m)');ylabel('temperature(c)')
                colorbar 
                disp(ii)
                legend(['Elapsed Time = ' num2str(0.3*ii) ' seconds'],'location','northoutside');
                pause(dt);
            end
        case 9
            for ii=1:TimeStep
                y_mesh=200;
                x=[0:mesh-1]*dz;
                y=[1:y_mesh];
                [X,Y]=meshgrid(x,y);
                inner_fluid_mesh=1:floor(D1_inner*y_mesh/(D2_outer));
                inner_wall_mesh=ceil(D1_inner*y_mesh/(D2_outer)):floor(D2_inner*y_mesh/(D2_outer));
                outer_fluid_mesh=ceil(D2_inner*y_mesh/(D2_outer)):floor(D1_outer*y_mesh/(D2_outer));
                outer_wall_mesh=ceil(D1_outer*y_mesh/(D2_outer)):floor(D2_outer*y_mesh/(D2_outer));
                for i=inner_fluid_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.outer.vl(1:end-1));
                end
                for i=inner_wall_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.outer.vl(1:end-1));
                end
                for i=outer_fluid_mesh
                    Z(i,1:mesh)=Temp_Time{ii}.outer.vl(1:end-1);
                end
                for i=outer_wall_mesh
                    Z(i,1:mesh)=mean(Temp_Time{2}.outer.vl(1:end-1));
                end
                surf(X,Y,Z);title('condensed water velocity (m/s) vs. length');
                xlim([0,(mesh-1)*dz]);ylim([ceil(D2_inner*y_mesh/(D2_outer)),floor(D1_outer*y_mesh/(D2_outer))]);
                view(0,90);xlabel('length of pipe (m)');ylabel('temperature(c)')
                colorbar 
                disp(ii)
                legend(['Elapsed Time = ' num2str(0.3*ii) ' seconds'],'location','northoutside');
                pause(dt);
            end
    end
end
% --- Executes on button press in Uniform.
function Uniform_Callback(hObject, eventdata, handles)
% hObject    handle to Uniform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Uniform,'value',1)
set(handles.Pulse,'value',0)
% Hint: get(hObject,'Value') returns toggle state of Uniform


% --- Executes on button press in Pulse.
function Pulse_Callback(hObject, eventdata, handles)
% hObject    handle to Pulse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Uniform,'value',0)
set(handles.Pulse,'value',1)
% Hint: get(hObject,'Value') returns toggle state of Pulse


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
