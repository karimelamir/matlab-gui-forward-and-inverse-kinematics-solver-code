function varargout = GUI(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 01-May-2024 10:38:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
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


% --- Executes just before GUI is made visible.
end
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
end
function varargout = GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



end
function Theta_1_Callback(hObject, eventdata, handles)
% hObject    handle to Theta_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Theta_1 as text
%        str2double(get(hObject,'String')) returns contents of Theta_1 as a double


% --- Executes during object creation, after setting all properties.
end
function Theta_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Theta_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


end
function Theta_2_Callback(hObject, eventdata, handles)
% hObject    handle to Theta_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Theta_2 as text
%        str2double(get(hObject,'String')) returns contents of Theta_2 as a double

end
% --- Executes during object creation, after setting all properties.
function Theta_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Theta_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


end
function Theta_3_Callback(hObject, eventdata, handles)
% hObject    handle to Theta_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Theta_3 as text
%        str2double(get(hObject,'String')) returns contents of Theta_3 as a double


end
% --- Executes during object creation, after setting all properties.
function Theta_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Theta_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


end
function Pos_x_Callback(hObject, eventdata, handles)
% hObject    handle to Pos_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pos_X as text
%        str2double(get(hObject,'String')) returns contents of Pos_X as a double

end
% 
% -
% -- Executes during object creation, after setting all properties.
function Pos_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pos_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


end
function Pos_y_Callback(hObject, eventdata, handles)
% hObject    handle to Pos_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pos_Y as text
%        str2double(get(hObject,'String')) returns contents of Pos_Y as a double

end
% -
% -- Executes during object creation, after setting all properties.
function Pos_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pos_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


end
function Pos_z_Callback(hObject, eventdata, handles)
% hObject    handle to Pos_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pos_Z as text
%        str2double(get(hObject,'String')) returns contents of Pos_Z as a double

end
% --- Executes during object creation, after setting all properties.
function Pos_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pos_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_forward.
% --- Executes on button press in btn_forward.

end
% Define your functions and scripts

function T_i = calculate_DH_transform(a, alpha, d, theta)
    % Calculate transformation matrix for the current joint using DH parameters
    T_i = [
        cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha), a*cos(theta);
        sin(theta), cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
        0, sin(alpha), cos(alpha), d;
        0, 0, 0, 1
    ];
end

% Define your GUI callback functions


% Other functions and scripts...

% GUI startup code

function btn_forward_Callback(hObject, eventdata, handles)
% hObject    handle to btn_forward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Theta_1 = str2double(get(handles.Theta_1, 'String'))*(pi/180) ;
Theta_2 = str2double(get(handles.Theta_2, 'String'))*(pi/180) ;
Theta_3 = str2double(get(handles.Theta_3, 'String')) *(pi/180);
L_1 = 20;
L_2 = 50;
L_3 = 40;

function T_i = calculate_DH_transform(theta, d, a, alpha)
    % Calculate transformation matrix for the current joint using DH parameters
    T_i = [
        cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha), a*cos(theta);
        sin(theta), cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
        0, sin(alpha), cos(alpha), d;
        0, 0, 0, 1
    ];
end

% Calculate forward kinematics
% Define DH parameters for each joint
a1 = 0;     % Link length
alpha1 = pi/2; % Link twist
d1 = 20;     % Link offset
theta1 = Theta_1; % Joint angle (updated to use Theta_1 input)

a2 = 50;     % Link length
alpha2 = 0; % Link twist
d2 = 0;     % Link offset
theta2 = Theta_2; % Joint angle (updated to use Theta_2 input)

a3 = 40;     % Link length
alpha3 = 0; % Link twist
d3 = 0;     % Link offset
theta3 = Theta_3; % Joint angle (updated to use Theta_3 input)

% Define DH parameters for each joint (replace these values with your actual parameters)
% a, alpha, d, and theta are DH parameters
DH_parameters = [
    theta1, d1, a1, alpha1;  % Joint 1 parameters
    theta2, d2, a2, alpha2;  % Joint 2 parameters
    theta3, d3, a3, alpha3;  % Joint 3 parameters
    % Add more rows if you have more joints
]

% Initialize identity transformation matrix
T = eye(4);

% Loop through each joint to calculate transformation matrices
for i = 1:size(DH_parameters, 1)
    % Extract DH parameters for the current joint
    a = DH_parameters(i, 3);
    alpha = DH_parameters(i, 4);
    d = DH_parameters(i, 2);
    theta = DH_parameters(i, 1);
    
    % Calculate transformation matrix for the current joint
    % (You need to implement this based on DH convention)
    % Example: 
    % T_i = [cos(theta) -sin(theta)*cos(alpha) sin(theta)*sin(alpha) a*cos(theta);
    %        sin(theta) cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta);
    %        0 sin(alpha) cos(alpha) d;
    %        0 0 0 1];
    T_i = calculate_DH_transform(theta, d, a, alpha);
    
    % Multiply with the previous transformation matrix
    T = T * T_i;
end

% Output the final transformation matrix
disp(T);



% Ensure T is a valid transformation matrix
if size(T, 1) == 4 && size(T, 2) == 4
    set(handles.Pos_x, 'String', num2str(floor(T(1,4))));
    set(handles.Pos_y, 'String', num2str(floor(T(2,4))));
    set(handles.Pos_z, 'String', num2str(floor(T(3,4))));
else
    errordlg('Error in calculating forward kinematics.', 'FK Error');
end

robot = SerialLink(DH_parameters, 'name', 'MyRobot');

% Calculate the forward kinematics transformation matrix
T2 = robot.fkine([theta1, theta2, theta3]);
robot.plot([Theta_1 Theta_2 Theta_3]);
% Display the transformation matrix
disp('Transformation Matrix (peter corke method):');
disp(T2);

end




function pushbutton2_Callback(hObject, eventdata, handles)
   Pos_x = str2double(get(handles.Pos_x, 'String'));
Pos_y = str2double(get(handles.Pos_y, 'String'));
Pos_z = str2double(get(handles.Pos_z, 'String'));

L_1 = 20;
L_2 = 50;
L_3 = 40;

L(1) = Link([0 L_1 0 pi/2]);
L(2) = Link([0 0 L_2 0]);
L(3) = Link([0 0 L_3 0]);

Robot = SerialLink(L);
Robot.name = 'robotics';

T = transl(Pos_x, Pos_y, Pos_z); % Create a translation matrix

% Inverse Kinematics - adjust the mask based on the DOF of the robot
J = Robot.ikine(T, 'mask', [1 1 1 0 0 0]) * 180 / pi;
set(handles.Theta_1, 'String', num2str(floor(J(1))));
set(handles.Theta_2, 'String', num2str(floor(J(2))));
set(handles.Theta_3, 'String', num2str(floor(J(3))));

Robot.plot(J * pi / 180);


%without using petercorke toolbox to calculate inverse kinematics
    x=Pos_x;
    y=Pos_y;
    z=Pos_z;
    % Link lengths
    L1 = 20; % Length of link 1
    L2 = 50; % Length of link 2
    L3 = 40; % Length of link 3

    % Calculate theta1
    theta1 = atan2(y, x);
    
    % Calculate distance from (x,y) to the end of link 1 projection on the xy-plane
    r = sqrt(x^2 + y^2);
    k=sqrt(r^2+(z-L1)^2);

    theta3=acos((k^2 - (L2)^2 - (L3)^2)/(2*L2*L3));
    d=cos(theta3);
    theta31=atan2(sqrt(1-d^2),d);
    theta32=atan2(-sqrt(1-d^2),d);
    beta=atan2((z-L1),r);
    h1=(L3*sin(theta31));
 
    t1=(L2+L3*cos(theta31));
   
    alpha=atan2(h1,t1);

    theta21=beta-alpha;
    beta=atan2((z-L1),r);
    h2=(L3*sin(theta32));
 
    t2=(L2+L3*cos(theta32));
   
    alpha=atan2(h2,t2);

    theta22=beta-alpha;
    th1=rad2deg(theta1);
    th2=rad2deg(theta21);
    th3=rad2deg(theta31);

    th1=rad2deg(theta1);
    th22=rad2deg(theta22);
    th32=rad2deg(theta32); 

    theta = [th1, th2, th3];
    theta2= [th1,th22,th32];
disp('inverse without petercorke(first solution)');
    disp(theta);
 disp('2nd solution');
 disp(theta2);
end

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
    % Get the user-entered joint angles
    t1 = str2double(get(handles.Theta_1, 'String'));
    t2 = str2double(get(handles.Theta_2, 'String'));
    t3 = str2double(get(handles.Theta_3, 'String'));

    % Define link parameters
    n = 3;  % Number of links
    L = ones(1, n);  % Lengths of each link (you can modify these values as needed)

    % Initialize variables to hold workspace points
    X = [];
    Y = [];

    % Generate joint angles
    resolution = 60;  % Resolution of the angle sampling
    theta = linspace(-pi, pi, resolution);

    % Compute the workspace
    for t1 = theta
        for t2 = theta
            for t3 = theta
                % Initialize the transformation matrix
                T = eye(4);
                angles = [t1, t2, t3];  % Array for three joint angles

                % Compute the overall transformation matrix
                for i = 1:n
                    Ti = [cos(angles(i)), -sin(angles(i)), 0, L(i) * cos(angles(i));
                          sin(angles(i)), cos(angles(i)), 0, L(i) * sin(angles(i));
                          0, 0, 1, 0;
                          0, 0, 0, 1];
                    T = T * Ti;
                end

                % Extract end-effector position
                x = T(1,4);
                y = T(2,4);

                % Store the position
                X = [X, x];
                Y = [Y, y];
            end
        end
    end

    % Plot the workspace
    figure;
    plot(X, Y, '.');
    title(sprintf('Workspace for %d-Link Robot', n));
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    axis equal;
end