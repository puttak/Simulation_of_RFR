function varargout = Simulation_of_RFR(varargin)
% Simulation_of_RFR MATLAB code for Simulation_of_RFR.fig
%      Simulation_of_RFR, by itself, creates a new Simulation_of_RFR or raises the existing
%      singleton*.
%
%      H = Simulation_of_RFR returns the handle to a new Simulation_of_RFR or the handle to
%      the existing singleton*.
%
%      Simulation_of_RFR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in Simulation_of_RFR.M with the given input arguments.
%
%      Simulation_of_RFR('Property','Value',...) creates a new Simulation_of_RFR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Simulation_of_RFR_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Simulation_of_RFR_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to about Simulation_of_RFR

% Last Modified by GUIDE v2.5 29-Dec-2018 20:32:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Simulation_of_RFR_OpeningFcn, ...
    'gui_OutputFcn',  @Simulation_of_RFR_OutputFcn, ...
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


% --- Executes just before Simulation_of_RFR is made visible.
function Simulation_of_RFR_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Simulation_of_RFR (see VARARGIN)


% Choose default command line output for Simulation_of_RFR
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Simulation_of_RFR wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Simulation_of_RFR_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listfun.
function listfun_Callback(hObject, eventdata, handles)
% hObject    handle to listfun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
getpara(handles);
% guidata(hObject,handles);

% Hints: contents = cellstr(get(hObject,'String')) returns listfun contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listfun


function handles = getpara(handles)
% Obtain the input to show the current Cin
handles.data.u = str2double(get(handles.u,'String'));
handles.data.To = str2double(get(handles.To,'String'));
handles.data.Tc = str2double(get(handles.Tc,'String'));
handles.data.t = str2double(get(handles.t,'String'));
handles.data.rev = str2double(get(handles.rev,'String'));
handles.data.A = str2double(get(handles.A,'String'));
handles.data.B = str2double(get(handles.B,'String'));
handles.data.w = str2double(get(handles.w,'String'));
handles.data.tmax=handles.data.t*handles.data.rev;
set(handles.tmax,'String',handles.data.tmax);
% Preset the thermal diffusivity of the medium chosen from list
listStrings = get(handles.listfun,'String');
typefun = listStrings{get(handles.listfun,'Value')};
if handles.data.w ==0
    typefun='Constant';
end
switch typefun
    case 'Periodic'
        set(handles.feedfun,'String',sprintf('Cin(t)=%g+%g*sin(%g*pi*t)',handles.data.A,handles.data.B,handles.data.w));
    case 'Square'
        set(handles.feedfun,'String',sprintf('Cin(t)=%g+%g*square(%g*pi*t)',...
            handles.data.A,handles.data.B,handles.data.w));
    case 'Constant'
        set(handles.w,'String','0');
        handles.data.w=0;
        set(handles.feedfun,'String',sprintf('Cin(t)=%g',handles.data.A+handles.data.B));
end
if isnan(handles.data.u)==1 || isnan(handles.data.To)==1 || ...
        isnan(handles.data.Tc)==1 || isnan(handles.data.t)==1 || ...
        isnan(handles.data.rev)==1 || isnan(handles.data.A)==1 || ...
        isnan(handles.data.B)==1 || isnan(handles.data.w)==1
    wrnHelp(handles)
    set(handles.stop,'UserData',1);
    return
elseif (handles.data.u)<0 || (handles.data.To)<523 || ...
        (handles.data.Tc)<23 || (handles.data.t)<=0 || ...
        (handles.data.rev)<=0 || (handles.data.A)<0 || ...
        (handles.data.B)<0 || handles.data.w<0 || (handles.data.A)<(handles.data.B) || ...
        uint16(handles.data.rev)~=(handles.data.rev) || uint16(handles.data.t)~=(handles.data.t)
    wrnHelp(handles)
    set(handles.stop,'UserData',1);
    return
elseif (handles.data.u)>0.5 || (handles.data.To)>923 || ...
        (handles.data.Tc)>523 || (handles.data.A)>1 || ...
        (handles.data.B)>1 || handles.data.w>10
    wrnHelp(handles)
    set(handles.stop,'UserData',1);
    return
end

function plotcin(handles)
set(handles.stop,'UserData',0);
listStrings = get(handles.listfun,'String');
typefun = listStrings{get(handles.listfun,'Value')};
if handles.data.w ==0
    typefun='Constant';
end
switch typefun
    case 'Periodic'
        Cin=@(t) handles.data.A+handles.data.B*sin(handles.data.w*pi*t);
        tit='Feed Gas Periodic Function (Every Half Reverse)';
    case 'Square'
        Cin=@(t) handles.data.A+handles.data.B*square(handles.data.w*pi*t);
        tit='Feed Gas Like Square Wave Function (Every Half Reverse)';
    case 'Constant'
        Cin=@(t) handles.data.A+handles.data.B;
        tit='Feed Gas Constant Function (Every Half Reverse)';
end
set(handles.load2,'String','Loading ...')
pause(0.01)
ti=[0:0.002:handles.data.t/2];
for i=1:length(ti)
    Cint(i)=Cin(ti(i));
end
set(handles.load2,'String','')
handles.fig.pl = plot(handles.cin,ti,Cint,'LineStyle','-','Marker','o','Color','red');
title(handles.cin,tit),
xlabel(handles.cin,'t (Time,second)'),ylabel(handles.cin,'C_i_n (Feed Concentration, mol/L)')
xlim(handles.cin,[0 handles.data.t/2])
ylim(handles.cin,[0 handles.data.A+handles.data.B+0.1])
pause(0.01)



function dif_rfr(handles)
set(handles.stop,'UserData',0);
u=handles.data.u;
To=handles.data.To;
Tc=handles.data.Tc;
t=handles.data.t/2;
rev=handles.data.rev;
A=handles.data.A;
B=handles.data.B;
w=handles.data.w;
tmax=handles.data.tmax;
listStrings = get(handles.listfun,'String');
typefun = listStrings{get(handles.listfun,'Value')};
if handles.data.w ==0
    typefun='Constant';
end
switch typefun
    case 'Periodic'
        Cin=@(t) handles.data.A+handles.data.B*sin(handles.data.w*pi*t);
    case 'Square'
        Cin=@(t) handles.data.A+handles.data.B*square(handles.data.w*pi*t);
    case 'Constant'
        Cin=@(t) handles.data.A+handles.data.B;
end
% Input parameter
eps=0.69;ros=1082*836;rog=0.452595*1122.921;lamda=5.64;rr=0.0147/2;
Uw=2;aw=2/rr;del_H=802000;D=6.52614*10^-6;
Ea=9.629*10^4;R=8.3145;k_inf=7.34*10^7;kc=0.115;av=2000;nu=1;
% Fungsi persamaan Arrhenius
g=@(x) nu*k_inf*kc*av*exp(-Ea/(R*x))/(kc*av+nu*k_inf*exp(-Ea/(R*x)));
% Konstanta beda hingga
K=(1-eps)*ros+eps*rog;
a1=lamda/K;b1=u*rog/K;c1=Uw*aw/K;d1=del_H/K;
a2=D;b2=u/eps;d2=1/eps;
e1=((u*rog)/lamda);e2=(u/(eps*D));
% Partisi beda hingga
Nz=53;L=0.26;
dz=L/(Nz-1);dt=0.002;
Nt=floor(t/dt)+1;
r=dt/dz^2;s=dt/dz;
z=[0:dz:L];tf=[0:dt:t];
% Temperatur awal gas umpan
Tin=323;
% Nilai awal
T(1:Nz,1)=To;C(1:Nz,1)=0;
times=dt;
handles.fig.Tplot = plot(handles.Tplot,z,T(:,1),'b','LineWidth',2);
title(handles.Tplot,'Dynamic of Temperature Reactor'),
xlabel(handles.Tplot,'z(Reactor Length, m)'),ylabel(handles.Tplot,'T(Temperatur, K)')
xlim(handles.Tplot,[0 L])
ylim(handles.Tplot,[0 1000])
handles.fig.Cplot = plot(handles.Cplot,z,C(:,1),'r','LineWidth',2);
title(handles.Cplot,'Dynamic of Gas Concentration')
xlabel(handles.Cplot,'z(Reactor Length, m)'),ylabel(handles.Cplot,'C(Methane Concentration, mol/L)')
xlim(handles.Cplot,[0 L])
ylim(handles.Cplot,[0 handles.data.A+handles.data.B+0.1])
set(handles.load,'String',sprintf('Time = %.2f s, Loading ...',times-dt))
pause(0.0001)
% Iterasi beda hingga
for k=1:(rev*2)
    if get(handles.stop,'UserData') == 1
        pause(0.00001)
        break
    end
    if mod(k,2)==0
        for j=1:Nz
            Tb(j,1) = T(Nz+1-j,Nt);
            Cb(j,1) = C(Nz+1-j,Nt);
        end
        T(:,1)=Tb(:,1);
        C(:,1)=Cb(:,1);
    elseif mod(k,2)==1 && k>2
        T(:,1)=Tb(:,Nt);
        C(:,1)=Cb(:,Nt);
    end
    for n=1:(Nt-1)
        T(2,n+1)=T(2,n)+a1*r*(T(3,n)-2*T(2,n)+(e1*dz*Tin+T(2,n))/(1+e1*dz))-b1*s*(T(2,n)-(e1*dz*Tin+T(2,n))/(1+e1*dz))-c1*dt*(T(2,n)-Tc)+d1*dt*g(T(2,n))*C(2,n);
        C(2,n+1)=C(2,n)+a2*r*(C(3,n)-2*C(2,n)+(e2*dz*Cin(tf(n))+C(2,n))/(1+e2*dz))-b2*s*(C(2,n)-(e2*dz*Cin(tf(n))+C(2,n))/(1+e2*dz))-d2*dt*g(T(2,n))*C(2,n);
        T(1,n+1)=(e1*dz*Tin+T(2,n+1))/(1+e1*dz);
        C(1,n+1)=(e2*dz*Cin(tf(n))+C(2,n+1))/(1+e2*dz);
        for j=3:Nz-1
            T(j,n+1)=T(j,n)+a1*r*(T(j+1,n)-2*T(j,n)+T(j-1,n))-b1*s*(T(j,n)-T(j-1,n))-c1*dt*(T(j,n)-Tc)+d1*dt*g(T(j,n))*C(j,n);
            C(j,n+1)=C(j,n)+a2*r*(C(j+1,n)-2*C(j,n)+C(j-1,n))-b2*s*(C(j,n)-C(j-1,n))-d2*dt*g(T(j,n))*C(j,n);
        end
        T(Nz,n+1)=T(Nz,n)+a1*r*(T(Nz-1,n)-2*T(Nz,n)+T(Nz-1,n))-b1*s*(T(Nz,n)-T(Nz-1,n))-c1*dt*(T(Nz,n)-Tc)+d1*dt*g(T(Nz,n))*C(Nz,n);
        C(Nz,n+1)=C(Nz,n)+a2*r*(C(Nz-1,n)-2*C(Nz,n)+C(Nz-1,n))-b2*s*(C(Nz,n)-C(Nz-1,n))-d2*dt*g(T(Nz,n))*C(Nz,n);
        if get(handles.stop,'UserData') == 1
            break
        end
    end
    if mod(k,2)==0
        for n=1:Nt
            for j=1:Nz
                Tb(j,n)=T(Nz+1-j,n);
                Cb(j,n)=C(Nz+1-j,n);
            end
            T(:,n)=Tb(:,n);
            C(:,n)=Cb(:,n);
        end
    end
    if get(handles.stop,'UserData') == 1
        break
    end
    for n=2:(Nt)
        times=times+dt;
        if get(handles.stop,'UserData') == 1
            break
        end
        if uint16(n/25)==n/25
            if max(T(:,n))>1000
                set(handles.stop,'UserData',1);
                wrnHelp(handles)
                break
            end
            set(handles.load,'String',sprintf('Time = %.2f s',times))
            handles.fig.Tplot = plot(handles.Tplot,z,T(:,n),'b','LineWidth',2);
            title(handles.Tplot,'Dynamic of Temperature Reactor'),
            xlabel(handles.Tplot,'z(Reactor Length, m)'),ylabel(handles.Tplot,'T(Temperatur, K)')
            xlim(handles.Tplot,[0 L])
            ylim(handles.Tplot,[0 1000])
            handles.fig.Cplot = plot(handles.Cplot,z,C(:,n),'r','LineWidth',2);
            title(handles.Cplot,'Dynamic of Gas Concentration')
            xlabel(handles.Cplot,'z(Reactor Length, m)'),ylabel(handles.Cplot,'C(Methane Concentration, mol/L)')
            xlim(handles.Cplot,[0 L])
            ylim(handles.Cplot,[0 handles.data.A+handles.data.B+0.1])
            pause(0.0001)
        end
    end
    set(handles.load,'String',sprintf('Time = %.2f s, Loading ...',times))
    pause(0.0001)
end
set(handles.load,'String',sprintf('Time = %.2f s',times))
pause(0.0001)

% --------------------------------------------------------------------
function about_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Create about message
dlgname = 'About Simulation_of_RFR';
txt = {'This app is reference to :';
    '';
    'Dynamics of Oxidation Reaction using Reverse Flow Reactor';
    'with Periodic Feed Gas Like Square Wave Function:';
    'A Numerical Approch';
    '';
    'By : Aang Nuryaman, Rahmad Riyanto, and Subian Saidi';
    '       Department of Mathematics';
    '       Faculty of Mathematics and Natural Science';
    '       Lampung University';
    'ICASMI 2018'
    '';
    '';
    'Created by :';
    'Rahmad Riyanto';
    ''};
helpdlg(txt,dlgname);

function wrnHelp(handles)
% Create help message
dlgname = 'Help';
txt = {'';
    '0 <= Gas Velocity <= 0.5';
    '523 <= Initial Temperature <= 923';
    '23 <= Cooled Temperature <= 523';
    'Time > 0, integer';
    'Reverse > 0, integer';
    '0 < B <= A <= 1';
    '0 <= w <= 10';
    'If temperature reactor > 1000, Reactor will blow up';
    'simulation stoped'};
warndlg(txt,dlgname);


% --- Executes on button press in about.
function help_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
wrnHelp(handles)

% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.stop,'UserData',1);
set(handles.load,'String','')
set(handles.u,'String','0.3')
set(handles.To,'String','723')
set(handles.Tc,'String','323')
set(handles.t,'String','60')
set(handles.rev,'String','3')
set(handles.tmax,'String','180')
set(handles.A,'String','0.25')
set(handles.B,'String','0.15')
set(handles.w,'String','0.2')
set(handles.feedfun,'String','Cin(t)=A+B*f(t)')
cla(handles.cin);
cla(handles.cin,'reset');
cla(handles.Tplot);
cla(handles.Tplot,'reset');
cla(handles.Cplot);
cla(handles.Cplot,'reset');



% --- Executes on button press in show.
function show_Callback(hObject, eventdata, handles)
% hObject    handle to show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.load,'String','')
set(handles.stop,'UserData',0);
cla(handles.cin);
cla(handles.cin,'reset');
cla(handles.Tplot);
cla(handles.Tplot,'reset');
cla(handles.Cplot);
cla(handles.Cplot,'reset');
handles = getpara(handles);
if get(handles.stop,'UserData') == 1
    return
end
plotcin(handles)
guidata(hObject,handles);

% --- Executes on button press in start.
function start_Callback(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% check if the Stop button is pressed: if not, proceed
set(handles.stop,'UserData',0);
cla(handles.Tplot);
cla(handles.Tplot,'reset');
cla(handles.Cplot);
cla(handles.Cplot,'reset');
cla(handles.cin);
cla(handles.cin,'reset');
set(handles.load,'String','')
handles = getpara(handles);
if get(handles.stop,'UserData') == 1
    return
end
plotcin(handles)
pause(0.0001)
dif_rfr(handles)
guidata(hObject,handles);

% --- Executes on button press in stop.
function stop_Callback(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.stop,'UserData',1);


function tmax_Callback(hObject, eventdata, handles)
% hObject    handle to tmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
getpara(handles);
% Hints: get(hObject,'String') returns contents of tmax as text
%        str2double(get(hObject,'String')) returns contents of tmax as a double



function rev_Callback(hObject, eventdata, handles)
% hObject    handle to rev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
getpara(handles);
% Hints: get(hObject,'String') returns contents of rev as text
%        str2double(get(hObject,'String')) returns contents of rev as a double



function t_Callback(hObject, eventdata, handles)
% hObject    handle to t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
getpara(handles);
% Hints: get(hObject,'String') returns contents of t as text
%        str2double(get(hObject,'String')) returns contents of t as a double



function A_Callback(hObject, eventdata, handles)
% hObject    handle to A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
getpara(handles);
% Hints: get(hObject,'String') returns contents of A as text
%        str2double(get(hObject,'String')) returns contents of A as a double



function B_Callback(hObject, eventdata, handles)
% hObject    handle to B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
getpara(handles);
% Hints: get(hObject,'String') returns contents of B as text
%        str2double(get(hObject,'String')) returns contents of B as a double



function w_Callback(hObject, eventdata, handles)
% hObject    handle to w (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
getpara(handles);
% Hints: get(hObject,'String') returns contents of w as text
%        str2double(get(hObject,'String')) returns contents of w as a double


% --- Executes on key press with focus on listfun and none of its controls.
function listfun_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to listfun (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
getpara(handles);


function Tc_Callback(hObject, eventdata, handles)
% hObject    handle to Tc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
getpara(handles);
% Hints: get(hObject,'String') returns contents of Tc as text
%        str2double(get(hObject,'String')) returns contents of Tc as a double



function To_Callback(hObject, eventdata, handles)
% hObject    handle to To (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
getpara(handles);
% Hints: get(hObject,'String') returns contents of To as text
%        str2double(get(hObject,'String')) returns contents of To as a double



function u_Callback(hObject, eventdata, handles)
% hObject    handle to u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
getpara(handles);
% Hints: get(hObject,'String') returns contents of u as text
%        str2double(get(hObject,'String')) returns contents of u as a double
