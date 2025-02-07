clear
%ddebiftool_path([pwd(),'/../../DDE_Biftool2025/']);
ddebiftool_path([getenv('HOME'),'/sourceforge-ddebiftool/releases/git-2023-06-27']);
%%
ntau=2;
parnames={'a','c_ext','delta','tau_s','tau_c'};
cind=[parnames;num2cell(1:length(parnames))];
in=struct(cind{:});
%% Define system using symbolic algebra
% define arbitrary variable names
x1=sym('x1',[1,ntau+1]); 
x2=sym('x2',[1,ntau+1]); 
x3=sym('x3',[1,ntau+1]); 
x4=sym('x4',[1,ntau+1]); 
syms(parnames{:});
par=sym(parnames);
%%
%
rhs=[delta-a*x1(2)+c_ext*(x2(3)-x1(3))+c_ext*(x3(3)-x1(3))+c_ext*(x4(3)-x1(3))-x1(1)^3;...
     delta-a*x2(2)+c_ext*(x1(3)-x2(3))+c_ext*(x3(3)-x2(3))+c_ext*(x4(3)-x2(3))-x2(1)^3;...
     delta-a*x3(2)+c_ext*(x1(3)-x3(3))+c_ext*(x2(3)-x3(3))+c_ext*(x4(3)-x3(3))-x3(1)^3;...
     delta-a*x4(2)+c_ext*(x1(3)-x4(3))+c_ext*(x2(3)-x4(3))+c_ext*(x3(3)-x4(3))-x4(1)^3];
%% Defferentiate and generate code, exporting it to sym_s4group ( creat file with the same name, and continue)
[fstr,erives]=dde_sym2funcs(rhs,[x1;x2;x3;x4],par,'filename','sym_oscillators_25',...
    'maxorder',2,'directional_derivative',true);