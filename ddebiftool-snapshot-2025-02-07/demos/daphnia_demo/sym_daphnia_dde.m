function varargout=sym_daphnia_dde(action,varargin)
%% Automatically generated with matlabFunction
% 
%#ok<*DEFNU,*INUSD,*INUSL>

switch action
  case 'ntau'
   varargout{1}=0;
   return
  case 'npar'
   varargout{1}=12;
   return
  case 'nf'
   varargout{1}=2;
   return
  case 'nx'
   varargout{1}=5;
   return
  case 'tp_del'
   varargout{1}=0;
   return
  case 'maxorder'
   varargout{1}=3;
   return
  case 'iscollected'
   varargout{1}=0;
   return
  case 'directional_derivative'
   varargout{1}=1;
   return
  case 'sys_tau_seq'
   varargout{1}={};
   return
  case 'sys_tau_in'
   varargout{1}=[];
   return
  case 'xpattern'
   varargout{1}=[1  1  2  1  4  1  5  1];
   return
end
ind=varargin{1};
order=varargin{2};
nout=varargin{3};
f=str2func(sprintf('sym_daphnia_dde_%s_%d_%d',action,ind,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{4:end});
end



function [out_1,out_2] = sym_daphnia_dde_rhs_1_0(in1,in2,in3,in4)
%SYM_DAPHNIA_DDE_RHS_1_0
%    [OUT_1,OUT_2] = SYM_DAPHNIA_DDE_RHS_1_0(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    07-Feb-2025 00:00:37

in_1_n1 = in1(:,1);
in_1_n2 = in1(:,2);
in_1_n4 = in1(:,4);
in_1_n5 = in1(:,5);
in_2_n1 = in2(:,1);
in_2_n2 = in2(:,2);
in_2_n4 = in2(:,4);
in_2_n5 = in2(:,5);
in_2_n6 = in2(:,6);
in_2_n7 = in2(:,7);
in_2_n9 = in2(:,9);
in_2_n10 = in2(:,10);
in_2_n11 = in2(:,11);
t2 = in_1_n1.*in_2_n5;
out_1 = -in_1_n1.*in_2_n10.*(in_1_n1./in_2_n9-1.0)-(in_1_n4.*in_2_n6.*t2)./(in_2_n7.*(t2+1.0));
if nargout > 1
    out_2 = in_1_n5-in_2_n2+in_2_n1.*exp(-in_1_n2.*in_2_n4.*in_2_n11);
end
end


function [out_1,out_2] = sym_daphnia_dde_rhs_1_1(in1,in2,in3,in4)
%SYM_DAPHNIA_DDE_RHS_1_1
%    [OUT_1,OUT_2] = SYM_DAPHNIA_DDE_RHS_1_1(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    07-Feb-2025 00:00:37

in_1_n1 = in1(:,1);
in_1_n2 = in1(:,2);
in_1_n4 = in1(:,4);
in_2_n1 = in2(:,1);
in_2_n4 = in2(:,4);
in_2_n5 = in2(:,5);
in_2_n6 = in2(:,6);
in_2_n7 = in2(:,7);
in_2_n9 = in2(:,9);
in_3_n1 = in3(:,1);
in_3_n2 = in3(:,2);
in_3_n4 = in3(:,4);
in_3_n5 = in3(:,5);
in_4_n1 = in4(:,1);
in_4_n2 = in4(:,2);
in_4_n4 = in4(:,4);
in_4_n5 = in4(:,5);
in_4_n6 = in4(:,6);
in_4_n7 = in4(:,7);
in_4_n9 = in4(:,9);
in_2_n10 = in2(:,10);
in_2_n11 = in2(:,11);
in_4_n10 = in4(:,10);
in_4_n11 = in4(:,11);
t2 = in_1_n1.*in_2_n5;
t3 = in_1_n2.*in_2_n4.*in_2_n11;
t4 = 1.0./in_2_n7;
t5 = 1.0./in_2_n9;
t6 = t2+1.0;
t7 = in_1_n1.*t5;
t8 = -t3;
t9 = exp(t8);
t10 = 1.0./t6;
t11 = t7-1.0;
out_1 = -in_1_n1.*in_2_n10.*(in_3_n1.*t5-in_4_n9.*t5.*t7)-in_3_n1.*in_2_n10.*t11-in_1_n1.*in_4_n10.*t11-in_1_n4.*in_4_n6.*t2.*t4.*t10-in_2_n6.*in_3_n4.*t2.*t4.*t10+in_1_n4.*in_2_n6.*in_4_n7.*t2.*t4.^2.*t10+in_1_n4.*in_2_n6.*t2.*t4.*t10.^2.*(in_1_n1.*in_4_n5+in_2_n5.*in_3_n1)-in_1_n1.*in_1_n4.*in_2_n6.*in_4_n5.*t4.*t10-in_1_n4.*in_2_n5.*in_2_n6.*in_3_n1.*t4.*t10;
if nargout > 1
    out_2 = in_3_n5-in_4_n2+in_4_n1.*t9-in_2_n1.*t9.*(in_1_n2.*in_4_n4.*in_2_n11+in_2_n4.*in_3_n2.*in_2_n11+in_1_n2.*in_2_n4.*in_4_n11);
end
end


function [out_1,out_2] = sym_daphnia_dde_rhs_1_2(in1,in2,in3,in4)
%SYM_DAPHNIA_DDE_RHS_1_2
%    [OUT_1,OUT_2] = SYM_DAPHNIA_DDE_RHS_1_2(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    07-Feb-2025 00:00:37

in_1_n1 = in1(:,1);
in_1_n2 = in1(:,2);
in_1_n4 = in1(:,4);
in_2_n1 = in2(:,1);
in_2_n4 = in2(:,4);
in_2_n5 = in2(:,5);
in_2_n6 = in2(:,6);
in_2_n7 = in2(:,7);
in_2_n9 = in2(:,9);
in_3_n1 = in3(:,1);
in_3_n2 = in3(:,2);
in_3_n4 = in3(:,4);
in_4_n1 = in4(:,1);
in_4_n4 = in4(:,4);
in_4_n5 = in4(:,5);
in_4_n6 = in4(:,6);
in_4_n7 = in4(:,7);
in_4_n9 = in4(:,9);
in_2_n10 = in2(:,10);
in_2_n11 = in2(:,11);
in_4_n10 = in4(:,10);
in_4_n11 = in4(:,11);
t2 = in_1_n1.*in_2_n5;
t3 = in_1_n1.*in_4_n5;
t4 = in_2_n5.*in_3_n1;
t5 = in_1_n2.*in_2_n4.*in_2_n11;
t6 = in_1_n2.*in_4_n4.*in_2_n11;
t7 = in_2_n4.*in_3_n2.*in_2_n11;
t8 = in_1_n2.*in_2_n4.*in_4_n11;
t9 = 1.0./in_2_n7;
t11 = 1.0./in_2_n9;
t10 = t9.^2;
t12 = t11.^2;
t13 = t2+1.0;
t14 = in_3_n1.*t11;
t15 = -t5;
t17 = t3+t4;
t22 = t6+t7+t8;
t16 = exp(t15);
t18 = in_1_n1.*in_4_n9.*t12;
t19 = 1.0./t13;
t20 = t19.^2;
t21 = -t18;
t23 = t14+t21;
et1 = -in_1_n1.*in_2_n10.*(in_1_n1.*in_4_n9.^2.*t11.^3.*2.0-in_3_n1.*in_4_n9.*t12.*2.0)-in_3_n1.*in_2_n10.*t23.*2.0-in_1_n1.*in_4_n10.*t23.*2.0-in_3_n1.*in_4_n10.*(in_1_n1.*t11-1.0).*2.0-in_1_n4.*in_4_n6.*t3.*t9.*t19.*2.0-in_2_n6.*in_3_n4.*t3.*t9.*t19.*2.0-in_1_n4.*in_4_n6.*t4.*t9.*t19.*2.0-in_2_n6.*in_3_n4.*t4.*t9.*t19.*2.0-in_3_n4.*in_4_n6.*t2.*t9.*t19.*2.0+in_1_n4.*in_2_n6.*t3.*t9.*t17.*t20.*2.0+in_1_n4.*in_2_n6.*t4.*t9.*t17.*t20.*2.0+in_1_n4.*in_4_n6.*t2.*t9.*t17.*t20.*2.0+in_2_n6.*in_3_n4.*t2.*t9.*t17.*t20.*2.0-in_1_n4.*in_2_n6.*in_4_n7.^2.*t2.*t9.^3.*t19.*2.0-in_1_n4.*in_2_n6.*t2.*t9.*t17.^2.*t19.^3.*2.0-in_1_n4.*in_2_n6.*in_3_n1.*in_4_n5.*t9.*t19.*2.0+in_1_n4.*in_2_n6.*in_4_n7.*t3.*t10.*t19.*2.0+in_1_n4.*in_2_n6.*in_4_n7.*t4.*t10.*t19.*2.0+in_1_n4.*in_4_n6.*in_4_n7.*t2.*t10.*t19.*2.0+in_2_n6.*in_3_n4.*in_4_n7.*t2.*t10.*t19.*2.0+in_1_n4.*in_2_n6.*in_3_n1.*in_4_n5.*t2.*t9.*t20.*2.0;
et2 = in_1_n4.*in_2_n6.*in_4_n7.*t2.*t10.*t17.*t20.*-2.0;
out_1 = et1+et2;
if nargout > 1
    out_2 = in_2_n1.*t16.*t22.^2-in_2_n1.*t16.*(in_3_n2.*in_4_n4.*in_2_n11.*2.0+in_1_n2.*in_4_n4.*in_4_n11.*2.0+in_2_n4.*in_3_n2.*in_4_n11.*2.0)-in_4_n1.*t16.*t22.*2.0;
end
end


function [out_1,out_2] = sym_daphnia_dde_rhs_1_3(in1,in2,in3,in4)
%SYM_DAPHNIA_DDE_RHS_1_3
%    [OUT_1,OUT_2] = SYM_DAPHNIA_DDE_RHS_1_3(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    07-Feb-2025 00:00:37

in_1_n1 = in1(:,1);
in_1_n2 = in1(:,2);
in_1_n4 = in1(:,4);
in_2_n1 = in2(:,1);
in_2_n4 = in2(:,4);
in_2_n5 = in2(:,5);
in_2_n6 = in2(:,6);
in_2_n7 = in2(:,7);
in_2_n9 = in2(:,9);
in_3_n1 = in3(:,1);
in_3_n2 = in3(:,2);
in_3_n4 = in3(:,4);
in_4_n1 = in4(:,1);
in_4_n4 = in4(:,4);
in_4_n5 = in4(:,5);
in_4_n6 = in4(:,6);
in_4_n7 = in4(:,7);
in_4_n9 = in4(:,9);
in_2_n10 = in2(:,10);
in_2_n11 = in2(:,11);
in_4_n10 = in4(:,10);
in_4_n11 = in4(:,11);
t2 = in_1_n1.*in_2_n5;
t3 = in_1_n1.*in_4_n5;
t4 = in_2_n5.*in_3_n1;
t5 = in_4_n7.^2;
t6 = in_4_n9.^2;
t7 = in_1_n2.*in_2_n4.*in_2_n11;
t8 = in_1_n2.*in_4_n4.*in_2_n11;
t9 = in_2_n4.*in_3_n2.*in_2_n11;
t10 = in_1_n2.*in_2_n4.*in_4_n11;
t11 = 1.0./in_2_n7;
t14 = 1.0./in_2_n9.^2;
t15 = 1.0./in_2_n9.^3;
t17 = in_3_n2.*in_4_n4.*in_2_n11.*2.0;
t18 = in_1_n2.*in_4_n4.*in_4_n11.*2.0;
t19 = in_2_n4.*in_3_n2.*in_4_n11.*2.0;
t12 = t11.^2;
t13 = t11.^3;
t16 = t2+1.0;
t20 = -t7;
t22 = t3+t4;
t26 = in_3_n1.*in_4_n9.*t14.*2.0;
t29 = in_1_n1.*t6.*t15.*2.0;
t30 = t8+t9+t10;
t31 = t17+t18+t19;
t21 = exp(t20);
t23 = 1.0./t16;
t27 = t22.^2;
t24 = t23.^2;
t25 = t23.^3;
et1 = in_1_n1.*in_2_n10.*(in_1_n1.*in_4_n9.^3.*t14.^2.*6.0-in_3_n1.*t6.*t15.*6.0)-in_3_n1.*in_4_n10.*(in_3_n1./in_2_n9-in_1_n1.*in_4_n9.*t14).*6.0+in_3_n1.*in_2_n10.*(t26-t29).*3.0+in_1_n1.*in_4_n10.*(t26-t29).*3.0-in_3_n4.*in_4_n6.*t3.*t11.*t23.*6.0-in_3_n4.*in_4_n6.*t4.*t11.*t23.*6.0-in_1_n4.*in_2_n6.*t3.*t5.*t13.*t23.*6.0-in_1_n4.*in_2_n6.*t4.*t5.*t13.*t23.*6.0-in_1_n4.*in_4_n6.*t2.*t5.*t13.*t23.*6.0-in_2_n6.*in_3_n4.*t2.*t5.*t13.*t23.*6.0-in_1_n4.*in_2_n6.*t3.*t11.*t25.*t27.*6.0-in_1_n4.*in_2_n6.*t4.*t11.*t25.*t27.*6.0+in_1_n4.*in_4_n6.*t3.*t11.*t22.*t24.*6.0+in_2_n6.*in_3_n4.*t3.*t11.*t22.*t24.*6.0+in_1_n4.*in_4_n6.*t4.*t11.*t22.*t24.*6.0+in_2_n6.*in_3_n4.*t4.*t11.*t22.*t24.*6.0-in_1_n4.*in_4_n6.*t2.*t11.*t25.*t27.*6.0-in_2_n6.*in_3_n4.*t2.*t11.*t25.*t27.*6.0+in_3_n4.*in_4_n6.*t2.*t11.*t22.*t24.*6.0+in_1_n4.*in_2_n6.*in_4_n7.^3.*t2.*t12.^2.*t23.*6.0;
et2 = in_1_n4.*in_2_n6.*t2.*t11.*t22.^3.*t24.^2.*6.0-in_1_n4.*in_3_n1.*in_4_n5.*in_4_n6.*t11.*t23.*6.0-in_2_n6.*in_3_n1.*in_3_n4.*in_4_n5.*t11.*t23.*6.0+in_1_n4.*in_4_n6.*in_4_n7.*t3.*t12.*t23.*6.0+in_2_n6.*in_3_n4.*in_4_n7.*t3.*t12.*t23.*6.0+in_1_n4.*in_4_n6.*in_4_n7.*t4.*t12.*t23.*6.0+in_2_n6.*in_3_n4.*in_4_n7.*t4.*t12.*t23.*6.0+in_3_n4.*in_4_n6.*in_4_n7.*t2.*t12.*t23.*6.0+in_1_n4.*in_2_n6.*in_3_n1.*in_4_n5.*in_4_n7.*t12.*t23.*6.0+in_1_n4.*in_2_n6.*in_3_n1.*in_4_n5.*t3.*t11.*t24.*6.0+in_1_n4.*in_2_n6.*in_3_n1.*in_4_n5.*t4.*t11.*t24.*6.0+in_1_n4.*in_2_n6.*in_3_n1.*in_4_n5.*t11.*t22.*t24.*6.0+in_1_n4.*in_3_n1.*in_4_n5.*in_4_n6.*t2.*t11.*t24.*6.0+in_2_n6.*in_3_n1.*in_3_n4.*in_4_n5.*t2.*t11.*t24.*6.0-in_1_n4.*in_2_n6.*in_4_n7.*t3.*t12.*t22.*t24.*6.0-in_1_n4.*in_2_n6.*in_4_n7.*t4.*t12.*t22.*t24.*6.0+in_1_n4.*in_2_n6.*in_4_n7.*t2.*t12.*t25.*t27.*6.0-in_1_n4.*in_4_n6.*in_4_n7.*t2.*t12.*t22.*t24.*6.0-in_2_n6.*in_3_n4.*in_4_n7.*t2.*t12.*t22.*t24.*6.0+in_1_n4.*in_2_n6.*t2.*t5.*t13.*t22.*t24.*6.0-in_1_n4.*in_2_n6.*in_3_n1.*in_4_n5.*in_4_n7.*t2.*t12.*t24.*6.0;
et3 = in_1_n4.*in_2_n6.*in_3_n1.*in_4_n5.*t2.*t11.*t22.*t25.*-1.2e+1;
out_1 = et1+et2+et3;
if nargout > 1
    out_2 = -in_2_n1.*t21.*t30.^3+in_4_n1.*t21.*t30.^2.*3.0-in_4_n1.*t21.*t31.*3.0+in_2_n1.*t21.*t30.*t31.*3.0-in_2_n1.*in_3_n2.*in_4_n4.*in_4_n11.*t21.*6.0;
end
end

