function stability=p_stabil(funcs,p,method,varargin)
%% compute stability information for point
% function stability=p_stabil(funcs,point,method)
% INPUT:
%   funcs problem functions
%	point solution point
%	method method parameters 
% OUTPUT:
%	stability stability information

% (c) DDE-BIFTOOL v. 2.00, 30/11/2001
% Update on 05/03/2007 ("flag_newhheur"  <=> (imag(method.lms_parameter_rho)~=0) )   
%
% $Id: p_stabil.m 308 2018-10-28 15:08:12Z jansieber $
%
%%
fullfuncs=funcs;
[kind,point,funcs]=dde_get_kind(p,fullfuncs,'solution_for_stability'); %#ok<ASGLU>
stability=dde_apply({'dde_',point.kind,'_eig'},funcs,point,method,varargin{:});
end