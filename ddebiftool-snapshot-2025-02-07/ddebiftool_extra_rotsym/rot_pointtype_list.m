function fcn=rot_pointtype_list(varargin)
%% list of possible point types and flags
%
% $Id: rot_pointtype_list.m 357 2019-06-30 23:59:31Z jansieber $
%%
defaults={'stabilityfield','l0','stability',{},'kind',''};
options=dde_set_options(defaults,varargin,'pass_on');
floquet=@(p)p.stability.mu;
ev=@(p)p.stability.(options.stabilityfield);
pstab=@(x)log(abs(x)+eps);
omega_i=@(p,i)p.nvec.omega(i);
estab=@(x)real(x);
edist=@(x,y)abs(x-y);
pdist=@dde_psol_mult_dist;
fcnarr=struct(...
    'stst',  {0, ev,      estab,edist, @(p)0},...
    'psol',  {0, floquet, pstab,pdist, @(p)[1;1]},...
    'RWFold',{1, ev,      estab,edist, @(p)[0;0]},...
    'fold',{1, ev,      estab,edist, @(p)[0;0]},...
    'hopf',  {1, ev,      estab,edist, @(p)[0;1i*get_omega(p)*[-1;1]]},...
    'POfold',{1, floquet, pstab,pdist, @(p)[1;1;1]},...
    'torus', {1, floquet, pstab,pdist, @(p)[1;1;exp(1i*p.omega*pi*[-1;1])]},...
    'PD',    {1, floquet, pstab,pdist, @(p)[1;1;-1]},...
    'hcli',  {1, @(p)0,   pstab,pdist, @(p)[1;1]},...
    'hoho',  {2, ev,      estab,edist, @(p)1i*[0;omega_i(p_tohoho(p),1)*[-1;1];omega_i(p_tohoho(p),2)*[-1;1]]},...
    'zeho',  {2, ev,      estab,edist, @(p)[0;0;1i*p.nvec.omega*[-1;1]]},...
    'genh',  {2, ev,      estab,edist, @(p)1i*p.omega*[0;-1;1]},...
    'CP',    {2, ev,      estab,edist, @(p)[0;0]},...
    'cusp',  {2, ev,      estab,edist, @(p)[0;0]},...
    'BT',    {2, ev,      estab,edist, @(p)[0;0;0]});
fields={'codim','getev','stab','dist','triv'};
for i=1:length(fields)
    fcn.(fields{i})=fcnarr(i);
end
for i=1:length(options.stability)
    for k=1:length(fields)
        fcn.(fields{k}).(options.stability{i,1})=options.stability{i,2}(k);
    end
end
if ~isempty(options.kind)
    for i=1:length(fields)
        fcn.(fields{i})=getfield(fcn.(fields{i}),options.kind); %#ok<GFLD>
    end
end
end
%% 
function oms=get_omega(p)
if isfield(p,'omega')
    oms=p.omega;
elseif isfield(p,'nvec') && isfield(p.nvec,'omega')
    oms=p.nvec.omega(1);
end
end
