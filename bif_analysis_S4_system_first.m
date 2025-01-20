clear;
base=[pwd(),'\..\..\DDE_Biftool2025\'];
addpath([base,'ddebiftool'],...
    [base,'ddebiftool_extra_psol'],...
    [base,'ddebiftool_utilities'],...
    [base,'ddebiftool_extra_rotsym'],...
    [base,'ddebiftool_extra_nmfm'],...
    [base,'ddebiftool_extra_symbolic'],...
    [base,'ddebiftool_coco']);
%%
parnames={'a','c_ext','delta','tau_s','tau_c'};
cind=[parnames;num2cell(1:length(parnames))];
in=struct(cind{:});
par([in.a,  in.c_ext, in.delta,  in.tau_s,  in.tau_c])=...
    [1,        1,        0,        0.025,     0.02];
x0=[0;0;0;0];
parbd={'min_bound',[in.c_ext,0;in.tau_s,0; in.tau_c,0],...
    'max_bound',[in.c_ext,3;in.tau_s,3; in.tau_c,5],...
    'max_step',[in.c_ext,0.05; in.tau_s,0.01;in.tau_c,0.01; 0,0.01]};
funcs=set_symfuncs(@sym_oscillators_25,'sys_tau',@()[in.tau_s,in.tau_c]);
%% compute and find trivial equilibri in tau_s
[branch0,suc0]=SetupStst(funcs,...
    'parameter',par,'x',x0,...
    'contpar',in.tau_s,'step',0.01,parbd{:},...
    'print_residual_info',true);
%
branch0.method.point.newton_max_iterations=10;
figure(1)
clf;
[branch0,ss1,fs1,rs1]=br_contn(funcs,branch0,250);
branch0=br_rvers(branch0);
[branch0,ss01,fs01,rs01]=br_contn(funcs,branch0,200);
%% compute stability
[branch0,unst_ind0,~,~]=br_stabl(funcs,branch0,0,0);
[branch0,bif1testfuncs,indbif1,bif1types]=LocateSpecialPoints(funcs,branch0,'debug',true);
[branch0,unst_inds,~,~]=br_stabl(funcs,branch0,0,0);
getpar=@(x,i)arrayfun(@(p)p.parameter(i),x.point);
getx=@(x,i)arrayfun(@(p)p.x(i),x.point);
par_ax=getpar(branch0,in.tau_s);
%
x0_ax=getx(branch0,1);
x0_ax2=getx(branch0,2);
indhopf=indbif1(strcmp(bif1types,'hopf'));
indfold=indbif1(strcmp(bif1types,'fold'));
%% Fix tau_s at 1.4 and continue in one-parameter continuation in tau_c.
[~,it]=min(abs(par_ax-1.4));
clrs=lines();
%%  Symmetry condition on Hopf bifurcation point 
nx=length(x0);
idnet=eye(4);
Mg_1234=[0,1,0,0; 0,0,1,0; 0,0,0,1; 1,0,0,0];
pmfix_1234=dde_stst_lincond('pmfix',nx,'v','trafo',Mg_1234,'rotation',[1,4]);

branch0=br_remove_extracolumns(branch0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now, we compute one-parameter bifurcation branch in $\delat$ 
parbd={'min_bound',[in.c_ext,0;in.tau_s,0; in.tau_c,0; in.delta,-1],...
    'max_bound',[in.c_ext,5;in.tau_s,4; in.tau_c,5; in.delta,1.1],...
    'max_step',[in.c_ext,0.01; in.tau_s,0.01;in.tau_c,0.01; in.delta,0.01; 0,0.02]};
% One-parameter continuation in $\delta$
[funcs_dt,branch0_dt,sucdt]=ChangeBranchParameters(funcs,branch0,it,'contpar',in.delta,'outputfuncs',true,parbd{:});
figure(6)
clf;
[branch0_dt,sdt,fdt,rdt]=br_contn(funcs_dt,branch0_dt,1000);
branch0_dt=br_rvers(branch0_dt);
[branch0_dt,sdt0,fdt0,rdt0]=br_contn(funcs_dt,branch0_dt,1000);
delta_x=getpar(branch0_dt,in.delta);
%% Then we compute one-paramter bifurcation in $tau_c$ at fixed $\delta=0.2$
[~,it2]=min(abs(delta_x-0.2));
[funcs_cdt,branch0_tauc_dt,succdt]=ChangeBranchParameters(funcs_dt,branch0_dt,it2,'contpar',in.tau_c,'outputfuncs',true,parbd{:});
figure(1)
hold on
[branch0_tauc_dt,scdt,fcdt,rcdt]=br_contn(funcs_cdt,branch0_tauc_dt,1000);
branch0_tauc_dt=br_rvers(branch0_tauc_dt);
[branch0_tauc_dt,scdt0,fcdt0,rcdt0]=br_contn(funcs_cdt,branch0_tauc_dt,1000);
%
[branch0_tauc_dt,unst_indscdt,~,~]=br_stabl(funcs_cdt,branch0_tauc_dt,0,0);
bifp_dt=find(diff(unst_indscdt));
% Bisection
[branch0_tauc_dt_bis,indbifcdt,indmapcdt,notcorrectedcdt]=br_bisection(funcs_cdt,branch0_tauc_dt,[bifp_dt(end),bifp_dt(end)+1],'hopf');
[branch0_tauc_dt_bis,unst_indscdt,~,~]=br_stabl(funcs_cdt,branch0_tauc_dt_bis,0,0);
par_axcdt=getpar(branch0_tauc_dt_bis,in.tau_c);
x0_axcdt=getx(branch0_tauc_dt_bis,1);
%% One-parameter bifurcation brnanch in $\tau_c$
parbd_c={'min_bound',[in.c_ext,0;in.tau_s,0; in.tau_c,0],...
    'max_bound',[in.c_ext,6;in.tau_s,5; in.tau_c,5],...
    'max_step',[in.c_ext,0.05; in.tau_s,0.05;in.tau_c,0.1; 0,0.1]};
% Symmetry condition on POs 
psfix_1234=dde_psol_lincond('psfix',nx,'profile',...
    'trafo',Mg_1234,'shift',[1,4],'condprojint',linspace(0.1,0.2,3)'*[1,1]);
[fpsol_dt,psol_dt,sucp_dt]=SetupPsol(funcs_cdt,branch0_tauc_dt_bis,indbifcdt,'contpar',in.tau_c,'extracolumns','auto','initcond',pmfix_1234,...
    'outputfuncs',true,'extra_condition',true,'usercond',psfix_1234,'intervals',60,'degree',4,...
    parbd_c{:},'print_residual_info',1,'matrix','sparse','remesh',false);
figure(1)
hold on
psol_dt=br_contn(fpsol_dt,psol_dt,300);
% Compute stability of POs
[psol_dt,unst_po_dt,dom_po_dt,triv_defect_po_dt]=br_stabl(fpsol_dt,psol_dt,0,1,'exclude_trivial',true);%,'locate_trivial',@(p)[1,1],'geteigenfuncs',true);
x_taus_po_dt=arrayfun(@(x)x.parameter(in.tau_c),psol_dt.point);
ymxs_po_dt=arrayfun(@(x)max(x.profile(1,:)),psol_dt.point);
ymins_po_dt=arrayfun(@(x)min(x.profile(1,:)),psol_dt.point);
branch0_tauc_dt=br_remove_extracolumns(branch0_tauc_dt);
branch0_tauc_dt_bis=br_remove_extracolumns(branch0_tauc_dt_bis);
psol_dt=br_remove_extracolumns(psol_dt);
%% %%%%%%%%%%%%%%%%%%%%%%%
branch0=br_remove_extracolumns(branch0);
%%
mk=15; % MarkerSise
lwidth={'LineWidth',6};
figure(10)
clf;
tiledlayout(2,4,'TileSpacing','compact')
nexttile([2,2])
hold on; grid on 
%%%%%%% plot equilibria %%%%%%%%%%
%%% 
plt_eqstb=plot(par_axcdt(unst_indscdt==0),x0_axcdt(unst_indscdt==0),'-','Color','k','LineWidth',4);
plt_equnstb=plot(par_axcdt(unst_indscdt==6),x0_axcdt(unst_indscdt==6),'k-.','LineWidth',4);
plot(par_axcdt(unst_indscdt>=7),x0_axcdt(unst_indscdt>=7),'.','Color',[0.7 0.7 0.7],'MarkerSize', 8,'LineWidth',4)
% %%%% Plotting POs stability  %%%%%%
bifloc_dt=find(diff(unst_po_dt));
plt_pod1stb=plot(x_taus_po_dt(bifloc_dt(2):bifloc_dt(3)),ymxs_po_dt(bifloc_dt(2):bifloc_dt(3)),'-','Color',clrs(7,:),lwidth{:},'markersize',mk);
plt_pod1unstb=plot(x_taus_po_dt(1:bifloc_dt(2)),ymxs_po_dt(1:bifloc_dt(2)),'-','Color',[0.7 0.7 0.7],lwidth{:},'markersize',mk);
plot(x_taus_po_dt(bifloc_dt(2):bifloc_dt(3)),ymins_po_dt(bifloc_dt(2):bifloc_dt(3)),'-','Color',clrs(7,:),lwidth{:},'markersize',mk)
plot(x_taus_po_dt(1:bifloc_dt(2)),ymins_po_dt(1:bifloc_dt(2)),'-','Color',[0.7 0.7 0.7],lwidth{:},'markersize',mk)
plot(x_taus_po_dt(bifloc_dt(3):end),ymxs_po_dt(bifloc_dt(3):end),'-','Color',[0.7 0.7 0.7],lwidth{:},'markersize',mk)
plot(x_taus_po_dt(bifloc_dt(3):end),ymins_po_dt(bifloc_dt(3):end),'-','Color',[0.7 0.7 0.7],lwidth{:},'markersize',mk)
%%%% Plot Hopf-Bifurcation points %%%%%
plt_hpf=plot(par_axcdt(indbifcdt),x0_axcdt(indbifcdt),'s','Marker', 's', 'MarkerSize', 12, ...
    'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black');
pltnan=plot(NaN,NaN,'Color',[1,1,1]);
legend_text={'POs: $\Pi\sim (1234)_{1/4}$','stst $\#$ unst=0','stst $\#$ unst$\geq 6$','equivariant Hopf-bifurcation','stable POs',...
   'unstable POs'};
vec_plt=[pltnan(1),plt_eqstb(1),plt_equnstb(1),plt_hpf(1),plt_pod1stb(1),plt_pod1unstb(1),pltnan(1)];%...
    %plot(NaN, NaN,'Color',[1,1,1])];
legend(vec_plt,legend_text,'Interpreter','latex','FontSize',24,'Location','northeastoutside')
ylabel('$x_{1}$','Interpreter','latex','FontName','Cambria','FontSize',26)
xlabel('$\tau_{c}$','Interpreter','latex','FontSize',26,'FontName','Cambria')
%xlim([0.3,1])
set(gca, 'FontSize',26,'FontWeight','bold','LineWidth',2)
nexttile([1,2])
[~,it2]=min(abs(x_taus_po_dt-3.5)); % % tau_c=3
po2=psol_dt.point(it2);
plot(po2.mesh*po2.period,po2.profile,'LineWidth',3)
title('\tau_c=3.5','FontSize',18)
hold on; grid on
%% Hopf bifurcation in $(\tau_c,\tau_s)$  with applying symmetry extension
pmfix_1234=dde_stst_lincond('pmfix',nx,'v',...
    'trafo',Mg_1234,'rotation',[1,4]);
[fhopf_dt,hopf_branch_dt,such]=SetupHopf(funcs,branch0_tauc_dt_bis,indbifcdt,'contpar',[in.tau_c,in.tau_s],...
    'dir',in.tau_s,'print_residual_info',1,'outputfuncs',true,'print_residual_info',1,...
    'extracolumns','auto','initcond',pmfix_1234,'extra_condition',true,'usercond',pmfix_1234,parbd{:});
figure(3)
clf
hold on%clf
hopf_branch_dt=br_contn(fhopf_dt,hopf_branch_dt,250);
hopf_branch_dt=br_rvers(hopf_branch_dt);
hopf_branch_dt=br_contn(fhopf_dt,hopf_branch_dt,250);
%
[hopf_branch_dt,unst_hopf_dt,dom_hopf,triv_defect_hopf]=br_stabl(fhopf_dt,hopf_branch_dt,0,0,'exclude_trivial',true,...
    'locate_trivial',@(p)1i*p.omega*[1,-1,1,-1,1,-1]);
x_taus_hf_dt=arrayfun(@(x)x.parameter(in.tau_s),hopf_branch_dt.point);
x_tauc_hf_dt=arrayfun(@(x)x.parameter(in.tau_c),hopf_branch_dt.point);
%%
figure(33)
clf;
tiledlayout(5,5,"TileSpacing","compact");
nexttile([5,3])
hold on; grid on
plt_eqstb_dt=plot(par_axcdt(unst_indscdt==0),x0_axcdt(unst_indscdt==0),'-','Color','k','LineWidth',3);
plt_equnstb_dt=plot(par_axcdt(unst_indscdt==6),x0_axcdt(unst_indscdt==6),'k--','LineWidth',3);
plot(par_axcdt(unst_indscdt>=7),x0_axcdt(unst_indscdt>=7),'.','Color',[0.7 0.7 0.7],'MarkerSize', 8,'LineWidth',2)
%%%%%% Plotting POs stability  %%%%%%
lwidth={'LineWidth',5};
% case delta=1
plt_pod1stb=plot(x_taus_po_dt(unst_po_dt==0),ymxs_po_dt(unst_po_dt==0),'-','Color',clrs(1,:),lwidth{:});
plt_pod1unstb=plot(x_taus_po_dt(unst_po_dt>=1),ymxs_po_dt(unst_po_dt>=1),'-.','Color',[0.7 0.7 0.7],lwidth{:});
plot(x_taus_po_dt(unst_po_dt==0),ymins_po_dt(unst_po_dt==0),'-','Color',clrs(1,:),lwidth{:})
plot(x_taus_po_dt(unst_po_dt>=1),ymins_po_dt(unst_po_dt>=1),'-.','Color',[0.7 0.7 0.7],lwidth{:})
%%%% Plot Hopf-Bifurcation points %%%%%
plt_hpf_dt=plot(par_axcdt(indbifcdt),x0_axcdt(indbifcdt),'s','Marker', 's', 'MarkerSize', 12, ...
    'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black');
legend_text_delta={'stable equilibria','unstable equilibria','equivariant Hopf-bifurcation',...
    'stable POs','unstable POs','POs: $\Pi\sim (1234)_{1/4}$'};
vec_plt_delta=[plt_eqstb_dt(1),plt_equnstb_dt(1),plt_hpf_dt(1),plt_pod1stb(1),plt_pod1unstb(1),...
     plot(NaN, NaN,'Color',[1,1,1])];
legend(vec_plt_delta,legend_text_delta,'Interpreter','latex','FontSize',24)%,'Location','northeastoutside')
ylabel('$x_{1}$','Interpreter','latex')%,'FontName','Cambria','FontSize',26)
xlabel('$\tau_{c}$','Interpreter','latex')%,'FontSize',26,'FontName','Cambria')
xlim([0.3,1])
title('(a)')
set(gca, 'FontSize',30,'FontWeight','bold','FontName','Courier','LineWidth',2,'Box', 'on')
%%%%% Hopf Bifurcation%%%%%%%%
nexttile([5,2])
hold on; grid on
plot(x_tauc_hf_dt(unst_hopf_dt==0),x_taus_hf_dt(unst_hopf_dt==0),'k.','LineWidth',4)
plot(x_tauc_hf_dt(unst_hopf_dt>=1),x_taus_hf_dt(unst_hopf_dt>=1),'r.','MarkerSize',10)
plt_hpf_dt=plot(par_axcdt(indbifcdt),branch0_tauc_dt_bis.point(10).parameter(in.tau_s),'s','Marker', 's', 'MarkerSize', 12, ...
    'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black');
legend('Hopf-bifurcation','starting point','Interpreter','latex')
ylabel('$\tau_{s}$','Interpreter','latex')
xlabel('$\tau_{c}$','Interpreter','latex')
title('(b)')
set(gca, 'FontSize',30,'FontWeight','bold','FontName','Courier','LineWidth',2,'Box', 'on')
hopf_branch_dt=br_remove_extracolumns(hopf_branch_dt);
psol_dt=br_remove_extracolumns(psol_dt);
%%
save('bif_analysis_S4_system_first.mat')
%%
% [~,it]=min(abs(x_taus_po_dt-1)); % % tau_c=3
% po=psol_dt.point(20);
% 
% fd23=@(t,x,Xd)funcs_dt.sys_rhs([x,Xd],po.parameter);
% his=@(t)dde_coll_eva(po.profile-0.1,po.mesh,1+t/po.period,po.degree);
% 
% sold_23=dde23(fd23, [po.parameter(in.tau_s),po.parameter(in.tau_c)] ,his,[0,30],ddeset('RelTol',1e-7,'AbsTol',1e-7));
% 
% figure(90)
% clf
% tiledlayout(1,2)   
% 
% nexttile
% hold on; grid on
% plot(sold_23.x,sold_23.y(1:2,:),'-',...
%     sold_23.x,sold_23.y(3:4,:),'-','LineWidth',3);
% title('simulation near PO')
% nexttile
% hold on; grid on
% plot(po.mesh*po.period,po.profile,'LineWidth',3)
% title('time profile for nera initial PO')
% %%
% 
% figure(7777)
% clf
% plot(real(dom_po_dt))
% yline(1 )
