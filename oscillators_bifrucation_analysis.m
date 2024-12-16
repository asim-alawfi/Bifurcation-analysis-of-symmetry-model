clear;
base=[pwd(),'\..\..\DDE_Biftool_Feb2024\'];
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
    'max_bound',[in.c_ext,3;in.tau_s,4; in.tau_c,4],...
    'max_step',[in.c_ext,0.05; in.tau_s,0.01;in.tau_c,0.01; 0,0.01]};
funcs=set_symfuncs(@sym_oscillators,'sys_tau',@()[in.tau_s,in.tau_c]);
%% compute and find trivial equilibri
[branch0,suc0]=SetupStst(funcs,...
    'parameter',par,'x',x0,...
    'contpar',in.tau_s,'step',0.01,parbd{:},...
    'print_residual_info',true);
%
branch0.method.point.newton_max_iterations=10;
figure(1)
clf;
[branch0,ss1,fs1,rs1]=br_contn(funcs,branch0,20000);
branch0=br_rvers(branch0);
[branch0,ss01,fs01,rs01]=br_contn(funcs,branch0,20000);
%%
[branch0,unst_inds,~,~]=br_stabl(funcs,branch0,0,0);
[branch0,bif1testfuncs,indbif1,bif1types]=LocateSpecialPoints(funcs,branch0,'debug',true);
%%
[branch0,unst_inds,~,~]=br_stabl(funcs,branch0,0,0);
getpar=@(x,i)arrayfun(@(p)p.parameter(i),x.point);
getx=@(x,i)arrayfun(@(p)p.x(i),x.point);
par_ax=getpar(branch0,in.tau_s);
%
x0_ax=getx(branch0,1);
x0_ax2=getx(branch0,2);
indhopf=indbif1(strcmp(bif1types,'hopf'));
indfold=indbif1(strcmp(bif1types,'fold'));
%%
figure(1)
clf;
hold on; grid on
plot(par_ax(unst_inds==0),x0_ax(unst_inds==0),'g.',par_ax(unst_inds==1),x0_ax(unst_inds==1),'k.','LineWidth',5);
plot(par_ax(unst_inds>=2),x0_ax(unst_inds>=2),'.','Color',[0.7 0.7 0.7],'LineWidth',5)
plot(par_ax(indhopf(1)),x0_ax(indhopf(1)),'s','Marker', 's', 'MarkerSize', 8, ...
    'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black');
%xlim([-0.5,3])
legend('stst $\#$ unst=0','stst $\#$ unst=1','hopf','Interpreter','latex','FontSize',12)
ylabel('$x_{1}$','Interpreter','latex','FontName','Cambria',FontSize=16)
xlabel('$\tau_{s}$','Interpreter','latex','FontSize',16,'FontName','Cambria')
%axis tight
set(gca, 'FontWeight','bold')
%save('computation_s4_group_p1.mat')
%% Fix tau_s at 1.4 and continue in one-parameter continuation in tau_c.
[~,it]=min(abs(par_ax-1.4));
[funcs_c,branch0_tauc,succ9]=ChangeBranchParameters(funcs,branch0,it,'contpar',in.tau_c,'outputfuncs',true,parbd{:});
figure(2)
clf;
[branch0_tauc,sc9,fc9,rc9]=br_contn(funcs_c,branch0_tauc,20000);
branch0_tauc=br_rvers(branch0_tauc);
[branch0_tauc,sc8,fc8,rc8]=br_contn(funcs_c,branch0_tauc,20000);
%%
[branch0_tauc,unst_indsc,~,~]=br_stabl(funcs_c,branch0_tauc,0,0);
bifp=find(diff(unst_indsc));
%% Bisection
[branch0_tauc_bis,indbifc,indmapc,notcorrectedc]=br_bisection(funcs_c,branch0_tauc,[bifp(end),bifp(end)+1],'hopf');
%%
[branch0_tauc_bis,unst_indsc2,~,~]=br_stabl(funcs_c,branch0_tauc_bis,0,0);
bifp2=find(diff(unst_indsc2));
par_axc=getpar(branch0_tauc_bis,in.tau_c);
x0_axc=getx(branch0_tauc_bis,1);
%%
clrs=lines();
figure(2)
clf;
hold on; %grid on
plot(par_axc(unst_indsc2==0),x0_axc(unst_indsc2==0),'-','Color','k','LineWidth',3)
plot(par_axc(unst_indsc2==6),x0_axc(unst_indsc2==6),'k--','LineWidth',3);
plot(par_axc(unst_indsc2>=7),x0_axc(unst_indsc2>=7),'.','Color',[0.7 0.7 0.7],'MarkerSize', 8,'LineWidth',2)
plot(par_axc(indbifc),x0_axc(indbifc),'s','Marker', 's', 'MarkerSize', 10, ...
    'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black');
xlim([0,1])
ylim([-0.1,0.1])
xticks([0,0.5,1])
yticks([-0.1,0,0.1])
legend('stst $\#$ unst=0','stst $\#$ unst=6','','Triple Hopf-Bif','Interpreter','latex','FontSize',12)
ylabel('$x_{1}$','Interpreter','latex','FontName','Cambria',FontSize=22)
xlabel('$\tau_{c}$','Interpreter','latex','FontSize',22,'FontName','Cambria')
set(gca, 'FontWeight','bold')
%%
%% Hopf bifurcation computation with applying symmetry extension
nx=length(x0);
idnet=eye(4);
Mg2=[0,1,0,0; 0,0,1,0; 0,0,0,1; 1,0,0,0];
pmfix=dde_stst_lincond('pmfix',nx,'v',...
    'trafo',Mg2,'rotation',[1,4]);
[fhopf,hopf_branch,such]=SetupHopf(funcs,branch0_tauc_bis,indbifc,'contpar',[in.tau_c,in.tau_s],...
    'dir',in.tau_s,'print_residual_info',1,'outputfuncs',true,'print_residual_info',1,...
    'extracolumns','auto','initcond',pmfix,'extra_condition',true,'usercond',pmfix,parbd{:});
figure(3)
clf
hold on%clf
hopf_branch=br_contn(fhopf,hopf_branch,300);
hopf_branch=br_rvers(hopf_branch);
hopf_branch=br_contn(fhopf,hopf_branch,300);
%%
[hopf_branch,unst_hopf,dom_hopf,triv_defect_hopf]=br_stabl(fhopf,hopf_branch,0,0,'exclude_trivial',true,...
    'locate_trivial',@(p)1i*p.omega*[1,-1,1,-1,1,-1]);
x_taus_hf=arrayfun(@(x)x.parameter(in.tau_s),hopf_branch.point);
x_tauc_hf=arrayfun(@(x)x.parameter(in.tau_c),hopf_branch.point);
%%
figure(3)
clf;
hold on; grid on
plot(x_tauc_hf(unst_hopf==0),x_taus_hf(unst_hopf==0),'k.','MarkerSize',10)
plot(x_tauc_hf(unst_hopf>=1),x_taus_hf(unst_hopf>=1),'.','MarkerSize',10)
legend('Hopf-bif $\#$ unst=0','Hopf-bif $\#$ unst$\geq 1$','Interpreter','latex','FontSize',16)
ylabel('$\tau_{s}$','Interpreter','latex','FontName','Cambria',FontSize=22)
xlabel('$\tau_{c}$','Interpreter','latex','FontSize',22,'FontName','Cambria')
set(gca, 'FontWeight','bold')
%%
hopf_branch=br_remove_extracolumns(hopf_branch);
branch0_tauc_bis=br_remove_extracolumns(branch0_tauc_bis);
branch0_tauc=br_remove_extracolumns(branch0_tauc);
branch0=br_remove_extracolumns(branch0);
save('part1_S4permutations.mat')
%%
parbd_c={'min_bound',[in.c_ext,0;in.tau_s,0; in.tau_c,0],...
    'max_bound',[in.c_ext,6;in.tau_s,5; in.tau_c,5],...
    'max_step',[in.c_ext,0.05; in.tau_s,0.05;in.tau_c,0.1; 0,0.1]};
psfix=dde_psol_lincond('psfix',nx,'profile',...
    'trafo',Mg2,'shift',[1,4],'condprojint',linspace(0.1,0.2,2)'*[1,1]);
[fpsol,psol,sucp]=SetupPsol(fhopf,branch0_tauc_bis,indbifc,'contpar',in.tau_c,'extracolumns','auto','initcond',pmfix,...
    'outputfuncs',true,'extra_condition',true,'usercond',psfix,'intervals',60,'degree',4,...
    parbd_c{:},'print_residual_info',1,'matrix','sparse','remesh',false);
%
figure(4)
clf
hold on
psol=br_contn(fpsol,psol,200);
% Compute stability of POs
[psol,unst_po,dom_po,triv_defect_po]=br_stabl(fpsol,psol,0,1,'exclude_trivial',true,'locate_trivial',@(p)[1,1],'geteigenfuncs',true);
x_taus_po=arrayfun(@(x)x.parameter(in.tau_c),psol.point);
%Plot stability of POs
ymxs_po=arrayfun(@(x)max(x.profile(1,:)),psol.point);
ymins_po=arrayfun(@(x)min(x.profile(1,:)),psol.point);
%%
psol=br_remove_extracolumns(psol);
figure(4)
clf;
hold on; grid on
plot(par_axc(unst_indsc2==0),x0_axc(unst_indsc2==0),'-','Color','k','LineWidth',3)
plot(par_axc(unst_indsc2==6),x0_axc(unst_indsc2==6),'k--','LineWidth',3);
plot(par_axc(unst_indsc2>=7),x0_axc(unst_indsc2>=7),'.','Color',[0.7 0.7 0.7],'MarkerSize', 8,'LineWidth',2)
plot(x_taus_po(unst_po==0),ymxs_po(unst_po==0),'ko',...
    x_taus_po(unst_po>=1),ymxs_po(unst_po>=1),'ro','LineWidth',2)
plot(x_taus_po(unst_po==0),ymins_po(unst_po==0),'ko',...
    x_taus_po(unst_po>=1),ymins_po(unst_po>=1),'ro','LineWidth',2)
plot(par_axc(indbifc),x0_axc(indbifc),'s','Marker', 's', 'MarkerSize', 12, ...
    'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black');
xlim([0,2])
%ylim([-0.1,0.1])
%xticks([0,0.5,1])
%yticks([-0.1,0,0.1])
legend('stst $\#$ unst=0','stst $\#$ unst=6','','stable POs','','triple Hopf-bif','Interpreter','latex','FontSize',12)
ylabel('$x_{1}$','Interpreter','latex','FontName','Cambria',FontSize=22)
xlabel('$\tau_{c}$','Interpreter','latex','FontSize',22,'FontName','Cambria')
set(gca, 'FontWeight','bold')
%

%%
Mg3=[0,0,1,0; 0,0,0,1; 1,0,0,0; 0,1,0,0];
pmfix3=dde_stst_lincond('pmfix',nx,'v','trafo',Mg3,'rotation',[1,2]);
psfix3=dde_psol_lincond('psfix',nx,'profile','trafo',Mg3,'shift',[1,2],'condprojint',linspace(0.1,0.2,2)'*[1,1]);
[fpsol2,psol2,sucp2]=SetupPsol(fhopf,branch0_tauc_bis,indbifc,'contpar',in.tau_c,'extracolumns','auto','initcond',pmfix3,...
    'outputfuncs',true,'extra_condition',true,'usercond',psfix3,'intervals',60,'degree',4,...
    parbd_c{:},'print_residual_info',1,'matrix','sparse','remesh',false);
%
figure(5)
hold on
psol2=br_contn(fpsol2,psol2,200);
%% Compute stability of POs
[psol2,unst_po2,dom_po2,triv_defect_po2]=br_stabl(fpsol2,psol2,0,1,'exclude_trivial',true,'locate_trivial',@(p)[1,1],'geteigenfuncs',true);
x_taus_po2=arrayfun(@(x)x.parameter(in.tau_c),psol2.point);
% Plot stability of POs
ymxs_po2=arrayfun(@(x)max(x.profile(1,:)),psol2.point);
ymins_po2=arrayfun(@(x)min(x.profile(1,:)),psol2.point);
figure(5)
%clf
hold on; grid on
plot(x_taus_po2(unst_po2==0),ymxs_po2(unst_po2==0),'bo',...
    x_taus_po2(unst_po2>=1),ymxs_po2(unst_po2>=1),'ro','LineWidth',2)
plot(x_taus_po2(unst_po2==0),ymins_po2(unst_po2==0),'bo',...
    x_taus_po2(unst_po2>=1),ymins_po2(unst_po2>=1),'ro','LineWidth',2)
%%
figure(12)
clf;
tiledlayout(6,2,"TileSpacing","compact");
nexttile([4,1])
hold on; grid on
plot(par_axc(unst_indsc2==0),x0_axc(unst_indsc2==0),'-','Color','k','LineWidth',3)
plot(par_axc(unst_indsc2==6),x0_axc(unst_indsc2==6),'k--','LineWidth',3);
plot(par_axc(unst_indsc2>=7),x0_axc(unst_indsc2>=7),'.','Color',[0.7 0.7 0.7],'MarkerSize', 8,'LineWidth',2)
plot(x_taus_po(unst_po==0),ymxs_po(unst_po==0),'ko',...
    x_taus_po(unst_po>=1),ymxs_po(unst_po>=1),'ro','LineWidth',2)
plot(x_taus_po(unst_po==0),ymins_po(unst_po==0),'ko',...
    x_taus_po(unst_po>=1),ymins_po(unst_po>=1),'ro','LineWidth',2)
plot(par_axc(indbifc),x0_axc(indbifc),'s','Marker', 's', 'MarkerSize', 12, ...
    'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black');
xlim([0,2])
%ylim([-0.1,0.1])
%xticks([0,0.5,1])
%yticks([-0.1,0,0.1])
legend('stst $\#$ unst=0','stst $\#$ unst=6','','stable POs','','triple Hopf-bif','Interpreter','latex','FontSize',18)
ylabel('$x_{1}$','Interpreter','latex','FontName','Cambria',FontSize=22)
xlabel('$\tau_{c}$','Interpreter','latex','FontSize',22,'FontName','Cambria')
set(gca, 'FontWeight','bold')
title('(a)','FontSize',16,'FontName','Cambria')
nexttile([4,1])
hold on; grid on
plot(par_axc(unst_indsc2==0),x0_axc(unst_indsc2==0),'-','Color','k','LineWidth',3)
plot(par_axc(unst_indsc2==6),x0_axc(unst_indsc2==6),'k--','LineWidth',3);
plot(par_axc(unst_indsc2>=7),x0_axc(unst_indsc2>=7),'.','Color',[0.7 0.7 0.7],'MarkerSize', 8,'LineWidth',2)
plot(x_taus_po2(unst_po2==0),ymxs_po2(unst_po2==0),'bo',...
    x_taus_po2(unst_po2>=1),ymxs_po2(unst_po2>=1),'ro','LineWidth',2)
plot(x_taus_po2(unst_po2==0),ymins_po2(unst_po2==0),'bo',...
    x_taus_po2(unst_po2>=1),ymins_po2(unst_po2>=1),'ro','LineWidth',2)
plot(par_axc(indbifc),x0_axc(indbifc),'s','Marker', 's', 'MarkerSize', 12, ...
    'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black');
xlim([0,2])
%ylim([-0.1,0.1])
%xticks([0,0.5,1])
%yticks([-0.1,0,0.1])
legend('stst $\#$ unst=0','stst $\#$ unst=6','','stable POs','','triple Hopf-bif','Interpreter','latex','FontSize',18)
ylabel('$x_{1}$','Interpreter','latex','FontName','Cambria',FontSize=22)
xlabel('$\tau_{c}$','Interpreter','latex','FontSize',22,'FontName','Cambria')
set(gca, 'FontWeight','bold')
title('(b)','FontSize',16,'FontName','Cambria')
nexttile([2,1])
po1=psol.point(20);
plot(po1.mesh*po1.period,po1.profile,'LineWidth',2)
xlabel('$t$','Interpreter','latex','FontSize',22,'FontName','Cambria')
xlim([0,po1.period])
ylim([-1,1])
set(gca, 'FontWeight','bold')
grid on
title('(c)','FontSize',16,'FontName','Cambria')
nexttile([2,1])
po2=psol2.point(20);
plot(po2.mesh*po2.period,po2.profile,'LineWidth',2)
legend('$x_{1}$','$x_{2}$','$x_{3}$','$x_{4}$','Interpreter','latex','FontSize',18,'Location','northeastoutside')
xlim([0,po2.period])
ylim([-1,1])
xlabel('$t$','Interpreter','latex','FontSize',22,'FontName','Cambria')
set(gca, 'FontWeight','bold')
grid on
title('(d)','FontSize',16,'FontName','Cambria')
%%
psol2=br_remove_extracolumns(psol2);
save('par3_S4permutations.mat')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%% No, we continue Hopf bifurcation in (tau_s,delta)-space%%%%%%%%%%
%% No, we continue Hopf bifurcation in (tau_s,delta)-space
% Then we fix delta and continue one-parameter continuation in tau_c for
% periodic solutions
parbd={'min_bound',[in.c_ext,0;in.tau_s,0; in.tau_c,0; in.delta,-2],...
    'max_bound',[in.c_ext,3;in.tau_s,4; in.tau_c,4; in.delta,2],...
    'max_step',[in.c_ext,0.05; in.tau_s,0.01;in.tau_c,0.01; in.delta,0.01; 0,0.01]};
[fhopf_delta,hopf_branch_delta,such_delta]=SetupHopf(funcs,branch0_tauc_bis,indbifc,'contpar',[in.tau_c,in.delta],...
    'dir',in.delta,'print_residual_info',1,'outputfuncs',true,'print_residual_info',1,...
    'extracolumns','auto','initcond',pmfix,'extra_condition',true,'usercond',pmfix,parbd{:});
figure(6)
clf
hold on
hopf_branch_delta=br_contn(fhopf_delta,hopf_branch_delta,300);
hopf_branch_delta=br_rvers(hopf_branch_delta);
hopf_branch_delta=br_contn(fhopf_delta,hopf_branch_delta,300);
%%
[hopf_branch_delta,unst_hopf_delta,dom_hopf_delta,triv_defect_hopf_delta]=br_stabl(fhopf_delta,hopf_branch_delta,0,0,'exclude_trivial',true,...
    'locate_trivial',@(p)1i*p.omega*[1,-1,1,-1,1,-1]);
%%
x_delt_hf_delta=arrayfun(@(x)x.parameter(in.delta),hopf_branch_delta.point);
x_tauc_hf_delta=arrayfun(@(x)x.parameter(in.tau_c),hopf_branch_delta.point);
figure(6)
clf;
hold on; grid on
plot(x_tauc_hf_delta(unst_hopf_delta==0),x_delt_hf_delta(unst_hopf_delta==0),'k.','MarkerSize',10)
plot(x_tauc_hf_delta(unst_hopf_delta>=1),x_delt_hf_delta(unst_hopf_delta>=1),'.','MarkerSize',10)
legend('Hopf-bif $\#$ unst=0','Hopf-bif $\#$ unst$\geq 1$','Interpreter','latex','FontSize',16)
ylabel('$\delta$','Interpreter','latex','FontName','Cambria',FontSize=22)
xlabel('$\tau_{c}$','Interpreter','latex','FontSize',22,'FontName','Cambria')
set(gca, 'FontWeight','bold')
%% We fixed delat=1.5 and branch off from hopf bifurcation tracking POs in tau_c
[~,ind_d1]=min(abs(x_delt_hf_delta-1.5));
hopf_branch=br_remove_extracolumns(hopf_branch);
%save('Hopf_br_sym_impose.mat')
%
[fpsol3,psol3,sucp3]=SetupPsol(fhopf_delta,hopf_branch_delta,ind_d1,'contpar',in.tau_c,'extracolumns','auto','initcond',pmfix,...
    'outputfuncs',true,'extra_condition',true,'usercond',psfix,'intervals',60,'degree',4,...
    parbd_c{:},'print_residual_info',1,'matrix','sparse','remesh',false);
figure(7)
clf
psol3=br_contn(fpsol3,psol3,200);
%% Compute stability of POs
[psol3,unst_po3,dom_po3,triv_defect_po3]=br_stabl(fpsol3,psol3,0,1,'exclude_trivial',true,'locate_trivial',@(p)[1,1],'geteigenfuncs',true);
x_taus_po3=arrayfun(@(x)x.parameter(in.tau_c),psol3.point);
%% Plot stability of POs
ymxs_po3=arrayfun(@(x)max(x.profile(1,:)),psol3.point);
ymins_po3=arrayfun(@(x)min(x.profile(1,:)),psol3.point);
figure(7)
clf
hold on; grid on
plot(x_taus_po3(unst_po3==0),ymxs_po3(unst_po3==0),'bo',...
    x_taus_po3(unst_po3>=1),ymxs_po3(unst_po3>=1),'ro','LineWidth',2)
plot(x_taus_po3(unst_po3==0),ymins_po3(unst_po3==0),'bo',...
    x_taus_po3(unst_po3>=1),ymins_po3(unst_po3>=1),'ro','LineWidth',2)
legend('stable PO: $\Pi:\sim (1234)_{1/4}$','stable PO: $\Pi:\sim (1234)_{1/4}$','','','Interpreter','latex','FontSize',16)
title('$\delta=1.5$','Interpreter','latex','FontSize',20)
set(gca,'FontWeight','bold')
%%
[~,ind_d2]=min(abs(x_delt_hf_delta+1.5));
[fpsol4,psol4,sucp4]=SetupPsol(fhopf_delta,hopf_branch_delta,ind_d2,'contpar',in.tau_c,'extracolumns','auto','initcond',pmfix,...
    'outputfuncs',true,'extra_condition',true,'usercond',psfix,'intervals',60,'degree',4,...
    parbd_c{:},'print_residual_info',1,'matrix','sparse','remesh',false);
%%
figure(8)
psol4=br_contn(fpsol4,psol4,200);
%% Compute stability of POs
[psol4,unst_po4,dom_po4,triv_defect_po4]=br_stabl(fpsol4,psol4,0,1,'exclude_trivial',true,'locate_trivial',@(p)[1,1],'geteigenfuncs',true);
x_taus_po4=arrayfun(@(x)x.parameter(in.tau_c),psol4.point);
% Plot stability of POs
ymxs_po4=arrayfun(@(x)max(x.profile(1,:)),psol4.point);
ymins_po4=arrayfun(@(x)min(x.profile(1,:)),psol4.point);
%%
figure(8)
clf
hold on; grid on
plot(x_taus_po4(unst_po4==0),ymxs_po4(unst_po4==0),'bo',...
    x_taus_po4(unst_po4>=1),ymxs_po4(unst_po4>=1),'ro','LineWidth',2)
plot(x_taus_po4(unst_po4==0),ymins_po4(unst_po4==0),'bo',...
    x_taus_po4(unst_po4>=1),ymins_po4(unst_po4>=1),'ro','LineWidth',2)
legend('stable PO: $\Pi:\sim (1234)_{1/4}$','stable PO: $\Pi:\sim (1234)_{1/4}$','','','Interpreter','latex','FontSize',16)
title('$\delta=-1.5$','Interpreter','latex','FontSize',20)
set(gca,'FontWeight','bold')
%%
[fpsoln,psoln,sucpn]=SetupPsol(fhopf_delta,hopf_branch_delta,ind_d1,'contpar',in.tau_c,'extracolumns','auto','initcond',pmfix3,...
    'outputfuncs',true,'extra_condition',true,'usercond',psfix3,'intervals',60,'degree',4,...
    parbd_c{:},'print_residual_info',1,'matrix','sparse','remesh',false);
figure(876)
psoln=br_contn(fpsoln,psoln,200);
%%
[psoln,unst_pon,dom_pon,triv_defect_pon]=br_stabl(fpsoln,psoln,0,1,'exclude_trivial',true,'locate_trivial',@(p)[1,1],'geteigenfuncs',true);
x_taus_pon=arrayfun(@(x)x.parameter(in.tau_c),psoln.point);
% Plot stability of POs
ymxs_pon=arrayfun(@(x)max(x.profile(1,:)),psoln.point);
ymins_pon=arrayfun(@(x)min(x.profile(1,:)),psoln.point);
%%
figure(435)
clf
hold on; grid on
plot(x_taus_pon(unst_pon==0),ymxs_pon(unst_pon==0),'bo',...
    x_taus_pon(unst_pon>=1),ymxs_pon(unst_pon>=1),'ro','LineWidth',2)
plot(x_taus_pon(unst_pon==0),ymins_pon(unst_pon==0),'bo',...
    x_taus_pon(unst_pon>=1),ymins_pon(unst_pon>=1),'ro','LineWidth',2)
%%
