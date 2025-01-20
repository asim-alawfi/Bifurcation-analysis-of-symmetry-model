clear;
base=[pwd(),'\..\..\DDE_Biftool2025\'];
addpath([base,'ddebiftool'],...
    [base,'ddebiftool_extra_psol'],...
    [base,'ddebiftool_utilities'],...
    [base,'ddebiftool_extra_rotsym'],...
    [base,'ddebiftool_extra_nmfm'],...
    [base,'ddebiftool_extra_symbolic'],...
    [base,'ddebiftool_coco']);

load('bif_analysis_S4_system_first.mat')
%% Permutation (13)(24) time shift=0
M_13_24=[0,0,1,0; 0,0,0,1; 1,0,0,0; 0,1,0,0]; %(13)(24)
pmfixs1=dde_stst_lincond('pmfix',nx,'v','trafo',M_13_24,'rotation',[0,1]);
psfixs1=dde_psol_lincond('psfix',nx,'profile','trafo',M_13_24,'shift',[0,1],'condprojint',linspace(0.0,0.2,2)'*[1,1]);
[fpsol_dts1,psol_dts1,sucp_dts1]=SetupPsol(funcs_cdt,branch0_tauc_dt_bis,indbifcdt,'contpar',in.tau_c,'extracolumns','auto',...
    'initcond',pmfixs1,'outputfuncs',true,'extra_condition',true,'usercond',psfixs1,'intervals',60,'degree',4,...
     parbd_c{:},'print_residual_info',1,'matrix','sparse','remesh',false);
%%
figure(200)
hold on
[psol_dts1,scdt,fcdt,rcdt]=br_contn(fpsol_dts1,psol_dts1,200);
%%
[psol_dts1,unst_po_dts1,dom_po_dts1,triv_defect_po_dts1]=br_stabl(fpsol_dts1,psol_dts1,0,1,'exclude_trivial',true);%,'locate_trivial',@(p)[1,1],'geteigenfuncs',true);
x_taus_po_dts1=arrayfun(@(x)x.parameter(in.tau_c),psol_dts1.point);
% Plot stability of POs
ymxs_po_dts1=arrayfun(@(x)max(x.profile(1,:)),psol_dts1.point);
ymins_po_dts1=arrayfun(@(x)min(x.profile(1,:)),psol_dts1.point);

% 
%psol_dts1=br_remove_extracolumns(psol_dts1);

%% Permutation (1432) time shift=0
M_1432=[0,0,0,1; 1,0,0,0; 0,1,0,0; 0,0,1,0]; %(12)
pmfixs2=dde_stst_lincond('pmfix',nx,'v','trafo',M_1432,'rotation',[1,4]);
psfixs2=dde_psol_lincond('psfix',nx,'profile','trafo',M_1432,'shift',[1,4],'condprojint',linspace(0.0,0.2,2)'*[1,1]);
[fpsol_dts2,psol_dts2,sucp_dts1]=SetupPsol(funcs_cdt,branch0_tauc_dt_bis,indbifcdt,'contpar',in.tau_c,'extracolumns','auto',...
    'initcond',pmfixs2,'outputfuncs',true,'extra_condition',true,'usercond',psfixs2,'intervals',60,'degree',4,...
     parbd_c{:},'print_residual_info',1,'matrix','sparse','remesh',false);
%
figure(200)
hold on
[psol_dts2,scdt,fcdt,rcdt]=br_contn(fpsol_dts2,psol_dts2,150);
%%
pp=psol_dts2.point(100)
figure(888)
clf
hold on
plot(pp.mesh*pp.period,pp.profile(1,:),'x')
    plot( pp.mesh*pp.period,pp.profile(2,:),'s')
    plot( pp.mesh*pp.period,pp.profile(3,:),'o')
   plot( pp.mesh*pp.period,pp.profile(4,:),'x','LineWidth',1)
%%
[psol_dts2,unst_po_dts2,dom_po_dts2,triv_defect_po_dts2]=br_stabl(fpsol_dts2,psol_dts2,0,1,'exclude_trivial',true);
x_taus_po_dts2=arrayfun(@(x)x.parameter(in.tau_c),psol_dts2.point);
% Plot stability of POs
ymxs_po_dts2=arrayfun(@(x)max(x.profile(1,:)),psol_dts2.point);
ymins_po_dts2=arrayfun(@(x)min(x.profile(1,:)),psol_dts2.point);


%%
pgray=0.5*[1 1 1];
clrscop=copper();
lwidth={'LineWidth',2,'MarkerSize',15};
figure(2201)
clf;
hold on; grid on 
%%%%%%% plot equilibria %%%%%%%%%%
%%%%%%% plot equilibria %%%%%%%%%%
%%% delat=0.2
% text(3.5,1,'$\delta=0.2$','Interpreter','latex','FontSize',20,'FontWeight','bold')
plot(par_axcdt(unst_indscdt==0),x0_axcdt(unst_indscdt==0),'-','Color','k','LineWidth',4)
plot(par_axcdt(unst_indscdt==6),x0_axcdt(unst_indscdt==6),'--','Color',0*[0.7 0.7 0.7],lwidth{:})
plot(par_axcdt(unst_indscdt>=7),x0_axcdt(unst_indscdt>=7),'--','Color',0*[0.7 0.7 0.7],lwidth{:})
%%%
%
hold on; grid on
plot(x_taus_po_dts1(unst_po_dts1==0),ymxs_po_dts1(unst_po_dts1==0),'o','Color',clrs(2,:),lwidth{:});
plot(x_taus_po_dts1(unst_po_dts1==0),ymins_po_dts1(unst_po_dts1==0),'o','Color',clrs(2,:),lwidth{:})
plot(x_taus_po_dts1(unst_po_dts1>=1),ymxs_po_dts1(unst_po_dts1>=1),'-','Color',pgray,lwidth{:})
plot(x_taus_po_dts1(unst_po_dts1>=1),ymins_po_dts1(unst_po_dts1>=1),'-','Color',pgray,lwidth{:})
% %%%
xlabel('\tau_c','FontSize',18)
%
plot(x_taus_po_dts2(unst_po_dts2==0),ymxs_po_dts2(unst_po_dts2==0),'s','Color','k',lwidth{:});
plot(x_taus_po_dts2(unst_po_dts2==0),ymins_po_dts2(unst_po_dts2==0),'s','Color','k',lwidth{:})
plot(x_taus_po_dts2(unst_po_dts2>=1),ymxs_po_dts2(unst_po_dts2>=1),'s','Color',clrs(3,:),lwidth{:})
plot(x_taus_po_dts2(unst_po_dts2>=1),ymins_po_dts2(unst_po_dts2>=1),'s','Color',clrs(3,:),lwidth{:})
%%%
bifloc_dt=find(diff(unst_po_dt));
plt_pod1stb=plot(x_taus_po_dt(bifloc_dt(2):bifloc_dt(3)),ymxs_po_dt(bifloc_dt(2):bifloc_dt(3)),'x','Color',pgray,lwidth{:});
            plot(x_taus_po_dt(bifloc_dt(2):bifloc_dt(3)),ymins_po_dt(bifloc_dt(2):bifloc_dt(3)),'x','Color',pgray,lwidth{:})
            plot(x_taus_po_dt(bifloc_dt(4):end),ymxs_po_dt(bifloc_dt(4):end),'x','Color',pgray,lwidth{:})
            plot(x_taus_po_dt(bifloc_dt(4):end),ymins_po_dt(bifloc_dt(4):end),'x','Color',pgray,lwidth{:})
plt_pod1unstb=plot(x_taus_po_dt(1:bifloc_dt(2)),ymxs_po_dt(1:bifloc_dt(2)),'x','Color',clrs(1,:),lwidth{:});
              plot(x_taus_po_dt(1:bifloc_dt(2)),ymins_po_dt(1:bifloc_dt(2)),'x','Color',clrs(1,:),lwidth{:})
              plot(x_taus_po_dt(bifloc_dt(3):bifloc_dt(4)),ymxs_po_dt(bifloc_dt(3):bifloc_dt(4)),'x','Color',clrs(1,:),lwidth{:});
              plot(x_taus_po_dt(bifloc_dt(3):bifloc_dt(4)),ymins_po_dt(bifloc_dt(3):bifloc_dt(4)),'x','Color',clrs(1,:),lwidth{:});
              %%%%%
% % % plot(x_taus_po_dtss,ymxs_po_dtss,'-','Color',clrs(4,:),lwidth{:});
% % % plot(x_taus_po_dtss,ymins_po_dtss,'-','Color',clrs(4,:),lwidth{:});
% plot(x_taus_po_dtss(unst_po_dtss==0),ymxs_po_dtss(unst_po_dtss==0),'-','Color',pgray,lwidth{:});
% plot(x_taus_po_dtss(unst_po_dtss==0),ymins_po_dtss(unst_po_dtss==0),'-','Color',pgray,lwidth{:})
% plot(x_taus_po_dtss(unst_po_dtss>=1),ymxs_po_dtss(unst_po_dtss>=1),'-','Color',clrs(3,:),lwidth{:})
% plot(x_taus_po_dtss(unst_po_dtss>=1),ymins_po_dtss(unst_po_dtss>=1),'-','Color',clrs(3,:),lwidth{:})
set(gca,'LineWidth',3,'Box','on'); 
% Zoom in
%lwidth={'Linewidth',4};
axes('Position', [0.3 0.35 0.1 0.15]); % Adjust position and size

plot(par_axcdt(unst_indscdt==0),x0_axcdt(unst_indscdt==0),'-','Color','k','LineWidth',4)
plot(par_axcdt(unst_indscdt==6),x0_axcdt(unst_indscdt==6),'--','Color',0*[0.7 0.7 0.7],lwidth{:})
plot(par_axcdt(unst_indscdt>=7),x0_axcdt(unst_indscdt>=7),'--','Color',0*[0.7 0.7 0.7],lwidth{:})
%%%
%
hold on; grid on
plot(x_taus_po_dts1(unst_po_dts1==0),ymxs_po_dts1(unst_po_dts1==0),'o','Color',clrs(2,:),lwidth{:});
plot(x_taus_po_dts1(unst_po_dts1==0),ymins_po_dts1(unst_po_dts1==0),'o','Color',clrs(2,:),lwidth{:})
plot(x_taus_po_dts1(unst_po_dts1>=1),ymxs_po_dts1(unst_po_dts1>=1),'-','Color',pgray,lwidth{:})
plot(x_taus_po_dts1(unst_po_dts1>=1),ymins_po_dts1(unst_po_dts1>=1),'-','Color',pgray,lwidth{:})
% %%%
xlabel('\tau_c','FontSize',18)
%
plot(x_taus_po_dts2(unst_po_dts2==0),ymxs_po_dts2(unst_po_dts2==0),'s','Color','k',lwidth{:});
plot(x_taus_po_dts2(unst_po_dts2==0),ymins_po_dts2(unst_po_dts2==0),'s','Color','k',lwidth{:})
plot(x_taus_po_dts2(unst_po_dts2>=1),ymxs_po_dts2(unst_po_dts2>=1),'s','Color',clrs(3,:),lwidth{:})
plot(x_taus_po_dts2(unst_po_dts2>=1),ymins_po_dts2(unst_po_dts2>=1),'s','Color',clrs(3,:),lwidth{:})
%%%
bifloc_dt=find(diff(unst_po_dt));
plt_pod1stb=plot(x_taus_po_dt(bifloc_dt(2):bifloc_dt(3)),ymxs_po_dt(bifloc_dt(2):bifloc_dt(3)),'x','Color',pgray,lwidth{:});
            plot(x_taus_po_dt(bifloc_dt(2):bifloc_dt(3)),ymins_po_dt(bifloc_dt(2):bifloc_dt(3)),'x','Color',pgray,lwidth{:})
            plot(x_taus_po_dt(bifloc_dt(4):end),ymxs_po_dt(bifloc_dt(4):end),'x','Color',pgray,lwidth{:})
            plot(x_taus_po_dt(bifloc_dt(4):end),ymins_po_dt(bifloc_dt(4):end),'x','Color',pgray,lwidth{:})
plt_pod1unstb=plot(x_taus_po_dt(1:bifloc_dt(2)),ymxs_po_dt(1:bifloc_dt(2)),'x','Color',clrs(1,:),lwidth{:});
              plot(x_taus_po_dt(1:bifloc_dt(2)),ymins_po_dt(1:bifloc_dt(2)),'x','Color',clrs(1,:),lwidth{:})
              plot(x_taus_po_dt(bifloc_dt(3):bifloc_dt(4)),ymxs_po_dt(bifloc_dt(3):bifloc_dt(4)),'x','Color',clrs(1,:),lwidth{:});
              plot(x_taus_po_dt(bifloc_dt(3):bifloc_dt(4)),ymins_po_dt(bifloc_dt(3):bifloc_dt(4)),'x','Color',clrs(1,:),lwidth{:});
              
xlim([0.486766 0.48678]); 
ylim([0.18 0.21]); 
set(gca, 'XTick', [], 'YTick', []); 
box on;
%%%%%%%
%lwidth={'Linewidth',4};
axes('Position', [0.7 0.45 0.2 0.2]); % Adjust position and size
hold on; grid on
plot(par_axcdt(unst_indscdt==0),x0_axcdt(unst_indscdt==0),'-','Color','k','LineWidth',4)
plot(par_axcdt(unst_indscdt==6),x0_axcdt(unst_indscdt==6),'--','Color',0*[0.7 0.7 0.7],lwidth{:})
plot(par_axcdt(unst_indscdt>=7),x0_axcdt(unst_indscdt>=7),'--','Color',0*[0.7 0.7 0.7],lwidth{:})
%%%
%
hold on; grid on
plot(x_taus_po_dts1(unst_po_dts1==0),ymxs_po_dts1(unst_po_dts1==0),'o','Color',clrs(2,:),lwidth{:});
plot(x_taus_po_dts1(unst_po_dts1==0),ymins_po_dts1(unst_po_dts1==0),'o','Color',clrs(2,:),lwidth{:})
plot(x_taus_po_dts1(unst_po_dts1>=1),ymxs_po_dts1(unst_po_dts1>=1),'-','Color',pgray,lwidth{:})
plot(x_taus_po_dts1(unst_po_dts1>=1),ymins_po_dts1(unst_po_dts1>=1),'-','Color',pgray,lwidth{:})
% %%%
xlabel('\tau_c','FontSize',18)
%
plot(x_taus_po_dts2(unst_po_dts2==0),ymxs_po_dts2(unst_po_dts2==0),'s','Color','k',lwidth{:});
plot(x_taus_po_dts2(unst_po_dts2==0),ymins_po_dts2(unst_po_dts2==0),'s','Color','k',lwidth{:})
plot(x_taus_po_dts2(unst_po_dts2>=1),ymxs_po_dts2(unst_po_dts2>=1),'s','Color',clrs(3,:),lwidth{:})
plot(x_taus_po_dts2(unst_po_dts2>=1),ymins_po_dts2(unst_po_dts2>=1),'s','Color',clrs(3,:),lwidth{:})
%%%
bifloc_dt=find(diff(unst_po_dt));
plt_pod1stb=plot(x_taus_po_dt(bifloc_dt(2):bifloc_dt(3)),ymxs_po_dt(bifloc_dt(2):bifloc_dt(3)),'x','Color',pgray,lwidth{:});
            plot(x_taus_po_dt(bifloc_dt(2):bifloc_dt(3)),ymins_po_dt(bifloc_dt(2):bifloc_dt(3)),'x','Color',pgray,lwidth{:})
            plot(x_taus_po_dt(bifloc_dt(4):end),ymxs_po_dt(bifloc_dt(4):end),'x','Color',pgray,lwidth{:})
            plot(x_taus_po_dt(bifloc_dt(4):end),ymins_po_dt(bifloc_dt(4):end),'x','Color',pgray,lwidth{:})
plt_pod1unstb=plot(x_taus_po_dt(1:bifloc_dt(2)),ymxs_po_dt(1:bifloc_dt(2)),'x','Color',clrs(1,:),lwidth{:});
              plot(x_taus_po_dt(1:bifloc_dt(2)),ymins_po_dt(1:bifloc_dt(2)),'x','Color',clrs(1,:),lwidth{:})
              plot(x_taus_po_dt(bifloc_dt(3):bifloc_dt(4)),ymxs_po_dt(bifloc_dt(3):bifloc_dt(4)),'x','Color',clrs(1,:),lwidth{:});
              plot(x_taus_po_dt(bifloc_dt(3):bifloc_dt(4)),ymins_po_dt(bifloc_dt(3):bifloc_dt(4)),'x','Color',clrs(1,:),lwidth{:});
              
xlim([2.8 4.5]); 
ylim([2.18 2.21]); 
set(gca, 'XTick', [], 'YTick', []); 
box on;
%%




%%
figure(888844)
clf;
tiledlayout(3,1)
nexttile
hold on; grid on
[~,it1]=min(abs(x_taus_po_dt-3.5));
p1=psol_dt.point(it1);
plot(p1.mesh*p1.period,p1.profile(1,:),'-','LineWidth',3)
plot(p1.mesh*p1.period,p1.profile(2,:),'-','LineWidth',3)
plot(p1.mesh*p1.period,p1.profile(3,:),'-','LineWidth',3)
plot(p1.mesh*p1.period,p1.profile(4,:),'-','LineWidth',3)
%
nexttile
hold on; grid on
[~,its1]=min(abs(x_taus_po_dts1-3.5));
p2=psol_dts1.point(its1);
plot(p2.mesh*p2.period,p2.profile(1,:),'s')
plot(p2.mesh*p2.period,p2.profile(2,:),'o')
plot(p2.mesh*p2.period,p2.profile(3,:),'x')
plot(p2.mesh*p2.period,p2.profile(4,:),'.','LineWidth',1)
nexttile
hold on; grid on
[~,its2]=min(abs(x_taus_po_dts2-3.5));
p3=psol_dts2.point(its2);
plot(p3.mesh*p3.period,p3.profile(1,:),'-','LineWidth',3)
plot(p3.mesh*p3.period,p3.profile(2,:),'-','LineWidth',3)
plot(p3.mesh*p3.period,p3.profile(3,:),'-','LineWidth',3)
plot(p3.mesh*p3.period,p3.profile(4,:),'-','LineWidth',3)
%%
pgray=0.5*[1 1 1];
clrscop=copper();
lwidth={'LineWidth',6};%,'MarkerSize',5};
figure(22001)
clf
tiledlayout(3,5,'TileSpacing','compact')
nexttile([1,3])
hold on;
%%%%%%% plot equilibria %%%%%%%%%%
%%% delat=0.2
% text(3.5,1,'$\delta=0.2$','Interpreter','latex','FontSize',20,'FontWeight','bold')
plot(par_axcdt(unst_indscdt==0),x0_axcdt(unst_indscdt==0),'-','Color','k','LineWidth',3)
plot(par_axcdt(unst_indscdt==6),x0_axcdt(unst_indscdt==6),'-.','Color',0*[1 1 1],'LineWidth',3)
plot(par_axcdt(unst_indscdt>=7),x0_axcdt(unst_indscdt>=7),'-.','Color',0*[1 1 1],'LineWidth',3)

bifloc_dt=find(diff(unst_po_dt));
plt_pod1stb=plot(x_taus_po_dt(bifloc_dt(2):bifloc_dt(3)),ymxs_po_dt(bifloc_dt(2):bifloc_dt(3)),'-','Color',clrs(1,:),lwidth{:});
            plot(x_taus_po_dt(bifloc_dt(2):bifloc_dt(3)),ymins_po_dt(bifloc_dt(2):bifloc_dt(3)),'-','Color',clrs(1,:),lwidth{:})
            plot(x_taus_po_dt(bifloc_dt(4):end),ymxs_po_dt(bifloc_dt(4):end),'-','Color',clrs(1,:),lwidth{:})
            plot(x_taus_po_dt(bifloc_dt(4):end),ymins_po_dt(bifloc_dt(4):end),'-','Color',clrs(1,:),lwidth{:})
plt_pod1unstb=plot(x_taus_po_dt(1:bifloc_dt(2)),ymxs_po_dt(1:bifloc_dt(2)),'-','Color',pgray,lwidth{:});
              plot(x_taus_po_dt(1:bifloc_dt(2)),ymins_po_dt(1:bifloc_dt(2)),'-','Color',pgray,lwidth{:})
              plot(x_taus_po_dt(bifloc_dt(3):bifloc_dt(4)),ymxs_po_dt(bifloc_dt(3):bifloc_dt(4)),'-','Color',pgray,lwidth{:});
              plot(x_taus_po_dt(bifloc_dt(3):bifloc_dt(4)),ymins_po_dt(bifloc_dt(3):bifloc_dt(4)),'-','Color',pgray,lwidth{:});
              % plot hopf
              plot(par_axcdt(indbifcdt),x0_axcdt(indbifcdt),'.','Color','r','MarkerSize',40,'LineWidth',3)

[~,it1]=min(abs(x_taus_po_dt-3.5));
plot(x_taus_po_dt(it1),ymxs_po_dt(it1),'k^','MarkerSize',12,'MarkerFaceColor','k')
legend('stable equilibria','unstable equilibria','','stable POs','','','','unstable POs','','','','equivariant Hopf',...
    'interpreter','latex','FontSize',20)
ylim([-3,3])
set(gca,'LineWidth',2,'Box','on','FontSize',14,'FontWeight','bold','XTick',[],'YTick',[-2,0,2])
%
nexttile([1,2])
hold on;
p1=psol_dt.point(it1);
plot(p1.mesh*p1.period,p1.profile(1,:),'-','LineWidth',3)
plot(p1.mesh*p1.period,p1.profile(2,:),'-','LineWidth',3)
plot(p1.mesh*p1.period,p1.profile(3,:),'-','LineWidth',3)
plot(p1.mesh*p1.period,p1.profile(4,:),'-','LineWidth',3)
plot(7.5,2.5,'k^','MarkerSize',12,'MarkerFaceColor','k')
legend('$x_1$', '$x_2$','$x_3$','$x_4$','$P_1$','Interpreter','latex','FontSize',20)
ylim([-3,3])
set(gca,'LineWidth',2,'Box','on','XTick',[],'YTick',[],'FontSize',14,'FontWeight','bold')
%
nexttile([1,3])
hold on
plot(par_axcdt(unst_indscdt==0),x0_axcdt(unst_indscdt==0),'-','Color','k','LineWidth',3)
plot(par_axcdt(unst_indscdt==6),x0_axcdt(unst_indscdt==6),'-.','Color',0*[0.7 0.7 0.7],'LineWidth',3)
plot(par_axcdt(unst_indscdt>=7),x0_axcdt(unst_indscdt>=7),'-.','Color',0*[0.7 0.7 0.7],'LineWidth',3)
plot(x_taus_po_dts1(unst_po_dts1>=0),ymxs_po_dts1(unst_po_dts1>=0),'-','Color',clrs(1,:),lwidth{:});
plot(x_taus_po_dts1(unst_po_dts1>=0),ymins_po_dts1(unst_po_dts1>=0),'-','Color',clrs(1,:),lwidth{:})
% plot(x_taus_po_dts1(unst_po_dts1>=1),ymxs_po_dts1(unst_po_dts1>=1),'x','Color',pgray,lwidth{:})
% plot(x_taus_po_dts1(unst_po_dts1>=1),ymins_po_dts1(unst_po_dts1>=1),'x','Color',pgray,lwidth{:})
            % plot hopf
plot(par_axcdt(indbifcdt),x0_axcdt(indbifcdt),'.','Color','r','MarkerSize',40,'LineWidth',3)
[~,its1]=min(abs(x_taus_po_dts1-3.5));
plot(x_taus_po_dts1(its1),ymxs_po_dts1(its1),'r^','MarkerSize',12,'MarkerFaceColor','r')
ylim([-3,3])
set(gca,'LineWidth',2,'Box','on','FontSize',14,'FontWeight','bold','XTick',[],'YTick',[-2,0,2])
%
% xlabel('\tau_c','FontSize',18)
nexttile([1,2])
hold on;
p2=psol_dts1.point(its1);
plot(p2.mesh*p2.period,p2.profile(1,:),'s','LineWidth',1)
plot(p2.mesh*p2.period,p2.profile(2,:),'o','LineWidth',1)
plot(p2.mesh*p2.period,p2.profile(3,:),'x','LineWidth',1)
plot(p2.mesh*p2.period,p2.profile(4,:),'.','LineWidth',1)
plot(7.5,2.5,'r^','MarkerSize',12,'MarkerFaceColor','r')
legend('$x_1$', '$x_2$','$x_3$','$x_4$','$P_2$','Interpreter','latex','FontSize',20)
ylim([-3,3])
set(gca,'LineWidth',2,'Box','on','XTick',[],'YTick',[],'FontSize',14,'FontWeight','bold')
%
nexttile([1,3])
hold on; 
biflocs2=find(diff(unst_po_dts2));
plot(par_axcdt(unst_indscdt==0),x0_axcdt(unst_indscdt==0),'-','Color','k','LineWidth',3)
plot(par_axcdt(unst_indscdt==6),x0_axcdt(unst_indscdt==6),'-.','Color',0*[0.7 0.7 0.7],'LineWidth',3)
plot(par_axcdt(unst_indscdt>=7),x0_axcdt(unst_indscdt>=7),'-.','Color',0*[0.7 0.7 0.7],'LineWidth',3)
plot(x_taus_po_dts2(biflocs2(2):biflocs2(3)),ymxs_po_dts2(biflocs2(2):biflocs2(3)),'-','Color',clrs(1,:),lwidth{:});
plot(x_taus_po_dts2(biflocs2(2):biflocs2(3)),ymins_po_dts2(biflocs2(2):biflocs2(3)),'-','Color',clrs(1,:),lwidth{:})
plot(x_taus_po_dts2(biflocs2(4):end),ymxs_po_dts2(biflocs2(4):end),'-','Color',clrs(1,:),lwidth{:});
plot(x_taus_po_dts2(biflocs2(4):end),ymins_po_dts2(biflocs2(4):end),'-','Color',clrs(1,:),lwidth{:})
plot(x_taus_po_dts2(1:biflocs2(2)),ymxs_po_dts2(1:biflocs2(2)),'-','Color',pgray,lwidth{:})
plot(x_taus_po_dts2(1:biflocs2(2)),ymins_po_dts2(1:biflocs2(2)),'-','Color',pgray,lwidth{:})
plot(x_taus_po_dts2(biflocs2(3):biflocs2(4)),ymxs_po_dts2(biflocs2(3):biflocs2(4)),'-','Color',pgray,lwidth{:})
plot(x_taus_po_dts2(biflocs2(3):biflocs2(4)),ymins_po_dts2(biflocs2(3):biflocs2(4)),'-','Color',pgray,lwidth{:})
% plot hopf
              plot(par_axcdt(indbifcdt),x0_axcdt(indbifcdt),'.','Color','r','MarkerSize',40,'LineWidth',3)
[~,its2]=min(abs(x_taus_po_dts2-3.5));
plot(x_taus_po_dts2(its2),ymxs_po_dts2(its2),'^','Color',[0 0.6 0],'MarkerSize',12,'MarkerFaceColor',[0 0.6 0])
ylim([-3,3])
set(gca,'LineWidth',2,'Box','on','FontSize',14,'FontWeight','bold','YTick',[-2,0,2])
nexttile([1,2])
hold on; 
p3=psol_dts2.point(its2);
plot(p3.mesh*p3.period,p3.profile(1,:),'-','LineWidth',3)
plot(p3.mesh*p3.period,p3.profile(2,:),'-','LineWidth',3)
plot(p3.mesh*p3.period,p3.profile(3,:),'-','LineWidth',3)
plot(p3.mesh*p3.period,p3.profile(4,:),'-','LineWidth',3)
set(gca,'LineWidth',2,'Box','on','FontSize',14,'FontWeight','bold','YTick',[])
plot(7.5,2.5,'^','Color',[0 0.6 0],'MarkerSize',12,'MarkerFaceColor',[0 0.6 0])
legend('$x_1$', '$x_2$','$x_3$','$x_4$','$P_3$','Interpreter','latex','FontSize',20)
ylim([-3,3])

%%
save('bif_analysis_S4_system_second.mat')
%%





% %%
% figure(220)
% clf;
% hold on; grid on 
% %%%%%%% plot equilibria %%%%%%%%%%
% %%% delat=0.2
%  text(3.5,1,'$\delta=0.2$','Interpreter','latex','FontSize',20,'FontWeight','bold')
% plot(par_axcdt(unst_indscdt==0),x0_axcdt(unst_indscdt==0),'-','Color','k','LineWidth',4)
% plot(par_axcdt(unst_indscdt==6),x0_axcdt(unst_indscdt==6),'--','Color',[0.7 0.7 0.7],'MarkerSize', 8,'LineWidth',3)
% plot(par_axcdt(unst_indscdt>=7),x0_axcdt(unst_indscdt>=7),'--','Color',[0.7 0.7 0.7],'MarkerSize', 8,'LineWidth',3)
% bifloc_dt=find(diff(unst_po_dt));
% plt_pod1stb=plot(x_taus_po_dt(bifloc_dt(2):bifloc_dt(3)),ymxs_po_dt(bifloc_dt(2):bifloc_dt(3)),'-','Color',clrs(7,:),lwidth{:},'markersize',mk);
% plt_pod1unstb=plot(x_taus_po_dt(1:bifloc_dt(2)),ymxs_po_dt(1:bifloc_dt(2)),'-','Color',[0.7 0.7 0.7],lwidth{:},'markersize',mk);
% plot(x_taus_po_dt(bifloc_dt(2):bifloc_dt(3)),ymins_po_dt(bifloc_dt(2):bifloc_dt(3)),'-','Color',clrs(7,:),lwidth{:},'markersize',mk)
% plot(x_taus_po_dt(1:bifloc_dt(2)),ymins_po_dt(1:bifloc_dt(2)),'-','Color',[0.7 0.7 0.7],lwidth{:},'markersize',mk)
% plot(x_taus_po_dt(bifloc_dt(3):end),ymxs_po_dt(bifloc_dt(3):end),'-','Color',[0.7 0.7 0.7],lwidth{:},'markersize',mk)
% plot(x_taus_po_dt(bifloc_dt(3):end),ymins_po_dt(bifloc_dt(3):end),'-','Color',[0.7 0.7 0.7],lwidth{:},'markersize',mk)
% xlabel('\tau_c','FontSize',18)
% hold on; grid on
% 
% plot(x_taus_po_dts1(unst_po_dts1==0),ymxs_po_dts1(unst_po_dts1==0),'o','Color',clrs(1,:),lwidth{:});
% plot(x_taus_po_dts1(unst_po_dts1==0),ymins_po_dts1(unst_po_dts1==0),'o','Color',clrs(1,:),lwidth{:})
% plot(x_taus_po_dts1(unst_po_dts1>=1),ymxs_po_dts1(unst_po_dts1>=1),'r-',lwidth{:})
% plot(x_taus_po_dts1(unst_po_dts1>=1),ymins_po_dts1(unst_po_dts1>=1),'r-',lwidth{:})

%%
% [~,it]=min(abs(x_taus_po_dtss-1)); % % tau_c=3
% pp=psol_dtss.point(it);
% % %
% % fd23=@(t,x,Xd)funcs_dt.sys_rhs([x,Xd],po.parameter);
% % his=@(t)dde_coll_eva(po.profile-0.1,po.mesh,1+t/po.period,po.degree);
% % 
% % sold2_23=dde23(fd23, [po.parameter(in.tau_s),po.parameter(in.tau_c)] ,his,[0,300],...
% %     ddeset('RelTol',1e-7,'AbsTol',1e-7,'Event',@(t,x,Xd)event_permutx(x)));
% % %%
% % save('dde23_solution.mat','sold_23')
% % %%
% figure(90)
% clf
% % tiledlayout(2,1)   
% % 
% % nexttile
% % hold on; grid on
% % plot(sold_23.x,sold_23.y(1:2,:),'-',...
% %     sold_23.x,sold_23.y(3:4,:),'-','LineWidth',3);
% % title('simulation near PO')
% % nexttile
% % hold on; grid on
% hold on
% plot(pp.mesh*pp.period,pp.profile(1,:),'x')
%     plot( pp.mesh*pp.period,pp.profile(2,:),'s')
%     plot( pp.mesh*pp.period,pp.profile(3,:),'o')
%    plot( pp.mesh*pp.period,pp.profile(4,:),'x','LineWidth',1)
% title('time profile for nera initial PO')


% 
% Mg_3=[0,1,0,0; 0,0,1,0; 1,0,0,0; 0,0,0,0]; %(12)
% S_3=[1,0,0,0; 0,1,0,0;0,0,1,0; 0,0,0,0]';
% %Mg_13=[0,0,1,0; 0,1,0,0; 1,0,0,0; 0,0,0,1];
% pmfixs3=dde_stst_lincond('pmfix',nx,'v','trafo',Mg_3,'rotation',[1,2]);
% psfixs3=dde_psol_lincond('psfix',nx,'profile','trafo',Mg_3,'shift',[1,2],'stateproj',S_3,'condprojint',linspace(0.0,0.2,2)'*[1,1]);
% [fpsol_dts3,psol_dts3,sucp_dts3]=SetupPsol(funcs_cdt,branch0_tauc_dt_bis,indbifcdt,'contpar',in.tau_c,'extracolumns','auto',...
%     'initcond',pmfixs3,'outputfuncs',true,'extra_condition',true,'usercond',psfixs3,'intervals',60,'degree',4,...
%      parbd_c{:},'print_residual_info',1,'matrix','sparse','remesh',false);
% %%
% figure(201)
% hold on
% [psol_dts3,scdt,fcdt,rcdt]=br_contn(fpsol_dts3,psol_dts3,5);
% %%
% Mg_13_24=[0,0,1,0; 0,0,0,1; 1,0,0,0; 0,1,0,0]; %(13)(24)
% S_13_24=[1,0,0,0; 0,1,0,0;0,0,1,0; 0,0,0,1]';
% pmfixs4=dde_stst_lincond('pmfix',nx,'v','trafo',Mg_13_24,'rotation',[1,2]);
% psfixs4=dde_psol_lincond('psfix',nx,'profile','trafo',Mg_13_24,'shift',[1,2],'stateproj',S_13_24,'condprojint',linspace(0.0,0.2,2)'*[1,1]);
% [fpsol_dts4,psol_dts4,sucp_dts4]=SetupPsol(funcs_cdt,branch0_tauc_dt_bis,indbifcdt,'contpar',in.tau_c,'extracolumns','auto',...
%     'initcond',pmfixs4,'outputfuncs',true,'extra_condition',true,'usercond',psfixs4,'intervals',60,'degree',4,...
%      parbd_c{:},'print_residual_info',1,'matrix','sparse','remesh',false);
% %% 
% figure(201)
% hold on
% [psol_dts4,scdt,fcdt,rcdt]=br_contn(fpsol_dts4,psol_dts4,5);