base=[pwd(),'\..\..\DDE_Biftool2025\'];
addpath([base,'ddebiftool'],...
    [base,'ddebiftool_extra_psol'],...
    [base,'ddebiftool_utilities'],...
    [base,'ddebiftool_extra_rotsym'],...
    [base,'ddebiftool_extra_nmfm'],...
    [base,'ddebiftool_extra_symbolic'],...
    [base,'ddebiftool_coco']);
load('bif_analysis_S4_system_second.mat')
%% We plot int_{0}^{1/2} u_A(t)-u_B(t+1/2) dt by setting:
Sint_A=dde_lincond_struct(size(psol_dt.point(1).profile,1),'profile','trafo',0,...
    'shift',[1,4],'condprojmat',-1,'stateproj',[1,0,0,0],'condprojint',[0,0.25]);
Sint_B=dde_lincond_struct(size(psol_dt.point(1).profile,1),'profile','trafo',0,...
    'shift',[1,4],'condprojmat',-1,'stateproj',[0,1,0,0],'condprojint',[0.25,0.75]);
yax_Sint_A=arrayfun(@(x)dde_psol_lincond(x,Sint_A),psol_dt.point);
yax_Sint_B=arrayfun(@(x)dde_psol_lincond(x,Sint_B),psol_dt.point);
yax_sym=yax_Sint_A-yax_Sint_B;

% Sint_A1=dde_lincond_struct(size(psol_dt.point(1).profile,1),'profile','trafo',0,...
%     'shift',[1,2],'condprojmat',-1,'stateproj',[1,0,0,0],'condprojint',[0,0.5]);
% Sint_B1=dde_lincond_struct(size(psol_dt.point(1).profile,1),'profile','trafo',0,...
%     'shift',[1,2],'condprojmat',-1,'stateproj',[0,1,0,0],'condprojint',[0.5,1]);
yax_Sint_A1=arrayfun(@(x)dde_psol_lincond(x,Sint_A),psol_dts1.point);
yax_Sint_B1=arrayfun(@(x)dde_psol_lincond(x,Sint_B),psol_dts1.point);
yax_sym1=yax_Sint_A1-yax_Sint_B1;

yax_Sint_A2=arrayfun(@(x)dde_psol_lincond(x,Sint_A),psol_dts2.point);
yax_Sint_B2=arrayfun(@(x)dde_psol_lincond(x,Sint_B),psol_dts2.point);
yax_sym2=yax_Sint_A2-yax_Sint_B2;

yax_Sint_A3=arrayfun(@(x)dde_psol_lincond(x,Sint_A),psol_dtss.point);
yax_Sint_B3=arrayfun(@(x)dde_psol_lincond(x,Sint_B),psol_dtss.point);
yax_sym3=yax_Sint_A3-yax_Sint_B3;
%%
figure(1)
clf
tiledlayout(8,6,'TileSpacing','compact')
nexttile([8,4])
hold on;
bifloc_dt=find(diff(unst_po_dt));
plt_pod1stb=plot(x_taus_po_dt(bifloc_dt(2):bifloc_dt(3)),yax_sym(bifloc_dt(2):bifloc_dt(3)),'-','Color',clrs(1,:),lwidth{:});
            plot(x_taus_po_dt(bifloc_dt(4):end),yax_sym(bifloc_dt(4):end),'-','Color',clrs(1,:),lwidth{:})
plt_pod1unstb=plot(x_taus_po_dt(1:bifloc_dt(2)),yax_sym(1:bifloc_dt(2)),'-','Color',pgray,lwidth{:});
             plot(x_taus_po_dt(bifloc_dt(3):bifloc_dt(4)),yax_sym(bifloc_dt(3):bifloc_dt(4)),'-','Color',pgray,lwidth{:});
sym_bif=plot(x_taus_po_dt(bifloc_dt(2:end)),yax_sym(bifloc_dt(2:end)),'k.','MarkerSize',55);
[~,it1]=min(abs(x_taus_po_dt-3.5));
p1legn=plot(x_taus_po_dt(it1),yax_sym(it1),'k^','MarkerSize',16,'MarkerFaceColor','k');

%%%%%%%
plt_pod2unstb=plot(x_taus_po_dts1(unst_po_dts1>=0),yax_sym1(unst_po_dts1>=0),'-','Color',clrs(2,:),lwidth{:});
[~,its1]=min(abs(x_taus_po_dts1-3.5));
p12egn=plot(x_taus_po_dts1(its1),yax_sym1(its1),'^','Color',[0 0.6 0],'MarkerSize',16,'MarkerFaceColor',[0 0.6 0]);
pltnan1=plot(NaN,NaN,'Color',[1,1,1]);
%%%%%%%
plot(x_taus_po_dts2(unst_po_dts2>=0),yax_sym2(unst_po_dts2>=0),'-','Color',clrs(3,:),lwidth{:})
[~,its2]=min(abs(x_taus_po_dts2-3.5));
p13egn=plot(x_taus_po_dts2(its2),yax_sym2(its2),'r^','MarkerSize',16,'MarkerFaceColor','r');
pltnan2=plot(NaN,NaN,'Color',[1,1,1]);

plt_pod3unst=plot(x_taus_po_dtss(unst_po_dtss>=1),yax_sym3(unst_po_dtss>=1),'-','Color',clrs(4,:),lwidth{:});
p14egn=plot(x_taus_po_dtss(itss),yax_sym3(itss),'^','Color',[1 0 1],'MarkerSize',16,'MarkerFaceColor',[1 0 1]);
pltnan3=plot(NaN,NaN,'Color',[1,1,1]);

plt_hpf=plot(x_taus_po_dt(1),yax_sym(1),'s','Marker', 's', 'MarkerSize', 18, ...
    'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black');
legend_text={'stable POs: $\Pi\sim (1234)_{1/4}$','unstable POs: $\Pi\sim (1234)_{1/4}$','stable POs: $\Pi\sim (13)_{0}(24)_{0}$  $\&$   $\Pi\sim (13)_{1/2}(24)_{1/2}$',...
    'unstable POs: $\Pi\sim (123)_{0}$', 'symmetry-breaking POs','equivariant Hopf-bifurcation'};
vec_plt=[plt_pod1stb(1),plt_pod1unstb(1),plt_pod2unstb(1),plt_pod3unst(1),sym_bif(1),plt_hpf(1)];%,p1legn(1),p12egn(1),p13egn(1),p14egn(1)];
legend(vec_plt,legend_text,'Interpreter','latex','FontSize',24,'Location','northwestoutside')
xlabel('$\tau_c$','FontSize',30,'Interpreter','latex')
set(gca,'LineWidth',2,'Box','on','FontWeight','bold','XTick',[],'YTick',[-2,0,2])
% Plot time profile PO~ (1234)_{1/4}
nexttile([2,2])
hold on;
p1=psol_dt.point(it1);
pl1=plot(NaN,NaN,'^','Color','k','MarkerSize',16,'MarkerFaceColor','k');
p1x=plot(p1.mesh*p1.period,p1.profile(1,:),'-','Color',clrs(1,:),'LineWidth',3);
p2x=plot(p1.mesh*p1.period,p1.profile(2,:),'-','Color',clrs(2,:),'LineWidth',3);
p3x=plot(p1.mesh*p1.period,p1.profile(3,:),'-','Color',clrs(3,:),'LineWidth',3);
p4x=plot(p1.mesh*p1.period,p1.profile(4,:),'-','Color',clrs(4,:),'LineWidth',3);
%plot(NaN,NaN,'Color',[1,1,1]);%plot(7.5,2.5,'k^','MarkerSize',12,'MarkerFaceColor','k')
xtext_leg1={'$P_1$','$x_1$', '$x_2$','$x_3$','$x_4$'};
xtext_vec1=[pl1(1),p1x,p2x,p3x,p4x];
legend(xtext_vec1,xtext_leg1','Interpreter','latex','FontSize',23,'location','northeastoutside')
ylim([-3,3])
xlim([0,p1.period])
set(gca,'LineWidth',2,'Box','on','XTick',[],'YTick',[],'FontSize',14,'FontWeight','bold')
%
nexttile([2,2])
hold on;
p2=psol_dts1.point(its1);
ypl0=plot(NaN,NaN,'^','Color',[0 0.6 0],'MarkerSize',16,'MarkerFaceColor',[0 0.6 0]);
ypl1=plot(p2.mesh*p2.period,p2.profile(1,:),'s','Color',clrs(1,:),'LineWidth',1);
ypl2=plot(p2.mesh*p2.period,p2.profile(2,:),'o','Color',clrs(2,:),'LineWidth',1);
ypl3=plot(p2.mesh*p2.period,p2.profile(3,:),'x','Color',clrs(3,:),'LineWidth',1);
ypl4=plot(p2.mesh*p2.period,p2.profile(4,:),'.','Color',clrs(4,:),'LineWidth',1);
xtext_leg2={'$P_2$','$x_1$', '$x_2$','$x_3$','$x_4$'};
xtext_vec2=[ypl0(1),ypl1,ypl2,ypl3,ypl4];
legend(xtext_vec2,xtext_leg2,'Interpreter','latex','FontSize',23,'location','northeastoutside')
ylim([-3,3])
xlim([0,p2.period])
set(gca,'LineWidth',2,'Box','on','XTick',[],'YTick',[],'FontSize',14,'FontWeight','bold')
%
nexttile([2,2])
hold on; 
p3=psol_dts2.point(its2);
plt30=plot(NaN,NaN,'r^','MarkerSize',16,'MarkerFaceColor','r');
plt31=plot(p3.mesh*p3.period,p3.profile(1,:),'x','Color',clrs(1,:),'LineWidth',1,'MarkerSize',10);
plt32=plot(p3.mesh*p3.period,p3.profile(2,:),'-','Color',clrs(2,:),'LineWidth',3);
plt33=plot(p3.mesh*p3.period,p3.profile(3,:),'-','Color',clrs(3,:),'LineWidth',3);
plt34=plot(p3.mesh*p3.period,p3.profile(4,:),'-','Color',clrs(4,:),'LineWidth',3);
xtext_leg3={'$P_3$','$x_1$', '$x_2$','$x_3$','$x_4$'};
xtext_vec3=[plt30(1),plt31,plt32,plt33,plt34];
legend(xtext_vec3,xtext_leg3,'Interpreter','latex','FontSize',23,'location','northeastoutside')
ylim([-3,3])
xlim([0,p3.period])
set(gca,'LineWidth',2,'Box','on','XTick',[],'YTick',[],'FontSize',14,'FontWeight','bold')
nexttile([2,2])
hold on; 
p4=psol_dtss.point(itss);
pl4=plot(NaN,NaN,'^','Color',[1,0,1],'MarkerSize',16,'MarkerFaceColor',[1 0 1]);
px2=plot(p4.mesh*p4.period,p4.profile(2,:),'o','Color',clrs(2,:),'LineWidth',1.5,'MarkerSize',10);
px1=plot(p4.mesh*p4.period,p4.profile(1,:),'-','Color',clrs(1,:),'LineWidth',5,'MarkerSize',12);
px3=plot(p4.mesh*p4.period,p4.profile(3,:),'-.','Color',clrs(3,:),'LineWidth',2)
px4=plot(p4.mesh*p4.period,p4.profile(4,:),'-','Color',clrs(4,:),'LineWidth',3)
set(gca,'LineWidth',2,'Box','on','FontSize',14,'FontWeight','bold','YTick',[])
xtext_leg={'$P_4$','$x_1$', '$x_2$','$x_3$','$x_4$'};
xtext_vec=[pl4(1),px1,px2,px3,px4];
legend(xtext_vec,xtext_leg,'Interpreter','latex','FontSize',23,'location','northeastoutside')
ylim([-3,3])
xlim([0,p4.period])
set(gca, 'FontSize',14,'FontWeight','bold','LineWidth',2,'Box', 'on')