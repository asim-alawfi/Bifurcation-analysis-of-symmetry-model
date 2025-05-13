clear;
%ddebiftool_path([pwd(),'/ddebiftool-snapshot-2025-02-07']);
% base=[pwd(),'\..\..\DDE_Biftool2025\'];
% addpath([base,'ddebiftool'],...
%     [base,'ddebiftool_extra_psol'],...
%     [base,'ddebiftool_utilities'],...
%     [base,'ddebiftool_extra_rotsym'],...
%     [base,'ddebiftool_extra_nmfm'],...
%     [base,'ddebiftool_extra_symbolic'],...
%     [base,'ddebiftool_coco']);
% try 
ddebiftool_path([pwd(),'/ddebiftool-snapshot-2025-02-07']);
format compact
format short g
%%
load("S4_demo_psol_results.mat")
%% Functions to extract parameter and mean of max values on a given branch 
parameter_array=@(branch,parind)arrayfun(@(x)x.parameter(parind),branch.point);
max_mean=@(branch)arrayfun(@(x)max(mean(x.profile,1)),branch.point);
%%
% 1) psol.p1234_1_4
% 2) psol.p123_0_1
% 3) psol.p12_0_1_34_1_2
% 4) psol.p12_0_1_34_0_1
x_p1234_1_4=parameter_array(psolref.p1234_1_4,ip.tau_c);
y_p1234_1_4=max_mean(psolref.p1234_1_4);

x_p123_0_1=parameter_array(psolref.p123_0_1,ip.tau_c);
y_p123_0_1=max_mean(psolref.p123_0_1);

x_p12_0_1_34_1_2=parameter_array(psolref.p12_0_1_34_1_2,ip.tau_c);
y_p12_0_1_34_1_2=max_mean(psolref.p12_0_1_34_1_2);


x_p12_0_1_34_0_1=parameter_array(psolref.p12_0_1_34_0_1,ip.tau_c);
y_p12_0_1_34_0_1=max_mean(psolref.p12_0_1_34_0_1);
%
ppsol_tests.p1234_1_4=psol_tests.p1234_1_4;
ppsol_tests.p123_0_1=psol_tests.p123_0_1;
ppsol_tests.p12_0_1_34_0_1=psol_tests.p12_0_1_34_0_1;
ppsol_tests.p12_0_1_34_1_2=psol_tests.p12_0_1_34_1_2;

%% Here change NaN value (branching point) 
ppsol_tests.p1234_1_4(1)=1;
ppsol_tests.p123_0_1(1)=1;
ppsol_tests.p12_0_1_34_0_1(1)=0;
ppsol_tests.p12_0_1_34_1_2(1)=3;
%
clrs=lines();
lwidth={'LineWidth',2};
figure(5)
clf
hold on; grid on
plot(x_p1234_1_4,y_p1234_1_4,'x',lwidth{:})
plot( x_p12_0_1_34_1_2,y_p12_0_1_34_1_2,'-',lwidth{:})
plot(x_p123_0_1,y_p123_0_1,'-',lwidth{:})
plot(x_p1234_1_4,y_p1234_1_4,'-',lwidth{:})
%%
% integral for x1 from 0 to 0.25
yfun1=dde_lincond_struct(size(psol.p123_0_1.point(1).profile,1),'profile','trafo',0,...
    'shift',[1,4],'condprojmat',-1,'stateproj',[1,0,0,0],'condprojint',[0,.25]);
% integral for x2 from 0.25 to 0.75
yfun2=dde_lincond_struct(size(psol.p123_0_1.point(1).profile,1),'profile','trafo',0,...
    'shift',[1,4],'condprojmat',-1,'stateproj',[0,1,0,0],'condprojint',[0.25,0.75]);
%% 
% branch psol.p1234_1_4
y1_p1234_1_4=arrayfun(@(x)dde_psol_lincond(x,yfun1),psolref.p1234_1_4.point);
y2_p1234_1_4=arrayfun(@(x)dde_psol_lincond(x,yfun2),psolref.p1234_1_4.point);
y_p1234_1_4=y1_p1234_1_4-y2_p1234_1_4;
%%%%%%%%%%

% branch psol.p123_0_1
y1_p123_0_1=arrayfun(@(x)dde_psol_lincond(x,yfun1),psolref.p123_0_1.point);
y2_p123_0_1=arrayfun(@(x)dde_psol_lincond(x,yfun2),psolref.p123_0_1.point);
y_p123_0_1=y1_p123_0_1-y2_p123_0_1;
%%%%%%%%%%%%

% psol.p12_0_1_34_1_2
y1_p12_0_1_34_1_2=arrayfun(@(x)dde_psol_lincond(x,yfun1),psolref.p12_0_1_34_1_2.point);
y2_p12_0_1_34_1_2=arrayfun(@(x)dde_psol_lincond(x,yfun2),psolref.p12_0_1_34_1_2.point);
y_p12_0_1_34_1_2=y1_p12_0_1_34_1_2-y2_p12_0_1_34_1_2;

%%%%%%%
% branch psol.p12_0_1_34_0_1
y1_p12_0_1_34_0_1=arrayfun(@(x)dde_psol_lincond(x,yfun1),psolref.p12_0_1_34_0_1.point);
y2_p12_0_1_34_0_1=arrayfun(@(x)dde_psol_lincond(x,yfun2),psolref.p12_0_1_34_0_1.point);
y_p12_0_1_34_0_1=y1_p12_0_1_34_0_1-y2_p12_0_1_34_0_1;
%%
lwidth={'LineWidth',4};
figure(8)
clf
tiledlayout(8,6,'TileSpacing','compact')
nexttile([8,4])
hold on;
plt_1234_14=plot(x_p1234_1_4(ppsol_tests.p1234_1_4==0),y_p1234_1_4(ppsol_tests.p1234_1_4==0),'color',clrs(1,:),lwidth{:})
plt_1234_14uns=plot(x_p1234_1_4(ppsol_tests.p1234_1_4>=1),y_p1234_1_4(ppsol_tests.p1234_1_4>=1),'color',0.1*[1 1 1],lwidth{:})
bifloc_1234_1_4=find(diff(ppsol_tests.p1234_1_4));
pltbif=plot(x_p1234_1_4(bifloc_1234_1_4),y_p1234_1_4(bifloc_1234_1_4),'k.','Markersize',40)
[~,it1]=min(abs(x_p1234_1_4-3.5));
p1t1=plot(x_p1234_1_4(it1),y_p1234_1_4(it1),'k^','MarkerSize',12,'MarkerFaceColor','k');
%%%%%%%%%
plt_123_01=plot(x_p123_0_1,y_p123_0_1,'color',clrs(2,:),lwidth{:});
[~,it2]=min(abs(x_p123_0_1-3.5));
pplt2=plot(x_p123_0_1(it2),y_p123_0_1(it2),'^','Color',[0 0.6 0],'MarkerSize',12,'MarkerFaceColor',[0 0.6 0]);
%%%%%%%%%%%%%%%%
[~,it3]=min(abs(x_p12_0_1_34_1_2-3.5));
plt_p12_0_1_34_1_2=plot(x_p12_0_1_34_1_2,y_p12_0_1_34_1_2,'color',clrs(3,:),lwidth{:});
plt3=plot(x_p12_0_1_34_1_2(it3),y_p12_0_1_34_1_2(it3),'r^','MarkerSize',12,'MarkerFaceColor','r');


%
[~,it4]=min(abs(x_p12_0_1_34_0_1-3.5));
plt_p12_0_1_34_0_1=plot(x_p12_0_1_34_0_1,y_p12_0_1_34_0_1,'color',clrs(4,:),lwidth{:});
plt4=plot(x_p12_0_1_34_0_1(it4),y_p12_0_1_34_0_1(it4),'m^','MarkerSize',12,'MarkerFaceColor','m');
%%%%%%%%%%%%%
legend_text={'stable POs: $\Pi\sim (1234)_{1/4}$','unstable POs: $\Pi\sim (1234)_{1/4}$','unstable POs: $\Pi\sim (123)_{0}$',...
    'unstable POs: $\Pi\sim (12)_{0}(34)_{1/2}$','unstable POs: $\Pi\sim (12)_{0}(34)_{0}$'}%,...
            % 'unstable POs: $\Pi\sim (12)_{0}(34)_{1/2}$','unstable POs: $\Pi\sim (123)_{0}$', 'symmetry-breaking POs','equivariant Hopf-bifurcation'};

vec_plt=[plt_1234_14(1),plt_1234_14uns(1),plt_123_01(1),plt_p12_0_1_34_1_2(1),plt_p12_0_1_34_0_1]%,plt_pod4unstb(1),plt_pod3unst(1),sym_bif(1),plt_hpf(1)];%,p1legn(1),p12egn(1),p13egn(1),p14egn(1)];
legend(vec_plt,legend_text,'Interpreter','latex','FontSize',24,'Location','northwestoutside')
set(gca,'LineWidth',2,'Box','on','FontWeight','bold','YTick',[-2,0,2],'FontSize',14)
xlabel('$\tau_c$','Interpreter','latex','FontSize',30)
nexttile([2,2])
hold on;
p1=psolref.p1234_1_4.point(it1);
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
p2=psolref.p123_0_1.point(it2);
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
p3=psolref.p12_0_1_34_1_2.point(it3);
pl4=plot(NaN,NaN,'r^','MarkerSize',12,'MarkerFaceColor','r');
px2=plot(p3.mesh*p3.period,p3.profile(2,:),'o','Color',clrs(2,:),'LineWidth',1,'MarkerSize',5);
px1=plot(p3.mesh*p3.period,p3.profile(1,:),'.','Color',clrs(1,:),'LineWidth',5,'MarkerSize',12);
px3=plot(p3.mesh*p3.period,p3.profile(3,:),'x','Color',clrs(3,:),'LineWidth',1);
px4=plot(p3.mesh*p3.period,p3.profile(4,:),'s','Color',clrs(4,:),'LineWidth',1);
set(gca,'LineWidth',2,'Box','on','FontSize',14,'FontWeight','bold','YTick',[])
xtext_leg={'$P_3$','$x_1$', '$x_2$','$x_3$','$x_4$'};
xtext_vec=[pl4(1),px1,px2,px3,px4];
legend(xtext_vec,xtext_leg,'Interpreter','latex','FontSize',23,'location','northeastoutside')
ylim([-3,3])
xlim([0,p3.period])
set(gca, 'FontSize',14,'FontWeight','bold','LineWidth',2,'Box', 'on')

nexttile([2,2])
hold on; 
p4=psolref.p12_0_1_34_0_1.point(it4);
pl4=plot(NaN,NaN,'m^','MarkerSize',12,'MarkerFaceColor','m');
px2=plot(p4.mesh*p4.period,p4.profile(2,:),'o','Color',clrs(2,:),'LineWidth',1.5,'MarkerSize',10);
px1=plot(p4.mesh*p4.period,p4.profile(1,:),'-','Color',clrs(1,:),'LineWidth',5,'MarkerSize',12);
px3=plot(p4.mesh*p4.period,p4.profile(3,:),'x','Color',clrs(3,:),'LineWidth',1);
px4=plot(p4.mesh*p4.period,p4.profile(4,:),'s','Color',clrs(4,:),'LineWidth',1);
set(gca,'LineWidth',2,'Box','on','FontSize',14,'FontWeight','bold','YTick',[])
xtext_leg={'$P_4$','$x_1$', '$x_2$','$x_3$','$x_4$'};
xtext_vec=[pl4(1),px1,px2,px3,px4];
legend(xtext_vec,xtext_leg,'Interpreter','latex','FontSize',23,'location','northeastoutside')
ylim([-3,3])
xlim([0,p4.period])
set(gca, 'FontSize',14,'FontWeight','bold','LineWidth',2,'Box', 'on')
%%
