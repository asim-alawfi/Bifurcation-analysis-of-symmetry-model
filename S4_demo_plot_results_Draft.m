clear;
%ddebiftool_path([pwd(),'/ddebiftool-snapshot-2025-02-07']);
base=[pwd(),'\..\..\DDE_Biftool2025\'];
addpath([base,'ddebiftool'],...
    [base,'ddebiftool_extra_psol'],...
    [base,'ddebiftool_utilities'],...
    [base,'ddebiftool_extra_rotsym'],...
    [base,'ddebiftool_extra_nmfm'],...
    [base,'ddebiftool_extra_symbolic'],...
    [base,'ddebiftool_coco']);
format compact
format short g
%%
load("S4_demo_psol_results.mat")
%% Functions to extract parameter and mean of max values on a given branch 
parameter_array=@(branch,parind)arrayfun(@(x)x.parameter(parind),branch.point);
max_mean=@(branch)arrayfun(@(x)max(mean(x.profile,1)),branch.point);
% 1) psol.p1234_1_4
% 2) psol.p123_0_1
% 3) psol.p12_0_1_34_1_2
% 4) psol.p12_0_1_34_0_1
x_p1234_1_4=parameter_array(psol.p1234_1_4,ip.tau_c);
y_p1234_1_4=max_mean(psol.p1234_1_4);

x_p123_0_1=parameter_array(psol.p123_0_1,ip.tau_c);
y_p123_0_1=max_mean(psol.p123_0_1);

x_p12_0_1_34_1_2=parameter_array(psol.p12_0_1_34_1_2,ip.tau_c);
y_p12_0_1_34_1_2=max_mean(psol.p12_0_1_34_1_2);


x_p12_0_1_34_0_1=parameter_array(psol.p12_0_1_34_0_1,ip.tau_c);
y_p12_0_1_34_0_1=max_mean(psol.p12_0_1_34_0_1);
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
%% branch psol.p1234_1_4
y1_p1234_1_4=arrayfun(@(x)dde_psol_lincond(x,yfun1),psol.p1234_1_4.point);
y2_p1234_1_4=arrayfun(@(x)dde_psol_lincond(x,yfun2),psol.p1234_1_4.point);
y_p1234_1_4=y1_p1234_1_4-y2_p1234_1_4;
% branch psol.p123_0_1
y1_p123_0_1=arrayfun(@(x)dde_psol_lincond(x,yfun1),psol.p123_0_1.point);
y2_p123_0_1=arrayfun(@(x)dde_psol_lincond(x,yfun2),psol.p123_0_1.point);
y_p123_0_1=y1_p123_0_1-y2_p123_0_1;
% psol.p12_0_1_34_1_2
y1_p12_0_1_34_1_2=arrayfun(@(x)dde_psol_lincond(x,yfun1),psol.p12_0_1_34_1_2.point);
y2_p12_0_1_34_1_2=arrayfun(@(x)dde_psol_lincond(x,yfun2),psol.p12_0_1_34_1_2.point);
y_p12_0_1_34_1_2=y1_p12_0_1_34_1_2-y2_p12_0_1_34_1_2;
% branch psol.p12_0_1_34_0_1
y1_p12_0_1_34_0_1=arrayfun(@(x)dde_psol_lincond(x,yfun1),psol.p12_0_1_34_0_1.point);
y2_p12_0_1_34_0_1=arrayfun(@(x)dde_psol_lincond(x,yfun2),psol.p12_0_1_34_0_1.point);
y_p12_0_1_34_0_1=y1_p12_0_1_34_0_1-y2_p12_0_1_34_0_1;
%
figure(6)
clf
hold on; grid on
plot(x_p1234_1_4,y_p1234_1_4,'color',clrs(1,:),lwidth{:})
plot(x_p123_0_1,y_p123_0_1,'color',clrs(2,:),lwidth{:});
plot(x_p12_0_1_34_1_2,y_p12_0_1_34_1_2,'color',clrs(3,:),lwidth{:});
plot(x_p12_0_1_34_0_1,y_p12_0_1_34_0_1,'color',clrs(4,:),lwidth{:})