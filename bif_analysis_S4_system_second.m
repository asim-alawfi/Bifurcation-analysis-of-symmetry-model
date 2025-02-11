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
% M_13_24=[0,0,1,0; 0,0,0,1; 1,0,0,0; 0,1,0,0]; %(13)(24)
%%
M_13_24=[0,0,1,0; 0,0,0,1; 1,0,0,0; 0,1,0,0]; %(13)(24)
M_13_24_shift=[0,1,0,0; 1,0,0,0; 0,0,0,1; 0,0,1,0]; %(12)(34)
pmfixs1=dde_stst_lincond('pmfix',nx,'v','trafo',M_13_24,'rotation',[0,1]);
psfixs1=dde_psol_lincond('psfix',nx,'profile','trafo',M_13_24,'shift',[0,1],'condprojint',linspace(0.0,0.2,2)'*[1,1]);
pmfixs_shift=dde_stst_lincond('pmfix',nx,'v','trafo',M_13_24_shift,'rotation',[1,2]);
psfixs_shift=dde_psol_lincond('psfix',nx,'profile','trafo',M_13_24_shift,'shift',[1,2],'condprojint',linspace(0.0,0.2,2)'*[1,1]);
[fpsol_dts1,psol_dts1,sucp_dts1]=SetupPsol(funcs_cdt,branch0_tauc_dt_bis,indbifcdt,'contpar',in.tau_c,'extracolumns','auto',...
    'initcond',[pmfixs1,pmfixs_shift],'outputfuncs',true,'extra_condition',true,'usercond',[psfixs1,psfixs_shift],'intervals',60,'degree',4,...
     parbd_c{:},'print_residual_info',1,'matrix','sparse','remesh',false);
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

%% Permutation 
M_flip=[0,1; 1,0];%0,0; 0,0,0,0; 0,0,0,0];%; 1,0,0,0; 0,0,0,0; 0,0,0,0]; %(12)
S12=[1,0,0,0; 0,1,0,0];
S34=[0,0,1,0; 0,0,0,1];
pmfixs12=dde_stst_lincond('pmfix',nx,'v','trafo',M_flip,'rotation',[0,1],'stateproj',S12);
psfixs12=dde_psol_lincond('psfix',nx,'profile','trafo',M_flip,'shift',[0,1],'condprojint',linspace(0.0,0.2,2)'*[1,1],'stateproj',S12);
pmfixs34=dde_stst_lincond('pmfix',nx,'v','trafo',M_flip,'rotation',[1,2],'stateproj',S34);
psfixs34=dde_psol_lincond('psfix',nx,'profile','trafo',M_flip,'shift',[1,2],'condprojint',linspace(0.0,0.2,2)'*[1,1],'stateproj',S34);
[fpsol_dts2,psol_dts2,sucp_dts2]=SetupPsol(funcs_cdt,branch0_tauc_dt_bis,indbifcdt,'contpar',in.tau_c,'extracolumns','auto',...
    'initcond',[pmfixs12,pmfixs34],'outputfuncs',true,'extra_condition',true,'usercond',[psfixs12,psfixs34],'intervals',60,'degree',4,...
     parbd_c{:},'print_residual_info',1,'matrix','sparse','remesh',false);
% [psol_dts2,scdt,fcdt,rcdt]=br_contn(fpsol_dts2,psol_dts2,9);
% [fpsol_dts2,psol_dts2,sucp_dts1]=SetupPsol(funcs_cdt,branch0_tauc_dt_bis,indbifcdt,'contpar',in.tau_c,'extracolumns','auto',...
%     'initcond',pmfixs2,'outputfuncs',true,'extra_condition',true,'usercond',psfixs2,'intervals',60,'degree',4,...
%      parbd_c{:},'print_residual_info',1,'matrix','sparse','remesh',false);
%
figure(200)
hold on
[psol_dts2,s2,fs2,rs2]=br_contn(fpsol_dts2,psol_dts2,150);
%%
pp=psol_dts2.point(end)
figure(888)
clf
hold on
plot(pp.mesh*pp.period,pp.profile(1,:),'x')
    plot( pp.mesh*pp.period,pp.profile(2,:),'s')
    plot( pp.mesh*pp.period,pp.profile(3,:),'o')
   plot( pp.mesh*pp.period,pp.profile(4,:),'x','LineWidth',1)
[psol_dts2,unst_po_dts2,dom_po_dts2,triv_defect_po_dts2]=br_stabl(fpsol_dts2,psol_dts2,0,1,'exclude_trivial',true);
x_taus_po_dts2=arrayfun(@(x)x.parameter(in.tau_c),psol_dts2.point);
% Plot stability of POs
ymxs_po_dts2=arrayfun(@(x)max(x.profile(1,:)),psol_dts2.point);
ymins_po_dts2=arrayfun(@(x)min(x.profile(1,:)),psol_dts2.point);
%%
M_123=[0,0,1,0; 0,0,1,0; 1,0,0,0; 0,0,0,1]; %(132)

pmfixss=dde_stst_lincond('pmfix',nx,'v','trafo',M_123,'rotation',[0,1]);
psfixss=dde_psol_lincond('psfix',nx,'profile','trafo',M_123,'shift',[0,1],'condprojint',linspace(0.0,0.2,2)'*[1,1]);
[fpsol_dtss,psol_dtss,sucp_dtss]=SetupPsol(funcs_cdt,branch0_tauc_dt_bis,indbifcdt,'contpar',in.tau_c,'extracolumns','auto',...
    'initcond',pmfixss,'outputfuncs',true,'extra_condition',true,'usercond',psfixss,'intervals',60,'degree',4,...
     parbd_c{:},'print_residual_info',1,'matrix','sparse','remesh',false);
%
figure(200)
clf
hold on
[psol_dtss,scdts,fcdts,rcdts]=br_contn(fpsol_dtss,psol_dtss,300);

[psol_dtss,unst_po_dtss,dom_po_dtss,triv_defect_po_dtss]=br_stabl(fpsol_dtss,psol_dtss,0,1,'exclude_trivial',true);%,'locate_trivial',@(p)[1,1],'geteigenfuncs',true);
x_taus_po_dtss=arrayfun(@(x)x.parameter(in.tau_c),psol_dtss.point);
% Plot stability of POs
ymxs_po_dtss=arrayfun(@(x)max(x.profile(1,:)),psol_dtss.point);
ymins_po_dtss=arrayfun(@(x)min(x.profile(1,:)),psol_dtss.point);
%%
save('bif_analysis_S4_system_second.mat')
%%
eigenvectors_values=branch0_tauc_dt_bis.point(indbifcdt).stability.v;

%%
fileID = fopen('eigenvectors_table.tex', 'w');

fprintf(fileID, '\\begin{table}[h]\n\\centering\n');
fprintf(fileID, '\\begin{tabular}{|c|c|c|c|}\n\\hline\n');
fprintf(fileID, 'v_1 & Eigenvector 2 & Eigenvector 3 & Eigenvector 4 \\\\\\hline\n');
M = rand(6,6) + 1i*rand(6,6); % Example complex eigenvector matrix (6 eigenvectors, each 1x6)
fileID = fopen('eigenvectors_table.tex', 'w');

fprintf(fileID, '\\begin{table}[h]\n\\centering\n');
fprintf(fileID, '\\begin{tabular}{|c|c|c|c|c|c|}\n\\hline\n');
fprintf(fileID, '$v_1$ & $v_2 $ & $v_3$ & $v_4$ & $v_5$ & $v_6$\\\\\\hline\n');

for i = 1:size(eigenvectors_values,1)
    
    fprintf(fileID, '$%.4f %+.4fi$ & $%.4f %+.4fi$ & $%.4f %+.4fi$ & $%.4f %+.4fi$ & $%.4f %+.4fi$ & $%.4f %+.4fi$ \\\\\\hline\n', ...
        real(eigenvectors_values(i,1)), imag(eigenvectors_values(i,1)), real(eigenvectors_values(i,2)), imag(eigenvectors_values(i,2)), ...
        real(eigenvectors_values(i,3)), imag(eigenvectors_values(i,3)), real(eigenvectors_values(i,4)), imag(eigenvectors_values(i,4)), ...
        real(eigenvectors_values(i,5)), imag(eigenvectors_values(i,5)), real(eigenvectors_values(i,6)), imag(eigenvectors_values(i,6)));
end

fprintf(fileID, '\\end{tabular}\n\\caption{...}\n\\end{table}\n');
fclose(fileID);

% 
% M_flip=[0,1; 1,0]%; 1,0,0,0; 0,0,0,0; 0,0,0,0]; %(12)
% S12=[1,0,0,0; 0,1,0,0];
% S34=[0,0,1,0; 0,0,0,1];
% pmfixs12=dde_stst_lincond('pmfix',nx,'v','trafo',M_flip,'rotation',[0,1],'stateproj',S12);
% psfixs12=dde_psol_lincond('psfix',nx,'profile','trafo',M_flip,'shift',[0,1],'stateproj',S12,'condprojint',linspace(0.0,0.2,2)'*[1,1]);
% pmfixs34=dde_stst_lincond('pmfix',nx,'v','trafo',M_flip,'rotation',[1,2],'stateproj',S34);
% psfixs34=dde_psol_lincond('psfix',nx,'profile','trafo',M_flip,'shift',[1,2],'stateproj',S34,'condprojint',linspace(0.0,0.2,2)'*[1,1]);
% [fpsol_dts2,psol_dts2,sucp_dts1]=SetupPsol(funcs_cdt,branch0_tauc_dt_bis,indbifcdt,'contpar',in.tau_c,'extracolumns','auto',...
%     'initcond',[pmfixs12,pmfixs34],'outputfuncs',true,'extra_condition',true,'usercond',[psfixs12,psfixs34],'intervals',60,'degree',4,...
%      parbd_c{:},'print_residual_info',1,'matrix','sparse','remesh',false);
% [psol_dts2,scdt,fcdt,rcdt]=br_contn(fpsol_dts2,psol_dts2,9);
% %%
%