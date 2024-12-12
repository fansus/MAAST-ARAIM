function [vpl, hpl, sig_acc, emt, subsets, pap_subset, p_not_monitored, rho] = mhss_raim_baseline_v4(G, sigpr2_int, sigpr2_acc, nom_bias_int, nom_bias_acc, p_sat, p_const,...
    opt_flag, rho_j, subsets, pap_subset, p_not_monitored)
%*************************************************************************
%*     Copyright c 2012 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%*************************************************************************

%Created 14 August 2012 by Juan Blanch
%Modified 17 March 2015 by Juan Blanch (refinement of exclusion list)
%Modified 24 August 2017 by Juan Blanch (allocation to exclusion, modified
%Q function, and subset selection algorithm) 
%MHSS_RAIM_BASELINE_V3 implements the vpl, hpl, sigma accuracy, and emt
%computation described in: Blanch, J., Lee, Y., Walter, T., Enge, P., Pervan, B., Belabbas, B., Spletter, A., Rippl, M.,
%"Advanced RAIM User Algorithm Description: Integrity Support Message Processing,
%Fault Detection, Exclusion, and Protection Level Calculation," Proceedings of the 25th 
%International Technical Meeting of The Satellite Division of the Institute of Navigation (ION GNSS 2012), Nashville, September 2012.

%G (Nsat X (3+Nconst)) geometry matrix in ENU. G(i,3+j)=1 if satellite i belongs to constellation j and zero otherwise
%sig2pr_int (Nsat X 1) nominal variance of the pseudorange error for integrity
%sig2pr_acc (Nsat X 1) nominal variance of the pseudorange error for accuracy
%nom_bias_int (Nsat X 1) nominal bias of the pseudorange error for integrity
%nom_bias_acc (Nsat X 1) nominal bias of the pseudorange error for accuracy
%p_sat (Nsat X 1) a priori probability of satellite fault
%p_const (Nsat X 1) a priori probability of constellation fault
%if p_sat or p_const are one scalar, it is assumed that it applies to all satellites or constellations
%rho_j is the fraction of the integrity budget given to exclusion mode j.
%For FD, rho_j =1

rho = [];

if nargin<9
    rho_j = 1;
end


global PHMI_VERT PHMI_HOR P_THRES PFA_VERT PFA_HOR P_EMT PL_TOL FC_THRES
global SIG_ACC_MAX_VERT SIG_ACC_MAX_HOR1 SIG_ACC_MAX_HOR2 



%%%%%%%%%%%% Determine subsets and associated probabilities %%%%%%%%%%%%%%%
if nargin<12
 [subsets, pap_subset, p_not_monitored]= determine_subsets_v4(G, p_sat, p_const, P_THRES, FC_THRES);
end
%%%%%%%%%%%%% Compute subset position solutions, sigmas, and biases %%%%%%%

[sigma, sigma_ss, bias, bias_ss, s1vec, s2vec, s3vec, x, chi2] = compute_subset_solutions(...
             G, sigpr2_int, sigpr2_acc, nom_bias_int, nom_bias_acc, subsets, zeros(size(G,1),1));
         
%%%%%%%%%%%%% Adjust all-in-view position coefficients %%%%%%%%%%%%%%%%%%%%
if opt_flag
nsets=size(subsets,1);


% s10 = s1vec(1,:) ;  %vec(1,:)+t1*(s1vec(idx_weak1,:)-s1vec(1,:));
% s20 = s2vec(1,:) ;  %vec(1,:)+t2*(s2vec(idx_weak2,:)-s2vec(1,:));
% s30 = s3vec(1,:) ;
% So(1,:)=s10;
% So(2,:)=s20;
% So(3,:)=s30;

[ s3opt ] = compute_adjusted_position_1D(pap_subset,sigma(:,3),sigma_ss(:,3), s3vec, sigpr2_acc,SIG_ACC_MAX_VERT);
[ s2opt ] = compute_adjusted_position_1D(pap_subset,sigma(:,2),sigma_ss(:,2), s2vec, sigpr2_acc,SIG_ACC_MAX_HOR2);
[ s1opt ] = compute_adjusted_position_1D(pap_subset,sigma(:,1),sigma_ss(:,1), s1vec, sigpr2_acc,SIG_ACC_MAX_HOR1);

s1vec(1,:) = s1opt;  %vec(1,:)+t1*(s1vec(idx_weak1,:)-s1vec(1,:));
s2vec(1,:) = s2opt;  %vec(1,:)+t2*(s2vec(idx_weak2,:)-s2vec(1,:));
s3vec(1,:) = s3opt;  %s3vec(1,:)+t3*(s3vec(idx_weak3,:)-s3vec(1,:));

sigma(1,1) = sqrt((s1vec(1,:).^2)* sigpr2_int);
sigma(1,2) = sqrt((s2vec(1,:).^2)* sigpr2_int);
sigma(1,3) = sqrt((s3vec(1,:).^2)* sigpr2_int);

bias(1,1)   = abs(s1vec(1,:))*nom_bias_int;
bias(1,2)   = abs(s2vec(1,:))*nom_bias_int;
bias(1,3)   = abs(s3vec(1,:))*nom_bias_int;

delta_s1vec = s1vec - ones(nsets,1)*s1vec(1,:);
delta_s2vec = s2vec - ones(nsets,1)*s2vec(1,:);
delta_s3vec = s3vec - ones(nsets,1)*s3vec(1,:);

sigma_ss(:,1) = sqrt((delta_s1vec.^2)* sigpr2_acc);
sigma_ss(:,2) = sqrt((delta_s2vec.^2)* sigpr2_acc);
sigma_ss(:,3) = sqrt((delta_s3vec.^2)* sigpr2_acc);

bias_ss(:,1)   = abs(delta_s1vec)*nom_bias_acc;
bias_ss(:,2)   = abs(delta_s2vec)*nom_bias_acc;
bias_ss(:,3)   = abs(delta_s3vec)*nom_bias_acc;
end        
%%%%%%%%%%%%% Filter out modes that cannot be monitored and adjust%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% integrity budget %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 [sigma, sigma_ss, bias, bias_ss, ~, ~, s3vec, ~,...
     pap_subset, p_not_monitored_2, ~] = filter_out_subsets(sigma, sigma_ss, bias,...
     bias_ss, s1vec, s2vec, s3vec, subsets, pap_subset, p_not_monitored, x, chi2);
 
 if nargin<12
     added_p_not_mon = 0;
     p_not_monitored = p_not_monitored_2;
 else
     added_p_not_mon = p_not_monitored_2 - p_not_monitored;
 end
 
 
 
 %Nsubsets = size(subsets,1);

%%%%%%%%%%%%%%%%%%%%%%% Compute test thresholds %%%%%%%%%%%%%%%%%%%%%%%%%%%

[T1, T2, T3] = compute_test_thresholds(sigma_ss, bias_ss, PFA_VERT, PFA_HOR);

%%%%%%%%%%%%%%%%%%%%% Compute Vertical Protection Level %%%%%%%%%%%%%%%%%%%

p_fault = pap_subset;
p_fault(1) = 2;
phmi_vert = rho_j*(PHMI_VERT - PHMI_VERT/(PHMI_VERT+PHMI_HOR)*p_not_monitored) -PHMI_VERT/(PHMI_VERT+PHMI_HOR)*added_p_not_mon ;
[vpl, alloc3]     = compute_protection_level_v4(sigma(:,3), bias(:,3) + T3, p_fault, phmi_vert, PL_TOL);

%%%%%%%%%%%%%%%%%%%%% Compute Horizontal Protection Level %%%%%%%%%%%%%%%%%
phmi_hor = rho_j*(PHMI_HOR/2 -.5*PHMI_HOR/(PHMI_VERT+PHMI_HOR)*p_not_monitored)-.5*PHMI_HOR/(PHMI_VERT+PHMI_HOR)*added_p_not_mon ;
[hpl1, alloc1] = compute_protection_level_v4(sigma(:,1), bias(:,1) + T1, p_fault,  phmi_hor, PL_TOL);
[hpl2, alloc2] = compute_protection_level_v4(sigma(:,2), bias(:,2) + T2, p_fault,  phmi_hor, PL_TOL);
hpl = sqrt(hpl1^2 + hpl2^2);

%Exclude modes that are double counted

idx = find((alloc3+alloc2+alloc1)>=1);

if ~isempty(idx)

p_not_monitored = p_not_monitored + sum(p_fault(idx));
idx = setdiff(1:length(p_fault),idx);
% s1vec = s1vec(idx,:);
% s2vec = s2vec(idx,:);
s3vec = s3vec(idx,:);
sigma = sigma(idx,:);
sigma_ss = sigma_ss(idx,:);      
bias = bias(idx,:);
bias_ss = bias_ss(idx,:);
subsets = subsets(idx,:);
%p_fault =p_fault(idx);
pap_subset= pap_subset(idx);
%Nsubsets =length(idx);

%%%%%%%%%%%%%%%%%%%%%%% Compute test thresholds %%%%%%%%%%%%%%%%%%%%%%%%%%%

[T1, T2, T3] = compute_test_thresholds(sigma_ss, bias_ss, PFA_VERT, PFA_HOR);

%%%%%%%%%%%%%%%%%%%%% Compute Vertical Protection Level %%%%%%%%%%%%%%%%%%%

p_fault = pap_subset;
p_fault(1) = 2;
phmi_vert = rho_j*(PHMI_VERT - PHMI_VERT/(PHMI_VERT+PHMI_HOR)*p_not_monitored) -PHMI_VERT/(PHMI_VERT+PHMI_HOR)*added_p_not_mon ;
[vpl, alloc3]     = compute_protection_level_v4(sigma(:,3), bias(:,3) + T3, p_fault,  phmi_vert, PL_TOL);

%%%%%%%%%%%%%%%%%%%%% Compute Horizontal Protection Level %%%%%%%%%%%%%%%%%
phmi_hor = rho_j*(PHMI_HOR/2 -.5*PHMI_HOR/(PHMI_VERT+PHMI_HOR)*p_not_monitored)-.5*PHMI_HOR/(PHMI_VERT+PHMI_HOR)*added_p_not_mon ;
[hpl1, alloc1] = compute_protection_level_v4(sigma(:,1), bias(:,1) + T1, p_fault,  phmi_hor, PL_TOL,1-alloc3);
[hpl2, ~] = compute_protection_level_v4(sigma(:,2), bias(:,2) + T2, p_fault,  phmi_hor, PL_TOL, 1-alloc1-alloc3);
hpl = sqrt(hpl1^2 + hpl2^2);


end


%%%%%%%%%%%%%%%%%%% Compute Effective Monitor threshold %%%%%%%%%%%%%%%%%%%
% idx = find(pap_subset>=P_EMT);
sigma3_acc = sqrt((s3vec.*s3vec)*sigpr2_acc);
% K_md_emt   = - norminv(.5*P_EMT./pap_subset);
%emt = max(T3(idx) + K_md_emt(idx).*sigma3_acc(idx));
emt = compute_emt(s3vec, sigpr2_acc, T3, pap_subset, P_EMT);

%%%%%%%%%%%%%%%%%%%% Compute sigma accuracy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sig_acc = sigma3_acc(1);
%End

