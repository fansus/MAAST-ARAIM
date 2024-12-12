function [sigma, sigma_ss, bias, bias_ss, s1vec, s2vec, s3vec, subsets,...
     pap_subsets, p_not_monitored, idx, x, chi2] = filter_out_subsets(sigma, sigma_ss, bias,...
     bias_ss, s1vec, s2vec, s3vec, subsets, pap_subsets, p_not_monitored, x, chi2)
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
 
 
 
 idx = find(min(sigma(:,1:3),[],2)<Inf);
 sigma = sigma(idx,:);
 sigma_ss = sigma_ss(idx,:);      
 bias = bias(idx,:);
 bias_ss = bias_ss(idx,:);
 s1vec = s1vec(idx,:);
 s2vec = s2vec(idx,:);
 s3vec = s3vec(idx,:);
 subsets =  subsets(idx,:);
 p_not_monitored = p_not_monitored +sum(pap_subsets)-sum(pap_subsets(idx,:));
 pap_subsets =pap_subsets(idx,:);
 x = x(idx,:);
 chi2 = chi2(idx,:);