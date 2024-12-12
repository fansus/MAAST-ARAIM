function S = compute_s_coefficients(G, W, GtW, invcov0, indb, cov0);
 
%*************************************************************************
%*     Copyright c 2012 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%************************************************************************* 

%Created by Juan Blanch 28/05/2012

    GtWindb = G(indb,:)'*W(indb,:);
    GtWindG = invcov0  - GtWindb*G;
    S       = Inf*ones(size(G,2),size(G,1));
        
    not_zero = find(sum(abs(GtWindG))>0);    
    GtWindGnew = GtWindG(not_zero,not_zero);
    n_unk = length(not_zero);
    GtWindbnew =  GtWindb(not_zero,:); 
    
    if (size(G,1)-length(indb))>= max(n_unk,4)  
       GtWindb = G(indb,:)'*W(indb,:);
       covd  = inv(GtWindGnew);
       Sred = covd*(GtW(not_zero,:) - GtWindbnew);
       S(not_zero,:) = Sred;
    else
       S = ones(size(G,2),size(G,1))*Inf;       
    end    
    
    
    