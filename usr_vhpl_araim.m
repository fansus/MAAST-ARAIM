function vhpl=usr_vhpl_araim(los_xyzb, usr_idx, sig2_i, prn, sig2acc_i, bnom_i, bcont_i, psat_i)

%*************************************************************************
%*     Copyright c 2015 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%*************************************************************************
%
%USR_VHPL_ARAIM calculates user vertical and horizontal protection levels (vpl & hpl)
%standard deviation of the 
%VHPL=USR_VHPL(LOS_XYZB, USR_IDX, SIG2_I, PRN, SIG2ACC_I, BNOM_I, BCONT_I, PSAT_I)
%   Given n_los of user lines of sight vectors in ECEF WGS-84 coordinates 
%   (X in first column, Y in second column ...) in LOS_XYZB(nlos,4), the
%   corresponding user identification numbers in USR_IDX (n_los,1), the
%   integrity variances SIG2_I (n_los,1), accuracy variances SIG2ACC_I(n_los,1), 
%   nominal biases BNOM_I(n_los,1), nominal biases for continuity BCONT_I(n_los,1),
%   prior probability of satellite fault PSAT_I(n_los,1)  for each los, this function will determine
%   the vertical protection limit (VPL), horizontal protection limit (HPL),
%   the standard deviation of the accuracy (SIG_ACC) and the Effective Monitor
%   Threshold in two formulations (EMT and EMT_NEW)using the function
%   MHSS_RAIM_BASELINE(G, SIGPR2_INT, SIGPR2_ACC, NOM_BIAS_INT,
%   NOM_BIAS_ACC, P_SAT, P_CONST).
%
%   See also: MHSS_RAIM_BASELINE_V3, MHSS_RAIM_BASELINE_ADJUSTED_POSITION_V3

%2013August13 Created by Juan Blanch
%2015March 18 Modified by Juan Blanch (includes 1 dimensional position
%optimization)

global MOPS_NOT_MONITORED 

global MOPS_MIN_GPSPRN MOPS_MAX_GPSPRN MOPS_MIN_GLOPRN MOPS_MAX_GLOPRN ...
       MOPS_MIN_GALPRN MOPS_MAX_GALPRN  ...
       MOPS_MIN_BDUPRN MOPS_MAX_BDUPRN
global ARAIM_PCONST_GPS ARAIM_PCONST_GAL ARAIM_PCONST_GLO ARAIM_PCONST_BDU  
global SIG_ACC_MAX_VERT SIG_ACC_MAX_HOR1 SIG_ACC_MAX_HOR2 VPLT HPLT EMTT FDE_FLAG
global ATTEMPT_OPT PL0_FDE FDE_WF_FLAG
%initialize return value
n_usr=max(usr_idx);
[n_los, ~]=size(los_xyzb);
vhpl=repmat(MOPS_NOT_MONITORED,n_usr,4);


for usr=1:n_usr
    
  if mod(usr,500)==0
      disp(usr);
  end
  if usr==3233
      disp(usr);
  end
  svidxgps = find((usr_idx==usr)&(prn >= MOPS_MIN_GPSPRN & prn <= MOPS_MAX_GPSPRN)&(sig2_i<Inf)&(psat_i<1));
  svidxgal = find((usr_idx==usr)&(prn >= MOPS_MIN_GALPRN & prn <= MOPS_MAX_GALPRN)&(sig2_i<Inf)&(psat_i<1));
  svidxglo = find((usr_idx==usr)&(prn >= MOPS_MIN_GLOPRN & prn <= MOPS_MAX_GLOPRN)&(sig2_i<Inf)&(psat_i<1));
  svidxbdu = find((usr_idx==usr)&(prn >= MOPS_MIN_BDUPRN & prn <= MOPS_MAX_BDUPRN)&(sig2_i<Inf)&(psat_i<1));
  
  ngps = length(svidxgps);
  ngal = length(svidxgal);
  nglo = length(svidxglo);
  nbdu = length(svidxbdu);
  
  svconst = [(~isempty(svidxgps)) (~isempty(svidxgal)) (~isempty(svidxglo)) (~isempty(svidxbdu))];
  
  %nconst = sum([(~isempty(svidxgps)) (~isempty(svidxgal)) (~isempty(svidxglo)) (~isempty(svidxbdu))]);
  nsat   = ngps + ngal + nglo + nbdu;
  svidx  = [svidxgps; svidxgal; svidxglo; svidxbdu];
     
  clk_gps = [ones(ngps,svconst(1)); zeros(nsat - ngps,svconst(1))];
  clk_gal = [zeros(ngps,svconst(2)); ones(ngal,svconst(2));zeros(nsat - ngps - ngal,svconst(2))];    
  clk_glo = [zeros(ngps+ngal,svconst(3)); ones(nglo,svconst(3));zeros(nbdu,svconst(3))];
  clk_bdu = [zeros(nsat-nbdu,svconst(4)); ones(nbdu,svconst(4))];
   
  G = [ los_xyzb(svidx,1:3) clk_gps  clk_gal clk_glo clk_bdu];

  p_const = [ones(svconst(1),1)*ARAIM_PCONST_GPS; ones(svconst(2),1)*ARAIM_PCONST_GAL;...
             ones(svconst(3),1)*ARAIM_PCONST_GLO; ones(svconst(4),1)*ARAIM_PCONST_BDU];
       
  n_view=length(svidx);
  
  if(n_view>3)  
    
    sigpr2_int = sig2_i(svidx);
    sigpr2_acc = sig2acc_i(svidx);
    nom_bias_int = bnom_i(svidx);
    nom_bias_acc = bcont_i(svidx);
    p_sat         = psat_i(svidx);
    
  if FDE_FLAG
      
     subsets_exc = ones(n_view,n_view) - eye(n_view);
     if FDE_WF_FLAG
        idwf = find(p_const>1e-7);
        subsets_wf = ones(length(idwf),n_view) - G(:,idwf+3)';
        subsets_exc = [subsets_exc;subsets_wf];    
     end
     
     Nsubsets = size(subsets_exc,1);
     Nexc = size(subsets_exc,1);
     hpl_exc = Inf*ones(Nsubsets,1);
     vpl_exc = Inf*ones(Nsubsets,1);
     sig_acc_exc = Inf*ones(Nsubsets,1);
     emt_exc     = Inf*ones(Nsubsets,1);
     rho_j = 1/(Nexc+1);
     
     [vpl, hpl, sig_acc, emt, subsets, pap_subset, p_not_monitored] = mhss_raim_baseline_v4(G, sigpr2_int, sigpr2_acc, nom_bias_int, nom_bias_acc, p_sat, p_const,0, rho_j);
          
     if (sig_acc<SIG_ACC_MAX_VERT)&&((vpl>VPLT)||(hpl>HPLT)||emt>EMTT)&&(ATTEMPT_OPT>0)
     %[vpl_adj, hpl_adj, sig_acc_adj, emt_adj] = mhss_raim_baseline_adjusted_position_v4(G, sigpr2_int, sigpr2_acc, nom_bias_int, nom_bias_acc, p_sat, p_const,rho_j);
     [vpl_adj, hpl_adj, sig_acc_adj, emt_adj,subsets, pap_subset, p_not_monitored] = mhss_raim_baseline_v4(G, sigpr2_int, sigpr2_acc,...
         nom_bias_int, nom_bias_acc, p_sat, p_const,ATTEMPT_OPT, rho_j);
         if (vpl_adj<vpl)||(hpl_adj<hpl)
            vpl = vpl_adj;
            emt = emt_adj;
            sig_acc = sig_acc_adj;
            hpl = hpl_adj;
         end 
     end 
      
     for i=1:Nsubsets    
         idx = logical(subsets_exc(i,:));
         subsets_i = subsets(:,idx);
         if sum(idx)>3
         [vpl_bl, hpl_bl, sig_acc_bl, emt_bl] = mhss_raim_baseline_v4(G(idx,:), sigpr2_int(idx), sigpr2_acc(idx), ...
             nom_bias_int(idx), nom_bias_acc(idx), p_sat(idx), p_const, 0, rho_j, subsets_i, pap_subset, p_not_monitored);
%          [vpl_bl, hpl_bl, sig_acc_bl, emt_bl] = mhss_raim_baseline_v4(G(idx,:), sigpr2_int(idx), sigpr2_acc(idx), ...
%              nom_bias_int(idx), nom_bias_acc(idx), p_sat(idx), p_const, 0, rho_j);
          vpl_exc(i) = vpl_bl;
          emt_exc(i) = emt_bl;
          sig_acc_exc(i) = sig_acc_bl;
          hpl_exc(i) = hpl_bl;
    
            if (sig_acc_exc(i)<SIG_ACC_MAX_VERT)&&((vpl_exc(i)>VPLT)||(hpl_exc(i)>HPLT)||emt_exc(i)>EMTT)&&(ATTEMPT_OPT>0)
              
              [vpl_adj, hpl_adj, sig_acc_adj, emt_adj] = mhss_raim_baseline_v4(G(idx,:), sigpr2_int(idx), sigpr2_acc(idx), ...
             nom_bias_int(idx), nom_bias_acc(idx), p_sat(idx), p_const,ATTEMPT_OPT,rho_j...
                  , subsets_i, pap_subset, p_not_monitored);
             
               if (vpl_adj<vpl_exc(i))
                       vpl_exc(i) = vpl_adj;
                       emt_exc(i) = emt_adj;
                       sig_acc_exc(i) = sig_acc_adj;
                   
               end
               
               if (hpl_adj<hpl_exc(i))           
                       hpl_exc(i) = hpl_adj;
               end
               
               
               
               
            end
         end
    end
    hpl = max(hpl,max(hpl_exc(1:Nsubsets))); 
    vpl = max(vpl,max(vpl_exc(1:Nsubsets)));  
    emt = max(emt,max(emt_exc(1:Nsubsets)));
    sig_acc = max(sig_acc, max(sig_acc_exc(1:Nsubsets)));
     
  else
    
      %% No detection FDE performance
      if PL0_FDE
          rho_j = 1/n_view;
      else
          rho_j = 1;
      end
      
      
   [vpl_bl, hpl_bl, sig_acc_bl, emt_bl] = mhss_raim_baseline_v4(G, sigpr2_int, sigpr2_acc, nom_bias_int, nom_bias_acc, p_sat, p_const,0, rho_j);
    
    vpl = vpl_bl;
    emt = emt_bl;
    sig_acc = sig_acc_bl;
    hpl = hpl_bl;
    
 
   
   
   if (sig_acc<SIG_ACC_MAX_VERT)&&((vpl>VPLT)||(hpl>HPLT)||emt>EMTT)&&(ATTEMPT_OPT>0)
     %[vpl_adj, hpl_adj, sig_acc_adj, emt_adj] = mhss_raim_baseline_adjusted_position_v4(G, sigpr2_int, sigpr2_acc, nom_bias_int, nom_bias_acc, p_sat, p_const,rho_j);
     [vpl_adj, hpl_adj, sig_acc_adj, emt_adj] = mhss_raim_baseline_v4(G, sigpr2_int, sigpr2_acc, nom_bias_int, nom_bias_acc, p_sat, p_const,ATTEMPT_OPT, rho_j);
     if (vpl_adj<vpl)||(hpl_adj<hpl)
        vpl = vpl_adj;
        emt = emt_adj;
        sig_acc = sig_acc_adj;
        hpl = hpl_adj;
     end
     
   end
  end
   vhpl(usr,1)= vpl;
   vhpl(usr,2)= emt;
   vhpl(usr,3)= sig_acc;
   vhpl(usr,4)= hpl;
   
  else                      
      vhpl(usr,1)=Inf;
      vhpl(usr,2)=Inf;                 
      vhpl(usr,3)=Inf;
      vhpl(usr,4)=Inf;
  end
  
   
end


%toc







