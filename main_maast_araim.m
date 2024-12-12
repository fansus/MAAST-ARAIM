%*************************************************************************
%*     Copyright c 2018 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%*************************************************************************
%
%MAIN_ARAIM_MAAST runs maast for ARAIM.
% all settings are edited in this file (RAIM ISM parameters, the user grid, the constellation
%configuration, and the receiver settings)
%allows for repeatable and recorded runs
%The parameters are described in:	Blanch, J., Walter, T., Enge, P., Lee, Y., Pervan, B., Rippl, M., Spletter, A., Kropp, V.,
%"Baseline Advanced RAIM User Algorithm and Possible Improvements," IEEE Transactions on Aerospace and Electronic Systems, 
%Volume 51,  No. 1, January 2015.

%In the position optimization option, the code attempts to lower the
%protection levels by applying the method described in: 
% Blanch, J., Walter, T., Enge, P., Kropp, V.,”A Simple Position Estimator that Improves Advanced RAIM Performance,” 
% IEEE Transactions on Aerospace and Electronic Systems Vol. 51, No. 3, July 2015.
%Some of the constants might be set to a different value than in the paper

%The protection levels are computed following the algorithm described in
%the ARAIM Airborne Algorithm Algorithm Description Document v3.0

%Created 2018 June 14 by Juan Blanch

% global TRUTH_FLAG
% 
% global BRAZPARAMS RTR_FLAG IPP_SPREAD_FLAG

global ARAIM_URA_GPS ARAIM_URA_GLO ARAIM_URA_GAL ARAIM_URA_BDU ...
       ARAIM_BIAS_GPS ARAIM_BIAS_GLO ARAIM_BIAS_GAL ARAIM_BIAS_BDU

global ARAIM_URE_GPS ARAIM_URE_GLO ARAIM_URE_GAL ARAIM_URE_BDU ...
       ARAIM_BIAS_CONT_GPS ARAIM_BIAS_CONT_GLO ARAIM_BIAS_CONT_GAL ARAIM_BIAS_CONT_BDU
   
global ARAIM_PSAT_GPS ARAIM_PSAT_GAL ARAIM_PSAT_GLO ARAIM_PSAT_BDU ...
       ARAIM_PCONST_GPS ARAIM_PCONST_GAL ARAIM_PCONST_GLO ARAIM_PCONST_BDU

global PHMI_VERT PHMI_HOR P_THRES PFA_VERT PFA_HOR P_EMT PL_TOL FC_THRES ...
    PL0_FDE FDE_FLAG FDE_WF_FLAG

%global P_EMT_NEW 
global SIG_ACC_MAX_VERT
global SIG_ACC_MAX_HOR1 SIG_ACC_MAX_HOR2
global VPLT HPLT EMTT ATTEMPT_OPT

global ARAIM_USRMASK_GPS ARAIM_USRMASK_GLO ARAIM_USRMASK_GAL ARAIM_USRMASK_BDU 
global ARAIM_SIN_USRMASK_GPS ARAIM_SIN_USRMASK_GLO ARAIM_SIN_USRMASK_GAL
global ARAIM_SIN_USRMASK_BDU PFAULT_EXC_THRES




%Integrity Support Message parameters
ARAIM_URA_GPS = 1;
ARAIM_URA_GAL = 1;
ARAIM_URA_GLO = 1;
ARAIM_URA_BDU = Inf;

ARAIM_BIAS_GPS = .75;
ARAIM_BIAS_GAL = .75;
ARAIM_BIAS_GLO = .75;
ARAIM_BIAS_BDU = Inf;

ARAIM_URE_GPS = 2/3*ARAIM_URA_GPS;
ARAIM_URE_GAL = 2/3*ARAIM_URA_GAL;
ARAIM_URE_GLO = 2/3*ARAIM_URA_GLO;
ARAIM_URE_BDU = Inf;

ARAIM_BIAS_CONT_GPS = 0;
ARAIM_BIAS_CONT_GAL = 0;
ARAIM_BIAS_CONT_GLO = 0;
ARAIM_BIAS_CONT_BDU = 0;

ARAIM_PSAT_GPS = 1e-5;
ARAIM_PSAT_GAL = 1e-5;
ARAIM_PSAT_GLO = 1e-3;
ARAIM_PSAT_BDU = 1;

ARAIM_PCONST_GPS = 1e-8; %For H-ARAIM, PCONST_GPS can be set to 10^-8
ARAIM_PCONST_GAL = 1e-4;
ARAIM_PCONST_GLO = 1e-4;
ARAIM_PCONST_BDU = 1;

FDE_FLAG = 0; %Computes worst case PLs,EMT, sig_acc under a fault or outage scenario. In this version, only single outages and faults are taken into account for the continuity assessment.
%Probability of failed exclusion is set at PFA

FDE_WF_FLAG = 0;

PL0_FDE = 0; %When FDE_FLAG is off,and this flag is on the PLS are computed assuming that some of the integrity budget is reserved for the exclusion function 


%Mask angles

ARAIM_USRMASK_GPS = 5;
ARAIM_USRMASK_GLO = 5;
ARAIM_USRMASK_GAL = 5;
ARAIM_USRMASK_BDU = 5;

ARAIM_SIN_USRMASK_GPS = sin(ARAIM_USRMASK_GPS*pi/180);
ARAIM_SIN_USRMASK_GLO = sin(ARAIM_USRMASK_GLO*pi/180);
ARAIM_SIN_USRMASK_GAL = sin(ARAIM_USRMASK_GAL*pi/180);
ARAIM_SIN_USRMASK_BDU = sin(ARAIM_USRMASK_BDU*pi/180);

%ARAIM receiver parameters

%Vertical-ARAIM
PHMI_VERT = 9.8e-8;
P_THRES = 6e-8;
PFA_VERT = 3.9e-6;
PFA_HOR  = 0.9e-6;
PL_TOL = 1e-2;
PHMI_HOR = 2e-9;
P_EMT = 1e-5;
FC_THRES = .01;

VPLT = 35;
HPLT = 40;
EMTT = 15;
% 
SIG_ACC_MAX_VERT = 1.86;
SIG_ACC_MAX_HOR1 = 3;  %Do not set to infinite
SIG_ACC_MAX_HOR2 = 3;  %Do not set to infinite


% Horizontal-ARAIM

% PHMI_VERT = 0.1e-8;
% P_THRES = 6e-8;
% PFA_VERT = .01e-7;
% PFA_HOR  = 1e-7;
% PL_TOL = 1e-2;
% PHMI_HOR = 9.9e-8;
% P_EMT = 1e-5;
% FC_THRES = .01;
% 
% % 
% VPLT = Inf;
% HPLT = 185;
% EMTT = Inf;
% 
% SIG_ACC_MAX_VERT = 10;
% SIG_ACC_MAX_HOR1 = 10;  %Do not set to infinite
% SIG_ACC_MAX_HOR2 = 10;  %Do not set to infinite



ATTEMPT_OPT = 1; %determines whether we attempt to optimize the position estimation:
% ATTEMPT_OPT = 0: no optimization
% ATTEMPT_OPT = 1: 1-dimensional optimization



%%%%%%%%
% SV Menu
%WGC scenarios

%svfile = {'almmops-1.txt','almanac Galileo 24-1 Week 703.alm.txt'};
%svfile = {'almmops.txt','almanac Galileo 24 Week 703.alm.txt','almglonass.txt'};
svfile = {'almmops.txt','almanac Galileo 24 Week 703.alm.txt'};
%svfile ={'almgps24+3.txt','almanac Galileo 24 + 3 Spare Week 703.alm.txt'};
%svfile ={'almgps24+3.txt','almanac Galileo 24 + 8 Spare Week 703.alm.txt'};
%svfile = {'almmops-1.txt'};
%svfile = {'almmops.txt'};
%svfile ={'almgps24+3.txt'};
%svfile = {'almmops_22.txt'};

%Other scenarios

%svfile = {'current.txt'};  
%svfile =  {'almmops.txt','almgalileo.txt','almglonass.txt'};
% svfile =  {'almmops.txt','almgalileo24.txt'};

init_const;      % global physical and gps constants
init_col_labels; % column indices 
init_mops;       % MOPS constants

close all;


%User signals:
% dual_freq = 0 : L1 only
% dual_freq = 1 : L1-L5
% dual_freq = 2 : L5 only


 dual_freq = 1;            
      

% USER CNMP Menu

      %select AAD-B model
      usrcnmpfun = 'af_cnmp_mops';      
      init_cnmp_mops;
      
      %select AAD-A model
%      usrcnmpfun = 'af_cnmpaad';
%      init_aada;
      
      %select AAD-B model
%      usrcnmpfun = 'af_cnmpaad';      
%      init_aadb;


       

% USER Menu
%select the world as the user area
      usrpolyfile = 'usrworld.dat';
      usrlatstep = 10;
      usrlonstep = 10;
      
       %Start time for simulation
      TStart = 0;
      
      %End time for simulation
      %TEnd = 862200;
      TEnd = 86400;
      
      % Size of time step
      TStep = 300;
      
%select CONUS as the user area
%      usrpolyfile = 'usrconus.dat';
      
      %select Alaska as the user area
%      usrpolyfile = 'usralaska.dat';
      
      %select Canada as the user area
%      usrpolyfile = 'usrcanada.dat';
      
      %select Mexico as the user area
%      usrpolyfile = 'usrmexico.dat';
      
      %select North America as the user area
      %usrpolyfile = 'usrn_america.dat';
      
      %select Europe as the user area
%      usrpolyfile = 'usreurope.dat';
      
      %select Japan as the user area
%      usrpolyfile = 'usrmsas.dat';
      
      %select Brazil as the user area
%      usrpolyfile = 'usrbrazil.dat';
      
      
    



        % check if file(s) exist
        i=1;
        while i<=size(svfile,2)
          if iscell(svfile)
            fid=fopen(svfile{i});
          else
            fid=fopen(svfile);
            i = size(svfile,2);
          end
          if fid==-1
              fprintf('Almanac file not found.  Please try again.\n');
              return;
          else
              fclose(fid);
          end 
          i=i+1;
        end
        
     
      


% Mode / Alert limit

      %choose PA mode vs NPA  
      pa_mode = 1;
      
      %choose VAL, HAL, EMT Threshold, and sigma accuracy threshold
      vhal = [35, 40, 15, 1.87];
     %vhal = [Inf, 1668];
      
% OUTPUT Menu

      %initialize histograms
      %init_hist;
        
      
 
      % turn on or off output options
      %1: Availability  2: V/HPL  3: EMT  4: sig_acc  
      
      outputs = [1 1 1 1];

      % Assign percentage
      percent = 0.995; % 1 = 100%
      
      % Coverage values computed for users with |lat|<lamax
      latmax = 90/90*pi/2;
      
      
% WRS Menu (not used in ARAIM)
 
      %wrsfile = 'wrs_foc.dat';
      
      
        
% RUN Simulation

      svmrun(usrcnmpfun,usrpolyfile, svfile, TStart, TEnd, ...
             TStep, usrlatstep, usrlonstep, outputs, percent, vhal, ...
             pa_mode, dual_freq, latmax);



