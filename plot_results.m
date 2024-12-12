%using plot_no_gui

clear; close; clc;
init_graph4();
init_col_labels();
load('outputs')


                   
%choose PA mode vs NPA  
pa_mode = 1;
      
%choose VAL, HAL, EMT Threshold, and sigma accuracy threshold
       vhal = [35, 40, 15, 1.8700];
      %vhal = [35, Inf,Inf,Inf];
%       vhal = [50, 40, Inf, Inf];
%vhal = [Inf, 185, Inf, Inf] ;     
% turn on or off output options
%1: Availability  2: V/HPL  3: EMT  4: sig_acc  
      
outputs = [1 1 1 1];

% Assign percentage
percent = 0.995; % 1 = 100%
      
% Coverage values computed for users with |lat|<latmax
latmax = 90/90*pi/2;
      
plot_no_gui(usrdata,vpl,hpl,usrlatgrid,usrlongrid,outputs,percent,vhal,pa_mode,...
            emt,sig_acc, latmax)                       