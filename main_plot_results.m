%using plot_no_gui

clear; close; clc;

load('outputs.mat')

VAL = 35;
HAL = 40;
EMT_th= 15;
ACC_th = 1.87;
% VAL = 50;
% HAL = 40;
% EMT_th= Inf;
% ACC_th = Inf;
%usrlatgrid = -85:5:85;

avail = .995;
latmax = 70/90*pi/2;

%Outputs = [availability, VPL/HPL, EMT, sig_acc]
outputs = [1 1 1 1];
percent = 0.995;
%vhal = [35, 40];
vhal = [VAL, HAL, EMT_th, ACC_th];
pa_mode = 1;

%emt = emt_new
% plot_no_gui(satdata,usrdata,wrsdata,igpdata,inv_igp_mask,...
%                        sat_xyz,udrei,givei,vpl,hpl,usrlatgrid,...
% 					   usrlongrid,outputs,percent,vhal,pa_mode,...
%                        udre_hist,give_hist, udrei_hist,givei_hist,...
%                        emt,sig_acc, EMT_th, ACC_th, latmax)

plot_no_gui(usrdata,vpl,hpl,usrlatgrid,usrlongrid,outputs,percent,vhal,pa_mode,...
            emt_new,sig_acc, latmax)
