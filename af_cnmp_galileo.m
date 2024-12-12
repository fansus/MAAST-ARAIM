function sig2=af_cnmp_galileo(el)
%*************************************************************************
%*     Copyright c 2011 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                              *
%*************************************************************************
%
%SIG2_AAD calculate airborne pseudorange confidence (variance) for
%L1 L5 Galileo per ToolCrossCheckScenarioDefinitionV1.4a
%

%created 11 May, 2011 by Juan Blanch

eldeg = max(el*180/pi,5);

bin_size = 5;

idx_bin_low = floor(eldeg/bin_size);
idx_bin_hi  = idx_bin_low+1;

el_table = [5:bin_size:90]';


sig_table = [0.4529 0.3553 0.3063  0.2638  0.2593  0.2555 0.2504 0.2438  0.2396...
    0.2359 0.2339 0.2302 0.2295 0.2278 0.2297 0.2310 0.2274 0.2277 0.2277]';

sig = sig_table(idx_bin_low) +( (eldeg - el_table(idx_bin_low))/...
bin_size).* (sig_table(idx_bin_hi)-sig_table(idx_bin_low));   

sig2 = sig.^2;
