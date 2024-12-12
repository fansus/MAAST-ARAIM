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
%USR_VHPL_ARAIM 计算用户的垂直保护等级（VPL）和水平保护等级（HPL）
%给定n_los个用户视线向量（在ECEF WGS-84坐标系中），并输入相关的用户标识号、完整性方差、精度方差、名义偏差、连续性偏差和卫星故障概率等，本函数将计算垂直保护等级（VPL）、水平保护等级（HPL）、精度的标准差（SIG_ACC）以及有效监控阈值（EMT 和 EMT_NEW）等。
%   调用MHSS_RAIM_BASELINE计算保护等级、精度等。

% 2013年8月13日 Juan Blanch 创建
% 2015年3月18日 Juan Blanch 修改（包括一维位置优化）

global MOPS_NOT_MONITORED 

global MOPS_MIN_GPSPRN MOPS_MAX_GPSPRN MOPS_MIN_GLOPRN MOPS_MAX_GLOPRN ...
       MOPS_MIN_GALPRN MOPS_MAX_GALPRN  ...
       MOPS_MIN_BDUPRN MOPS_MAX_BDUPRN
global ARAIM_PCONST_GPS ARAIM_PCONST_GAL ARAIM_PCONST_GLO ARAIM_PCONST_BDU  
global SIG_ACC_MAX_VERT SIG_ACC_MAX_HOR1 SIG_ACC_MAX_HOR2 VPLT HPLT EMTT FDE_FLAG
global ATTEMPT_OPT PL0_FDE FDE_WF_FLAG

% 初始化返回值
n_usr = max(usr_idx);  % 获取用户数量
[n_los, ~] = size(los_xyzb);  % 获取视线数量
vhpl = repmat(MOPS_NOT_MONITORED, n_usr, 4);  % 初始化 vhpl 矩阵

% 对每个用户进行处理
for usr = 1:n_usr
    
    if mod(usr,500) == 0
        disp(usr);  % 每处理500个用户打印一次
    end
    
    if usr == 3233
        disp(usr);  % 特定用户的调试信息
    end
    
    % 查找各卫星系统的索引
    svidxgps = find((usr_idx == usr) & (prn >= MOPS_MIN_GPSPRN & prn <= MOPS_MAX_GPSPRN) & (sig2_i < Inf) & (psat_i < 1));
    svidxgal = find((usr_idx == usr) & (prn >= MOPS_MIN_GALPRN & prn <= MOPS_MAX_GALPRN) & (sig2_i < Inf) & (psat_i < 1));
    svidxglo = find((usr_idx == usr) & (prn >= MOPS_MIN_GLOPRN & prn <= MOPS_MAX_GLOPRN) & (sig2_i < Inf) & (psat_i < 1));
    svidxbdu = find((usr_idx == usr) & (prn >= MOPS_MIN_BDUPRN & prn <= MOPS_MAX_BDUPRN) & (sig2_i < Inf) & (psat_i < 1));
    
    % 各卫星系统的数量
    ngps = length(svidxgps);
    ngal = length(svidxgal);
    nglo = length(svidxglo);
    nbdu = length(svidxbdu);
    
    % 记录卫星系统的可用性
    svconst = [(~isempty(svidxgps)) (~isempty(svidxgal)) (~isempty(svidxglo)) (~isempty(svidxbdu))];
    
    % 总卫星数目和有效卫星的索引
    nsat = ngps + ngal + nglo + nbdu;
    svidx = [svidxgps; svidxgal; svidxglo; svidxbdu];
    
    % 创建时钟标识向量
    clk_gps = [ones(ngps, svconst(1)); zeros(nsat - ngps, svconst(1))];
    clk_gal = [zeros(ngps, svconst(2)); ones(ngal, svconst(2)); zeros(nsat - ngps - ngal, svconst(2))];    
    clk_glo = [zeros(ngps + ngal, svconst(3)); ones(nglo, svconst(3)); zeros(nbdu, svconst(3))];
    clk_bdu = [zeros(nsat - nbdu, svconst(4)); ones(nbdu, svconst(4))];
    
    % 创建 G 矩阵，包含了所有卫星的位置信息和时钟标识
    G = [los_xyzb(svidx, 1:3) clk_gps clk_gal clk_glo clk_bdu];

    % 创建常数矩阵
    p_const = [ones(svconst(1), 1) * ARAIM_PCONST_GPS; ones(svconst(2), 1) * ARAIM_PCONST_GAL; ...
               ones(svconst(3), 1) * ARAIM_PCONST_GLO; ones(svconst(4), 1) * ARAIM_PCONST_BDU];
    
    n_view = length(svidx);
    
    if(n_view > 3)
        % 提取视线对应的完整性方差、精度方差、名义偏差和卫星故障概率
        sigpr2_int = sig2_i(svidx);
        sigpr2_acc = sig2acc_i(svidx);
        nom_bias_int = bnom_i(svidx);
        nom_bias_acc = bcont_i(svidx);
        p_sat = psat_i(svidx);
        
        if FDE_FLAG
            % 如果启用了 FDE（故障检测与隔离），则构建故障子集
            subsets_exc = ones(n_view, n_view) - eye(n_view);
            if FDE_WF_FLAG
                idwf = find(p_const > 1e-7);
                subsets_wf = ones(length(idwf), n_view) - G(:, idwf + 3)';
                subsets_exc = [subsets_exc; subsets_wf];    
            end
            
            Nsubsets = size(subsets_exc, 1);
            Nexc = size(subsets_exc, 1);
            hpl_exc = Inf * ones(Nsubsets, 1);
            vpl_exc = Inf * ones(Nsubsets, 1);
            sig_acc_exc = Inf * ones(Nsubsets, 1);
            emt_exc = Inf * ones(Nsubsets, 1);
            rho_j = 1 / (Nexc + 1);
            
            % 调用 RAIM 基线方法计算保护等级
            [vpl, hpl, sig_acc, emt, subsets, pap_subset, p_not_monitored] = mhss_raim_baseline_v4(G, sigpr2_int, sigpr2_acc, nom_bias_int, nom_bias_acc, p_sat, p_const, 0, rho_j);
            
            % 如果精度不超过最大垂直精度，并且 VPL 或 HPL 超过阈值，则进行位置调整
            if (sig_acc < SIG_ACC_MAX_VERT) && ((vpl > VPLT) || (hpl > HPLT) || emt > EMTT) && (ATTEMPT_OPT > 0)
                [vpl_adj, hpl_adj, sig_acc_adj, emt_adj, subsets, pap_subset, p_not_monitored] = mhss_raim_baseline_v4(G, sigpr2_int, sigpr2_acc, nom_bias_int, nom_bias_acc, p_sat, p_const, ATTEMPT_OPT, rho_j);
                if (vpl_adj < vpl) || (hpl_adj < hpl)
                    vpl = vpl_adj;
                    emt = emt_adj;
                    sig_acc = sig_acc_adj;
                    hpl = hpl_adj;
                end
            end
            
            % 计算不同子集的保护等级
            for i = 1:Nsubsets    
                idx = logical(subsets_exc(i, :));
                subsets_i = subsets(:, idx);
                if sum(idx) > 3
                    [vpl_bl, hpl_bl, sig_acc_bl, emt_bl] = mhss_raim_baseline_v4(G(idx, :), sigpr2_int(idx), sigpr2_acc(idx), ...
                        nom_bias_int(idx), nom_bias_acc(idx), p_sat(idx), p_const, 0, rho_j, subsets_i, pap_subset, p_not_monitored);
                    vpl_exc(i) = vpl_bl;
                    emt_exc(i) = emt_bl;
                    sig_acc_exc(i) = sig_acc_bl;
                    hpl_exc(i) = hpl_bl;
                    
                    % 如果精度不超过最大垂直精度，并且 VPL 或 HPL 超过阈值，则进行位置调整
                    if (sig_acc_exc(i) < SIG_ACC_MAX_VERT) && ((vpl_exc(i) > VPLT) || (hpl_exc(i) > HPLT) || emt_exc(i) > EMTT) && (ATTEMPT_OPT > 0)
                        [vpl_adj, hpl_adj, sig_acc_adj, emt_adj] = mhss_raim_baseline_v4(G(idx, :), sigpr2_int(idx), sigpr2_acc(idx), ...
                            nom_bias_int(idx), nom_bias_acc(idx), p_sat(idx), p_const, ATTEMPT_OPT, rho_j, subsets_i, pap_subset, p_not_monitored);
                        
                        if (vpl_adj < vpl_exc(i))
                            vpl_exc(i) = vpl_adj;
                            emt_exc(i) = emt_adj;
                            sig_acc_exc(i) = sig_acc_adj;
                        end
                        
                        if (hpl_adj < hpl_exc(i))           
                            hpl_exc(i) = hpl_adj;
                        end
                    end
                end
            end
            
            % 计算最大保护等级
            hpl = max(hpl, max(hpl_exc(1:Nsubsets))); 
            vpl = max(vpl, max(vpl_exc(1:Nsubsets)));  
            emt = max(emt, max(emt_exc(1:Nsubsets)));
            sig_acc = max(sig_acc, max(sig_acc_exc(1:Nsubsets)));
            
        else
            % 没有 FDE 时的计算
            if PL0_FDE
                rho_j = 1 / n_view;
            else
                rho_j = 1;
            end
            
            [vpl_bl, hpl_bl, sig_acc_bl, emt_bl] = mhss_raim_baseline_v4(G, sigpr2_int, sigpr2_acc, nom_bias_int, nom_bias_acc, p_sat, p_const, 0, rho_j);
            
            vpl = vpl_bl;
            emt = emt_bl;
            sig_acc = sig_acc_bl;
            hpl = hpl_bl;
            
            % 如果精度不超过最大垂直精度，并且 VPL 或 HPL 超过阈值，则进行位置调整
            if (sig_acc < SIG_ACC_MAX_VERT) && ((vpl > VPLT) || (hpl > HPLT) || emt > EMTT) && (ATTEMPT_OPT > 0)
                [vpl_adj, hpl_adj, sig_acc_adj, emt_adj] = mhss_raim_baseline_v4(G, sigpr2_int, sigpr2_acc, nom_bias_int, nom_bias_acc, p_sat, p_const, ATTEMPT_OPT, rho_j);
                if (vpl_adj < vpl) || (hpl_adj < hpl)
                    vpl = vpl_adj;
                    emt = emt_adj;
                    sig_acc = sig_acc_adj;
                    hpl = hpl_adj;
                end
            end
        end
        % 保存计算结果
        vhpl(usr, 1) = vpl;
        vhpl(usr, 2) = emt;
        vhpl(usr, 3) = sig_acc;
        vhpl(usr, 4) = hpl;
    else
        % 如果视线数小于等于3，设置为不监控
        vhpl(usr, 1) = Inf;
        vhpl(usr, 2) = Inf;                 
        vhpl(usr, 3) = Inf;
        vhpl(usr, 4) = Inf;
    end
end

%toc
