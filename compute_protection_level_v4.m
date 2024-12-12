function [vpl, alloc] = compute_protection_level_v2mod(sigma, threshold_plus_bias, pfault, phmi, pl_tol, alloc_max)
%*************************************************************************
%*     Copyright c 2015 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%*************************************************************************

%Created 14 August 2012 by Juan Blanch
%Updated 24 February 2015 by Juan Blanch
%COMPUTE_PROTECTION_LEVEL_V2 solves the equation:   pfault.*normcdf((threshold_plus_bias - vpl)./sigma) = phmi
%alloc is the vector of corresponding allocations normcdf((threshold_plus_bias - vpl)./sigma) 
%alloc max is the maximum necessary allocation

% 检查参数 alloc_max 是否传入，如果没有传入，默认为每个元素的最大分配为 1
if nargin<6
    alloc_max = ones(length(sigma),1);  % 如果没有传入 alloc_max，则将其初始化为全 1 向量
end

%%%% Exclude sigmas that are infinite and evaluate their integrity contribution
% 排除 sigma 为无穷大的元素，并计算其对应的完整性贡献
index_Inf = find(sigma == Inf);  % 找出所有 sigma 为无穷大的索引
index_fin = setdiff(1:length(sigma),index_Inf);  % 获取有效的索引，去除无穷大值
p_not_monitorable = sum(pfault(index_Inf));  % 计算不可监控的卫星故障概率

% 如果不可监控的概率大于等于完整性要求，则返回无穷大保护等级
if p_not_monitorable >= phmi    
    vpl = Inf;  % 设置保护等级为无穷大
    alloc = zeros(length(sigma),1);  % 分配为全 0 向量
else
    % 去除 sigma 为无穷大的元素，重新定义有效数据
    sigma = sigma(index_fin);
    threshold_plus_bias = threshold_plus_bias(index_fin);
    pfault = pfault(index_fin);
    phmi = phmi - p_not_monitorable;  % 更新完整性预算
    maxCount = 10;  % 最大迭代次数

    % 计算保护等级的下界
    Klow = -norminv(min(1, phmi ./ (pfault .* alloc_max)));  % 计算保护等级下界因子
    vpl_low = max(threshold_plus_bias + Klow .* sigma);  % 计算保护等级下界

    % 计算保护等级的上界
    Khigh = max(0, -norminv(phmi ./ (length(sigma) * pfault)));  % 计算保护等级上界因子
    vpl_high = max(threshold_plus_bias + Khigh .* sigma);  % 计算保护等级上界

    % 计算完整性预算的对数
    log10phmi = log10(phmi);  

    count = 0;
    % 使用二分法迭代逼近保护等级，直到满足精度要求或达到最大迭代次数
    while ((vpl_high - vpl_low > pl_tol) && (count < maxCount))
        count = count + 1;
        vpl_half = (vpl_low + vpl_high) / 2;  % 计算保护等级的中值
        cdfhalf = log10(sum(pfault .* min(modnormcdf((threshold_plus_bias - vpl_half) ./ sigma), alloc_max)));  % 计算当前中值的累计分布函数

        % 根据当前中值与目标完整性预算的对数比较，调整上下界
        if cdfhalf > log10phmi
            vpl_low = vpl_half;  % 如果当前中值的 CDF 大于目标完整性预算，则更新下界
        else
            vpl_high = vpl_half;  % 否则更新上界
        end
    end

    vpl = vpl_high;  % 最终保护等级为上界值
    alloc = min(modnormcdf((threshold_plus_bias - vpl) ./ sigma), alloc_max);  % 计算对应的分配，限制不超过最大分配值
end
%End
