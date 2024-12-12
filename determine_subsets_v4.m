function  [subsets, pap_subset, p_not_monitored] = determine_subsets_v4(G, p_sat, p_const, p_thres, fc_thres)
%*************************************************************************
%*     Copyright c 2015 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%*************************************************************************

%Created september 2014 by Juan Blanch
%Updated april 2016 by Juan Blanch 
%DETERMINE_SUBSETS_NEW integrates changes proposed in WGC ARAIM briefing
%from April 2016



Nsat = size(G,1);
Nconst = size(G,2)-3;
if max(size(p_sat))==1              %if p_sat is a scalar, apply to all satellites
   p_sat = ones(Nsat,1)*p_sat;
end
if max(size(p_const))==1
   p_const = ones(Nconst,1)*p_const;
end
%Remove constellations for which there are no satellites
idc = find(sum(abs(G(:,4:end))>0));
G = G(:,[1:3 3+idc]);
Nconst = size(G,2)-3;
p_const = p_const(idc);


%Compute probability of no fault

pnofault = prod(1-p_sat)*prod(1-p_const);

%Initialize pnotmonitored

p_not_monitored = 1-pnofault;


%Determine upper bound of subset size
[nsatconstmax, ~] = determine_nmax([p_sat;p_const],p_thres);
subsetsize =0;
for j=0:nsatconstmax
    subsetsize = subsetsize + nchoosek(Nconst+Nsat,j);
end


subsets_ex = zeros(subsetsize,Nsat+Nconst);
subsets_ex(1,:) = zeros(1,Nsat+Nconst); %subset corresponding to the all-in-view
pap_subset = zeros(subsetsize,1);
pap_subset(1) = 1;
%Initialize k (number of simultaneous faults) and subset index j
k =0;
j = 1;

while (k<=Nsat)&&(p_not_monitored>p_thres)
    k=k+1;
    subsets_k = determine_k_subsets(Nsat+Nconst,k);
    pap_subsets_k = prod((subsets_k*diag([p_sat./(1-p_sat);p_const./(1-p_const)])...
        + (1 - subsets_k)),2)*pnofault;
    
    %sort subsets by decreasing probability
    [p_subsets_k_s, index_k_s]=sort(pap_subsets_k,1,'descend');
    subsets_k_s = subsets_k(index_k_s,:);
    
    n_k = size(subsets_k,1);
    h=0;
    while (h<n_k)&&(p_not_monitored>p_thres)
        h=h+1;
        if p_subsets_k_s(h)>0
         j=j+1;
         subsets_ex(j,:) = subsets_k_s(h,:);
         pap_subset(j) = p_subsets_k_s(h);
         p_not_monitored = p_not_monitored - pap_subset(j);
        end
    end
end
subsets_ex = subsets_ex(1:j,:);
pap_subset = pap_subset(1:j);
idconst = find(sum(subsets_ex(:,end-length(idc)+1:end))>0);

%find subsets included in one constellation

nconst_m = length(idconst);



%Transform Nsat+Nconst subset vector in Nsat subset vector
G = G(:,[1 2 3 3+idc]);
subsets_const = (G(:,4:(3+Nconst))*subsets_ex(:,Nsat+1:Nsat+Nconst)')';
subsets_sat  = min(subsets_ex(:,1:Nsat) + subsets_const ,1);

   


%subset consolidation

for  jj=1:nconst_m
    %find subset corresponding to constellation wide fault first
    id_const_c=find((sum(subsets_ex(:,1:Nsat),2)==0).*...
        (sum(subsets_ex,2)==1).*(subsets_ex(:,Nsat+idconst(jj))==1));
    %subsets_sat(id_const_c,:)
    
    %find subsets that include a constellation fault and satellites whithin
    %that constellation
    
    index_cs =find((((subsets_sat*subsets_sat(id_const_c,:)'))>0).*...
    (((subsets_sat*(1-subsets_sat(id_const_c,:))'))<=0));
    %subsets_sat(index_cs,:)
    
    idremove = find(pap_subset(index_cs)<fc_thres*pap_subset(id_const_c));
    idrmv    = index_cs(idremove);
    %remove from the list and add probability to constellation wide fault
    
    pap_subset(id_const_c) =pap_subset(id_const_c) + sum(pap_subset(idrmv));
    
    idnew = setdiff(1:length(pap_subset),idrmv);
    
    pap_subset = pap_subset(idnew);
    subsets_sat = subsets_sat(idnew,:);
%     idoutconst = setdiff(1:Nconst,idconst(jj));
%     find(G(:,3+idconst(jj))==1);
%     
%     id_const_cs = 

    
end
%subsets_sat

subsets =1 - subsets_sat;




%End