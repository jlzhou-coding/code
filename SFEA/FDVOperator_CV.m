function CV = FDVOperator_CV(CV,i)
%normalization
if i==0%stage 1
    [m,n]=size(CV);
    CV=zeros(m,n);
else
    v_min=min(CV);
    v_max=max(CV);
    v_t=(CV-v_min)./(v_max-v_min+1e-8);
    lower=0;
    upper=1;
    R=upper-lower;
    gamma_a = R*10^-i.*floor(10^i*R.^-1.*(v_t-lower)) + lower;
    gamma_b = R*10^-i.*ceil(10^i*R.^-1.*(v_t-lower)) + lower;
    miu1    = 1./(v_t-gamma_a);
    miu2    = 1./(gamma_b-v_t);
    logical = miu1-miu2>0;
    CV  = gamma_b;%CV
    CV(find(logical)) = gamma_a(find(logical));%CV 
end
end