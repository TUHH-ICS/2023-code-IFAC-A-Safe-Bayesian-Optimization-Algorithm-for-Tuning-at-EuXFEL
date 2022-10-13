function [L_f,L_g] = updateLipschitz(x,y,c,L_fold,L_gold)
    y_n = y(end);
    x_n = x(end,:);
    y_diff = abs(y_n-y(1:end-1));
    x_norm = vecnorm(x_n-x(1:end-1,:),2,2);
    L_f = max(L_fold,max(y_diff./x_norm));
    c_n = c(end,:);
    c_diff = abs(c_n-c(1:end-1,:));
    L_g = max(L_gold,max(c_diff./x_norm,[],1));
end