function [f_app,f_unc,f_u,f_l] = buildObjectiveBounds(y,L_f,xk,x,mu)
    xk = permute(xk,[3,2,1]);
    x_norm=squeeze(vecnorm(x-xk,2,2))';
    f_u = min(y+mu*L_f*x_norm,[],1)';
    f_l = max(y-mu*L_f*x_norm,[],1)';
    f_app=0.5*(f_u+f_l);
    f_unc = f_u-f_l;
end

