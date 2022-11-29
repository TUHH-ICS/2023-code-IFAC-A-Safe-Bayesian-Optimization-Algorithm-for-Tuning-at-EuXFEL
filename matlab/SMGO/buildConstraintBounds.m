function [g_app,g_unc] = buildConstraintBounds(c,L_g,xk,x,mu)
    S = size(c,2);
    xk = permute(xk,[3,2,1]);
    x_norm=squeeze(vecnorm(x-xk,2,2))';

    g_u = zeros(size(x,1),S);
    g_l = zeros(size(x,1),S);
    
    for i=1:S
        g_u(:,i) = min(c(:,i)+mu*L_g(i)*x_norm,[],1)';
        g_l(:,i) = max(c(:,i)-mu*L_g(i)*x_norm,[],1)';
    end
    
    g_app=0.5*(g_u+g_l);
    g_unc = g_u-g_l;
end

