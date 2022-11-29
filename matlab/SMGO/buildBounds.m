function [f_app,g_app,f_unc,g_unc,f_u,f_l] = buildBounds(y,c,L_f,L_g,xk,x,mu)
    [f_app,f_unc,f_u,f_l]=buildObjectiveBounds(y,L_f,xk,x,mu);
    [g_app,g_unc] = buildConstraintBounds(c,L_g,xk,x,mu);
end