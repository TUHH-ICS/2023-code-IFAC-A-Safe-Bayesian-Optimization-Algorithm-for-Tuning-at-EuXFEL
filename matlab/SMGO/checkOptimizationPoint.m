function b=checkOptimizationPoint(f_est,yopt,alpha,L_f)
    b = f_est <= yopt-alpha*L_f;
end

