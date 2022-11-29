function [E] = buildSamplePoints(E_old,x,B,cond)
    D = size(x,2);
    x_n = x(end,:);
    b_p = cond(:,2) - x_n';
    b_m = x_n' - cond(:,1);
    
    k = (1:B-1)'/B;
    K = zeros((B-1)*D,D);
    for i = 1 : D
        K(1+(i-1)*(B-1):(B-1)+(i-1)*(B-1),i) = k;
    end
    Yd_p = x_n+(K.*b_p');
    Yd_m = x_n-(K.*b_m');
    Yd = union(Yd_m,Yd_p,'rows');
    
    ad = x(1:end-1,:)-x_n;
    if ~isempty(ad)
        norm_ = vecnorm(ad,2,2);
        ad = ad./norm_;
        ad = unique(ad,"rows");
    
        D2 = size(ad,1);
        K = zeros((B-1)*D2,D);
        cp1 = (cond(:,2)'-x_n)./ad;
        cpp1 = cp1;
        cpp1(cp1 < 0) = max(cp1,[],'all')+1; 
        cp2 = (cond(:,1)'-x_n)./ad;
        cpp2 = cp2;
        cpp2(cp2 < 0) = max(cp2,[],'all')+1;

        b_pp = min(cpp1,cpp2);
        b_pp = min(b_pp,[],2);
        %tempb_pp(~nonz) = max(tempb_pp(nonz),[],'ComparisonMethod','abs')+1;
        b_pp = repelem(b_pp,B-1,1);
        
        cm1 = -cp1;
        cm2 = -cp2;
        cmm1 = cm1;
        cmm2 = cm2;
        cmm1(cm1 < 0) = max(cm1,[],'all')+1;
        cmm2(cm2 < 0) = max(cm2,[],'all')+1;

        b_mm = min(cmm1,cmm2);
        b_mm = min(b_mm,[],2);
        %tempb_mm(~nonz) = max(tempb_mm(nonz),[],'ComparisonMethod','abs')+1;
        b_mm = repelem(b_mm,B-1,1);
    
        for i = 1:D2
            K(1+(i-1)*(B-1):(B-1)+(i-1)*(B-1),:)=ad(i,:).*k;
        end
    
        Yd_pp = x_n+(K.*b_pp);
        Yd_mm = x_n-(K.*b_mm);

        Ydd = union(Yd_pp,Yd_mm,"rows");
        Y = union(Yd,Ydd,'rows');
    else
        Y = Yd;
    end
    E = uniquetol(union(Y,E_old,"rows"),1e-2,'ByRow',true);
end

