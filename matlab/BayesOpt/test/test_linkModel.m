kp = linspace(0,0.04,1000);
Fr.u = 'w(1)';
Fr.y = 'r';
Glaser.u = {'w(2)','u(1)'};
Glaser.y = 'y(1)';
sum{1} = sumblk('e(1) = r - y(1)');
C1 = pid(30,0)/sys.k_phi;
C1.u = 'e(1)';
C1.y = 'u(1)';
Gg1 = connect(Glaser,C1,Fr,sum{:},'w','y');
for i = 1:length(kp)
    C = pid(kp(i),0)/sys_link.k_phi;
    C.u = 'e';
    C.y = 'u';
    Glink.u = {'y(1)','w(3)','u'};
    Glink.y = {'l','y(2)'};
    sums{1} = sumblk('e = y(1) - l');
    sums{2} = sumblk('z = r-y(2)');
    Gg = connect(Gg1,C,Fr,Glink,sums{:},'w','z');
%     poles = pole(Gg(1));
%     figure(1)
%     hold on
%     plot(real(poles),imag(poles),'b.')
%     hold off
    
    figure(2)
    hold on
    h2Norm = norm(Gg(2),2);
    if h2Norm ~= inf
        plot(kp(i),norm(Gg(2),2),'r.')
    end
    hold off
end
%rltool(Glink(2,3),C)
%%
Glink_pade = sys_link.Gpade;
Glink = sys_link.G;
C = pid(0,0);
C.u = 'e';
C.y = 'u';
sum = sumblk('e = r - l');
Gg = connect(C,Glink,sum,{'r'},{'y'});

step(Gg,6e-3)
