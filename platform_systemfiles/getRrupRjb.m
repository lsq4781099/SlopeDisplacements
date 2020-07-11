function[SC]=getRrupRjb(SC)

Ztor    = SC.Ztor;
dip     = SC.dip;
W       = SC.W;
x       = SC.Rx;
SC.Rrup = x*0;
SC.Rjb  = x*0;

for i=1:length(x)
    xsource = [0 W*cosd(-dip)];
    ysource = [0 W*sind(-dip)]-Ztor(i);
    Wp = Ztor(i)/sind(dip);
    xp = Ztor(i)/tand(dip);
    xq = Ztor(i)*tand(dip);
    xr = (Wp+W)/cosd(dip)-xp;
    
    xrb1 = 0;
    xrb2 = W*cosd(dip);
    
    % Rjb
    if x(i)<0 || x(i)> xrb2
        SC.Rjb(i)=min(abs(x(i)-[xrb1,xrb2]));
    else
        SC.Rjb(i) =0;
    end
    
    % Rrup
    if x(i)<=xq
        prup   = [xsource(1),ysource(1)];
    elseif and(xq<x(i),x(i)<xr)
        wx = cosd(dip)*(x(i)+xp);
        prup = [wx*cosd(dip)-xp,-wx*sind(dip)];
    else
        prup = [xsource(2),ysource(2)];
    end
    SC.Rrup(i) = norm(prup-[x(i),0]);
end

SC.Rhyp = SC.Rrup;
SC.Zhyp = Ztor+W/2*sind(dip);


