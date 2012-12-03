
Hmax = 20;

bPositive = cell(Hmax,1);
bNegative = cell(Hmax,1);

for H = 10:Hmax
    
    t = 1 : 1 : H ;
    b = zeros(1,H);
    d = zeros(1,H);
    
    for h=1:H,
        a = 4*H - 3*h + 1;
        c = (H -h + 1)^2;
        m = 4 * ((a+h)^2);
        p = 12 * ((2*c) + (h*a));
        r = 2 * (a+h);
        b(h) = ( r + ((m - p)^0.5)) / 6;
        d(h) = ( r - ((m - p)^0.5)) / 6;
        
    end
    
    bPositive{H,1} = (b);
    bNegative{H,1} = (d);
end

% hold all;
for H = 1:Hmax
    %     plot(bPositive{H,1},'*-');
%     plot(bNegative{H,1},'o-');
end

pkmin = 0;
pkmax = 20;
x = randi([-15 20],10,1);

tprev = 0;

for i = 1:100
    
    pk = 0.5 * (pkmin + pkmax);
    
    cvx_begin
    
    variable t;
    dual variables y1 y2 y3
    variable w(10,1);
    
    maximize(t)
    
    subject to
    
    x' * (w) + pk * sum((w)) >= t : y1
    
    w <= 1 : y2;
    
    x .* w >= x .* (1 - w) : y3;
    
    
    cvx_end
    
   
    if tprev < t
        tprev = t;
        pkmin = pk;
    else
        pkmax = pk;
    end
    
    if abs(pkmax - pkmin) < 0.002
        break;
    end
    
end



