function q = invKin(p,R,model)

l1 = model.l1;
l2 = model.l2;

maxNm = l1+l2;
norm(p);
percent = norm(p) / maxNm;
clampLen = getIntVal(norm(p), maxNm*(1-1e-6), .99);
p = p/norm(p)*clampLen ;

nmp = norm(p)^2;
%nmp

%other = (nmp - l1^2 - l2^2) / (2*l1*l2)


t3 = acos( (nmp - l1^2 - l2^2) / (2 * l1 * l2));
t1 = atan2(p(2), -p(3));
a = l1 + l2*cos(t3);
b = -l2*sin(t3);
c = p(1);


t2 = 2*atan2(sqrt(a^2 + b^2-c^2) + b, a+c);

q = [t1 t2 t3]';



end


function y = getIntVal(x,mx,cutoff)
    if x < cutoff*mx
        y = x;
    elseif x < mx
        left = x - mx*cutoff;
        norm = left/ ( (1-cutoff)*mx);  
        y = cutoff*mx + (1-cutoff)*mx*smoothRamp( norm);
    else
        y = mx;    
    end
end

function y = smoothRamp(x)
    y = -126*x^11 + 686*x^10 - 1505*x^9 + 1665*x^8 - 930*x^7 + 210*x^6 + x;
end
