function [x] = pinghengjinjiexuexi(x,worse,pop,fit,fobj)
m=rand;
n=rand;
bestX=[];
for j=1:pop
    [ ~, b] = min( fit );
    bestX = x( b, : );
end
x_avg=(1./pop).*(sum(bestX));
for i=1:pop
    if rand<0.5
    mm=x(i,:)+(m.*(worse-x_avg)+n.*(bestX-x_avg)).*log(1./rand);
    else
     mm=x(i,:)-(m.*(worse-x_avg)+n.*(bestX-x_avg)).*log(1./rand);
    end
if fobj(mm)<fobj(x(i,:))
    x(i,:)=mm;
end
end

end


