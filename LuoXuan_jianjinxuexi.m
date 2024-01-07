function [x] = LuoXuan_jianjinxuexi(x,fit,sortIndex,pop,fobj)
l=2*rand-1;
b=1;
%   此处显示详细说明
[A,B]=sort(fit);
for i=1:pop
    if (B(i)>sortIndex(i))
        q=find(B==(sortIndex(i)));
           zz=x(i,:)+x(q(1),:).*exp(-b*l*(B(i)-sortIndex(i))).*cos(2*pi*l);
     if(fobj(zz)<fobj(x(B(i),:)))
         x(B(i),:)=zz;
     end
    end
end
end

