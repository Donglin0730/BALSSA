
function [lb,ub,dim,fobj] = Get_Functions_details(F)


switch F
    case 'F1'
        fobj = @F1;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F2'
        fobj = @F2;
        lb=-10;
        ub=10;
        dim=30;
        
    case 'F3'
        fobj = @F3;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F4'
        fobj = @F4;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F5'
        fobj = @F5;
        lb=-30;
        ub=30;
        dim=30;
        
    case 'F6'
        fobj = @F6;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F7'
        fobj = @F7;
        lb=-1.28;
        ub=1.28;
        dim=30;
        
    case 'F8'
        fobj = @F8;
        lb=-500;
        ub=500;
        dim=30;
        
    case 'F9'
        fobj = @F9;
        lb=-5.12;
        ub=5.12;
        dim=30;
        
    case 'F10'
        fobj = @F10;
        lb=-32;
        ub=32;
        dim=30;
        
    case 'F11'
        fobj = @F11;
        lb=-600;
        ub=600;
        dim=30;
        
    case 'F12'
        fobj = @F12;
        lb=-50;
        ub=50;
        dim=30;
        
    case 'F13'
        fobj = @F13;
        lb=-50;
        ub=50;
        dim=30;
        
    case 'F14'
        fobj = @F14;
        lb=-65.536;
        ub=65.536;
        dim=2;
        
    case 'F15'
        fobj = @F15;
        lb=-5;
        ub=5;
        dim=4;
        
    case 'F16'
        fobj = @F16;
        lb=-5;
        ub=5;
        dim=2;
        
    case 'F17'
        fobj = @F17;
        lb=[-5,0];
        ub=[10,15];
        dim=2;
        
    case 'F18'
        fobj = @F18;
        lb=-2;
        ub=2;
        dim=2;
        
    case 'F19'
        fobj = @F19;
        lb=0;
        ub=1;
        dim=3;
        
    case 'F20'
        fobj = @F20;
        lb=0;
        ub=1;
        dim=6;     
        
    case 'F21'
        fobj = @F21;
        lb=0;
        ub=10;
        dim=4;    
        
    case 'F22'
        fobj = @F22;
        lb=0;
        ub=10;
        dim=4;    
        
    case 'F23'
        fobj = @F23;
        lb=0;
        ub=10;
        dim=4;    
        
       case 'F24'
        fobj = @bukin6;
        lb=[-15,-5];
        ub=[-3,3];
        dim=2;    
        case 'F25'%-1
        fobj = @drop;
        lb=-5.12;
        ub=5.12;
        dim=2;    
        case 'F26'%-959.6407
        fobj = @egg;
        lb=-5.12;
        ub=5.12;
        dim=2; 
        case 'F27'%-19.2085
        fobj = @holder;
        lb=-10;
        ub=10;
        dim=2; 
         case 'F28'%0
        fobj = @levy13;
        lb=-10;
        ub=10;
        dim=2; 
        case 'F29'%0
        fobj = @schwef;
        lb=-500;
        ub=500;
        dim=30; 
        case 'F30'%-186.7309
        fobj = @shubert;
        lb=-10;
        ub=10;
        dim=2; 
         case 'F31'%0
        fobj = @zakharov;
        lb=-5;
        ub=10;
        dim=30; 
         case 'F32'%0
        fobj = @camel3;
        lb=-5;
        ub=5;
        dim=2; 
         case 'F33'%0
        fobj = @ dixonpr;
        lb=-10;
        ub=10;
        dim=30; 
         case 'F34'%-1
        fobj = @ easom;
        lb=-100;
        ub=100;
        dim=2;  
         case 'F35'%-9.66015
        fobj = @ michal;
        lb=0;
        ub=pi;
        dim=10;  
        case 'F36'%0
        fobj = @ beale;
        lb=-4.5;
        ub=4.5;
        dim=2;  
         case 'F37'%0
        fobj = @ colville;
        lb=-10;
        ub=10;
        dim=4;  
         case 'F38'%0
        fobj = @ powell;
        lb=-4;
        ub=5;
        dim=30;  
         case 'F39'%0
        fobj = @ schaffer2;
        lb=-100;
        ub=100;
        dim=2;  
         case 'F40'%0
        fobj = @ dejong5;
        lb=-65.536;
        ub=65.536;
        dim=2;
         case 'F41'%0
        fobj = @ goldpr;
        lb=-2;
        ub=2;
        dim=2;
         case'F42'%Bent Cigar Function
        fobj = @ bcf;
        lb=-10^10;
        ub=10^10;
        dim=10;
       case'F43'%Sum of Different Power Function
        fobj = @ sdp;
        lb=-10*30;
        ub=10*30;
        dim=30;  
        case'F44'% Zakharov Function
        fobj = @ zf;
        lb=-5^10;
        ub=10^10;
        dim=10;
        case'F45'% Non-continuous Rotated Rastrigin¡¯s Function
        fobj = @ nrrf;
        lb=-10^10;
        ub=10^10;
        dim=10;
         case'F46'% Levy Function 
        fobj = @ lf;
        lb=-10^30;
        ub=10^30;
        dim=30;
         case'F47'% Ackley¡¯s Function
        fobj = @ af;
        lb=-32^30;
        ub=32^30;
        dim=30;
        case'F48'% Griewank¡¯s Function
        fobj = @ gf;
        lb=-600;
        ub=600;
        dim=30;
         case'F49'%  Rosenbrock¡¯s Function
        fobj = @ rf;
        lb=-5.12^30;
        ub=5.12^30;
        dim=30;
end

end

% F1

function o = F1(x)
o=sum(x.^2);

end

% F2

function o = F2(x)
o=sum(abs(x))+prod(abs(x));
end

% F3

function o = F3(x)
dim=size(x,2);
o=0;
for i=1:dim
    o=o+sum(x(1:i))^2;
end
end

% F4

function o = F4(x)
o=max(abs(x));
end

% F5

function o = F5(x)
dim=size(x,2);
o=sum(100*(x(2:dim)-(x(1:dim-1).^2)).^2+(x(1:dim-1)-1).^2);
end

% F6

function o = F6(x)
o=sum(abs((x+.5)).^2);
end

% F7

function o = F7(x)
dim=size(x,2);
o=sum([1:dim].*(x.^4))+rand;
end

% F8

function o = F8(x)
o=sum(-x.*sin(sqrt(abs(x))));
end

% F9

function o = F9(x)
dim=size(x,2);
o=sum(x.^2-10*cos(2*pi.*x))+10*dim;
end

% F10

function o = F10(x)
dim=size(x,2);
o=-20*exp(-.2*sqrt(sum(x.^2)/dim))-exp(sum(cos(2*pi.*x))/dim)+20+exp(1);
end

% F11

function o = F11(x)
dim=size(x,2);
o=sum(x.^2)/4000-prod(cos(x./sqrt([1:dim])))+1;
end

% F12

function o = F12(x)
dim=size(x,2);
o=(pi/dim)*(10*((sin(pi*(1+(x(1)+1)/4)))^2)+sum((((x(1:dim-1)+1)./4).^2).*...
(1+10.*((sin(pi.*(1+(x(2:dim)+1)./4)))).^2))+((x(dim)+1)/4)^2)+sum(Ufun(x,10,100,4));
end

% F13

function o = F13(x)
dim=size(x,2);
o=.1*((sin(3*pi*x(1)))^2+sum((x(1:dim-1)-1).^2.*(1+(sin(3.*pi.*x(2:dim))).^2))+...
((x(dim)-1)^2)*(1+(sin(2*pi*x(dim)))^2))+sum(Ufun(x,5,100,4));
end

% F14

function o = F14(x)
aS=[-32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32;,...
-32 -32 -32 -32 -32 -16 -16 -16 -16 -16 0 0 0 0 0 16 16 16 16 16 32 32 32 32 32];

for j=1:25
    bS(j)=sum((x'-aS(:,j)).^6);
end
o=(1/500+sum(1./([1:25]+bS))).^(-1);
end

% F15

function o = F15(x)
aK=[.1957 .1947 .1735 .16 .0844 .0627 .0456 .0342 .0323 .0235 .0246];
bK=[.25 .5 1 2 4 6 8 10 12 14 16];bK=1./bK;
o=sum((aK-((x(1).*(bK.^2+x(2).*bK))./(bK.^2+x(3).*bK+x(4)))).^2);
end

% F16

function o = F16(x)
o=4*(x(1)^2)-2.1*(x(1)^4)+(x(1)^6)/3+x(1)*x(2)-4*(x(2)^2)+4*(x(2)^4);
end

% F17

function o = F17(x)
o=(x(2)-(x(1)^2)*5.1/(4*(pi^2))+5/pi*x(1)-6)^2+10*(1-1/(8*pi))*cos(x(1))+10;
end

% F18

function o = F18(x)
o=(1+(x(1)+x(2)+1)^2*(19-14*x(1)+3*(x(1)^2)-14*x(2)+6*x(1)*x(2)+3*x(2)^2))*...
    (30+(2*x(1)-3*x(2))^2*(18-32*x(1)+12*(x(1)^2)+48*x(2)-36*x(1)*x(2)+27*(x(2)^2)));
end

% F19

function o = F19(x)
aH=[3 10 30;.1 10 35;3 10 30;.1 10 35];cH=[1 1.2 3 3.2];
pH=[.3689 .117 .2673;.4699 .4387 .747;.1091 .8732 .5547;.03815 .5743 .8828];
o=0;
for i=1:4
    o=o-cH(i)*exp(-(sum(aH(i,:).*((x-pH(i,:)).^2))));
end
end

% F20

function o = F20(x)
aH=[10 3 17 3.5 1.7 8;.05 10 17 .1 8 14;3 3.5 1.7 10 17 8;17 8 .05 10 .1 14];
cH=[1 1.2 3 3.2];
pH=[.1312 .1696 .5569 .0124 .8283 .5886;.2329 .4135 .8307 .3736 .1004 .9991;...
.2348 .1415 .3522 .2883 .3047 .6650;.4047 .8828 .8732 .5743 .1091 .0381];
o=0;
for i=1:4
    o=o-cH(i)*exp(-(sum(aH(i,:).*((x-pH(i,:)).^2))));
end
end

% F21

function o = F21(x)
aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];

o=0;
for i=1:5
    o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
end
end

% F22

function o = F22(x)
aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];

o=0;
for i=1:7
    o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
end
end

% F23

function o = F23(x)
aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];

o=0;
for i=1:10
    o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
end
end

function o=Ufun(x,a,k,m)
o=k.*((x-a).^m).*(x>a)+k.*((-x-a).^m).*(x<(-a));   
end
function [y] = bukin6(xx)
x1 = xx(1);
x2 = xx(2);

term1 = 100 * sqrt(abs(x2 - 0.01*x1^2));
term2 = 0.01 * abs(x1+10);

y = term1 + term2;

end
function [y] = drop(xx)
x1 = xx(1);
x2 = xx(2);


frac1 = 1 + cos(12*sqrt(x1^2+x2^2));
frac2 = 0.5*(x1^2+x2^2) + 2;

y = -frac1/frac2;

end
function [y] = egg(xx)
x1 = xx(1);
x2 = xx(2);

term1 = -(x2+47) * sin(sqrt(abs(x2+x1/2+47)));
term2 = -x1 * sin(sqrt(abs(x1-(x2+47))));

y = term1 + term2;

end
function [y] = holder(xx)
x1 = xx(1);
x2 = xx(2);

fact1 = sin(x1)*cos(x2);
fact2 = exp(abs(1 - sqrt(x1^2+x2^2)/pi));

y = -abs(fact1*fact2);

end
function [y] = levy13(xx)
x1 = xx(1);
x2 = xx(2);

term1 = (sin(3*pi*x1))^2;
term2 = (x1-1)^2 * (1+(sin(3*pi*x2))^2);
term3 = (x2-1)^2 * (1+(sin(2*pi*x2))^2);

y = term1 + term2 + term3;

end
function [y] = schwef(xx)
dim = length(xx);
sum = 0;
for ii = 1:dim 
	xi = xx(ii);
	sum = sum + xi*sin(sqrt(abs(xi)));
end

y = 418.9829*dim - sum;

end
function [y] = shubert(xx)
x1 = xx(1);
x2 = xx(2);
sum1 = 0;
sum2 = 0;

for ii = 1:5
	new1 = ii * cos((ii+1)*x1+ii);
	new2 = ii * cos((ii+1)*x2+ii);
	sum1 = sum1 + new1;
	sum2 = sum2 + new2;
end

y = sum1 * sum2;

end
function [y] = zakharov(xx)
dim = length(xx);
sum1 = 0;
sum2 = 0;

for ii = 1:dim
	xi = xx(ii);
	sum1 = sum1 + xi^2;
	sum2 = sum2 + 0.5*ii*xi;
end

y = sum1 + sum2^2 + sum2^4;

end
function [y] = camel3(xx)
x1 = xx(1);
x2 = xx(2);

term1 = 2*x1^2;
term2 = -1.05*x1^4;
term3 = x1^6 / 6;
term4 = x1*x2;
term5 = x2^2;

y = term1 + term2 + term3 + term4 + term5;

end
function [y] = dixonpr(xx)
x1 = xx(1);
dim = length(xx);
term1 = (x1-1)^2;

sum = 0;
for ii = 2:dim
	xi = xx(ii);
	xold = xx(ii-1);
	new = ii * (2*xi^2 - xold)^2;
	sum = sum + new;
end

y = term1 + sum;

end
function [y] = easom(xx)
x1 = xx(1);
x2 = xx(2);

fact1 = -cos(x1)*cos(x2);
fact2 = exp(-(x1-pi)^2-(x2-pi)^2);

y = fact1*fact2;

end
function [y] = michal(xx, m)
if (nargin == 1)
    m = 10;
end

dim = length(xx);
sum = 0;

for ii = 1:dim
	xi = xx(ii);
	new = sin(xi) * (sin(ii*xi^2/pi))^(2*m);
	sum  = sum + new;
end

y = -sum;

end
function [y] = beale(xx)
x1 = xx(1);
x2 = xx(2);

term1 = (1.5 - x1 + x1*x2)^2;
term2 = (2.25 - x1 + x1*x2^2)^2;
term3 = (2.625 - x1 + x1*x2^3)^2;

y = term1 + term2 + term3;

end
function [y] = colville(xx)
x1 = xx(1);
x2 = xx(2);
x3 = xx(3);
x4 = xx(4);

term1 = 100 * (x1^2-x2)^2;
term2 = (x1-1)^2;
term3 = (x3-1)^2;
term4 = 90 * (x3^2-x4)^2;
term5 = 10.1 * ((x2-1)^2 + (x4-1)^2);
term6 = 19.8*(x2-1)*(x4-1);

y = term1 + term2 + term3 + term4 + term5 + term6;

end
function [y] = powell(xx)
dim = length(xx);
sum = 0;

for ii = 1:(dim/4)
	term1 = (xx(4*ii-3) + 10*xx(4*ii-2))^2;
	term2 = 5 * (xx(4*ii-1) - xx(4*ii))^2;
	term3 = (xx(4*ii-2) - 2*xx(4*ii-1))^4;
	term4 = 10 * (xx(4*ii-3) - xx(4*ii))^4;
	sum = sum + term1 + term2 + term3 + term4;
end

y = sum;

end
function [y] = schaffer2(xx)
x1 = xx(1);
x2 = xx(2);

fact1 = (sin(x1^2-x2^2))^2 - 0.5;
fact2 = (1 + 0.001*(x1^2+x2^2))^2;

y = 0.5 + fact1/fact2;

end
function [y] = dejong5(xx)
x1 = xx(1);
x2 = xx(2);
sum = 0;

A = zeros(2, 25);
a = [-32, -16, 0, 16, 32];
A(1, :) = repmat(a, 1, 5);
ar = repmat(a, 5, 1);
ar = ar(:)';
A(2, :) = ar;

for ii = 1:25
    a1i = A(1, ii);
    a2i = A(2, ii);
    term1 = ii;
    term2 = (x1 - a1i)^6;
    term3 = (x2 - a2i)^6;
    new = 1 / (term1+term2+term3);
    sum = sum + new;
end

y = 1 / (0.002 + sum);

end
function [y] = goldpr(xx)
x1 = xx(1);
x2 = xx(2);

fact1a = (x1 + x2 + 1)^2;
fact1b = 19 - 14*x1 + 3*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2;
fact1 = 1 + fact1a*fact1b;

fact2a = (2*x1 - 3*x2)^2;
fact2b = 18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2;
fact2 = 30 + fact2a*fact2b;

y = fact1*fact2;

end
function [z] = bcf(u)
dim=size(u,2);
z=u(1)^2+10^6*(sum(u(2:dim).^2));
end
function [z] = sdp(u)
dim=size(u);
sum_1=0;
for i=1:dim
  sum_1=sum_1+abs(u(i))^(i+1);
end
z=sum_1;
end
function [z] = zf(u)
dim=size(u,2);
sum_1=0;
for i=1:dim
    sum_1=sum_1+0.5*i*u(i);
end
z=sum(u.^2)+sum_1^2+sum_1^4;
end
function [z] = nrrf(u)
dim=size(u,2);
sum_1=0;
for i=1:dim
    if abs(u(i))<1/2
        y(i)=u(i);
    else 
        y(i)=(round(2*u(i)))/2;
    end
    sum_1=sum_1+y(i)^2-10*cos(2*pi*y(i))+10;
end
z=sum_1;
end
function [z] = lf(u)
dim=size(u,2);
w=1+(u-1)/4;
sum_1=sum(((w(1:dim-1)-1).^2).*(1+(10*((sin(pi*w(1:dim-1)+1)).^2))));
sum_2=(sin(pi*w(1)))^2;
sum_3=((w(dim)-1)^2)*(1+10*((sin(pi*w(dim)))^2));
z=sum_1+sum_2+sum_3;
end
function [z] = af(u)
numv=size(u,2);
x=u(1:numv); 
z=-20*exp(-0.2*sqrt(sum(x.^2)/numv))-exp(sum(cos(2*pi*x))/numv)+20+exp(1);
end
function [z] = gf(u)
dim=size(u);
sum_1=0;
for i=1:dim
    sum_1=sum_1+u(i)^2;
end
z1=1/4000*sum_1;
z2=1;
for i=1:dim
    z2=z2*cos(u(i)/sqrt(i));
end
z=z1-z2+1;
end
function [z] = rf(u)
dim=size(u,2);
z=sum(100*(u(2:dim)-(u(1:dim-1).^2)).^2+(u(1:dim-1)-1).^2);
end