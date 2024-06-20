function raiz = secante(f,x0,tol)
k=0;
h=0.000001;
x1=f(x0,h);
while abs(x1-x0)>=tol
  fprintf("|      %4.8f      |      %4.8f     \n",k,x1);
  k+=1;
  x0=x1;
  x1=f(x0,h);
endwhile
raiz=x1
endfunction




function p=lagranger(x,y)
  n=length(x);
  f=poly(x(2:n));
  p=y(1)*f/polyval(f,x(1));
  for k=2:n
    f=poly([x(1:k-1),x(k+1:n)]);
    p=p+(y(k)*f/polyval(f,x(k)));

  endfor


endfunction







function r=newtondos(x,y)
r=[y(1)];
n=length(x);
for k=2:n
  p=poly(x(1:k-1));
  h= (y(k)-polyval(r,x(k)))/polyval(p,x(k));
  r=[0 r]+ h*p;

endfor

endfunction

function I=trapecios(f,a,b,n)
h=(b-a)/n;
I=0;
for k=1:n
    I=I+f(a+(k-1)*h)+f(a+k*h);
endfor
I=I*h/2;
endfunction


function A = romberg(f,a,b,error)
  r(1,1) = trapecios(f,a,b,1);
  r(2,1) = trapecios(f,a,b,2);
  r(2,2) = (4*r(2,1)-r(1,1))/3;
  k = 2;
  while abs(r(k,k) - r(k,k-1))>=error
    k = k+1;
    r(k,1) = trapecios(f,a,b,2^(k-1));
    for j = 2:k
      r(k,j) = (4^(j-1)*r(k,j-1)-r(k-1,j-1))/(4^(j-1)-1);
    endfor
  endwhile
  A = r(k,k);
endfunction




function [t,y]=EULER(f,t0,y0,T,p)
h=T/p;
t=zeros(p+1,1);
y=zeros(p+1,1);
t(1)=t0;
y(1)=y0;
for k=2:1:p+1
  t(k)=t(k-1)+h;
  y(k)=y(k-1)+h*f(t(k-1),y(k-1));
endfor
plot(t,y,'*');


function [t,y]=TAYLOR(f,ft,fy,t0,y0,T,p)
h=T/p;
t=zeros(p+1,1);
y=zeros(p+1,1);
t(1)=t0;
y(1)=y0;
for k=2:1:p+1
  t(k)=t(k-1)+h;
  y(k)=y(k-1)+h*f(t(k-1),y(k-1))+(ft(t(k-1),y(k-1))+fy(t(k-1),y(k-1))*f(t(k-1),y(k-1)))*h^2/2;
endfor;
plot(t,y,'*');




function [t,y]=RK4(f,t0,y0,T,p)
h=T/p;
t=zeros(p+1,1);
y=zeros(p+1,1);
t(1)=t0;
y(1)=y0;
for k=2:1:p+1
  t(k)=t(k-1)+h;
  k1=h*f(t(k-1),y(k-1));
  k2=h*f(t(k-1)+h/2,y(k-1)+k1/2);
  k3=h*f(t(k-1)+h/2,y(k-1)+k2/2);
  k4=h*f(t(k-1)+h,y(k-1)+k3);
  y(k)=y(k-1)+k1/6+2*k2/6+2*k3/6+k4/6;
endfor
plot(t,y,'*');
