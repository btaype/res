function m=biseccion(f,a,b,t)
  m1=(a + b)/2;
  fprintf('m = %.6f\n',m1 );
  while abs(b-a)>t
    if f(a)*f(m1)<0;
      b=m1;
    else
      a=m1;
    endif
    m1=(a + b)/2;
    fprintf('m = %.6f\n',m1 );
  endwhile
  m=m1;
endfunction



function m = puntofalso(f, a, b, t)
  m1 = (a*f(b) - b*f(a)) / (f(b) - f(a));
  fprintf('m = %.6f\n', m1);
  while abs(b - a) > t
    if f(a) * f(m1) < 0
      b = m1;
    else
      a = m1;
    endif
    m1 = (a*f(b) - b*f(a)) / (f(b) - f(a));
    fprintf('m = %.6f\n', m1);
  endwhile
  m = m1;
endfunction




function m = newton(x,g,t)
  x1=g(x)
  fprintf('m = %.6f\n',x1);
  while abs(x1-x)>=t
    x=x1;
    x1=g(x);
    fprintf('m = %.6f\n',x1);
  endwhile
  m=x1;
endfunction







function m=secante(g,x,t)
  h=t;
  x1=g(x,h);
  fprintf('m = %.8f\n',x1);
  while abs(x1-x)>=t
    x=x1;
    x1=g(x,h);
    fprintf('m = %.8f\n',x1);
  endwhile
  m=x1;
endfunction




function inversa1=inversa_matriz(a)
  n=size(a,1);
  ext=eye(n);
  a=[a,ext];
  for i=1:n
    a(i,:)= a(i,:)/a(i,i);
    for j=1:n
       if i~=j
        a(j,:) = a(j,:) - a(i,:)*a(j,i);
       endif
    endfor
  endfor
  inversa1=a(:,n+1:end)
endfunction


function inversa1=inversa_matriz(a,amplia)
  n=size(a,1);
  ext=amplia;
  a=[a,ext];
  for i=1:n
    a(i,:)= a(i,:)/a(i,i);
    for j=1:n
       if i~=j
        a(j,:) = a(j,:) - a(i,:)*a(j,i);
       endif
    endfor
  endfor

  inversa1=a;
endfunction


function [u, j] = halla(a, b, n)
    for i = 1:n-1
        for j = i+1:n
            m = a(j,i) / a(i,i);
            a(j,i) = 0;
            for k = i+1:n
                a(j,k) = a(j,k) - m * a(i,k);
            endfor
            b(j) = b(j) - m * b(i);
        endfor
    endfor
    u = a;
    j = b;
endfunction

function x=Rtri(A,b)
  n=size(A,1);
  x=zeros(n,1);
  x(n)=b(n)/A(n,n);
  for k=(n-1):-1:1
    s=0;
    for j=(k+1):n
      s=s+A(k,j)*x(j);
    endfor
    x(k)=(b(k)-s)/A(k,k);
  endfor
endfunction
