function y=eerfcs(x,scal)
%the complementary error function

t=1./(1.+0.5.*abs(x));
pa=[0.17087277,-0.82215223,1.48851587,-1.13520398,0.27886807,...
    -0.18628806,0.09678418,0.37409196,1.00002368,-1.2655123];
ans1=polyval(pa,t);


y=zeros(size(x));
jj=numel(x);
for k=1:jj
    if(x(k)>=0.0)
        %erf=1-tau, erfc=tau
        y(k)=t(k).*exp(-x(k)^2+ans1(k)+scal(k));
    else
        %erf=tau-1, erfc=2-tau
        y(k)=2.*exp(scal(k))-t(k).*exp(-x(k)^2+ans1(k)+scal(k));
    end
    
end
end