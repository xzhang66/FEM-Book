function [ex,S,M,w]=exact_m
L=8; a1=4; p=10; m=10; s=-10; 
E=1e4; I=1; 

i=[0:.01:L]; 
ex=i; 

ind=1; 
for x=i     
    if x<a1
        w(ind)=(-x^4/24+14*x^3/3-71*x^2)/(E*I); 
        M(ind)=-x^2/2+28*x-142; 
        S(ind)=x-28;
    else 
        w(ind)=(-x^4/24+3*x^3-51*x^2-80*x+106.6667)/(E*I); 
        M(ind)=-x^2/2+18*x-102; 
        S(ind)=x-18;
    end
    
    ind=ind+1;                         
end 