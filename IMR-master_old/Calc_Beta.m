ofunction [ B ] = Calc_Beta( chi,w )
%Pass chi number and returns appropiate Beta factor 

X=sqrt(w*1i/chi);

B=real(((X*coth(X)-1)^-1-3/X^2)^-1);


end

