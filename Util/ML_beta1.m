function [beta,  beta_prime] = ML_beta1(U,par)
% Function for the opening and closing rates in the Morris-Lecar equation.
% Alpha and beta are written to have separate parameters. 


 %alpha = 0.5.*( cosh( (U - par.u3a)./(2.*par.u4a)) + tanh( (U - par.u3a)./par.u4a ).*cosh( (U - par.u3a)./(2.*par.u4a) ) );
 beta = 0.5.*( 1 - tanh( (U - par.u3b)./par.u4b)).*cosh( (U - par.u3b)./(2.*par.u4b) );

if nargout > 2
    
    %alpha_prime = ( sinh((U - par.u3a)./(2.*par.u4a)).*( 0.25 - 0.25.*tanh((par.u3a - U)./(par.u4a))) + 0.5.*cosh((par.u3a - U)./(2.*par.u4a)).*(sech((par.u3a - U)./(par.u4a))).^2 )./par.u4a;
   beta_prime = ( sinh((U - par.u3b)./(2.*par.u4b)).*( 0.25.*tanh((par.u3b - U)./par.u4b) + 0.25) - 0.5.*cosh((par.u3b - U)./(2.*par.u4b)).*(sech((par.u3b - U)./par.u4b)).^2 )./par.u4b;
    
end