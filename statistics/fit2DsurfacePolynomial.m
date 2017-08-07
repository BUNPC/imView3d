function [cAll] = fit2DsurfacePolynomial(y,x,z)

Zprofile = @(c,y,x) (c(1).*y.*y + c(2).*x.*x + c(3).*y.*x + c(4).*y + c(5).*x +c(6) );

Zprofile_error = @(c,y,x,z) (z - (c(1).*y.*y + c(2).*x.*x + c(3).*y.*x + c(4).*y + c(5).*x +c(6) ) ) ;

MyOptions = optimset('Display','off','LevenbergMarquardt','on');

c0 = [1 1 1 1 1 0];

[cAll,resnorm,resid,exitflag,output,lambda,jacobian_c] = lsqnonlin(Zprofile_error, c0, [],[],MyOptions,y,x,z);
    