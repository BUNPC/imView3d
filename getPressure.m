function pres = getPressure( diam, vesselType )

% polyval values from polyfit of 4th order
% from grabLipowsky
% note that pCV is for capillaries and veins >=8 um
% note that pCA is for arterioles with diam > 8um
%      and that you must solve with (16-diam)
diam = max(diam,8);
if vesselType==1
    pCA = [-2.0559e-7 1.6279e-4 -8.9432e-4 -1.4090 45.0691];
    pres = polyval( pCA, (16-diam) );
else
    pCV = [6.7517e-6 -0.001 0.0562 -1.5670 43.8772];
    pres = polyval( pCV, diam );
end

