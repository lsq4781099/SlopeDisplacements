function [lny,sigma,tau,phi] = KM06(T,Mw,Rrup,mechanism)

% Steven L. Kramer and Robert A. Mitchell (2006) Ground Motion Intensity 
% Measures for Liquefaction Hazard Evaluation. Earthquake Spectra: 
% May 2006, Vol. 22, No. 2, pp. 413-438.
% DOI: https://doi.org/10.1193/1.2194970

lny   = nan(size(Mw));
sigma = nan(size(Mw));
tau   = nan(size(Mw));
phi   = nan(size(Mw));

if T~=-4.1 % CAV5
    return
end

% model coefficients, 
%El valor absoluto de los coeficientes son los mismos del paper, pero deben
%haber algunos que lleven signos negativos puesto que si no ocurre esto, el
%CAV ser�a creciente en relaci�n con la distancia a la ruptura.

c1      =3.495;     %No puede ser negativo ya que sino parten muy abajo las curvas
c2      =2.764;
c3      =-8.539;    %Esta negativo, para que se parezca lo m�s posible al gr�fico del paper
c4      =-1.008;    %Debe ser negativo para que la pendiente de la recta sea negativa
h       =6.155;
f1      =-0.464;    %Lo puse negativo debido a que en la funcion de BTA03 es negativo y ademas adquiere similitud al gr�fico del paper si es negativo
f2      =0.165;


switch mechanism % SGS site category 
    case 'strike-slip',     FN = 0; FR = 0; 
    case 'normal',          FN = 1; FR = 0; 
    case 'reverse',         FN = 0; FR = 1; 
    case 'reverse-oblique', FN = 0; FR = 1;
end


lny   = c1 + c2*(Mw-6) + c3*log(Mw/6) + c4*log(sqrt(Rrup.^2+h.^2)) + f1*FN + f2*FR;
sigma = 0.708;


