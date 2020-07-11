function [A]=rupRelation(M,epsilon,model)

switch model
    case 1 %'null'
        A = M*0;
        
    case 2 %'wc1994'
        a=1.00; b =4.00; sigma = 0.25;
        logA   = (M-b+sigma*epsilon)/a;
        A      = 10.^logA;
        
    case 3 %'ellsworth'
        a=1.00; b =4.20; sigma = 0.12;
        logA   = (M-b+sigma*epsilon)/a;
        A      = 10.^logA;
        
    case 4% 'hanksbakun2001'
        a=1.00; b =3.98; sigma = 0.12; a2=4/3; b2= 3.09;
        log10A  = (M-b +sigma*epsilon)/a;
        log10A2 = (M-b2+sigma*epsilon)/a2; A2  = 10.^log10A2;
        A      = 10.^log10A;
        ind    = (A>468);
        A(ind) = A2(ind);
        
    case 5% 'somerville1999'
        a=1.00; b =3.95; sigma = 0.12;
        logA   = (a*M-b+sigma*epsilon)/a;
        A      = 10.^logA;
        
    case 6% 'wellscoppersmithr1994'
        a=0.90; b =4.33; sigma = 0.25;
        logA   = (M-b+sigma*epsilon)/a;
        A      = 10.^logA;
        
    case 7% 'wellscoppersmithss1994'
        a=1.02; b =3.94; sigma = 0.23;
        logA   = (M-b+sigma*epsilon)/a;
        A      = 10.^logA;
        
    case 8% 'strasser2010'
        a=0.846; b =4.441; sigma = 0.286;
        logA   = (M-b+sigma*epsilon)/a;
        A      = 10.^logA;
end
