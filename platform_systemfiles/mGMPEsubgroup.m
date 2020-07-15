function[str,C]=mGMPEsubgroup(indx)

methods  = pshatoolbox_methods(1);
gmpetype = {methods.type}';
C0       = strcmp(gmpetype,'regular');
methods  = methods(C0);
%% Copy and past from GMM_summary.xls
IMlist = [
1	0	0	0	0	0	0	1	0	0	0      % Youngs et al. 1997           
1	0	0	0	0	0	0	1	0	0	0      % Atkinson & Boore, 2003       
1	0	0	0	0	0	0	1	0	0	0      % Zhao et al. 2006             
1	0	0	0	0	0	0	1	0	0	0      % McVerry et al. 2006          
1	0	0	0	0	0	0	1	0	0	0      % Boroschek et al. 2012        
1	0	0	0	0	0	0	1	0	0	0      % Abrahamson et al. 2016       
1	0	0	0	0	0	0	1	0	0	0      % Abrahamson et al. 2018  
1	1	0	0	0	0	0	1	0	0	0      % Kuehn et al. 2020
1	1	0	0	0	0	0	1	0	0	0      % Parker et al. 2020
1	0	0	0	0	0	0	1	0	0	0      % Arteta et al. 2018           
1	0	0	0	0	0	0	1	0	0	0      % Idini et al. 2016            
1	0	0	0	0	0	0	1	0	0	0      % Montalva et al. 2017         
1	0	0	0	0	0	0	1	0	0	0      % Montalva et al. 2017 (HQ) 
1	0	0	0	0	0	0	1	0	0	0      % Montalva et al. 2018
1	1	0	0	0	0	0	1	0	0	0      % SIBER-RISK 2019              
1	1	0	0	0	0	0	1	0	0	0      % Garcia et al. 2005           
1	0	0	0	0	0	0	1	0	0	0      % Jaimes et al. 2006           
1	1	0	0	0	0	0	1	0	0	0      % Jaimes et al. 2015           
1	1	0	0	0	0	0	1	0	0	0      % Jaimes et al. 2016           
1	1	0	0	0	0	0	1	0	0	0      % Garcia-Soto Jaimes 2017      
0	0	0	0	0	0	0	0	0	0	1      % Garcia-Soto Jaimes 2017 (HV) 
1	0	0	0	0	0	0	1	0	0	0      % Bernal et al. 2014           
1	0	0	0	0	0	0	1	0	0	0      % Sadigh et al. 1997           
1	0	0	0	0	0	0	1	0	0	0      % Idriss 2008 - NGA            
1	1	0	0	0	0	0	1	0	0	0      % Chiou Youngs 2008 - NGA      
1	1	0	0	0	0	0	1	0	0	0      % Boore Atkinson 2008 - NGA    
1	1	1	0	0	0	0	1	0	0	0      % Campbell Bozorgnia 2008 - NGA
1	1	0	0	0	0	0	1	0	0	0      % Abrahamson Silva 2008 - NGA  
1	0	0	0	0	0	0	1	0	0	0      % Abrahamson Silva 1997 (Horz) 
1	0	0	0	0	0	0	1	0	0	0      % Idriss 2014 - NGAW2          
1	1	0	0	0	0	0	1	0	0	0      % CY 2014 - NGAW2              
1	1	0	0	0	0	0	1	0	0	0      % CB 2014 - NGAW2              
1	1	0	0	0	0	0	1	0	0	0      % BSSA 2014 - NGAW2            
1	1	0	0	0	0	0	1	0	0	0      % ASK 2014 - NGAW2             
1	0	0	0	0	0	0	0	0	1	0      % Akkar & Boomer 2007          
1	1	0	0	0	0	0	1	0	0	0      % Akkar & Boomer 2010          
1	0	0	0	0	0	0	1	0	0	0      % Arroyo et al. 2010           
1	1	0	0	0	0	0	1	0	0	0      % Bindi et al. 2011            
1	1	0	0	0	0	0	1	0	0	0      % Kanno et al. 2006            
1	1	0	0	0	0	0	1	0	1	0      % Cauzzi et al., 2015          
0	0	0	0	1	0	0	0	0	0	0      % Du & Wang, 2012              
0	0	0	0	1	1	0	0	0	0	0      % Foulser-Piggott, Goda 2015   
0	0	0	0	0	1	0	0	0	0	0      % Travasarou, Bray, Abra  2003 
0	0	0	0	1	1	1	0	0	0	0      % Bullock et al, 2017          
0	0	0	0	1	0	0	0	0	0	0      % Campbell,Bozorgnia 2010      
0	0	0	0	1	0	0	0	0	0	0      % Campbell,Bozorgnia 2011      
0	0	0	0	1	1	0	0	0	0	0      % Campbell,Bozorgnia 2019      
0	0	0	0	1	0	0	0	0	0	0      % Kramer & Mitchell, 2006  
1	0	0	0	0	0	0	1	0	0	0      % Median PCE BCHydro
1	0	0	0	0	0	0	1	0	0	0];    % Median PCE NGA
                                                                             
mech = [
1	1	0      % Youngs et al. 1997           
1	1	0      % Atkinson & Boore, 2003       
1	1	1      % Zhao et al. 2006             
1	1	1      % McVerry et al. 2006          
1	0	0      % Boroschek et al. 2012        
1	1	0      % Abrahamson et al. 2016       
1	1	0      % Abrahamson et al. 2018     
1   1   0      % Kuehn et al. 2020
1   1   0      % Parker et al. 2020
0	1	0      % Arteta et al. 2018           
1	1	0      % Idini et al. 2016            
1	1	0      % Montalva et al. 2017         
1	1	0      % Montalva et al. 2017 (HQ)   
1	1	0      % Montalva et al. 2018
1	1	0      % SIBER-RISK 2019              
0	1	0      % Garcia et al. 2005           
1	0	0      % Jaimes et al. 2006           
0	1	0      % Jaimes et al. 2015           
0	0	1      % Jaimes et al. 2016           
1	0	0      % Garcia-Soto Jaimes 2017      
1	0	0      % Garcia-Soto Jaimes 2017 (HV) 
1	1	0      % Bernal et al. 2014           
0	0	1      % Sadigh et al. 1997           
0	0	1      % Idriss 2008 - NGA            
0	0	1      % Chiou Youngs 2008 - NGA      
0	0	1      % Boore Atkinson 2008 - NGA    
0	0	1      % Campbell Bozorgnia 2008 - NGA
0	0	1      % Abrahamson Silva 2008 - NGA  
0	0	1      % Abrahamson Silva 1997 (Horz) 
0	0	1      % Idriss 2014 - NGAW2          
0	0	1      % CY 2014 - NGAW2              
0	0	1      % CB 2014 - NGAW2              
0	0	1      % BSSA 2014 - NGAW2            
0	0	1      % ASK 2014 - NGAW2             
0	0	1      % Akkar & Boomer 2007          
0	0	1      % Akkar & Boomer 2010          
1	0	0      % Arroyo et al. 2010           
0	0	1      % Bindi et al. 2011            
1	1	0      % Kanno et al. 2006            
0	0	1      % Cauzzi et al., 2015          
0	0	1      % Du & Wang, 2012              
0	0	1      % Foulser-Piggott, Goda 2015   
0	0	1      % Travasarou, Bray, Abra  2003 
1	1	1      % Bullock et al, 2017          
0	0	1      % Campbell,Bozorgnia 2010      
0	0	1      % Campbell,Bozorgnia 2011      
0	0	1      % Campbell,Bozorgnia 2019      
0	0	1      % Kramer & Mitchell, 2006 
1	1	0      % Median PCE BCHydro
0	0   1];    % Median PCE NGA
                                              

%% Selection
%       1     2     3     4     5          6     7    8     9    10   11   12    13        14           15           16  
%fn = {'All','PGA','PGV','PGD','Duration','CAV','AI','VGI','SA','SV','SD','H/V','NGA-2008','NGA-West2','Subduction','Crustal'};
N = size(IMlist,1);
switch indx
    case 1 , C = ones(N,1);
    case 2 , C = IMlist(:,1);
    case 3 , C = IMlist(:,2);
    case 4 , C = IMlist(:,3);
    case 5 , C = IMlist(:,5);
    case 6 , C = IMlist(:,6);
    case 7 , C = IMlist(:,7);
    case 8 , C = IMlist(:,8);
    case 9 , C = IMlist(:,9);
    case 10, C = IMlist(:,10);
    case 11, C = IMlist(:,11);
    case 12, C = zeros(N,1);[~,B]=intersect({methods.str},{'I2008','CY2008','BA2008','CB2008','AS2008'},'stable'); C(B)=1; 
    case 13, C = zeros(N,1);[~,B]=intersect({methods.str},{'I2014','CY2014','BSSA2014','CB2014','ASK2014'},'stable'); C(B)=1; 
    case 14, C = or(mech(:,1),mech(:,2));
    case 15, C = mech(:,3);
end
str={methods(C==1).label}';
