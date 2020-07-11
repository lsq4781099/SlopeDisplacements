function[ind]=selectsource(MaxDistance,xyz,source)
% SELECT SOURCES THAT ARE CLOSER TO sys.MaxDistance KM FROM THE SITE with
% coordinates xyz

% Analysis Parameters
Nsources = length(source);
ind     = true(1,Nsources);
for i=1:Nsources
    delta = source(i).hypm-xyz;
    dist = min(delta.^2*[1;1;1]);
    ind(i)=(dist<MaxDistance^2);
end


