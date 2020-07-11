function [sys] =process_model(sys,opt)

n=max(sys.branch);

% process pointsources
for i=1:n(1)
    sys.point1(i)  = process_point1  (sys.src1(i) , opt.ellipsoid);
    sys.line1(i)   = process_line1   (sys.src2(i) , opt.ellipsoid);
    sys.area1(i)   = process_area1   (sys.src3(i) , opt.ellipsoid);
    sys.area2(i)   = process_area2   (sys.src4(i) , opt.ellipsoid);
    sys.area3(i)   = process_area3   (sys.src5(i) , opt.ellipsoid);
    sys.area4(i)   = process_area4   (sys.src6(i) , opt.ellipsoid);
    sys.volume1(i) = process_volume1 (sys.src7(i) , opt.ellipsoid);
end

% process magnitude recurrence relations
for i=1:n(3)
    sys.delta(i)     = process_delta     (sys.delta(i));
    sys.magtable(i)  = process_magtable  (sys.magtable(i));
    sys.truncexp(i)  = process_truncexp  (sys.truncexp(i)  , opt.MagDiscrete);
    sys.truncnorm(i) = process_truncnorm (sys.truncnorm(i) , opt.MagDiscrete);
    sys.yc1985(i)    = process_yc1985    (sys.yc1985(i)    , opt.MagDiscrete);
end