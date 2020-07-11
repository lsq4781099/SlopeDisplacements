function[point1]=process_obj1(point1,ellipsoid)

% num = mechanism gmmptr strike dip style rupArea
np     = size(point1.txt,1);
normal = zeros(np,3);
for j=1:np
    strike = point1.num(j,6);
    dip = point1.num(j,7);
    dsv = dsv_point(point1.vert(j,:),strike,dip,ellipsoid);
    normal(j,:)= dsv(:,3)';
end

point1.normal = normal;

