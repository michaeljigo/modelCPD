% calculate the eucledian distance between two points in 2D space.

function dist = calcDistance(v1, v2)

dist = sqrt((v2(1)-v1(1))^2+(v2(2)-v1(2))^2);

if isnan(dist)
    dist = inf;
end