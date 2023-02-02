function [Grid, new_order] = GridLebedev(Q)

degree = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, ...
    350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, ...
    3470, 3890, 4334, 4802, 5294, 5810];

order = [3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,35,41,47,53,59,65,71,77, ...
    83,89,95,101,107,113,119,125,131];

new_order = order((order - Q) > 0);
Q = degree((order - Q) > 0);
new_order = new_order(1);
Q = Q(1);

leb_tmp = getLebedevSphere(Q);

[Grid.Phi,Grid.Theta] = cart2sph(leb_tmp.x,leb_tmp.y,leb_tmp.z);

Grid.Theta = pi/2 - Grid.Theta;

Grid.Phi = pi + Grid.Phi;

Grid.Weight = leb_tmp.w;

end