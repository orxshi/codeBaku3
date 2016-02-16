//=INPUTS=====================================================

H = 20; // 30
D = 1;
//HH = 10;
mx = 200; // 150 // 100 // 125
my = 200; // 150 // 100 // 125
mz = 1;
bump = 200; // 200 // 1000

//============================================================

//HHH = (H - HH) / 2;

pLL = newp; Point(pLL) = {-H,-H,0}; // lower left
pUL = newp; Point(pUL) = {-H,H,0};  // upper left
pLR = newp; Point(pLR) = {H,-H,0};  // lower right
pUR = newp; Point(pUR) = {H,H,0};   // upper right

//pILL = newp; Point(pILL) = {-HHH,-HHH,0}; // inner lower left
//pIUL = newp; Point(pIUL) = {-HHH, HHH,0}; // inner upper left
//pILR = newp; Point(pILR) = { HHH,-HHH,0}; // inner lower right
//pIUR = newp; Point(pIUR) = { HHH, HHH,0}; // inner upper right

lL = newl; Line(lL) = {pLL,pUL};
lB = newl; Line(lB) = {pLR,pLL};
lR = newl; Line(lR) = {pUR,pLR};
lT = newl; Line(lT) = {pUL,pUR};

Transfinite Line{lL} = my Using Bump bump;
Transfinite Line{lR} = my Using Bump bump;
Transfinite Line{lB} = mx Using Bump bump;
Transfinite Line{lT} = mx Using Bump bump;

Line Loop(1) = {lL, lT, lR, lB};
sBase = news; Plane Surface(sBase) = {1};

Transfinite Surface{sBase}; Recombine Surface{sBase};

outB[] = Extrude {0,0,D} {Surface{sBase}; Layers{mz}; Recombine;};

Physical Surface(1) = {outB[0], sBase}; // laterals
Physical Surface(2) = {outB[2], outB[4], outB[3], outB[5]}; // sides
Physical Volume(1) = {outB[1]};


