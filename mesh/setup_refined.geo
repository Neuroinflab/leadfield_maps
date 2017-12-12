Mesh.Algorithm = 5;
Mesh.CharacteristicLengthExtendFromBoundary = 1;
Mesh.CharacteristicLengthMax = 3;
Mesh.CharacteristicLengthMin = 0.2;

Macro Cuboid
      p1 = newp; Point(p1) = {xmin,ymin,zmin,lc};
      p2 = newp; Point(p2) = {xmax,ymin,zmin,lc};
      p3 = newp; Point(p3) = {xmax,ymax,zmin,lc};
      p4 = newp; Point(p4) = {xmin,ymax,zmin,lc};
      l1 = newl; Line(l1) = {p4,p3};
      l2 = newl; Line(l2) = {p3,p2};
      l3 = newl; Line(l3) = {p2,p1};
      l4 = newl; Line(l4) = {p1,p4};
      l5 = newreg; Line Loop(l5) = {l2,l3,l4,l1};
      s6 = 1000; Plane Surface(s6) = {l5};
      tmp1[] = Extrude {0,0.0,zmax-zmin} {
             Surface{s6};
      };     
      theloops[t] = newreg;
      Surface Loop(theloops[t]) = {tmp1[0], tmp1[2], tmp1[3], tmp1[4], tmp1[5], s6};

      // psurfloops[t] = newreg;
      // Physical Surface(psurfloops[t]) = {tmp[0], tmp[2], tmp[3], tmp[4], tmp[5], s6};
      // pvolloops[t] = newreg;
      // Physical Volume(pvolloops[t]) = tmp1[1];
Return

// Large external volume
xmin = 0.0; xmax = 20.0;
ymin = -1.0; ymax = 35.0;
zmin = -20.0; zmax = 0.0;
// xmin = -15.0; xmax = 35.0;
// ymin = -16.0; ymax = 55.0;
// zmin = -40.0; zmax = 20.0;
lc = 3.; t = 1;
Call Cuboid;

extrnairVol = newreg;
Volume(extrnairVol) = theloops[1];

extrnsurfPhy = newreg;
Physical Surface (extrnsurfPhy) = {theloops[1]};


Merge "case3.msh";
// SEE t13.geo

brainsurfLoop = newreg;
Surface Loop(brainsurfLoop) = {2 : 19};
brainvol = newreg;
Volume(brainvol) = {brainsurfLoop};

brainsurfPhy = newreg;
Physical Surface(brainsurfPhy) = {2 : 19};

brainvolPhy = newreg;
Physical Volume(brainvolPhy) = brainvol;
groundsurfPhy = newreg;
Physical Surface(groundsurfPhy) = {2};

enclosingbrain = newreg;
Volume(enclosingbrain) = {theloops[1], brainsurfLoop};
enclosingbrainvol = newreg;
Physical Volume (enclosingbrainvol) = enclosingbrain;

// Refinement!!![5.196, 22.913, -4.9957]
// Merge "case6.msh";
// Surface{1} In Volume {extrnairVol};

// lc2 = 0.001;
// pp = newp;
// Point(pp) = {5.196, 22.913, -4.9957, lc2};
// Point{pp} In Volume {brainvol};
