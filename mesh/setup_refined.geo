Mesh.Algorithm = 1;
// Mesh.CharacteristicLengthExtendFromBoundary = 1;
// Mesh.CharacteristicLengthMax = 3;
// Mesh.CharacteristicLengthMin = 0.2;

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

Function CheeseHole 
  // taken from t5.geo
  p1 = newp; Point(p1) = {x,  y,  z,  lcar3} ;
  p2 = newp; Point(p2) = {x+r,y,  z,  lcar3} ;
  p3 = newp; Point(p3) = {x,  y+r,z,  lcar3} ;
  p4 = newp; Point(p4) = {x,  y,  z+r,lcar3} ;
  p5 = newp; Point(p5) = {x-r,y,  z,  lcar3} ;
  p6 = newp; Point(p6) = {x,  y-r,z,  lcar3} ;
  p7 = newp; Point(p7) = {x,  y,  z-r,lcar3} ;

  c1 = newreg; Circle(c1) = {p2,p1,p7};
  c2 = newreg; Circle(c2) = {p7,p1,p5};
  c3 = newreg; Circle(c3) = {p5,p1,p4};
  c4 = newreg; Circle(c4) = {p4,p1,p2};
  c5 = newreg; Circle(c5) = {p2,p1,p3};
  c6 = newreg; Circle(c6) = {p3,p1,p5};
  c7 = newreg; Circle(c7) = {p5,p1,p6};
  c8 = newreg; Circle(c8) = {p6,p1,p2};
  c9 = newreg; Circle(c9) = {p7,p1,p3};
  c10 = newreg; Circle(c10) = {p3,p1,p4};
  c11 = newreg; Circle(c11) = {p4,p1,p6};
  c12 = newreg; Circle(c12) = {p6,p1,p7};

  l1 = newreg; Line Loop(l1) = {c5,c10,c4};   Ruled Surface(newreg) = {l1};
  l2 = newreg; Line Loop(l2) = {c9,-c5,c1};   Ruled Surface(newreg) = {l2};
  l3 = newreg; Line Loop(l3) = {c12,-c8,-c1}; Ruled Surface(newreg) = {l3};
  l4 = newreg; Line Loop(l4) = {c8,-c4,c11};  Ruled Surface(newreg) = {l4};
  l5 = newreg; Line Loop(l5) = {-c10,c6,c3};  Ruled Surface(newreg) = {l5};
  l6 = newreg; Line Loop(l6) = {-c11,-c3,c7}; Ruled Surface(newreg) = {l6};
  l7 = newreg; Line Loop(l7) = {-c2,-c7,-c12};Ruled Surface(newreg) = {l7};
  l8 = newreg; Line Loop(l8) = {-c6,-c9,c2};  Ruled Surface(newreg) = {l8};

  theloopss[t] = newreg ; 
  Surface Loop(theloopss[t]) = {l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1};
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
RefineMesh;

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


lcar3 = 0.01;
x = 15.374994; y= 15.199986; z=-4.125000; r=0.17; t=1;
Call CheeseHole;
Surface{l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1} In Volume {enclosingbrain};

// pp = newp;
// Point(pp) = {19.0, 32.0, -2.0, lc2};
// Point{pp} In Volume {enclosingbrain};
// pp2 = newp;
// Point(pp2) = {19.0, 32.0, -1.0, lc2};
// lx = newl;
// Line(lx) = {pp,pp2};	
// Line{lx} In Volume {enclosingbrain};

// Refinement!!![5.196, 22.913, -4.9957]
// Merge "case6_post.msh";
// Surface{1} In Volume {enclosingbrain};

