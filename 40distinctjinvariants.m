load "jInvariant.m";
/*
Use parametrisation of moduli space of abelian surfaces over Q
with full level-3 structure
*/

F := Rationals();
R<y0,y1,y2,y3,y4> := PolynomialRing(F,5);
f := y0*(y0^3+y1^3+y2^3+y3^3+y4^3) + 3*y1*y2*y3*y4;
B := Scheme(ProjectiveSpace(R), [f]);


t1 := 2;
t2 := -1;
t3 := 3;

curve:=[t1^3-3*t1^2*t3-3*t1*t2^2-3*t1*t2*t3-t2^3-1,
-t1^3+3*t1^2*t3-3*t1*t3^2+t2^3+1,
-t1^4+t1^3*t2+3*t1^3*t3-3*t1^2*t2*t3-3*t1^2*t3^2-2*t1*t2^3-3*t1*t2^2*t3+t1-t2^4-t2,
-t1^4+4*t1^3*t3+3*t1^2*t2^2+3*t1^2*t2*t3-3*t1^2*t3^2+t1*t2^3-3*t1*t2^2*t3-3*t1*t2*t3^2+t1-t2^3*t3-t3,
-t1^4-t1^3*t2+2*t1^3*t3+3*t1^2*t2*t3+t1*t2^3+3*t1*t2^2*t3+t1+t2^4+t2^3*t3+t2+t3];





a1 := curve[2]/curve[1];
a2 := curve[3]/curve[1];
a3 := curve[4]/curve[1];
a4 := curve[5]/curve[1];



//Ka<a1,a2,a3,a4>:=FunctionField(k,4);


KaX<X>:=PolynomialRing(F);
//the defining homogeneous equation of the model of the Burkhardt Quartic we use:
B:=y0^4+y0*y1^3+y0*y2^3+y0*y3^3+3*y1*y2*y3*y4+y0*y4^3;

//Below we give 4 lists [H,lambda,G], where lambda is a rational function in (a1,...,a4)
//and H,G are polynomials in X with coefficients that are rational functions in (a1,...,a4).
//The expression G^2+4*lambda*H^3 yields the same sextic in X for each triple [H,lambda,G].
//When (1:a1:a2:a3:a4) is a point on B that does not lie in a j-plane and has a4 != 0
//then this corresponds to the 3-torsion point corresponding to the j-plane y[0]=y[i]=0
//(for i=1,...,4, in the order given)
//(The article describes H,G as homogeneous forms in (x,z), of degrees 2 and 3 respectively.
//Here we set (x,z)=(X,1) to get a more compact representation)
HLGs:=[[(a1^3*a3*a4^3+a1^2*a2^2*a4^5+a1^2*a2^2*a4^2+2*a1*a2*a3^2*a4^4-a1*a2*a3^2*a4-a2^3*a3+
a3^4*a4^3)*X^2+(a1^4*a4^4+2*a1^2*a2*a3*a4^5-a1^2*a2*a3*a4^2+2*a1*a2^3*a4^4+a1*a2^3*a4+
2*a1*a3^3*a4^4+a1*a3^3*a4+2*a2^2*a3^2*a4^3+2*a2^2*a3^2)*X+a1^3*a2*a4^3+a1^2*a3^2*a4^5+
a1^2*a3^2*a4^2+2*a1*a2^2*a3*a4^4-a1*a2^2*a3*a4+a2^4*a4^3-a2*a3^3,
(-a4^3-1)/(a1^6*a4^6-6*a1^4*a2*a3*a4^4-2*a1^3*a2^3*a4^3-2*a1^3*a3^3*a4^3+9*a1^2*a2^2*a3^2*a4^2+
6*a1*a2^4*a3*a4+6*a1*a2*a3^4*a4+a2^6+2*a2^3*a3^3+a3^6),
(a1^6*a4^6+3*a1^4*a2*a3*a4^7-3*a1^4*a2*a3*a4^4+2*a1^3*a2^3*a4^9+4*a1^3*a2^3*a4^6+3*a1^3*a3^3*a4^6
+a1^3*a3^3*a4^3+6*a1^2*a2^2*a3^2*a4^8+3*a1^2*a2^2*a3^2*a4^5+6*a1^2*a2^2*a3^2*a4^2-
3*a1*a2^4*a3*a4^4+3*a1*a2^4*a3*a4+6*a1*a2*a3^4*a4^7+a2^6-3*a2^3*a3^3*a4^3-a2^3*a3^3+
2*a3^6*a4^6+a3^6*a4^3)/(a1^3*a4^3-3*a1*a2*a3*a4-a2^3-a3^3)*X^3+(3*a1^5*a2*a4^8+
3*a1^5*a2*a4^5+6*a1^4*a3^2*a4^7+6*a1^4*a3^2*a4^4+6*a1^3*a2^2*a3*a4^9+6*a1^3*a2^2*a3*a4^6+
6*a1^2*a2^4*a4^8+9*a1^2*a2^4*a4^5+3*a1^2*a2^4*a4^2+12*a1^2*a2*a3^3*a4^8+3*a1^2*a2*a3^3*a4^5-
9*a1^2*a2*a3^3*a4^2+12*a1*a2^3*a3^2*a4^7+9*a1*a2^3*a3^2*a4^4-3*a1*a2^3*a3^2*a4+6*a1*a3^5*a4^7+
6*a1*a3^5*a4^4-3*a2^5*a3*a4^3-3*a2^5*a3+6*a2^2*a3^4*a4^6+9*a2^2*a3^4*a4^3+
3*a2^2*a3^4)/(a1^3*a4^3-3*a1*a2*a3*a4-a2^3-a3^3)*X^2+(3*a1^5*a3*a4^8+3*a1^5*a3*a4^5+
6*a1^4*a2^2*a4^7+6*a1^4*a2^2*a4^4+6*a1^3*a2*a3^2*a4^9+6*a1^3*a2*a3^2*a4^6+12*a1^2*a2^3*a3*a4^8
+3*a1^2*a2^3*a3*a4^5-9*a1^2*a2^3*a3*a4^2+6*a1^2*a3^4*a4^8+9*a1^2*a3^4*a4^5+3*a1^2*a3^4*a4^2+
6*a1*a2^5*a4^7+6*a1*a2^5*a4^4+12*a1*a2^2*a3^3*a4^7+9*a1*a2^2*a3^3*a4^4-3*a1*a2^2*a3^3*a4+
6*a2^4*a3^2*a4^6+9*a2^4*a3^2*a4^3+3*a2^4*a3^2-3*a2*a3^5*a4^3-3*a2*a3^5)/(a1^3*a4^3-
3*a1*a2*a3*a4-a2^3-a3^3)*X+(a1^6*a4^6+3*a1^4*a2*a3*a4^7-3*a1^4*a2*a3*a4^4+3*a1^3*a2^3*a4^6
+a1^3*a2^3*a4^3+2*a1^3*a3^3*a4^9+4*a1^3*a3^3*a4^6+6*a1^2*a2^2*a3^2*a4^8+3*a1^2*a2^2*a3^2*a4^5
+6*a1^2*a2^2*a3^2*a4^2+6*a1*a2^4*a3*a4^7-3*a1*a2*a3^4*a4^4+3*a1*a2*a3^4*a4+2*a2^6*a4^6+
a2^6*a4^3-3*a2^3*a3^3*a4^3-a2^3*a3^3+a3^6)/(a1^3*a4^3-3*a1*a2*a3*a4-a2^3-a3^3)
],[
a1*a4*X^2+a2*X-a3,
-a1^3*a4^9-a1^3*a4^6+3*a1*a2*a3*a4^7+3*a1*a2*a3*a4^4+a2^3*a4^6+a2^3*a4^3+a3^3*a4^6+a3^3*a4^3,
(2*a1^3*a4^6+a1^3*a4^3-3*a1*a2*a3*a4^4+a2^3-a3^3*a4^3)*X^3+(3*a1^2*a2*a4^5+3*a1^2*a2*a4^2-
3*a2^2*a3*a4^3-3*a2^2*a3)*X^2+(-3*a1^2*a3*a4^5-3*a1^2*a3*a4^2+3*a2*a3^2*a4^3+3*a2*a3^2)*X-
a1^3*a4^3-3*a1*a2*a3*a4^4-a2^3*a4^3-2*a3^3*a4^3-a3^3
],[
a2*X^2-a3*X-a1*a4,
a1^3*a4^9+a1^3*a4^6-3*a1*a2*a3*a4^7-3*a1*a2*a3*a4^4-a2^3*a4^6-a2^3*a4^3-a3^3*a4^6-a3^3*a4^3,
(a1^3*a4^3+3*a1*a2*a3*a4^4+2*a2^3*a4^3+a2^3+a3^3*a4^3)*X^3+(3*a1^2*a2*a4^5+3*a1^2*a2*a4^2-
3*a2^2*a3*a4^3-3*a2^2*a3)*X^2+(-3*a1^2*a3*a4^5-3*a1^2*a3*a4^2+3*a2*a3^2*a4^3+3*a2*a3^2)*X-
2*a1^3*a4^6-a1^3*a4^3+3*a1*a2*a3*a4^4+a2^3*a4^3-a3^3
],[
(a1*a2^2*a4^3+a1*a2^2)*X^2+(a1^3*a4^2+a1*a2*a3*a4^3-2*a1*a2*a3+a2^3*a4^2+a3^3*a4^2)*X+
a1*a3^2*a4^3+a1*a3^2,
(-a1^3*a4^3+3*a1*a2*a3*a4+a2^3+a3^3)/(a1^6+6*a1^4*a2*a3*a4+2*a1^3*a2^3+2*a1^3*a3^3+
9*a1^2*a2^2*a3^2*a4^2+6*a1*a2^4*a3*a4+6*a1*a2*a3^4*a4+a2^6+2*a2^3*a3^3+a3^6),
(a1^6*a4^3+6*a1^4*a2*a3*a4^4+2*a1^3*a2^3*a4^6+5*a1^3*a2^3*a4^3+a1^3*a2^3+2*a1^3*a3^3*a4^3+
9*a1^2*a2^2*a3^2*a4^5+3*a1*a2^4*a3*a4^4-3*a1*a2^4*a3*a4+6*a1*a2*a3^4*a4^4-a2^6+a2^3*a3^3*a4^3
-a2^3*a3^3+a3^6*a4^3)/(a1^3+3*a1*a2*a3*a4+a2^3+a3^3)*X^3+(3*a1^5*a2*a4^5+3*a1^5*a2*a4^2+
3*a1^3*a2^2*a3*a4^6-3*a1^3*a2^2*a3+3*a1^2*a2^4*a4^5+3*a1^2*a2^4*a4^2+3*a1^2*a2*a3^3*a4^5+
3*a1^2*a2*a3^3*a4^2+9*a1*a2^3*a3^2*a4^4+9*a1*a2^3*a3^2*a4+3*a2^5*a3*a4^3+3*a2^5*a3+
3*a2^2*a3^4*a4^3+3*a2^2*a3^4)/(a1^3+3*a1*a2*a3*a4+a2^3+a3^3)*X^2+(-3*a1^5*a3*a4^5-
3*a1^5*a3*a4^2-3*a1^3*a2*a3^2*a4^6+3*a1^3*a2*a3^2-3*a1^2*a2^3*a3*a4^5-3*a1^2*a2^3*a3*a4^2-
3*a1^2*a3^4*a4^5-3*a1^2*a3^4*a4^2-9*a1*a2^2*a3^3*a4^4-9*a1*a2^2*a3^3*a4-3*a2^4*a3^2*a4^3-
3*a2^4*a3^2-3*a2*a3^5*a4^3-3*a2*a3^5)/(a1^3+3*a1*a2*a3*a4+a2^3+a3^3)*X+(-a1^6*a4^3-
6*a1^4*a2*a3*a4^4-2*a1^3*a2^3*a4^3-2*a1^3*a3^3*a4^6-5*a1^3*a3^3*a4^3-a1^3*a3^3-
9*a1^2*a2^2*a3^2*a4^5-6*a1*a2^4*a3*a4^4-3*a1*a2*a3^4*a4^4+3*a1*a2*a3^4*a4-a2^6*a4^3-
a2^3*a3^3*a4^3+a2^3*a3^3+a3^6)/(a1^3+3*a1*a2*a3*a4+a2^3+a3^3)]];

//Below is a triple (H,lambda,G) such that G^2+4*lambda*H^3 is -3*F, where F is the sextic
//defined by HLGs above. This triple corresponds to the 3-torsion points that the j-plane
//y0+...+y4=y0+y4=0 marks if (1:a1:a2:a3:a4) is a point on the Burkhardt quartic. Note that
//the identity of the sextics holds regardless of whether (1:a1:a2:a3:a4) satisfy the Burkhardt
//relation. We do need the Burkhardt relation to get the other cyclic order 3 subgroups defined
//over the base field.
HLGdual:=[(a1^2*a4^2+a1*a2*a4^3-a1*a2*a4^2-a1*a2*a4-a1*a3*a4^2+a2^2+a2*a3*a4+a3^2*a4^2)*X^2+(-a1^2*a4^3+
a1^2*a4^2+a1*a2*a4^3+a1*a2*a4+a1*a3*a4^3+a1*a3*a4+a2^2*a4^2-a2^2*a4-2*a2*a3+a3^2*a4^2-
a3^2*a4)*X+a1^2*a4^2-a1*a2*a4^2+a1*a3*a4^3-a1*a3*a4^2-a1*a3*a4+a2^2*a4^2+a2*a3*a4+a3^2,
(-3*a1^2*a4^4+3*a1^2*a4^3-3*a1^2*a4^2-3*a1*a2*a4^3+3*a1*a2*a4^2-3*a1*a2*a4-3*a1*a3*a4^3+
3*a1*a3*a4^2-3*a1*a3*a4-3*a2^2*a4^2+3*a2^2*a4-3*a2^2+3*a2*a3*a4^2-3*a2*a3*a4+3*a2*a3-
3*a3^2*a4^2+3*a3^2*a4-3*a3^2)/(a1^2*a4^4+2*a1^2*a4^3+a1^2*a4^2-2*a1*a2*a4^3-4*a1*a2*a4^2-
2*a1*a2*a4-2*a1*a3*a4^3-4*a1*a3*a4^2-2*a1*a3*a4+a2^2*a4^2+2*a2^2*a4+a2^2+2*a2*a3*a4^2+
4*a2*a3*a4+2*a2*a3+a3^2*a4^2+2*a3^2*a4+a3^2),
(-3*a1^4*a4^5+3*a1^4*a4^4-6*a1^3*a2*a4^6+6*a1^3*a2*a4^5-3*a1^3*a2*a4^4-3*a1^3*a2*a4^3+
6*a1^3*a3*a4^5-3*a1^3*a3*a4^4+3*a1^3*a3*a4^3+6*a1^2*a2^2*a4^6+6*a1^2*a2^2*a4^4+6*a1^2*a2^2*a4^2+
3*a1^2*a2*a3*a4^6-3*a1^2*a2*a3*a4^5+6*a1^2*a2*a3*a4^4+6*a1^2*a2*a3*a4^3-6*a1^2*a2*a3*a4^2-
6*a1^2*a3^2*a4^5+6*a1^2*a3^2*a4^4-6*a1^2*a3^2*a4^3-6*a1*a2^3*a4^4+6*a1*a2^3*a4^3-3*a1*a2^3*a4^2-
3*a1*a2^3*a4-3*a1*a2^2*a3*a4^5+3*a1*a2^2*a3*a4^4-6*a1*a2^2*a3*a4^3-6*a1*a2^2*a3*a4^2+
6*a1*a2^2*a3*a4+9*a1*a2*a3^2*a4^5-3*a1*a2*a3^2*a4^4-6*a1*a2*a3^2*a4^3+6*a1*a2*a3^2*a4^2+
3*a1*a3^3*a4^5-3*a1*a3^3*a4^4+6*a1*a3^3*a4^3-3*a2^4*a4+3*a2^4-6*a2^3*a3*a4^2+3*a2^3*a3*a4-
3*a2^3*a3-6*a2^2*a3^2*a4^3+6*a2^2*a3^2*a4^2-6*a2^2*a3^2*a4-3*a2*a3^3*a4^4+3*a2*a3^3*a4^3-
6*a2*a3^3*a4^2+3*a3^4*a4^4-3*a3^4*a4^3)/(a1*a4^2+a1*a4-a2*a4-a2-a3*a4-a3)*X^3+(6*a1^4*a4^6
-6*a1^4*a4^5+6*a1^4*a4^4+3*a1^3*a2*a4^7-9*a1^3*a2*a4^6+12*a1^3*a2*a4^5-9*a1^3*a2*a4^4+
3*a1^3*a2*a4^3-6*a1^3*a3*a4^6+12*a1^3*a3*a4^5-12*a1^3*a3*a4^4+6*a1^3*a3*a4^3+3*a1^2*a2^2*a4^6-
15*a1^2*a2^2*a4^5+18*a1^2*a2^2*a4^4-15*a1^2*a2^2*a4^3+3*a1^2*a2^2*a4^2+9*a1^2*a2*a3*a4^6-
9*a1^2*a2*a3*a4^5+9*a1^2*a2*a3*a4^3-9*a1^2*a2*a3*a4^2+6*a1^2*a3^2*a4^6-12*a1^2*a3^2*a4^5+
18*a1^2*a3^2*a4^4-12*a1^2*a3^2*a4^3+6*a1^2*a3^2*a4^2+12*a1*a2^3*a4^5-6*a1*a2^3*a4^4+
12*a1*a2^3*a4^3+6*a1*a2^3*a4+3*a1*a2^2*a3*a4^5+3*a1*a2^2*a3*a4^4+3*a1*a2^2*a3*a4^2+
3*a1*a2^2*a3*a4+6*a1*a2*a3^2*a4^5-12*a1*a2*a3^2*a4^4+6*a1*a2*a3^2*a4^2-12*a1*a2*a3^2*a4+
6*a1*a3^3*a4^5-12*a1*a3^3*a4^4+12*a1*a3^3*a4^3-6*a1*a3^3*a4^2-6*a2^4*a4^3+6*a2^4*a4^2-6*a2^4*a4
-3*a2^3*a3*a4^4+3*a2^3*a3*a4^3-12*a2^3*a3*a4^2+9*a2^3*a3*a4-9*a2^3*a3+9*a2^2*a3^2*a4^4-
9*a2^2*a3^2*a4^3+18*a2^2*a3^2*a4^2-9*a2^2*a3^2*a4+9*a2^2*a3^2+12*a2*a3^3*a4^3-12*a2*a3^3*a4^2+
12*a2*a3^3*a4+6*a3^4*a4^4-6*a3^4*a4^3+6*a3^4*a4^2)/(a1*a4^2+a1*a4-a2*a4-a2-a3*a4-a3)*X^2+
(6*a1^4*a4^6-6*a1^4*a4^5+6*a1^4*a4^4-6*a1^3*a2*a4^6+12*a1^3*a2*a4^5-12*a1^3*a2*a4^4+
6*a1^3*a2*a4^3+3*a1^3*a3*a4^7-9*a1^3*a3*a4^6+12*a1^3*a3*a4^5-9*a1^3*a3*a4^4+3*a1^3*a3*a4^3+
6*a1^2*a2^2*a4^6-12*a1^2*a2^2*a4^5+18*a1^2*a2^2*a4^4-12*a1^2*a2^2*a4^3+6*a1^2*a2^2*a4^2+
9*a1^2*a2*a3*a4^6-9*a1^2*a2*a3*a4^5+9*a1^2*a2*a3*a4^3-9*a1^2*a2*a3*a4^2+3*a1^2*a3^2*a4^6-
15*a1^2*a3^2*a4^5+18*a1^2*a3^2*a4^4-15*a1^2*a3^2*a4^3+3*a1^2*a3^2*a4^2+6*a1*a2^3*a4^5-
12*a1*a2^3*a4^4+12*a1*a2^3*a4^3-6*a1*a2^3*a4^2+6*a1*a2^2*a3*a4^5-12*a1*a2^2*a3*a4^4+
6*a1*a2^2*a3*a4^2-12*a1*a2^2*a3*a4+3*a1*a2*a3^2*a4^5+3*a1*a2*a3^2*a4^4+3*a1*a2*a3^2*a4^2+
3*a1*a2*a3^2*a4+12*a1*a3^3*a4^5-6*a1*a3^3*a4^4+12*a1*a3^3*a4^3+6*a1*a3^3*a4+6*a2^4*a4^4-
6*a2^4*a4^3+6*a2^4*a4^2+12*a2^3*a3*a4^3-12*a2^3*a3*a4^2+12*a2^3*a3*a4+9*a2^2*a3^2*a4^4-
9*a2^2*a3^2*a4^3+18*a2^2*a3^2*a4^2-9*a2^2*a3^2*a4+9*a2^2*a3^2-3*a2*a3^3*a4^4+3*a2*a3^3*a4^3-
12*a2*a3^3*a4^2+9*a2*a3^3*a4-9*a2*a3^3-6*a3^4*a4^3+6*a3^4*a4^2-6*a3^4*a4)/(a1*a4^2+a1*a4-
a2*a4-a2-a3*a4-a3)*X+(-3*a1^4*a4^5+3*a1^4*a4^4+6*a1^3*a2*a4^5-3*a1^3*a2*a4^4+3*a1^3*a2*a4^3
-6*a1^3*a3*a4^6+6*a1^3*a3*a4^5-3*a1^3*a3*a4^4-3*a1^3*a3*a4^3-6*a1^2*a2^2*a4^5+6*a1^2*a2^2*a4^4-
6*a1^2*a2^2*a4^3+3*a1^2*a2*a3*a4^6-3*a1^2*a2*a3*a4^5+6*a1^2*a2*a3*a4^4+6*a1^2*a2*a3*a4^3-
6*a1^2*a2*a3*a4^2+6*a1^2*a3^2*a4^6+6*a1^2*a3^2*a4^4+6*a1^2*a3^2*a4^2+3*a1*a2^3*a4^5-
3*a1*a2^3*a4^4+6*a1*a2^3*a4^3+9*a1*a2^2*a3*a4^5-3*a1*a2^2*a3*a4^4-6*a1*a2^2*a3*a4^3+
6*a1*a2^2*a3*a4^2-3*a1*a2*a3^2*a4^5+3*a1*a2*a3^2*a4^4-6*a1*a2*a3^2*a4^3-6*a1*a2*a3^2*a4^2+
6*a1*a2*a3^2*a4-6*a1*a3^3*a4^4+6*a1*a3^3*a4^3-3*a1*a3^3*a4^2-3*a1*a3^3*a4+3*a2^4*a4^4-
3*a2^4*a4^3-3*a2^3*a3*a4^4+3*a2^3*a3*a4^3-6*a2^3*a3*a4^2-6*a2^2*a3^2*a4^3+6*a2^2*a3^2*a4^2-
6*a2^2*a3^2*a4-6*a2*a3^3*a4^2+3*a2*a3^3*a4-3*a2*a3^3-3*a3^4*a4+3*a3^4)/(a1*a4^2+a1*a4-a2*a4-
a2-a3*a4-a3)];

V:={c[3]^2+4*c[2]*c[1]^3:c in HLGs};
/*
We check that all 4 expressions of the form P^2 + 4 \lambda * Q^3 give the same polynomial
*/
assert #V eq 1;


f := Random(V);
C := HyperellipticCurve(f);

"Genus 2 curve:", C;

/*
We now compute all possible decompositions f = P^2-Q^3
to use the formulas in our paper
*/
T<b0,b1,b2,c0,c1,c2,c3> := PolynomialRing(F,7);
S<t> := PolynomialRing(T);
coefs := Coefficients((c0+c1*t+c2*t^2+c3*t^3)^2-(b0+b1*t+b2*t^2)^3);
coefs2 := Coefficients(f);

SchemeThreeTorsionPoints := Scheme(AffineSpace(T), [coefs[i]-coefs2[i] : i in [1..7]]);

ThreeTorsionPoints3, F3 := PointsOverSplittingField(SchemeThreeTorsionPoints);

F3 := NumberField(AbsolutePolynomial(F3));

SchemeThreeTorsionPoints3 := ChangeRing(SchemeThreeTorsionPoints, F3);
ThreeTorsionPoints3 := RationalPoints(SchemeThreeTorsionPoints3);

/*
It turns out that all the decompositions P^2-Q^3 are defined over the degree-6
number field F3. Next we compute the j-invariants of the corresponding elliptic curves.
*/
"Number field of definition", F3;

/*
We check that we have indeed found the desired decompositions
*/
R<x> := PolynomialRing(F3);
for tp in ThreeTorsionPoints3 do
	Q := R![tp[1], tp[2], tp[3]];
	P := R![tp[4], tp[5], tp[6], tp[7]];
	assert P^2-Q^3 eq f;
end for;

/*
Given a decomposition P^2-Q^3, the corresponding elliptic curve is
given by w^3-3Q(x)-w*2P(x)=0. We compute these elliptic curves and
their j-invariants.
*/


jInvs := {};


R<x,w> := PolynomialRing(F3, 2);
for tp in ThreeTorsionPoints3 do
	b0, b1, b2, c0, c1, c2, c3 := Explode([tp[i] : i in [1..7]]);
	Q := b0 + b1*x + b2*x^2;
	P := c0 + c1*x + c2*x^2 + c3*x^3;
	E := w^3 - 3*Q*w - 2*P;
	jInvs := jInvs join {jInvariantCubicInhomogeneous(E)};
end for;

/*
First of all, we check that we have obtained 40 distinct j-invariants.
*/
assert #jInvs eq 40;

/*
Finally, we prove that they correspond to geometrically non-isogenous
elliptic curves. We first need a helper function.
*/

/*
Given two j-invariants lying in a commong algebraic number field K,
returns true if the elliptic curves with j-invariants j1, j2 can
be shown to be geometrically non-isogenous
*/
function ProveGeometricallyNonisogenous(j1, j2, K)
	j1 := K!j1;
	j2 := K!j2;
	O := RingOfIntegers(K);
	for p in [5..100] do
		I := ideal<O | [p]>;
		P := Factorisation(I)[1][1];
		FP, proj := ResidueClassField(P);
		n1 := Numerator(j1);
		d1 := Denominator(j1);
		n2 := Numerator(j2);
		d2 := Denominator(j2);
		if proj(d1) ne 0 and proj(d2) ne 0 then
			redj1 := proj(n1) / proj(d1);
			redj2 := proj(n2) / proj(d2);
			E1 := EllipticCurveFromjInvariant(redj1);
			E2 := EllipticCurveFromjInvariant(redj2);
			/*
			If exactly one of the two elliptic curves is supersingular,
			then they cannot be geometrically isogenous
			*/
			if IsOrdinary(E1) xor IsOrdinary(E2) then
				return true;
			end if;
			/*
			If both are ordinary, we check whether they have the same
			endomorphism field, which is a (geometric) isogeny invariant
			*/
			if IsOrdinary(E1) and IsOrdinary(E2) then
				EndoField1 := QuadraticField( Trace(E1)^2 - 4*#FP );
				EndoField2 := QuadraticField( Trace(E2)^2 - 4*#FP );
				if not IsIsomorphic(EndoField1, EndoField2) then
					return true;
				end if;
			end if;
		end if;
	end for;
	return false;
end function;

/*
Test that the 40 elliptic curves are pairwise geometrically non-isogenous
*/
jInvs := [j : j in jInvs];
F3 := OptimisedRepresentation(F3);
for i in [1..40] do
for k in [(i+1)..40] do
	j1 := jInvs[i];
	j2 := jInvs[k];
	assert ProveGeometricallyNonisogenous(j1, j2, F3);
end for;
end for;

