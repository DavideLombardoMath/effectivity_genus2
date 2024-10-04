/*
=========== Two identically dependent sections ===========
*/
R<b3,c2> := PolynomialRing(Rationals(),2);
F<b3,c2> := FieldOfFractions(R);
b1 := 0;

A<x, w> := PolynomialRing(F, 2);
Q := c2 * x^2 + 3;
P := b3 * x^3 + b1 * x;

Factorisation(P^2-Q^3);

// f := P^2 - Q^3;
// C := HyperellipticCurve(f);

E := Curve( AffineSpace(A), [ w^3 - 3*Q*w - 2*P ] );
E;
Support(Flexes(E));



/*
=========== The sections are generically independent ===========
*/
A<x, w> := PolynomialRing(Rationals(), 2);
b1 := 0;
b3 := 1;
c2 := -1;
Q := c2 * x^2 + 7;
P := b3 * x^3 + b1 * x + 10;

E := Curve( AffineSpace(A), [ w^3 - 3*Q*w - 2*P ] );
p1 := E![0,-4];
p2 := E![0,-1];
p3 := E![0, 5];

EElliptic, map := EllipticCurve(ProjectiveClosure(E), p1);
S1 := map(p1);
S2 := map(p2);
S3 := map(p3);

sigma1 := S1-S2;
sigma2 := S2-S3;

R := RealField(30);
PairingMatrix := Matrix(R, 2,2, [HeightPairing(sigma1, sigma1), HeightPairing(sigma1, sigma2), HeightPairing(sigma2, sigma1), HeightPairing(sigma2, sigma2)]);
Determinant(PairingMatrix);





/*
=========== One section is 2-torsion, the other has infinite order ===========
*/

/*
Finding the equation
*/
R<b1,b3,c2> := PolynomialRing(Rationals(),3);
F<b1,b3,c2> := FieldOfFractions(R);

A<x, w> := PolynomialRing(F, 2);
Q := c2 * x^2 + 7;
P := b3 * x^3 + b1 * x + 10;

E := Curve( AffineSpace(A), [ w^3 - 3*Q*w - 2*P ] );
p1 := E![0,-4];
p2 := E![0,-1];
p3 := E![0, 5];

EElliptic, map := EllipticCurve(ProjectiveClosure(E), p1);
S1 := map(p1);
S2 := map(p2);
S3 := map(p3);

sigma1 := S1-S2;
sigma2 := S2-S3;
y1 := (2*sigma1)[2];
Factorisation(Denominator(y1));
y2 := (2*sigma2)[2];
Factorisation(Denominator(y2));

/*
Resulting equation: b1^3 - 567/4*b1*c2 + 2187/4*b3 = 0
*/


/*
Test
*/
R<b1,c2> := PolynomialRing(Rationals(),2);
F<b1,c2> := FieldOfFractions(R);

b3 := -4/2187*(b1^3 - 567/4*b1*c2);

A<x, w> := PolynomialRing(F, 2);
Q := c2 * x^2 + 7;
P := b3 * x^3 + b1 * x + 10;

R<t> := PolynomialRing(F);
h := hom<A -> R | [t, 0]>;
C := HyperellipticCurve(h(P^2-Q^3));
C;
G2Invariants(C);

E := Curve( AffineSpace(A), [ w^3 - 3*Q*w - 2*P ] );
p1 := E![0,-4];
p2 := E![0,-1];
p3 := E![0, 5];

EElliptic, map := EllipticCurve(ProjectiveClosure(E), p1);
S1 := map(p1);
S2 := map(p2);
S3 := map(p3);

sigma1 := S1-S2;
sigma2 := S2-S3;

2*sigma1;

/*
Second test: we take a further specialisation to check that sigma2
is generically non-torsion
*/

b1 := 1;
c2 := 1;

b3 := -4/2187*(b1^3 - 567/4*b1*c2);

A<x, w> := PolynomialRing(Rationals(), 2);
Q := c2 * x^2 + 7;
P := b3 * x^3 + b1 * x + 10;

R<t> := PolynomialRing(Rationals());
h := hom<A -> R | [t, 0]>;
C := HyperellipticCurve(h(P^2-Q^3));
C;
G2Invariants(C);

E := Curve( AffineSpace(A), [ w^3 - 3*Q*w - 2*P ] );
p1 := E![0,-4];
p2 := E![0,-1];
p3 := E![0, 5];

EElliptic, map := EllipticCurve(ProjectiveClosure(E), p1);
S1 := map(p1);
S2 := map(p2);
S3 := map(p3);

sigma1 := S1-S2;
sigma2 := S2-S3;

Order(sigma1);
Order(sigma2);











/*
=========== One section is 3-torsion, the other has infinite order ===========
*/
R<b1,b3,c2> := PolynomialRing(Rationals(),3);
F<b1,b3,c2> := FieldOfFractions(R);

A<x, w> := PolynomialRing(F, 2);
Q := c2 * x^2 + 7;
P := b3 * x^3 + b1 * x + 10;

E := Curve( AffineSpace(A), [ w^3 - 3*Q*w - 2*P ] );
p1 := E![0,-4];
p2 := E![0,-1];
p3 := E![0, 5];

EElliptic, map := EllipticCurve(ProjectiveClosure(E), p1);
S1 := map(p1);
S2 := map(p2);
S3 := map(p3);

sigma1 := S1-S2;
sigma2 := S2-S3;
y := (3*sigma1)[2];
Factorisation(Denominator(y));

/*
We find the equation
b1^8 - 837/2*b1^6*c2 + 3645/2*b1^5*b3 + 951345/16*b1^4*c2^2 - 4113747/8*b1^3*b3*c2 +
        10451673/16*b1^2*b3^2 - 42338133/16*b1^2*c2^3 + 301327047/8*b1*b3*c2^2 -
        1420541793/16*b3^2*c2 - 129140163/4*c2^4
*/

/*
Now we look for a rational point of small height
*/

for b1n in [-10..10] do
for b1d in [1..10] do
for b3n in [-10..10] do
for b3d in [1..10] do

S := Scheme( AffineSpace(R), [b1^8 - 837/2*b1^6*c2 + 3645/2*b1^5*b3 + 951345/16*b1^4*c2^2 - 4113747/8*b1^3*b3*c2 +
        10451673/16*b1^2*b3^2 - 42338133/16*b1^2*c2^3 + 301327047/8*b1*b3*c2^2 -
        1420541793/16*b3^2*c2 - 129140163/4*c2^4, b3-b3n/b3d, b1-b1n/b1d] );
	rps := RationalPoints(S);
	if #rps gt 0 then
		rps;
	end if;
end for;
end for;
end for;
end for;

/*
We find the rational point (-9, 2, -11/4)
*/

/*
Test
*/
b1 := -9;
b3 := 2;
c2 := -11/4;


A<x, w> := PolynomialRing(Rationals(), 2);
Q := c2 * x^2 + 7;
P := b3 * x^3 + b1 * x + 10;

R<t> := PolynomialRing(Rationals());
h := hom<A -> R | [t, 0]>;
C := HyperellipticCurve(h(P^2-Q^3));
C;
G2Invariants(C);

E := Curve( AffineSpace(A), [ w^3 - 3*Q*w - 2*P ] );
p1 := E![0,-4];
p2 := E![0,-1];
p3 := E![0, 5];

EElliptic, map := EllipticCurve(ProjectiveClosure(E), p1);
S1 := map(p1);
S2 := map(p2);
S3 := map(p3);

sigma1 := S1-S2;
sigma2 := S2-S3;

Order(sigma1);
Order(sigma2);




/*
=========== Both sections are 2-torsion ===========
*/

R<b1> := PolynomialRing(Rationals());
F<b1> := FieldOfFractions(R);

b3 := 0;

A<x, w> := PolynomialRing(F, 2);
Q := x^2 + 7;
P := b1 * x^2 + 10;

R<t> := PolynomialRing(F);
h := hom<A -> R | [t, 0]>;
C := HyperellipticCurve(h(P^2-Q^3));
C;
G2Invariants(C);


E := Curve( AffineSpace(A), [ w^3 - 3*Q*w - 2*P ] );
p1 := E![0,-4];
p2 := E![0,-1];
p3 := E![0, 5];

d1 := Divisor(p1)-Divisor(p2);
d2 := Divisor(p2)-Divisor(p3);


EElliptic, map := EllipticCurve(ProjectiveClosure(E), p1);
S1 := map(p1);
S2 := map(p2);
S3 := map(p3);

sigma1 := S1-S2;
sigma2 := S2-S3;

2*sigma1;
2*sigma2;


/*
Explicit divisors
*/

R<c2> := PolynomialRing(Rationals());
F<c2> := FieldOfFractions(R);

b1 := 0;
b3 := 7/27* b1 * c2;

A<x, w> := PolynomialRing(F, 2);
Q := c2 * x^2 + 7;
P := b3*x^3 + b1 * x + 10;

E := Curve( AffineSpace(A), [ w^3 - 3*Q*w - 2*P ] );
p1 := E![0,-4];
p2 := E![0,-1];
p3 := E![0, 5];

sigma1 := Divisor(p1) - Divisor(p2);
sigma2 := Divisor(p2) - Divisor(p3);

IsPrincipal(2*sigma1);
IsPrincipal(2*sigma2);










