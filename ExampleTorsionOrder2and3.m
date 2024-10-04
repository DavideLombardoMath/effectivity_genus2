/*
Finding the equation
*/
K<u> := QuadraticField(-3);
R<b1,b2,b3> := PolynomialRing(K,3);
F<b1,b2,b3> := FieldOfFractions(R);


A<x, w> := PolynomialRing(F, 2);
Q := x^2 -9;
P := b3 * x^3 + b2 * x^2 + b1*x;

E := Curve( AffineSpace(A), [ w^3 - 3*Q*w - 2*P ] );
p1 := E![0,0];
p2 := E![0,3*u];
p3 := E![0,-3*u];

EElliptic, map := EllipticCurve(ProjectiveClosure(E), p1);
S1 := map(p1);
S2 := map(p2);
S3 := map(p3);

sigma1 := S1-S2;
sigma2 := S2-S3;

y1 := (2*sigma1)[2];
e1 := Factorisation(Denominator(y1))[2][1];

y2 := (3*sigma2)[2];
e2 := Factorisation(Denominator(y2))[2][1];

/*
e1 = e2 = 0 is the equation of the locus where sigma_1 is 2-torsion and
sigma_3 is 3-torsion. The scheme they define is not irreducible; we take
an irreducible component which is a curve L of genus 0.
*/

S := Scheme(AffineSpace(R), [e1, e2]);
ic := IrreducibleComponents(S);
ic1 := ic[1];

conic, map := Conic(ProjectiveClosure(Curve(ic1)));
conicQ := ChangeRing(conic, Rationals());

/*
We list some points on L
*/
ratpts := RationalPoints(conicQ : Bound := 1000);
for pt in ratpts do
	Inverse(map)(conic![pt[1], pt[2], pt[3]]);
end for;


/*
Test: the coefficients b1, b2, b3 are among the smallest that we obtain from
the above construction.
*/
K<u> := QuadraticField(-3);

b1 := 27;
b2 := 54*u; 
b3 := -53;

A<x, w> := PolynomialRing(F, 2);
Q := x^2 -9;
P := b3 * x^3 + b2 * x^2 + b1*x;
P^2 - Q^3;

E := Curve( AffineSpace(A), [ w^3 - 3*Q*w - 2*P ] );
p1 := E![0,0];
p2 := E![0,3*u];
p3 := E![0,-3*u];

EElliptic, map := EllipticCurve(ProjectiveClosure(E), p1);
S1 := map(p1);
S2 := map(p2);
S3 := map(p3);

sigma1 := S1-S2;
sigma2 := S2-S3;

/*
Let us check that the torsion orders are correct
*/
2*sigma1;
3*sigma2;

R<x> := PolynomialRing(K);

f := 2808*x^6 - 5724*u*x^5 - 11583*x^4 + 2916*u*x^3 + 486*x^2 + 729; // P^2-Q^3

C := HyperellipticCurve(f);
G2Invariants(C);

/*
We find and simplify the planar quartic model
*/

qInfty := C![0,27];
D4 := 4*Divisor(qInfty);

CPlanar := Image(DivisorMap(D4));
fPlanar := DefiningPolynomial(CPlanar);
AmbientPolyRing := Parent(fPlanar);
AssignNames(~AmbientPolyRing, ["z", "w", "U"]);
z := AmbientPolyRing.1;
w := AmbientPolyRing.2;

Evaluate(fPlanar, [2*3^3*z, 2*3^3*w, 1]) / (2^4*3^9);

A<z, w> := PolynomialRing(K, 2);
FinalResult := 27*z^4 + 54*u*z^3 - 9*z^2*w - 27*w^3 - 107*z^2 + 44*u*z*w + 108*w^2 - 52*u*z -
    159*w + 78;
X := Curve(AffineSpace(A), [FinalResult]);
SingularPoints(X);

/*
Translate the singularity to the origin
*/

f2 := Evaluate(FinalResult, [z-u,w]);
/*
27*z^4 - 54*u*z^3 - 9*z^2*w - 107*z^2 + 62*u*z*w - 27*w^3 + 108*w^2
*/

/*
Since our method is geometric, we are allowed to twist the result: we replace z by u*z to get an equation with rational coefficients
*/
1/3 * Evaluate(f2, [u*z,w]);

/*
81*z^4 - 162*z^3 + 9*z^2*w + 107*z^2 - 62*z*w - 9*w^3 + 36*w^2
*/

/*
Finally, we test for the triviality of the geometric endomorphism ring of this curve
*/

O := MaximalOrder(K);
R<x> := PolynomialRing(O);
f := 2808*x^6 - 5724*O!u*x^5 - 11583*x^4 + 2916*O!u*x^3 + 486*x^2 + 729; // P^2-Q^3
C := HyperellipticCurve(f);

test := 0;
for p in [q : q in [5..100] | IsPrime(q) and not (Integers()!Norm(Discriminant(C)) mod q eq 0) and (q mod 3) eq 1 ] do
	I := ideal<O | p>;
	I := Factorisation(I)[1][1];
	FF := ResidueClassField(I);
	_, u := IsSquare(FF!(-3));
	
	RFF<x> := PolynomialRing(FF);
	fFF := 2808*x^6 - 5724*u*x^5 - 11583*x^4 + 2916*u*x^3 + 486*x^2 + 729; // P^2-Q^3

	lpoly := LPolynomial(HyperellipticCurve(fFF));
	if IsIrreducible(lpoly) then
		F<alpha> := NumberField(lpoly);
		if IsIrreducible(MinimalPolynomial(alpha^12)) then	// test for geometric irreducibility
			test := GCD(test, Discriminant(MaximalOrder(F)));
			if test lt 5 then
				print "Trivial geometric endomorphism ring";
			end if;
		end if;
	end if;
end for;
