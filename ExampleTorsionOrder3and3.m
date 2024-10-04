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
y1 := (3*sigma1)[2];
Factorisation(Denominator(y1));


e1 := b1^8 + 162*u*b1^6*b2 - 3645/2*b1^6 - 19683/2*b1^5*b3 - 98415/4*b1^4*b2^2 -
        177147*u*b1^4*b2 + 17537553/16*b1^4 - 1594323*u*b1^3*b2*b3 +
        100442349/8*b1^3*b3 - 1594323/8*u*b1^2*b2^3 + 129140163/16*b1^2*b2^2 +
        387420489/8*u*b1^2*b2 + 387420489/16*b1^2*b3^2 - 3486784401/16*b1^2 +
        387420489/2*b1*b2^2*b3 + 3486784401/4*u*b1*b2*b3 - 31381059609/8*b1*b3 -
        129140163/4*b2^4 - 1162261467/4*u*b2^3 + 31381059609/16*b2^2 +
        31381059609/8*u*b2*b3^2 - 282429536481/16*b3^2;
e2 := b1^8 + 162*(-u)*b1^6*b2 - 3645/2*b1^6 - 19683/2*b1^5*b3 - 98415/4*b1^4*b2^2 -
        177147*(-u)*b1^4*b2 + 17537553/16*b1^4 - 1594323*(-u)*b1^3*b2*b3 +
        100442349/8*b1^3*b3 - 1594323/8*(-u)*b1^2*b2^3 + 129140163/16*b1^2*b2^2 +
        387420489/8*(-u)*b1^2*b2 + 387420489/16*b1^2*b3^2 - 3486784401/16*b1^2 +
        387420489/2*b1*b2^2*b3 + 3486784401/4*(-u)*b1*b2*b3 - 31381059609/8*b1*b3 -
        129140163/4*b2^4 - 1162261467/4*(-u)*b2^3 + 31381059609/16*b2^2 +
        31381059609/8*(-u)*b2*b3^2 - 282429536481/16*b3^2;

/*
e1, e2 are the equations expressing that sigma1, sigma2 are 3-torsion points. We take
their real and imaginary parts f1, f2 below.
*/


R<b1,b2,b3> := PolynomialRing(Rationals(),3);
f1 := 2*b1^8 - 3645*b1^6 - 19683*b1^5*b3 - 98415/2*b1^4*b2^2 + 17537553/8*b1^4 +
    100442349/4*b1^3*b3 + 129140163/8*b1^2*b2^2 + 387420489/8*b1^2*b3^2 -
    3486784401/8*b1^2 + 387420489*b1*b2^2*b3 - 31381059609/4*b1*b3 -
    129140163/2*b2^4 + 31381059609/8*b2^2 - 282429536481/8*b3^2;
f2 := 162*b1^6*b2 - 177147*b1^4*b2 - 1594323*b1^3*b2*b3 - 1594323/8*b1^2*b2^3 +
    387420489/8*b1^2*b2 + 3486784401/4*b1*b2*b3 - 1162261467/4*b2^3 +
    31381059609/8*b2*b3^2;


/*
Together, f1 and f2 give the torsion conditions we're looking for. We look for
points that satisfy both of them. The scheme defined by f1=f2=0 is reducible,
and we only consider an irreducible component that happens to be a curve L of
genus 0 with a rational point.
*/
S := Scheme(AffineSpace(R), [f1, f2]);
S := ReducedSubscheme(S);
ic := IrreducibleComponents(S);

conic, map := Conic(ProjectiveClosure(Curve(ic[2])));

/*
We list some points on L
*/
ratpts := RationalPoints(conic : Bound := 1000);
for pt in ratpts do
	Inverse(map)(pt);
end for;



/*
Test: the coefficients b1, b2, b3 are among the smallest that we obtain from
the above construction.
*/
K<u> := QuadraticField(-3);


b1 := 27;
b2 := 54; 
b3 := 19;

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
Test the torsion orders
*/
3*sigma1;
3*sigma2;


R<x> := PolynomialRing(Rationals());

f := 728*x^6 + 2916*x^5 + 3969*x^4 + 2052*x^3 + 118*x^2 + 729; // P^2-Q^3

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

Evaluate(fPlanar, [2^2*3^3*z, 4*9*w, 1]) / (2^5*3^8);

A<z, w> := PolynomialRing(Rationals(), 2);
FinalResult := 648*z^4 + 648*z^3 - 36*z^2*w - 12*w^3 + 142*z^2 - 56*z*w - 24*w^2 - 20*z - 19*w
    - 5;
X := Curve(AffineSpace(A), [FinalResult]);
SingularPoints(X);

/*
Translate the singularity to the origin
*/

1/2*Evaluate(FinalResult, [z-1/2,w]);
/*
324*z^4 - 324*z^3 - 18*z^2*w + 71*z^2 - 10*z*w - 6*w^3 - 12*w^2
*/




/*
Finally, we test for the triviality of the geometric endomorphism ring of this curve
*/

R<x> := PolynomialRing(Rationals());
f := 728*x^6 + 2916*x^5 + 3969*x^4 + 2052*x^3 + 118*x^2 + 729; // P^2-Q^3
C := HyperellipticCurve(f);

test := 0;
for p in [q : q in [5..100] | IsPrime(q) and not (Integers()!Discriminant(C) mod q eq 0) ] do
	lpoly := LPolynomial(ChangeRing(C, GF(p)));
	if IsIrreducible(lpoly) then
		K<alpha> := NumberField(lpoly);
		if IsIrreducible(MinimalPolynomial(alpha^12)) then	// test for geometric irreducibility
			test := GCD(test, Discriminant(MaximalOrder(K)));
			if test eq 1 then
				print "Trivial geometric endomorphism ring";
			end if;
		end if;
	end if;
end for;
