load "jInvariant.m";

/*
Construct (a version of) the universal curve
with two marked 3-torsion points
*/

R<r, s1, t1, u> := PolynomialRing(Rationals(),4);
F<r, s1, t1, u> := FieldOfFractions(R);
S<x, w> := PolynomialRing(F, 2);

s := -1/4*s1^3;
t := t1^3;

Q1 := s1*(x^2 + r*x + t);
P1 := (s-s*t-1)*x^3 + 3*s*(r-t)*x^2 + 3*s*r*(r-t)*x-s*t^2+s*r^3+t;

Q2 := s1*t1*(x^2 + x + r);
P2 := (s-s*t+1)*x^3 + 3*s*(r-t)*x^2 + 3*s*r*(r-t)*x-s*t^2+s*r^3-t;

H3 := s*x^2 + (2*s*r -s*t - 1)*x + s*r^2;
l3 := 4*t / (s*t+1)^2;

G3 := ((s^2 * t^2 - s^2 * t + 2 * s * t + s + 1) * x^3 
+ (3 * s^2 * t^2 - 3 * s^2 * t * r + 3 * s * t + 3 * s * r) * x^2
+ (3 * s^2 * t^2 * r - 3 * s^2 * t * r^2 + 3 * s * t * r + 3 * s * r^2) * x 
+ s^2 * t^3 - s^2 * t * r^3 + 2 * s * t^2 + s * r^3 + t) / (s*t+1);

// The polynomial which is equal to P1^2-Q1^3 and P2^2-Q2^3 is G3^2 + l3*H3^3
// So we formally define a new variable u with u^3 = -l3 and take
// Q3 = u*H3, P3 = G3

Q3 := u * H3;
P3 := G3;

/*
Auxiliary functions to replace every instance of u^3 with -l3
*/
function ReplaceuInPolynomial(p)
	coeffs := Coefficients(p);
	mons := [Exponents(m) : m in Monomials(p)];
	ExponentsU := [e[3] : e in mons];
	assert GCD(ExponentsU) mod 3 eq 0;
	mons := [ [m[1], m[2], m[3], Integers()!(m[4]/3)] : m in mons ];
	return &+[ coeffs[i] * R.1^mons[i][1] * R.2^mons[i][2] * R.3^mons[i][3] * (-l3)^mons[i][4] : i in [1..#mons] ];
end function;

function ReplaceuInRationalFunction(f)
	n := Numerator(f);
	d := Denominator(f);
	return ReplaceuInPolynomial(n)/ReplaceuInPolynomial(d);
end function;

/*
Compute the j-invariant of the corresponding elliptic curves and
check, using the Jacobian criterion, that they are algebraically
independent
*/

E1 := Curve(AffineSpace(S), [w^3 - 3*Q1*w - P1]);
E2 := Curve(AffineSpace(S), [w^3 - 3*Q2*w - P2]);
E3 := Curve(AffineSpace(S), [w^3 - 3*Q3*w - P3]);

E1 := ProjectiveClosure(E1);
j1 := jInvariantCubic(DefiningPolynomial(E1));

E2 := ProjectiveClosure(E2);
j2 := jInvariantCubic(DefiningPolynomial(E2));

E3 := ProjectiveClosure(E3);
j3 := jInvariantCubic(DefiningPolynomial(E3));
j3 := ReplaceuInRationalFunction(j3);




/*
See https://mathoverflow.net/questions/41535/how-to-show-a-set-of-polynomials-is-algebraically-independent for this version of the Jacobian criterion
*/
d1 := Denominator(j1);
d2 := Denominator(j2);
d3 := Denominator(j3);
n1 := Numerator(j1);
n2 := Numerator(j2);
n3 := Numerator(j3);

ds := [d1, d2, d3];
ns := [n1, n2, n3];

function JacobianMatrixEntry(i, j)
	return (Derivative(ns[i], Parent(ns[i]).j) * ds[i] - ns[i] * Derivative(ds[i], Parent(ds[i]).j)) / ds[i]^2;
end function;

M := Matrix( Parent(j1), 3, 3, [ JacobianMatrixEntry(i,j) : i in [1..3], j in [1..3] ] );

assert Determinant(M) ne 0;
