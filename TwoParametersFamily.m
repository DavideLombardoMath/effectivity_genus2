load "jInvariant.m";
/*
We will work over a very large finite field that contains the
coordinates of (most of) the points we are interested in.
*/
q := 103^60;

/*
Construct a finite cover of the affine line A^2(a,b) over which
a three-torsion point of Jac(X_{a,b}) is defined, together with
a determination of a trivialisation.
*/
A<c3,c2,c1,c0,b2,b1,b0,a,b> := PolynomialRing(GF(q), 9);
S<m> := PolynomialRing(A);
FF := FieldOfFractions(A);
AssignNames(~FF, ["c3", "c2", "c1", "c0", "b2", "b1", "b0", "a", "b"]);
T<x,w> := PolynomialRing(FF, 2);
h := hom<S->T | [x]>;

P := c3*m^3+c2*m^2+c1*m+c0;
Q := b2*m^2+b1*m+b0;

f1 := P^2 - Q^3;
f2 := 1-4*m^4*( b+a*m^2-m );	// this is the polynomial giving the hyperelliptic model

coef1 := Coefficients(f1);
coef2 := Coefficients(f2);

/*
We define the scheme S given by equating the coefficients of f2 to the
coefficients of a formal expression P^2-Q^3. This is the cover of A^2(a,b)
over which a 3-torsion point + trivialisation is naturally defined. We also
define the universal elliptic curve depending on the pair (P, Q).
*/
S := Scheme(AffineSpace(A), [ coef1[i] - coef2[i] : i in [1..#coef1] ]);

E := Curve(AffineSpace(T), [ w^3-3*h(Q)*w-2*h(P) ]);

/*
Now we restrict to the subscheme where j takes a prescribed value, for example j=1.
To do so, we compute the j-invariant of E as a curve over the appropriate function field,
and then impose the equation j_E = constant (=1).
*/

/*
We start by computing the j-invariant as a function on the moduli space.
*/
fE := DefiningPolynomial(E);
jE := jInvariantCubicInhomogeneous(fE);

/*
We also fix an elliptic curve E_3 over F_q with j-invariant 1. All the curves we consider
will be at least geometrically isomorphic to E_3. Note that we work with a short Weierstrass model of E_3, given by y^2 = g(x) for a suitable g(x). This will be important below.
*/
E3 := EllipticCurveFromjInvariant(GF(q)!1);
E3 := WeierstrassModel(E3);

/*
The principle of the proof is now the following. We explore our moduli space S_1 (which is 1-dimensional: it is the subvariety of S, which is 2-dimensional, where the elliptic curve E has j-invariant 1) by trying all values of a \in F_{103}. For a given a=a_0, there will be many values of b, not necessarily defined over F_{103} (but certainly over a finite extension), and for each of those -- since a decomposition f_{a,b}=P^2-Q^3 is defined over S -- there are corresponding polynomials P, Q. What we are trying to show is that the restriction of the section to S_1 is non-constant.
To check this, for each s_1 in S, we write down an isomorphism from E_{s_1} to E_3 and consider the image of the section. Since E_3 does not have non-trivial automorphisms apart from [-1], the image of the section in E_3 is well-defined up to [-1]. In particular, in short Weierstrass form, its x-coordinate is well-defined. 

If we knew that S_1 is geometrically connected, it would suffice to find two different x-coordinates in this way. Unfortunately, since we do not know that S_1 is geometrically connected (and it seems that checking this computationally is too hard), we need to find many different specialisations. We bound the number of connected components with the degree of the projective closure of the scheme S_1; the bound is 540.

There is one more detail to take into account: since we do not want to label the three points on E_{s_1} with given x-coordinate (doing this would require passing to a further cover of S_1), there is also an ambiguity in the definition of the differences p_{1, s_1}-p_{2, s_1}, p_{1, s_2}-p_{2, s_3}. Notice however that, if the section were to be constant, these differences (up to the action of [-1]) could take on at most 3 * 540 values. Thus, if we find more than this number of possible x-coordinates for the differences p_{i, s_1}-p_{j, s_1} as s_1 varies, we are done. This is what we check below.
*/
xCoords := {GF(q)!0};

jeq := Numerator(jE-1);
Sconstantj := Scheme(AffineSpace(A), [ coef1[i] - coef2[i] : i in [1..#coef1] ] cat [jeq]);
"Bound on the number of connected components", Degree(ProjectiveClosure(Sconstantj));


tot := 0;
for a0 in GF(103) do

Sconstantj := Scheme(AffineSpace(A), [ coef1[i] - coef2[i] : i in [1..#coef1] ] cat [jeq] cat [a-a0] );

rp:=RationalPoints(Sconstantj);
tot := tot + #rp;

/*
Here s1 is the point s_1 \in S_1 of the above discussion
*/
for s1 in rp do
	try
	/*
	MAGMA sometimes fails in applying the isomorphism between the different models of elliptic curves involved. We simply skip these points.
	*/
	/*
	Given a point of S_1, we compute the corresponding polynomials P, Q and elliptic curve E
	*/
	c3p,c2p,c1p,c0p,b2p,b1p,b0p,ap,bp := Explode( [s1[i] : i in [1..9]] );
	K := GF(q);
	T<xp,wp> := PolynomialRing(K, 2);
	P := c3p*xp^3+c2p*xp^2+c1p*xp+c0p;
	Q := b2p*xp^2+b1p*xp+b0p;
	E := Curve(AffineSpace(T), [ wp^3-3*Q*wp-2*P ]);
	EProj := ProjectiveClosure(E);
	E2, map2 := EllipticCurve(EProj);
	/*
	It is possible that E_{s_1} (here called E_2) and E_3 become isomorphic
	only over a quadratic extension. For simplicity, we only work with those
	values of s_1 that make E_{s_1}, E_3 isomorphic over the ground field.
	*/
	test := IsIsomorphic(E2, E3);
	if test then
		_, map3 := IsIsomorphic(E2, E3);
	end if;

	/*
	Next we compute the coordinates of the points p_{i, s_1} for i=1,2,3
	*/
	PP := AmbientSpace(EProj);
	SectPoints := Scheme(PP, [ PP.1, DefiningPolynomials(EProj)[1] ]);

	if #RationalPoints(SectPoints) eq 3 then
		SectRationalPoints := RationalPoints(SectPoints);
		Q1 := SectRationalPoints[1];
		Q2 := SectRationalPoints[2];
		Q3 := SectRationalPoints[3];
		// we map the section to points on E3
		R1 := map3(map2(Q1));
		R2 := map3(map2(Q2));
		R3 := map3(map2(Q3));
		xCoords := xCoords join { (R1-R2)[1] } ;
		xCoords := xCoords join { (R1-R3)[1] } ;
		xCoords := xCoords join { (R2-R3)[1] } ;
		"Distinct values", #xCoords;
		if #xCoords gt 3*(631) then
			"Found sufficiently many distinct specialisations, the proof is finished.";
		end if;
	end if;
	catch e
		;
	end try;
end for;

end for;

