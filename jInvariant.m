/*
Functions to compute the j-invariant of a general plane cubic
*/
function myCoefficient(polynomial, monomial)
	coeffs := Coefficients(polynomial);
	monomials := Monomials(polynomial);
	for i in [1..#monomials] do
		if monomials[i] eq monomial then
			return coeffs[i];
		end if;
	end for;
	return 0;
end function;

function SandT(a,b,c,d,e,f,g,h,i,j)
	S := a*g*e*c - a*g*h^2 - a*j*b*c + a*j*e*h + a*f*b*h - a*f*e^2 - d^2*e*c + d^2*h^2 + d*i*b*c - d*i*e*h + d*g*j*c - d*g*f*h - 2*d*j^2*h + 3*d*j*f*e - d*f^2*b - i^2*b*h + i^2*e^2 - i*g^2*c + 3*i*g*j*h - i*g*f*e - 2*i*j^2*e + i*j*f*b + g^2*f^2 - 2*g*j^2*f + j^4;

	T := a^2*b^2*c^2 - 3*a^2*e^2*h^2 - 6*a^2*b*e*h*c + 4*a^2*b*h^3 + 4*a^2*e^3*c - 6*a*d*g*b*c^2 + 18*a*d*g*e*h*c - 12*a*d*g*h^3 + 12*a*d*j*b*h*c - 24*a*d*j*e^2*c + 12*a*d*j*e*h^2 - 12*a*d*f*b*h^2 + 6*a*d*f*b*e*c + 6*a*d*f*e^2*h + 6*a*i*g*b*h*c - 12*a*i*g*e^2*c + 6*a*i*g*e*h^2 + 12*a*i*j*b*e*c + 12*a*i*j*e^2*h - 6*a*i*f*b^2*c + 18*a*i*f*b*e*h - 24*a*g^2*j*h*c - 24*a*i*j*b*h^2 - 12*a*i*f*e^3 + 4*a*g^3*c^2 - 12*a*g^2*f*e*c + 24*a*g^2*f*h^2 + 36*a*g*j^2*e*c + 12*a*g*j^2*h^2 + 12*a*g*j*f*b*c - 60*a*g*j*f*e*h - 12*a*g*f^2*b*h + 24*a*g*f^2*e^2 - 20*a*j^3*b*c - 12*a*j^3*e*h + 36*a*j^2*f*b*h + 12*a*j^2*f*e^2 - 24*a*j*f^2*b*e + 4*a*f^3*b^2 + 4*d^3*b*c^2 - 12*d^3*e*h*c + 8*d^3*h^3 + 24*d^2*i*e^2*c - 12*d^2*i*e*h^2 + 12*d^2*g*j*h*c + 6*d^2*g*f*e*c - 24*d^2*j^2*h^2 - 12*d^2*i*b*h*c - 3*d^2*g^2*c^2 - 24*g^2*j^2*f^2 + 24*g*j^4*f - 12*d^2*g*f*h^2 + 12*d^2*j^2*e*c - 24*d^2*j*f*b*c - 27*d^2*f^2*e^2 + 36*d^2*j*f*e*h + 24*d^2*f^2*b*h + 24*d*i^2*b*h^2 - 12*d*i^2*b*e*c - 12*d*i^2*e^2*h + 6*d*i*g^2*h*c - 60*d*i*g*j*e*c + 36*d*i*g*j*h^2 + 18*d*i*g*f*b*c - 6*d*i*g*f*e*h + 36*d*i*j^2*b*c - 12*d*i*j^2*e*h - 60*d*i*j*f*b*h + 36*d*i*j*f*e^2 + 6*d*i*f^2*b*e + 12*d*g^2*j*f*c - 12*d*g*j^3*c - 12*d*g*j^2*f*h + 36*d*g*j*f^2*e - 12*d*g*f^3*b + 24*d*j^4*h + 12*d*j^2*f^2*b + 4*i^3*b^2*c + 24*i^2*g^2*e*c - 27*i^2*g^2*h^2 - 36*d*j^3*f*e - 12*i^3*b*e*h + 8*i^3*e^3 - 24*i^2*g*j*b*c + 36*i^2*g*j*e*h + 6*i^2*g*f*b*h + 12*i^2*j^2*b*h - 3*i^2*f^2*b^2 - 12*d*g^2*f^2*h - 12*i^2*g*f*e^2 - 24*i^2*j^2*e^2 + 12*i^2*j*f*b*e - 12*i*g^3*f*c + 12*i*g^2*j^2*c + 36*i*g^2*j*f*h - 12*i*g^2*f^2*e - 36*i*g*j^3*h - 12*i*g*j^2*f*e + 12*i*g*j*f^2*b + 24*i*j^4*e - 12*i*j^3*f*b + 8*g^3*f^3 - 8*j^6;


	return S,T;

end function;

function SandTHomogeneous(cubic)

	R := Parent(cubic);
	x := R.1;
	y := R.2;
	z := R.3;


	a := myCoefficient(cubic, x^3);
	b := myCoefficient(cubic, y^3);
	c := myCoefficient(cubic, z^3);
	d := 1/3 * myCoefficient(cubic, x^2*y);
	e := 1/3 * myCoefficient(cubic, y^2*z);
	f := 1/3 * myCoefficient(cubic, z^2*x);

	g := 1/3 * myCoefficient(cubic, x*y^2);
	h := 1/3 * myCoefficient(cubic, y*z^2);
	i := 1/3 * myCoefficient(cubic, z*x^2);

	j := 1/6 * myCoefficient(cubic, x*y*z);



	return SandT(a,b,c,d,e,f,g,h,i,j);

end function;


function jInvariantCubicHomogeneous(cubic)
	S, T:= SandTHomogeneous(cubic);
	return 1728*(4*S)^3/((4*S)^3 - T^2);
end function;

function SandTInhomogeneous(cubic)

	R := Parent(cubic);
	x := R.1;
	y := R.2;


	a := myCoefficient(cubic, x^3);
	b := myCoefficient(cubic, y^3);
	c := myCoefficient(cubic, 1);
	d := 1/3 * myCoefficient(cubic, x^2*y);
	e := 1/3 * myCoefficient(cubic, y^2*1);
	f := 1/3 * myCoefficient(cubic, 1^2*x);

	g := 1/3 * myCoefficient(cubic, x*y^2);
	h := 1/3 * myCoefficient(cubic, y*1^2);
	i := 1/3 * myCoefficient(cubic, 1*x^2);

	j := 1/6 * myCoefficient(cubic, x*y*1);

	return SandT(a,b,c,d,e,f,g,h,i,j);

end function;


function jInvariantCubicHomogeneous(cubic)
	S, T:= SandTHomogeneous(cubic);
	return 1728*(4*S)^3/((4*S)^3 - T^2);
end function;

function jInvariantCubicInhomogeneous(cubic)
	S, T:= SandTInhomogeneous(cubic);
	return 1728*(4*S)^3/((4*S)^3 - T^2);
end function;
