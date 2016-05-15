# Algorithms
	Find a quadratic residue (mod p) of 'a'. p must be an odd prime.
	Solve the congruence of the form: x^2 = a (mod p) and returns x. 
	Note that p - x is also a root.

	0 is returned if no square root exists for these a and p.

	The Tonelli-Shanks algorithm is used (except for some simple cases in which the solution
	is known from an identity). This algorithm runs in polynomial time
	(unless the generalized Riemann hypothesis is false).

#License
	Do what ever you want