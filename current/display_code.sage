from sage.all import *
from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_codomain_kohel


def apply_ideal(E, n, v=-1):
    """Returns (n, pi-1)E, (n)E, or (n, pi+1)E depending on v=-1,0,+1

    Arguments:
        - E, an elliptic curve of the form y^2 = f(x) over F_p for a prime p
                satisfying p>3 and p=3 mod 8.
        - n, an odd prime dividing p+1.
        - v, an integer between -1 and 1, specifying what ideal to apply to
                the elliptic curve.

    Output:
        - E', an elliptic curve in the F_p-isomorphism class of
                (n, pi-1)E, (n)E, or (n, pi+1)E, depending on `v`.
              The resulting elliptic curve is found using Kohel's algorithm.
    """
    assert E.a1() == E.a3() == 0
    assert is_prime(n) and n & 1 and (p+1) % n == 0
    assert -1 <= v <= 1 and int(v) == v

    # Do the case v=0 separately, as this is much faster.
    if v == 0:
        Esw = E.short_weierstrass_model()
        return EllipticCurve([0, 0, 0, n^4 * Esw.a4(), n^6 * Esw.a6()])

    # We find the function f(x) such that E: y^2 = f(x).
    f = E.hyperelliptic_polynomials()[0]

    # We aspire to find a generator Q of the n-torsion points on the elliptic curve, i.e.,
    # a generator of the set {P in E: [n]P = O}
    if v == 1:
        # We want to apply the ideal (n, pi+1) to E.
        # ker phi = {P in E: [n]P = O} ∩ {P in E: pi(P)+P = O}
        # If pi(P) + P = O, then pi(P) = -P, and if P=(x,y), then x^p=x and y^p=-y is required.
        # So x in F_p and y^(p-1) = -1, so P in E(F_{p^2}).
        # We go through random points on the curve an check whether they are generators of the
        # n-torsion points of E(F_{p^2}).
        F = GF(p)
        EF2 = E.base_extend(GF(p^2))
        while True:
            # We still have that x in F_p.
            x_coor = F.random_element()
            # Check whether this point satisfies y^(p-1)=-1, so that since y^2=f(x) we find that
            # f(x)^((p-1)/2) = -1, giving that f(x) is not a quadratic residue.
            if not f(x_coor).is_square():
                # This point is in E(F_{p^2}), find the corresponding coordinates.
                G = EF2.lift_x(x_coor)
                # Assume that G is a generator of the n-torsion points of E(F_{p^2}), so that
                # G has order equal to p+1, then we find that [n]*[(p+1)/n]*G = O as
                # n divides p+1, therefore [(p+1)/n]*G has order n (or a divisor of n).
                Q = ((p+1)//n)*G
                # The only point of order a divisor of n is the point at infinity, check that
                # we did not find this point.
                if Q != EF2.point(0):
                    break
    elif v == -1:
        # The commented code
        # while True:
        #     x_coor = GF(p).random_element()
        #     if f(x_coor).is_square():
        #         G = E.lift_x(x_coor)
        #         Q = ((p+1)//n)*G
        #         if Q != E.point(0):
        #             break
        # is actually the same as the following built-in code.
        G = E.gens()[0]
        Q = ((p+1)//n)*G

    # The multiples [k]*Q = [k*(p+1)/n]*G still satisfy [n]*[k*(p+1)/n]*G = O, such that
    # each of [k]*Q is an n-torsion point, these are in fact all the n-torsion points that
    # we need to check.
    # We only need to know the distinct x-coordinates of [k]*Q.
    # We will repetitively add Q to itself, stopping at [n]*Q (exclusive).
    # Since [n-1]*Q = [-1]*Q = -Q has the same x-coordinate as [1]Q = Q and -Q has the same
    # x-coordinate as Q, we can even stop our search at [n//2]*Q (inclusive).
    P = Q
    xs = {Q.x()}
    for _ in range(n//2-1):
        P += Q
        xs.add(P.x())

    # We are ready to compute the monic kernel polynomial of the isogeny phi, remember that
    # ker phi = {P in E: [n]P = O} ∩ {P in E: pi(P) \pm P = O}
    # and `xs` are the x-coordinates of all points in {P in E: [n]P = O}.
    # A monic kernel polynomial of the resulting isogeny will be the product
    # of (x-P_x), where P_x is the x-coordinate of the point in the kernel.
    x = polygen(GF(p))
    if v == -1:
        # If (pi-1)(P) = O, then pi(P) = P, so that if P=(x,y) we require that x^p=x and y^p=y.
        # This gives that both x in F_p and y in F_p. Then, we need y^(p-1)=1. But, y^2=f(x),
        # so that the equality becomes y^(p-1)=f(x)^((p-1)/2)=1, implying that we need f(x)
        # to be a quadratic residue, i.e., a square in F_p.
        kernel_pol = prod(x - px for px in xs if f(px).is_square())
    elif v == 1:
        # If (pi+1)(P) = O, then pi(P) = -P, so that if P=(x,y) we require that x^p=x and y^p=-y.
        # We again get that x in F_p, but require that y^(p-1)=f(x)^((p-1)/2)=-1, implying that
        # f(x) may not be a square in F_p.
        # First, map each x-coordinate in F_{p^2} to a point in F_p (we already know that these
        # are in F_p so that we can apply this mapping).
        F = GF(p)
        xs = [F(px) for px in xs]
        kernel_pol = prod(x - px for px in xs if not f(px).is_square())

    # Specifying the algorithm that SageMath uses to compute the codomain of the isogeny
    # corresponding to the kernel speeds up the calculation significantly over using
    # `E.isogeny`, I presume that this is because SageMath does not need to compute the
    # isogeny in this case.
    return compute_codomain_kohel(E, kernel_pol)

def apply_ideals(E, ps, mults):
    """Applies the ideal (ps[i], pi-1)**mults[i] to E for all i
    Note: if mults[i] < 0, then we apply (ps[i], pi+1)**(-mults[i])."""
    e = E
    for n, mult in zip(ps, mults):
        if mult < 0:
            for _ in range(-mult):
                e = apply_ideal(e, n, 1)
        else:
            for _ in range(mult):
                e = apply_ideal(e, n, -1)
    return e



### Define your prime `p` here as the product of small primes in `ps`
ps = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
p = 4*prod(ps)-1

# Make sure that p>3 is prime and p = 3 (mod 8).
# Then, make sure that p is -1+4*[the product of distinct odd primes]
assert is_prime(p)
assert p % 8 == 3 and p > 3
assert factor((p+1)/4).radical_value() == (p+1)/4

# Create the starting curve, E0
E0 = EllipticCurve(GF(p), [1, 0])

# Display some information
print()
print(f'p = {p}, with starting curve {str(E0)[26:-27-len(str(p))]}')
K.<pi> = NumberField(E0.frobenius_polynomial())
print(f'The number of nodes in the isogeny graph is {K.order(pi).class_number()}')
print()

### Applying the ideals to the elliptic curves according to the example.
print()
print("Alice's resulting curve after applying her private key is")
E_A = apply_ideals(E0, (3, 13, 17, 29, 31), (1, -2, 2, -2, 1)).montgomery_model()
print(E_A)

print()
print("Bob's resulting curve after applying his private key is")
E_B = apply_ideals(E0, (5, 7, 13, 17, 19, 23), (2, 1, -2, 1, -1, 2)).montgomery_model()
print(E_B)

print()
print()
print("Alice applies her private key to E_B and gets")
print(apply_ideals(E_B, (3, 13, 17, 29, 31), (1, -2, 2, -2, 1)).montgomery_model())
print()
print("Bob applies his private key to E_A and gets")
print(apply_ideals(E_A, (5, 7, 13, 17, 19, 23), (2, 1, -2, 1, -1, 2)).montgomery_model())
print()
print("They thus get the same curve E_S.")
