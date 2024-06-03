from sage.all import *


def find_kernel_pol(E, gen):
    """Returns the kernel polynomial of a generator of an ideal"""
    assert E.a1() == E.a3() == 0
    eq_y = x^3 + E.a2()*x^2 + E.a4()*x + E.a6()

    pol = gen.polynomial()
    coeff = pol.padded_list()
    # Note that the generator of the ideal can be thought of as a morphism
    # of the elliptic curve.
    if len(coeff) == 1:
        # This morphism is only a multiplication of the form `[coeff[0]]`
        return E.scalar_multiplication(coeff[0]).kernel_polynomial()
    elif len(coeff) == 2:
        c0, c1 = coeff
        # In this case the morphism is of the form `[c0] + [c1]pi`
        if c0 == 0:
            # In this case the morphism is of the form `[c1]pi`
            if coeff[1] == 1:
                # In this case the morphism is of the form `pi`
                # The kernel_polynomial of the map `pi` is just 1
                return 1
            else:
                # We know that `pi` has rational_maps equal to (x^p, y^p).
                # We also know that there is no `y` in the x-coordinate
                # of Sage's scalar multiplication.
                scalar = E.scalar_multiplication(c1).rational_maps()
                return scalar[0].subs(x=x^p).denominator()
        elif c1 == 1:
            # We need to find the kernel polynomial of the morphism `[c0] + pi`.
            scalar = E.scalar_multiplication(c0).rational_maps()
            # We know that `pi` has rational_maps equal to (x^p, y^p).
            # Scalar multiplication in Sage has only one factor of `y` in the
            # y-coordinate, we make use of this fact by removing this factor of `y`
            # from the scalar multiplication as well as from the `pi` rational map.
            # We add the factor back in afterwards.
            x1, y1 = (x^p, (eq_y)^((p-1)//2))
            x2, y2 = (scalar[0], scalar[1]/y)
            l = (y2 - y1) / (x2 - x1)
            x3 = eq_y * l**2 - E.a2() - x1 - x2
            return x3.denominator()
        else:
            # We need to find the kernel polynomial of the morphism `[c0] + [c1]pi`
            # First, we find the rational_maps of `[c1]pi`.
            # Similar to above, we divide out a factor of `y` from the y-coordinate
            # of `[c1]` as well as from `pi`, we will add this factor at the end.
            x1, y1 = (x^p, (eq_y)^((p-1)//2))
            scalar = E.scalar_multiplication(c1).rational_maps()
            maps = (scalar[0].subs(x=x1, y=y1), scalar[1].subs(x=x1, y=y1))
            # `maps` represents the rational_maps of `[c1]pi` where the
            # y-coordinate is divided by `y`.

            # Here we find the rational_maps of `[c0]`, again, we divide the
            # y-coordinate by `y`.
            scalar = E.scalar_multiplication(c0).rational_maps()
            x2, y2 = (scalar[0], scalar[1]/y)

            # At this point we add the maps `[c0]` and `[c1]pi`, assuming that
            # the rational maps are unequal (as polynomials).
            l = (y2 - maps[1]) / (x2 - maps[0])
            # We add a factor of `y^2` back by multiplying with `eq_y`.
            x3 = eq_y * l**2 - E.a2() - x1 - x2
            return x3.denominator()

    raise BaseException('Unknown morhpism type')

def apply_ideal(E, ideal):
    """Perform the action of `ideal` on the elliptic curve `E`"""
    assert p & 1
    assert E.a1() == E.a3() == 0
    eq_y = x^3 + E.a2()*x^2 + E.a4()*x + E.a6()

    kernel_pols = []
    for b in ideal.gens():  # Typically, `gens_reduced` is desired here.
        # We need to find the kernel polynomial of b.
        # E.g. if b = s+16 then we need to find the kernel polynomial of the
        # morhpism `pi + [16]` where + signifies elliptic curve addition.
        kernel_pol = find_kernel_pol(E, b)
        kernel_pols.append(kernel_pol)

    if len(kernel_pols) > 1:
        kernel_pol = gcd(*kernel_pols)
    else:
        kernel_pol = kernel_pols[0]

    # We need to make an isogeny from k with the kernel polynomial 'kernel_pol'
    # To this end we define an ordinary polynomial (not over any finite field)
    # 'equal' to kernel_pol.
    _.<t> = ZZ[]
    pol = kernel_pol.subs(x=t).monic()
    # We map all the multiple roots to single roots
    if any(multiplicity >= 2 for _, multiplicity in pol.roots()):
        pol = prod(t-i for i, _multiplicity in pol.roots())
    isog = E.isogeny(pol)

    return isog.codomain()

def apply_ideals(E, ideals, mults):
    """Applies multiple ideals to an elliptic curve E"""
    e = E
    for f, mult in zip(ideals, mults):
        for _ in range(mult):
            e = apply_ideal(e, f).montgomery_model()
    return e


### Define your prime `p` here
p = 4*3*5*7*43-1
assert is_prime(p)

# We create `F`, which equals $F_p$ and `E` which equals $E_0$.
F = GF(p)
E = EllipticCurve(F, [1, 0])

# Redefine x and y
_.<x, y> = F[]

poly = E.frobenius_polynomial()
K.<pi> = NumberField(poly)

print()
print(f'p = {p}, with starting curve {str(E)[26:-27-len(str(p))]}')
print(f'The Frobenius polynomial is {poly}')
print(f'The class group has class number {K.order(pi).class_number()}')

# Print disclaimer
print()
print('The below result can take about a minute, this implementation is not efficient over large finite fields.')
print('However, this implementation is very versatile for smaller examples')
print()

### Choose your ideals here.
I3 = K.ideal(3, pi-1)
I3inv = K.ideal(3, pi+1)
I5 = K.ideal(5, pi-1)
I5inv = K.ideal(5, pi+1)
I7inv = K.ideal(7, pi+1)
I43 = K.ideal(43, pi-1)

# Applying the ideals to the elliptic curves according to the example.
print("Alice's resulting curve after applying her ideals is")
E_A = apply_ideals(E, (I3, I5inv, I43), (1, 1, 2)).montgomery_model()
print(E_A)

print()
print("Bob's resulting curve after applying his ideals is")
E_B = apply_ideals(E, (I3inv, I5, I7inv, I43), (2, 2, 1, 1)).montgomery_model()
print(E_B)

print()
print("Alice applies her ideals to E_B and gets")
print(apply_ideals(E_B, (I3, I5inv, I43), (1, 1, 2)).montgomery_model())
print("Bob applies his ideals to E_A and gets")
print(apply_ideals(E_A, (I3inv, I5, I7inv, I43), (2, 2, 1, 1)).montgomery_model())
print("They thus get the same curve E_S.")
