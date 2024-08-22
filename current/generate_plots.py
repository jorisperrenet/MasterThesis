import numpy as np
import matplotlib.pyplot as plt
import sympy
import os

def EC(x, y, a=1., b=-1.):
    return pow(y, 2) - pow(x, 3) - x * a - b

def elliptic_curve_addition():
    a, b = -1.5, 2

    xl, xr = -1.9, 3
    yl, yr = -5, 5
    plt.xlim(xl, xr)
    plt.ylim(yl, yr)

    y, x = np.ogrid[yl:yr:100j, xl:xr:100j]

    plt.contour(x.ravel(), y.ravel(), EC(x, y, a, b), [0])

    p1, p2 = (-9, 9), (4, -5)
    plt.plot((p1[0], p2[0]), (p1[1], p2[1]))

    # Find intersection points
    # equation of the line: y = (p2[1]-p1[1])/(p2[0]-p1[0]) x + ex
    dx = (p2[1]-p1[1])/(p2[0]-p1[0])
    ex = p1[1] - dx*p1[0]
    # equation of the curve y^2 = x^3 + ax + b
    # Plugging in
    # equation of the curve (dx*x + ex)^2 = x^3 + ax + b
    x, y = sympy.symbols('x y')
    xs = sympy.solve(sympy.Eq(x*x*x + a*x + b, (dx*x + ex)**2))
    xs = [i.as_real_imag()[0] for i in xs]
    ys = [dx*x + ex for x in xs]
    plt.scatter(xs, ys, color='k')

    P, Q, PQ = zip(xs, ys)
    kwargs = {'fontsize': 20}
    dx, dy = -.1, .3
    plt.annotate(r'$P$', xy=P, xytext=(P[0]+dx, P[1]+dy), **kwargs)
    plt.annotate(r'$Q$', xy=Q, xytext=(Q[0]+dx, Q[1]+dy), **kwargs)
    plt.annotate(r'$P*Q$', xy=PQ, xytext=(PQ[0]+dx, PQ[1]+dy), **kwargs)

def elliptic_curve_doubling():
    a, b = -1.5, 2

    xl, xr = -1.9, 3
    yl, yr = -5, 5
    plt.xlim(xl, xr)
    plt.ylim(yl, yr)

    y, x = np.ogrid[yl:yr:100j, xl:xr:100j]

    plt.contour(x.ravel(), y.ravel(), EC(x, y, a, b), [0])

    # Find a randon point on the curve
    x = 1.5
    y = np.sqrt(x*x*x + a*x + b)

    # Find intersection points of tangent
    dx = (3*x*x + a) / (2*y)
    ex = y - dx*x

    plt.plot((xl, xr), (dx*xl+ex, dx*xr+ex))
    # Plugging in
    # equation of the curve (dx*x + ex)^2 = x^3 + ax + b
    x, y = sympy.symbols('x y')
    xs = sympy.solve(sympy.Eq(x*x*x + a*x + b, (dx*x + ex)**2))
    xs = [i.as_real_imag()[0] for i in xs]
    ys = [dx*x + ex for x in xs]
    plt.scatter(xs, ys, color='k')

    PP, P = zip(xs, ys)
    kwargs = {'fontsize': 20}
    plt.annotate(r'$P$', xy=P, xytext=(P[0]+.05, P[1]-.6), **kwargs)
    plt.annotate(r'$P*P$', xy=PP, xytext=(PP[0]+.05, PP[1]-.7), **kwargs)


def elliptic_curve_composition_1():
    a, b = -2.5, 2

    xl, xr = -2.3, 3
    yl, yr = -5, 5
    plt.xlim(xl, xr)
    plt.ylim(yl, yr)

    y, x = np.ogrid[yl:yr:100j, xl:xr:100j]

    plt.contour(x.ravel(), y.ravel(), EC(x, y, a, b), [0])

    p1, p2 = (-8, 5), (4, -3)
    rx = -.2
    ry = -np.sqrt(rx*rx*rx + a*rx + b)

    # Find intersection points of P*Q
    dx = (p2[1]-p1[1])/(p2[0]-p1[0])
    ex = p1[1] - dx*p1[0]
    plt.plot((xl, xr), (dx*xl+ex, dx*xr+ex))
    x, y = sympy.symbols('x y')
    xs = sympy.solve(sympy.Eq(x*x*x + a*x + b, (dx*x + ex)**2))
    xs = [i.as_real_imag()[0] for i in xs]
    ys = [dx*x + ex for x in xs]
    plt.scatter(xs, ys, zorder=2, color='k')
    P, Q, PQ = zip(xs, ys)
    kwargs = {'fontsize': 20}
    dx, dy = -.1, .3
    plt.annotate(r'$P$', xy=P, xytext=(P[0]+dx, P[1]+dy), **kwargs)
    plt.annotate(r'$Q$', xy=Q, xytext=(Q[0]+dx, Q[1]+dy), **kwargs)
    plt.annotate(r'$P*Q$', xy=PQ, xytext=(PQ[0]+dx, PQ[1]+dy), **kwargs)

    # Find intersection points of (P*Q)*R
    dx = (ry-PQ[1])/(rx-PQ[0])
    ex = PQ[1] - dx*PQ[0]
    plt.plot((xl, xr), (dx*xl+ex, dx*xr+ex), zorder=1)
    x, y = sympy.symbols('x y')
    xs = sympy.solve(sympy.Eq(x*x*x + a*x + b, (dx*x + ex)**2))
    xs = [i.as_real_imag()[0] for i in xs]
    ys = [dx*x + ex for x in xs]
    plt.scatter(xs[:2], ys[:2], zorder=3, color='k')
    PQR, R, _ = zip(xs, ys)
    kwargs = {'fontsize': 20}
    plt.annotate(r'$R$', xy=R, xytext=(R[0]-.1, R[1]+.3), **kwargs)
    # plt.annotate(r'$T$', xy=T, xytext=(T[0]+dx, T[1]+dy), **kwargs)
    plt.annotate(r'$(P*Q)*R$', xy=PQR, xytext=(PQR[0]-.8, PQR[1]-.8), **kwargs)


def elliptic_curve_composition_2():
    a, b = -2.5, 2

    xl, xr = -2.3, 3
    yl, yr = -5, 5
    plt.xlim(xl, xr)
    plt.ylim(yl, yr)

    y, x = np.ogrid[yl:yr:100j, xl:xr:100j]

    plt.contour(x.ravel(), y.ravel(), EC(x, y, a, b), [0])

    p1, p2 = (-8, 5), (4, -3)
    rx = -.2
    ry = -np.sqrt(rx*rx*rx + a*rx + b)
    R = (rx, ry)

    # Find intersection points of P*Q
    dx = (p2[1]-p1[1])/(p2[0]-p1[0])
    ex = p1[1] - dx*p1[0]
    x, y = sympy.symbols('x y')
    xs = sympy.solve(sympy.Eq(x*x*x + a*x + b, (dx*x + ex)**2))
    xs = [i.as_real_imag()[0] for i in xs]
    ys = [dx*x + ex for x in xs]
    P, Q, PQ = zip(xs, ys)


    # Find intersection points of Q*R
    dx = (Q[1] - R[1]) / (Q[0] - R[0])
    ex = R[1] - dx*R[0]
    plt.plot((xl, xr), (dx*xl+ex, dx*xr+ex))
    x, y = sympy.symbols('x y')
    xs = sympy.solve(sympy.Eq(x*x*x + a*x + b, (dx*x + ex)**2))
    xs = [i.as_real_imag()[0] for i in xs]
    ys = [dx*x + ex for x in xs]
    plt.scatter(xs, ys, zorder=2, color='k')
    R, QR, Q = zip(xs, ys)
    kwargs = {'fontsize': 20}
    dx, dy = -.1, .3
    plt.annotate(r'$Q$', xy=Q, xytext=(Q[0]+dx, Q[1]+dy), **kwargs)
    plt.annotate(r'$R$', xy=R, xytext=(R[0]-.15, R[1]+dy), **kwargs)
    plt.annotate(r'$Q*R$', xy=QR, xytext=(QR[0]-.35, QR[1]-.9), **kwargs)

    # Find intersection points of P*(Q*R)
    dx = (P[1]-QR[1])/(P[0]-QR[0])
    ex = QR[1] - dx*QR[0]
    plt.plot((xl, xr), (dx*xl+ex, dx*xr+ex), zorder=1)
    x, y = sympy.symbols('x y')
    xs = sympy.solve(sympy.Eq(x*x*x + a*x + b, (dx*x + ex)**2))
    xs = [i.as_real_imag()[0] for i in xs]
    ys = [dx*x + ex for x in xs]
    plt.scatter(xs[:2], ys[:2], zorder=3, color='k')
    P, QR, PQR = zip(xs, ys)
    kwargs = {'fontsize': 20}
    plt.annotate(r'$P$', xy=P, xytext=(P[0]-.1, P[1]+.3), **kwargs)
    plt.annotate(r'$P*(Q*R)$', xy=PQR, xytext=(PQR[0]-1.5, PQR[1]-.8), **kwargs)
    plt.scatter(PQR[0], PQR[1], zorder=2, color='k')


def elliptic_curve_example_addition():
    a, b = -2, 2

    xl, xr = -1.9, 3
    yl, yr = -5, 5
    plt.xlim(xl, xr)
    plt.ylim(yl, yr)

    y, x = np.ogrid[yl:yr:100j, xl:xr:100j]

    plt.contour(x.ravel(), y.ravel(), EC(x, y, a, b), [0])

    def f(x): return np.sqrt(x*x*x + a*x + b)
    P_x = .5
    P = (P_x, f(P_x))
    Q_x = -1.5
    Q = (Q_x, -f(Q_x))

    # Find intersection points
    # equation of the line: y = (p2[1]-p1[1])/(p2[0]-p1[0]) x + ex
    dx = (Q[1]-P[1])/(Q[0]-P[0])
    ex = P[1] - dx*P[0]
    plt.plot((xl, xr), (dx*xl+ex, dx*xr+ex))
    # equation of the curve y^2 = x^3 + ax + b
    # Plugging in
    # equation of the curve (dx*x + ex)^2 = x^3 + ax + b
    x, y = sympy.symbols('x y')
    xs = sympy.solve(sympy.Eq(x*x*x + a*x + b, (dx*x + ex)**2))
    xs = [i.as_real_imag()[0] for i in xs]
    ys = [dx*x + ex for x in xs]
    plt.scatter(xs, ys, color='k')

    Q, P, PQ = zip(xs, ys)
    kwargs = {'fontsize': 20}
    dx, dy = -.1, .3
    plt.annotate(r'$P$', xy=P, xytext=(P[0]+dx, P[1]+dy), **kwargs)
    plt.annotate(r'$Q$', xy=Q, xytext=(Q[0]+dx, Q[1]+dy), **kwargs)
    plt.annotate(r'$P*Q$', xy=PQ, xytext=(PQ[0]-0.9, PQ[1]+0.1), **kwargs)
    plt.scatter(PQ[0], -PQ[1], color='k')
    plt.plot((PQ[0], PQ[0]), (yl, yr), color='k', zorder=1, linestyle='--')
    plt.annotate(r'$P+Q$', xy=PQ, xytext=(PQ[0]-1.0, -PQ[1]-.6), **kwargs)

def elliptic_curve_example_doubling():
    a, b = -2, 2

    xl, xr = -1.9, 3
    yl, yr = -5, 5
    plt.xlim(xl, xr)
    plt.ylim(yl, yr)

    y, x = np.ogrid[yl:yr:100j, xl:xr:100j]

    plt.contour(x.ravel(), y.ravel(), EC(x, y, a, b), [0])

    # Find a randon point on the curve
    x = 1.5
    y = np.sqrt(x*x*x + a*x + b)

    # Find intersection points of tangent
    dx = (3*x*x + a) / (2*y)
    ex = y - dx*x

    plt.plot((xl, xr), (dx*xl+ex, dx*xr+ex))
    # Plugging in
    # equation of the curve (dx*x + ex)^2 = x^3 + ax + b
    x, y = sympy.symbols('x y')
    xs = sympy.solve(sympy.Eq(x*x*x + a*x + b, (dx*x + ex)**2))
    xs = [i.as_real_imag()[0] for i in xs]
    ys = [dx*x + ex for x in xs]
    plt.scatter(xs, ys, color='k')

    PP, P = zip(xs, ys)
    kwargs = {'fontsize': 20}
    plt.annotate(r'$P$', xy=P, xytext=(P[0]+.05, P[1]-.6), **kwargs)
    plt.annotate(r'$P*P$', xy=PP, xytext=(PP[0]+.1, PP[1]-.7), **kwargs)
    plt.scatter(PP[0], -PP[1], color='k')
    plt.plot((PP[0], PP[0]), (yl, yr), color='k', zorder=1, linestyle='--')
    plt.annotate(r'$P+P$', xy=PP, xytext=(PP[0]+.1, -PP[1]+.2), **kwargs)


if __name__ == '__main__':
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "Helvetica"
    })
    plt.style.use('grayscale')

    plots = [
        (elliptic_curve_addition, False, 7/16),
        (elliptic_curve_doubling, False, 7/16),
        (elliptic_curve_composition_1, False, 7/16),
        (elliptic_curve_composition_2, False, 7/16),
        (elliptic_curve_example_addition, True, 7/16),
        (elliptic_curve_example_doubling, True, 7/16),
    ]
    for (plot, include_axis, aspect_ratio) in plots:
        # Set the kwargs
        plt.axis('on' if include_axis else 'off')
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.gca().set_aspect(aspect_ratio)


        # Generate the plot
        plot()
        # Do the file name structuring
        fn = f'../build_plots/{plot.__name__}'
        os.system(f'mkdir -p ../build_plots')
        os.system(f'rm {fn}.svg')
        os.system(f'rm {fn}.pdf')
        os.system(f'rm {fn}.pdf_tex')
        # Save the figure
        plt.savefig(f'{fn}.svg', bbox_inches='tight', pad_inches=0)
        plt.cla()
        plt.clf()
        # Convert the figure to LaTeX material
        os.system(f'inkscape {fn}.svg -o {fn}.pdf --export-latex')

    os.system(f'inkscape leiden_logo.svg -o ../build_plots/leiden_logo.pdf --export-latex')
