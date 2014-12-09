import numpy
import dadi

def modelOne(params, ns, pts):
    s,nu1,nu2,T = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: s * (nu1/s)**(t/T)
    nu2_func = lambda t: (1-s) * (nu2/(1-s)) ** (t/T)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func,
            m12 = 0, m21 = 0)

    fs = dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs


def modelTwo(params, ns, pts):
    s,nu1,nu2,T,m12 = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: s * (nu1/s)**(t/T)
    nu2_func = lambda t: (1-s) * (nu2/(1-s)) ** (t/T)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func,
            m12 = m12, m21 = 0)

    fs = dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs


