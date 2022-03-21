#!/usr/bin/env python3
'''
The fast fourier transform (FFT) algorithm applied to polynomial evaluation and interpolation.
'''

import argparse
import cmath
import math
import typing

Coefficient = typing.NewType('Coefficient', complex)

# an N-item list of coefficients representing a degree-N polynomial
# e.g. [3, 2, 1] would represent the polynomial 3x² + 2x¹ + 1
Polynomial = typing.NewType('Polynomial', list[Coefficient])

# an N-item list of each value of a polynomial when x is set to each Nth root
# of unity
ValuesAtNthRootsOfUnity = typing.NewType('ValuesAtNthRootsOfUnity', list[complex])


def _fft(p: [complex], w) -> complex:
    degree = len(p)
    if degree == 1:
        return p

    q_even = _fft(p[::2], w)
    q_odd = _fft(p[1::2], w)

    q = [0] * degree
    half_degree = int(degree // 2)
    for i in range(half_degree):
        odd_term = (w ** i) * q_odd[i]
        q[i] = q_even[i] + odd_term
        q[i + half_degree] = q_even[i] - odd_term
    return q


def is_power_2(x):
    return x & (x - 1) == 0


def fft(cs: Polynomial) -> ValuesAtNthRootsOfUnity:
    degree = len(cs)
    while not is_power_2(degree):
        cs.append(0)
        degree += 1
    w = math.e ** (2 * math.pi * 1j / degree)
    return _fft(cs, w)


def ifft(ys: ValuesAtNthRootsOfUnity) -> Polynomial:
    degree = len(ys)
    assert is_power_2(degree), 'number of points must be a power of 2'
    w = math.e ** (-2 * math.pi * 1j / degree)
    return [y / degree for y in _fft(ys, w)]


################################
# end of the interesting part; #
# the rest just handles i/o    #
################################


def fmt_values(ys: ValuesAtNthRootsOfUnity, polar=False):
    degree = len(ys)
    w = math.e ** (2 * math.pi * 1j / degree)
    for i, y in enumerate(ys):
        yield f'(ω{superscript(i)} = {fmtc(w ** i, polar=polar)}, {fmtc(y, polar=polar)})'


def superscript(s):
    return str(s).translate(str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹"))


def fmt_polynomial(cs):
    yield 'P(𝑥) ='
    first = True
    for i, c in list(enumerate(cs))[::-1]:
        c = normalize(c)
        if c == 0:
            continue
        if first:
            first = False
            op = ''
        elif c.real < 0 or (c.real == 0 and c.imag < 0):
            op = '- '
            c = -c
        else:
            op = '+ '
        factor = f'𝑥{superscript(i)}' if i else ''
        if i and c == 1:
            c_fmt = ''
        else:
            c_fmt = fmtc(c, parens_on_complex=(i != 0))
        yield f'{op}{c_fmt}{factor}'


def normalize(x, thresh=1e-15):
    if type(x) is complex:
        return complex(normalize(x.real), normalize(x.imag))
    n = round(x)
    if abs(x - n) < thresh:
        return n
    return x


def fmtc(x, polar=False, precision=2, parens_on_complex=False) -> complex or float:
    if type(x) is complex:
        x = normalize(x)

        if polar:
            r, phase = cmath.polar(x)
            return f'({fmtc(r, precision=precision)} ∠ {fmtc(math.degrees(phase), precision=precision)}°)'

        if x.imag == 0:
            return fmtc(x.real, precision=precision)

        if x.real == 0:
            if x.imag == 1:
                return '𝑖'
            elif x.imag == -1:
                return '-𝑖'
            elif x.imag < 0:
                return f'-{fmtc(-x.imag, precision=precision)}𝑖'
            else:
                return f'{fmtc(x.imag, precision=precision)}𝑖'

        if x.imag < 0:
            s = f'{fmtc(x.real, precision=precision)}-{fmtc(-x.imag, precision=precision)}𝑖'
        else:
            s = f'{fmtc(x.real, precision=precision)}+{fmtc(x.imag, precision=precision)}𝑖'

        return f'({s})' if parens_on_complex else s

    if type(x) is float:
        if x.is_integer():
            return fmtc(int(x))
        x_round = round(x, precision)
        if x_round == x:
            return str(x)
        trunc = '{:.' + str(precision) + 'f}…'
        return trunc.format(x_round)

    return str(x)


def cmd_evaluate(args):
    cs = [complex(c_raw.replace('i', 'j')) for c_raw in [args.A, *args.B]][::-1]
    if args.verbose:
        print("evaluating polynomial:", *fmt_polynomial(cs))
        print('points:')
    ys = fft(cs)
    print(*fmt_values(ys, polar=args.polar), sep='\n')


def cmd_interpolate(args):
    ys = [complex(y.replace('i', 'j')) for y in [args.Y0, *args.Y1]]
    if args.verbose:
        print("interpolating points:", *fmt_values(ys, polar=args.polar), sep='\n')
        print("polynomial: ", end='')
    cs = ifft(ys)
    print(*fmt_polynomial(cs))


def gen_parser():
    parse = {'': argparse.ArgumentParser(description=__doc__)}
    parse[''].set_defaults(func=lambda x: parse[''].print_help())
    parse[''].add_argument('--verbose', action='store_true')
    parse[''].add_argument('--polar', action='store_true', help='show points in polar format')
    parse['_'] = parse[''].add_subparsers()

    msg = (
        'FFT: evaluate a polynomial of degree N (given as a list of coefficients) when x is set to each of the Nth'
        ' roots of unity'
    )
    parse['evaluate'] = parse['_'].add_parser('evaluate', help=msg, description=msg)
    parse['evaluate'].set_defaults(func=cmd_evaluate)
    parse['evaluate'].add_argument('A', help='coefficient of the x^N term')
    parse['evaluate'].add_argument('B', help='coefficient of the x^(N-1) term, and so on', nargs=argparse.REMAINDER)

    msg = (
        'inverse FFT: calculate the coefficients of a polynomial of degree N given each of the values of the'
        ' polynomial when x is set to each of the Nth roots of unity'
    )
    parse['interpolate'] = parse['_'].add_parser('interpolate', help=msg, description=msg)
    parse['interpolate'].set_defaults(func=cmd_interpolate)
    parse['interpolate'].add_argument('Y0', help='value of the polynomial at x=ω⁰=1, the first Nth root of unity')
    parse['interpolate'].add_argument(
        'Y1',
        help='value of the polynomial at x=ω¹, the second Nth root of unity, etc.',
        nargs=argparse.REMAINDER,
    )

    return parse['']


if __name__ == '__main__':
    args = gen_parser().parse_args()
    args.func(args)
