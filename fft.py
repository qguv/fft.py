#!/usr/bin/env python3
'''
fast fourier transform in a polynomial context
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

    q_odd = _fft(p[::2], w)
    q_even = _fft(p[1::2], w)

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
        cs.insert(0, 0)
        degree += 1
    w = math.e ** (2 * math.pi * 1j / degree)
    return _fft(cs, w)


def ifft(ys: ValuesAtNthRootsOfUnity) -> Polynomial:
    degree = len(ys)
    assert is_power_2(degree), 'number of points must be a power of 2'
    w = (math.e ** (-2 * math.pi * 1j / degree)) / degree
    return _fft(ys, w)


################################
# end of the interesting part; #
# the rest just handles i/o    #
################################


def fmt_values(ys: ValuesAtNthRootsOfUnity):
    degree = len(ys)
    w = math.e ** (2 * math.pi * 1j / degree)
    for i, y in enumerate(ys):
        yield f'(ω{superscript(i)} = {fmtc(w ** i)}, {fmtc(y)})'


def superscript(s):
    return str(s).translate(str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹"))


def fmt_polynomial(cs):
    max_power = len(cs) - 1
    yield 'P(x) ='
    first = True
    for i, c in enumerate(cs):
        if c == 0:
            continue
        if first:
            first = False
            op = ''
        elif c.real <= 0:
            op = '- '
            c = -c
        else:
            op = '+ '
        power = max_power - i
        factor = f'x{superscript(power)}' if power else ''
        c_fmt = fmtc(c, polar=False, parens_on_complex=(power != 0))
        yield f'{op}{c_fmt}{factor}'


def fmtc(x, polar=True, precision=2, parens_on_complex=False) -> complex or float:
    if type(x) is complex:
        if x.imag == 0:
            return fmtc(x.real, precision=precision)
        if polar:
            r, phase = cmath.polar(x)
            s = f'{fmtc(r, precision=precision)} ∠ {fmtc(math.degrees(phase), precision=precision)}°'
        elif x.imag < 0:
            s = f'{fmtc(x.real, precision=precision)}-{fmtc(-x.imag, precision=precision)}i'
        else:
            s = f'{fmtc(x.real, precision=precision)}+{fmtc(x.imag, precision=precision)}i'
        if parens_on_complex:
            s = f'({s})'
        return s

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
    cs = [float(c_raw) for c_raw in [args.A, *args.B]]
    if args.verbose:
        print("evaluating polynomial:", *fmt_polynomial(cs))
        print('points:')
    ys = fft(cs)
    print(*fmt_values(ys), sep='\n')


def cmd_interpolate(args):
    ys = [float(y) for y in [args.Y0, *args.Y1]]
    if args.verbose:
        print("interpolating points:", *fmt_values(ys), sep='\n')
        print("polynomial: ", end='')
    cs = ifft(ys)
    print(*fmt_polynomial(cs))


if __name__ == '__main__':
    parse = {'': argparse.ArgumentParser(description=__doc__)}
    parse[''].set_defaults(func=lambda x: parse[''].print_help())
    parse[''].add_argument('--verbose', action='store_true')
    parse['_'] = parse[''].add_subparsers()

    msg = (
        'FFT: evaluate a polynomial of degree N (given as a list of coefficients) when x is set to each of the Nth'
        ' roots of unity'
    )
    parse['evaluate'] = parse['_'].add_parser('evaluate', help=msg, description=msg)
    parse['evaluate'].set_defaults(func=cmd_evaluate)
    parse['evaluate'].add_argument('A', help='coefficient of the x^N term')
    parse['evaluate'].add_argument('B', help='coefficient of the x^(N-1) term, and so on', nargs='*')

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
        nargs='*',
    )

    args = parse[''].parse_args()
    args.func(args)
