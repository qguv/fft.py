#!/usr/bin/env python3

from fft import __doc__ as fft_doc

import subprocess


def example(cmd):
    return f'''```
$ ./fft.py {' '.join(cmd)}
{subprocess.check_output(['./fft.py', *cmd]).decode('utf-8').rstrip()}
```'''


print(f'''# fast fourier transform

{fft_doc.strip()}

## evaluate

{example(['evaluate', '1', '1', '-1', '-1'])}
{example(['evaluate', 'i', 'i', 'i', 'i'])}
{example(['--polar', '--verbose', 'evaluate', '1', '1', '-1', '-1'])}

## interpolate

{example(['interpolate', '0', '-2-2i', '0', '-2+2i'])}
{example(['--verbose', 'interpolate', '4i', '0', '0', '0'])}
{example(['--polar', '--verbose', 'interpolate', '4i', '0', '0', '0'])}

## references

- https://faculty.sites.iastate.edu/jia/files/inline-files/polymultiply.pdf
- https://www.cs.cmu.edu/afs/cs/academic/class/15451-s10/www/lectures/lect0423.txt
- https://www.youtube.com/watch?v=h7apO7q16V0
- https://www.geeksforgeeks.org/fast-fourier-transformation-poynomial-multiplication/
''')
