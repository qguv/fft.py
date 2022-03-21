# fast fourier transform

The fast fourier transform (FFT) algorithm applied to polynomial evaluation and interpolation.

## evaluate

```
$ ./fft.py evaluate 1 1 -1 -1
(ω⁰ = 1, 0)
(ω¹ = 𝑖, -2-2𝑖)
(ω² = -1, 0)
(ω³ = -𝑖, -2+2𝑖)
```
```
$ ./fft.py evaluate i i i i
(ω⁰ = 1, 4𝑖)
(ω¹ = 𝑖, 0)
(ω² = -1, 0)
(ω³ = -𝑖, 0)
```
```
$ ./fft.py --polar --verbose evaluate 1 1 -1 -1
evaluating polynomial: P(𝑥) = 𝑥³ + 𝑥² - 𝑥¹ - 1
points:
(ω⁰ = (1 ∠ 0°), (0 ∠ 0°))
(ω¹ = (1 ∠ 90°), (2.83… ∠ -135°))
(ω² = (1 ∠ 180°), (0 ∠ 0°))
(ω³ = (1 ∠ -90°), (2.83… ∠ 135°))
```

## interpolate

```
$ ./fft.py interpolate 0 -2-2i 0 -2+2i
P(𝑥) = 𝑥³ + 𝑥² - 𝑥¹ - 1
```
```
$ ./fft.py --verbose interpolate 4i 0 0 0
interpolating points:
(ω⁰ = 1, 4𝑖)
(ω¹ = 𝑖, 0)
(ω² = -1, 0)
(ω³ = -𝑖, 0)
polynomial: P(𝑥) = 𝑖𝑥³ + 𝑖𝑥² + 𝑖𝑥¹ + 𝑖
```
```
$ ./fft.py --polar --verbose interpolate 4i 0 0 0
interpolating points:
(ω⁰ = (1 ∠ 0°), (4 ∠ 90°))
(ω¹ = (1 ∠ 90°), (0 ∠ 0°))
(ω² = (1 ∠ 180°), (0 ∠ 0°))
(ω³ = (1 ∠ -90°), (0 ∠ 0°))
polynomial: P(𝑥) = 𝑖𝑥³ + 𝑖𝑥² + 𝑖𝑥¹ + 𝑖
```

## references

- https://faculty.sites.iastate.edu/jia/files/inline-files/polymultiply.pdf
- https://www.cs.cmu.edu/afs/cs/academic/class/15451-s10/www/lectures/lect0423.txt
- https://www.youtube.com/watch?v=h7apO7q16V0
- https://www.geeksforgeeks.org/fast-fourier-transformation-poynomial-multiplication/

