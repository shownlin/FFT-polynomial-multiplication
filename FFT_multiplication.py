import cmath


def recursive_FFT(a):
    n = len(a)
    if n == 1:
        return a
    w_n = cmath.exp((2 * cmath.pi * 1j) / n)
    w = 1
    a0 = a[::2]
    a1 = a[1::2]
    y0 = recursive_FFT(a0)
    y1 = recursive_FFT(a1)
    y = [0] * n   # y is a list that initialize to n zeros
    for k in range(n // 2):
        y[k] = y0[k] + w * y1[k]
        y[k + n // 2] = y0[k] - w * y1[k]
        w = w * w_n
    return y


def inv_FFT(y):
    def inv_FFT_subroutine(y):
        n = len(y)
        if n == 1:
            return y
        w_n = cmath.exp((2 * cmath.pi * 1j) / n)
        w = 1
        y0 = y[::2]
        y1 = y[1::2]
        a0 = inv_FFT_subroutine(y0)
        a1 = inv_FFT_subroutine(y1)
        a = [0] * n
        for k in range(n // 2):
            a[k] = a0[k] + a1[k] / w
            a[k + n // 2] = a0[k] - a1[k] / w
            w = w * w_n
        return a
    n = len(y)
    a_nTimes = inv_FFT_subroutine(y)
    a = list()
    for i in range(len(a_nTimes)):
        a += [round(a_nTimes[i].real / n)]
    return a


def multiply(a, b):
    degree = len(a) + len(b) - 1
    n = max(len(a), len(b))
    if n == 0:
        return 0
    padding = 2 ** (n - 1).bit_length()
    a = a + [0] * (2 * padding - len(a))
    b = b + [0] * (2 * padding - len(b))
    c = list()
    for ai, bi in zip(recursive_FFT(a), recursive_FFT(b)):
        c = c + [ai * bi]
    return inv_FFT(c)[:degree]


if __name__ == "__main__":
    a = [9, -10, 7, 6, 21]  # 9 - 10x + 7x^2 + 6x^3 + 21x^4
    b = [-5, 4, 5, -2]  # -5 + 4x + 5x^2 - 2x^3
    c = multiply(a, b)
    print(c)    # output: [-45, 86, -30, -70, -26, 100, 93, -42]
