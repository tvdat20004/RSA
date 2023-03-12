import sys
from sympy.solvers import solve
from sympy import Symbol

# A reimplementation of pablocelayes rsa-wiener-attack for this purpose
# https://github.com/pablocelayes/rsa-wiener-attack/


class WienerAttack(object):
    def rational_to_contfrac(self, x, y):
        a = x // y
        if a * y == x:
            return [a]
        else:
            pquotients = self.rational_to_contfrac(y, x - a * y)
            pquotients.insert(0, a)
            return pquotients

    def convergents_from_contfrac(self, frac):
        convs = []
        for i in range(len(frac)):
            convs.append(self.contfrac_to_rational(frac[0:i]))
        return convs

    def contfrac_to_rational(self, frac):
        if len(frac) == 0:
            return (0, 1)
        elif len(frac) == 1:
            return (frac[0], 1)
        else:
            remainder = frac[1:len(frac)]
            (num, denom) = self.contfrac_to_rational(remainder)
            return (frac[0] * num + denom, num)

    def is_perfect_square(self, n):
        h = n & 0xF
        if h > 9:
            return -1

        if (h != 2 and h != 3 and h != 5 and h != 6 and h != 7 and h != 8):
            t = self.isqrt(n)
            if t*t == n:
                return t
            else:
                return -1

        return -1

    def isqrt(self, n):
        if n == 0:
            return 0
        a, b = divmod(n.bit_length(), 2)
        x = 2**(a+b)
        while True:
            y = (x + n//x)//2
            if y >= x:
                return x
            x = y

    def __init__(self, n, e):
        self.d = None
        self.p = None
        self.q = None
        sys.setrecursionlimit(100000)
        frac = self.rational_to_contfrac(e, n)
        convergents = self.convergents_from_contfrac(frac)

        for (k, d) in convergents:
            if k != 0 and (e * d - 1) % k == 0:
                phi = (e * d - 1) // k
                s = n - phi + 1
                discr = s*s - 4*n
                if(discr >= 0):
                    t = self.is_perfect_square(discr)
                    if t != -1 and (s + t) % 2 == 0:
                        self.d = d
                        x = Symbol('x')
                        roots = solve(x**2 - s * x + n, x)
                        if len(roots) == 2:
                            self.p = roots[0]
                            self.q = roots[1]
                        break

def power(a, b, mod):
    if b == 0:
        return 1 % mod
    elif b % 2 == 0:
        p = power(a, b // 2, mod)
        return (p * p) % mod
    else:
        p = power(a, (b - 1) // 2, mod)
        return (a * p * p) % mod

def int_to_ascii(m):
    # Decode to ascii (from https://crypto.stackexchange.com/a/80346)
    m_hex = hex(int(m))[2:-1]  # Number to hex
    m_ascii = "".join(
        chr(int(m_hex[i : i + 2], 16)) for i in range(0, len(m_hex), 2)
    )  # Hex to Ascii
    return m_ascii

def main():
    n = 126195321438426325889716627116088211916808669792474407440301779806981435873626284541866710664632373089344843582167312068425281144271868756149171245603482112363524007097935192482405928918584664493142854348593840207198282355637714006653469894958550331395150279236267780381821979466426457998575877602789118135451
    e = 113721594484292381101826547460907850536435973795939197203008793985549861693675655856072986559461351931681920051422877617859944776259898462135910449234020538104928819334754254754461756250472605284429998036747119417402737943057270490132499150033715124734786185529743187845603122183408394056376075188510084527461
    c = 79807729753785550021039570684651802991350606986676965714016962649149395811147460108821090484456996575067504095930131563014131211831122454321692125962185192201466793775632638270804580588643092468950724240817804726029151604950481062839759281331028272334399986680526884946477089349267392158703235810567870373951
    Wiener = WienerAttack(n,e)
    d = Wiener.d
    p = Wiener.p
    q = Wiener.q
    m = power(c,d,n)
    #print(m)
    m_ascii = int_to_ascii(m)
    print("Flag: %s" % (m_ascii+'}').strip())


if __name__ == '__main__':
    main()