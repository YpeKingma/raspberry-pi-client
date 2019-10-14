# -*- coding: utf-8 -*-

# Copyright 2019 Ype Kingma
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#     http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function # for python2 development only

"""  Homomorphic encryption as defined in this paper by Pascal Paillier (in bibtex format):
@inproceedings{paillier1999public,
  title={Public-key cryptosystems based on composite degree residuosity classes},
  author={Paillier, Pascal},
  booktitle={International Conference on the Theory and Applications of Cryptographic Techniques},
  pages={223--238},
  year={1999},
  organization={Springer}
}
Paillier 1999 is used to refer to this paper here.

The is an implementation of Scheme 1 of Paillier 1999.
The implemention was developed for use with python2 and python3.

Other references:

Paillier and Pointcheval 1999:
@inproceedings{paillier1999efficient,
  title={Efficient public-key cryptosystems provably secure against active adversaries},
  author={Paillier, Pascal and Pointcheval, David},
  booktitle={International Conference on the Theory and Application of Cryptology and Information Security},
  pages={165--179},
  year={1999},
  organization={Springer}
}

Koblitz:
Neal Koblitz, A course in number theory and cryptography, Springer, 1994, 2nd ed.

@inproceedings{rivest1999arestrong,
  title={AreStrong'Primes Needed for RSA?},
  author={Rivest, Ronald L and Silverman, Robert D},
  booktitle={IN THE 1997 RSA LABORATORIES SEMINAR SERIES, SEMINARS PROCEEDINGS},
  year={1999}
}
Referred to as Rivest 1999.

@inproceedings{vsvenda2016million,
  title={The Million-Key Question—Investigating the Origins of $\{$RSA$\}$ Public Keys},
  author={{\v{S}}venda, Petr and Nemec, Mat{\'u}{\v{s}} and Sekan, Peter and Kva{\v{s}}{\v{n}}ovsk{\`y}, Rudolf and Form{\'a}nek, David and Kom{\'a}rek, David and Maty{\'a}{\v{s}}, Vashek},
  booktitle={25th $\{$USENIX$\}$ Security Symposium ($\{$USENIX$\}$ Security 16)},
  pages={893--910},
  year={2016}
}
Referred to as Svenda 2016.

List of TBD:

Implement the fast variant of Paillier, Scheme 3, and follow section 7
on efficiency and implementation aspects. Scheme 3 uses a shorter prime (160 bit) alpha,
1 <= alpha <= lmbda (Pailler, top of p. 232, end of section 6, p. 233).
See also the table on p. 235 for the expected speeds of the different schemes.

Make this a pip installable library.

More possible performance optimizations are in the comments.
"""

import random
import numbers

try:
    import secrets
    secureRandom = secrets.SystemRandom()
except ImportError: # before python 3.6
    secureRandom = random.SystemRandom()

try:
    # raise ImportError # simulate gmpy2 not available
    import gmpy2

    from gmpy2 import mpz
    from gmpy2 import powmod
    from gmpy2 import gcd

    from gmpy2 import random_state, mpz_random, mpz_urandomb

    gmpy_random_state = random_state(random.randrange(2 ** 63))

    def randomIntBitSize(bitSize):
        p2 = mpz(2) ** (bitSize - 1)
        res = p2 + mpz_urandomb(gmpy_random_state, bitSize - 1)
        assert res.bit_length() == bitSize
        return res

    def secureRandomIntBitSize(bitSize): # avoid gmpy2 mpz_urandomb
        p2 = 2 ** (bitSize - 1)
        res = secureRandom.randrange(p2, 2 * p2)
        assert res.bit_length() == bitSize
        return mpz(res)

    from gmpy2 import is_prime
    def isProbablePrime(n, trials=25):
        return is_prime(n, trials)

    from gmpy2 import next_prime
    nextPrime = next_prime #  """ For given q, return the smallest p > q that is probable prime """

    from gmpy2 import lcm

    from gmpy2 import gcdext
    # gcdext(a, b) returns a 3-element tuple (g, s, t) such that
    # g == gcd(a, b) and g == a * s + b * t

except ImportError: # gmpy2 not available

    def mpz(n): # to avoid duplicate code
        return n

    def _randIntBitSize(bitSize, rnd):
        p2 = 2 ** (bitSize - 1)
        res = rnd.randrange(p2, 2 * p2)
        assert res.bit_length() == bitSize
        return res

    def randomIntBitSize(bitSize):
        return _randIntBitSize(bitSize, random)

    def secureRandomIntBitSize(bitSize):
        return _randIntBitSize(bitSize, secureRandom)

    powmod = pow # only with non negative powers when modulo argument is present.

    try:
        from math import gcd
    except ImportError: # before python 3.5
        from fractions import gcd
    except ImportError: # before python 2.6
        def gcd(a, b): # Warning: untested.
            """Calculate the Greatest Common Divisor of a and b.
            Unless b==0, the result will have the same sign as b (so that when
            b is divided by it, the result comes out positive).
            """
            while b:
                a, b = b, (a % b)
            return a

    _smallPrimesList =  [
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
        59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127,
        131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
        197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269,
        271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
        353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431,
        433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503,
        509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599,
        601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673,
        677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761,
        769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857,
        859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947,
        953, 967, 971, 977, 983, 991, 997]

    _smallPrimesDict = {}
    for p in _smallPrimesList:
        _smallPrimesDict[p] = 1

    def isProbablePrime(n, trials=8, rnd=random):
        """ Return value indicates that n is probably a prime.
        The probability that True is returned for a non prime n is less than 4**(-trials).
        See also https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test .
        """
        if n < 2:
            return False
        if n in _smallPrimesDict:
            return True
        for p in _smallPrimesList:
            if n % p == 0:
                return False
        nMinus1 = n - 1
        r = 0
        d = nMinus1
        while d % 2 == 0:
            d >>= 1
            r += 1

        def checkWitness(a): # for n composite
            x = powmod(a, d, n)
            if x == 1 or x == nMinus1:
                return False # continue witness loop
            for _ in range(r - 1):
                x = powmod(x, 2, n)
                if x == nMinus1:
                    return False # continue witness loop
            return True

        for _ in range(trials): # witness loop
            if checkWitness(rnd.randrange(2, nMinus1)):
                return False
        return True


    def nextPrime(q):
        """ Return the smallest p > q that is probable prime """
        p = q | 1
        if p == q:
            p += 2
        while True:
            if isProbablePrime(p):
                return p
            p += 2

    def lcm(a, b):
        return a * b // gcd(a, b)



def L(u, n): # See Paillier, p. 227, and https://en.wikipedia.org/wiki/Paillier_cryptosystem
    return (u - 1) // n # floor division
    # TBD: Paillier, p. 233, Decryption, recommends to speed this up by precomputing
    # n^-1 mod 2^bitSize(n) # or abs(n) ? Paillier writes |n| here.


class PaillierScheme1PublicKey(object): # Paillier 1999, p. 229
    def __init__(self, n, nSquared, g, rnd=secureRandom):
        self.n = mpz(n)
        self.nSquared = mpz(nSquared)
        self.g = mpz(g)
        r = mpz(rnd.randrange(5, self.n)) # avoid r <= 4
        self.rpown = powmod(r, n, self.nSquared)

    def encrypt(self, m):
        """ Encrypt plaintext m using Scheme 1. """
        assert m >= 0 and m < self.n
        return (powmod(self.g, m, self.nSquared) * self.rpown) % self.nSquared


class PaillierScheme1PrivateKey(object):
    # def __init__(self, n, nSquared, g, lmbda, phiN):
    def __init__(self, n, nSquared, g, lmbda, phiN, p, q):
        self.n = mpz(n)
        self.nSquared = mpz(nSquared)
        self.lmbda = mpz(lmbda)
        gl = powmod(mpz(g), self.lmbda, self.nSquared)
        lgl = L(gl, self.n)
        # Inverse modulo not yet available in python, i.e. pow(, -1, n).
        # Since pow(i, phi(m), m) == 1, for i != 0
        # (see https://www.algorithmist.com/index.php/Modular_inverse )
        # self.mu = inverseLgl = pow(lgl, phi(n) - 1, n)
        self.mu = powmod(lgl, mpz(phiN - 1), self.n) # precomputed constant, see Paillier 1999, top of p. 234.
        assert (self.mu * lgl) % n == 1

        # Precomputations for CRT, Paillier p. 234.
        self.p = mpz(p)
        self.pSquared = p ** 2
        gpm1p2 = powmod(g, p - 1, self.pSquared)
        lgp = L(gpm1p2, p)
        self.hp = powmod(lgp, p - 2, p)
        assert (self.hp * lgp) % p == 1

        self.q = mpz(q)
        self.qSquared = q ** 2
        gqm1q2 = powmod(g, q - 1, self.qSquared)
        lgq = L(gqm1q2, q)
        self.hq = powmod(lgq, q - 2, q)
        assert (self.hq * lgq) % q == 1


    def decrypt(self, c):
        """ Decrypt an encrypted number, using Scheme 1. """
        assert c >= 0 and c < self.nSquared

        # Paillier 1999, p. 234, Decryption using Chinese-remaindering
        # to speed decryption up by splitting up over p and q :
        cpm1p2 = powmod(c, self.p - 1, self.pSquared)
        lcp = L(cpm1p2, self.p)
        mp = (lcp * self.hp) % self.p

        cqm1q2 = powmod(c, self.q - 1, self.qSquared)
        lcq = L(cqm1q2, self.q)
        mq = (lcq * self.hq) % self.q

        (gcdmpmq, s, t) = gcdext(mp, mq)
        return (mp * s + mq * t) % self.n # CRT chinese remainder


def generateKeysPaillierScheme1(nBitSize, useSecureRandom=True):
    """ Return a tuple of (PaillierPublicKey, PaillierPrivateKey) instances with n of the given size.
        When useSecureRandom True, the python SystemRandom module is used for generating p, q and g.
        Otherwise the python random module is used.
    """
    assert nBitSize > 210 # Allow p and q at least 105 bits, p and q at least 104, see Fermat factorization below.
    assert nBitSize <= 8192 # Depends on available processing speed.

    # Choose a random prime number p of nBitSize//2, and q of the remaining bitsize for nBitSize.

    # See Koblitz chapter IV.2 RSA on choosing p and q: differ by a few decimal digits,
    # a weaker check is done below to avoid easy Fermat factorization,
    # see https://en.wikipedia.org/wiki/Fermat%27s_factorization_method

    # See also Rivest 1999: safe (Germain) primes do not hurt, but do not bring better protection
    # than large enough primes.
    # Strong primes (with a large factor in p+1) are not considered here,
    # see also https://gmpy2.readthedocs.io/en/latest/advmpz.html for corresponding checks.

    bitSizeP = nBitSize // 2

    while True:
        if useSecureRandom:
            nMin = secureRandomIntBitSize(nBitSize) # to start looking for q after generating p; uniform distribution
            pMin = secureRandomIntBitSize(bitSizeP) # uniform distribution.
        else:
            nMin = randomIntBitSize(nBitSize) # to start looking for q after generating p; uniform distribution
            pMin = randomIntBitSize(bitSizeP) # uniform distribution.
        p = nextPrime(pMin)
        qMin = nMin // p # uses only the upper half of the bits of nMin; the distribution of qMin is not uniform
        q = nextPrime(qMin)

        # See https://en.wikipedia.org/wiki/Fermat%27s_factorization_method
        # and https://crypto.stackexchange.com/questions/5262/rsa-and-prime-difference
        # and https://crypto.stackexchange.com/questions/5698/ansi-x9-31-standards-for-generating-random-numbers
        # p and q should differ in their first 100 bits: abs(p-q) > 2^(nBitSize/2−100)
        # Here p.bit_length() = k//2
        minBitLengthPQ = min(p.bit_length(), q.bit_length())
        assert minBitLengthPQ > 100
        if (abs(p-q) >> (minBitLengthPQ - 100)) == 0: # guard against Fermat factorization, very very unlikely for nBitSize > 250
            continue # retry
        n = p * q
        if n.bit_length() == nBitSize: # could be too large, but very very unlikely.
            break

    lmbda = lcm(p - 1, q - 1)
    phiN = (p - 1) * (q - 1)
    nSquared = n ** 2

    # Paillier uses g as an element of Z*nSquared (see Paillier, p. 225,  under 3).
    # randomly select a base g from B by verifying eq (4) on p. 229:
    # gcd(L(g^lmbda mod n^2, n), n) = 1
    while True:
        # TBD: Paillier 1999, p. 233 under Encryption, use a small g for encryption efficiency.
        if useSecureRandom:
            g = mpz(secureRandom.randrange(4, nSquared))
        else:
            g = mpz(random.randrange(4, nSquared))
        if gcd(L(powmod(g, lmbda, nSquared), n), n) == 1:
            break

    publicKey = PaillierScheme1PublicKey(n, nSquared, g)
    privateKey = PaillierScheme1PrivateKey(n, nSquared, g, lmbda, phiN, p, q)
    return (publicKey, privateKey)


if __name__ == "__main__":

    def testProbablePrime():
        numBits = 500
        startRange = 2 ** numBits
        rangeSize = 10000
        numFound = 0
        print("startRange", startRange)
        print("rangeSize", rangeSize)
        for i in range(startRange, startRange + rangeSize):
            if isProbablePrime(i):
                numFound += 1
        print("numFound", numFound)

    def testPaillierKeySize(nBitSize, mes, useSecureRandom=True):
        pub, prv = generateKeysPaillierScheme1(nBitSize=nBitSize, useSecureRandom=useSecureRandom)
        print("generated keys, nBitSize", nBitSize)
        print("n", pub.n)
        enc = pub.encrypt(mes)
        dec = prv.decrypt(enc)
        assert mes == dec
        print("testPaillierKeySize passed")

    def testHomomorphic1AddProdPow(m1, m2, pub, prv):
        enc1 = pub.encrypt(m1 % pub.n)
        enc2 = pub.encrypt(m2 % pub.n)
        decProd = prv.decrypt((enc1 * enc2) % pub.nSquared)
        assert decProd == (m1 + m2) % pub.n
        decPow12 = prv.decrypt(powmod(enc1, m2, pub.nSquared))
        decPow21 = prv.decrypt(powmod(enc2, m1, pub.nSquared))
        prod = (m1 * m2) % pub.n
        assert decPow12 == prod
        assert decPow21 == prod
        print("testHomomorphic1AddProdPow passed")

    def test1Blinding(mes, pub, prv):
        r = random.randrange(1, pub.n)
        enc = pub.encrypt(mes % pub.n)
        blinded = (enc * (powmod(r, pub.n, pub.nSquared))) % pub.nSquared
        decBlinded = prv.decrypt(blinded)
        assert decBlinded == mes
        print("test1Blinding passed")


    random.seed(65539)

    testProbablePrime()

    mes = 301
    testPaillierKeySize(211, mes)
    testPaillierKeySize(256, mes)
    testPaillierKeySize(256, mes, useSecureRandom=False)
    testPaillierKeySize(512, mes)
    testPaillierKeySize(1024, mes)
    #testPaillierKeySize(2048, mes)
    #testPaillierKeySize(4096, mes)
    #testPaillierKeySize(8192, mes)

    pub, prv = generateKeysPaillierScheme1(nBitSize=2048)
    testHomomorphic1AddProdPow(mes, mes + 987, pub, prv)
    test1Blinding(mes, pub, prv)
