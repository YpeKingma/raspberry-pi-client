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

""" Chinese remainder for constant moduli. """

from functools import reduce
from operator import __add__, __mul__

def mult_inv(a, m):
    """ Return b such that (a * b) % m == 1 """
    # pow(a, phi(m), m) == 1, for a != 0
    # so: b = pow(a, phi(m)-1, m)
    # for the moment assume m is a prime with phi(m) = m - 1
    b = pow(a, m-2, m)
    assert (a * b) % m == 1
    return b

class CRT(object):
    def __init__(self, moduli):
        self.moduli = moduli
        print("CRT moduli", moduli)
        self.P = reduce(__mul__, moduli)
        def invpp(md):
            pp = self.P // md
            inv = mult_inv(pp, md)
            return inv * pp
        self.invpps = map(invpp, moduli)

    def remainder(self, remainders):
        assert len(remainders) == len(self.moduli)
        s = reduce(__add__, [(rem * invpp) for rem, invpp in zip(remainders, self.invpps)])
        print("remainder sum", s)
        res = s % self.P
        print("remainder", res)
        return res


if __name__ == "__main__":

    def testCRT(moduli, remainders, expected):
        crt1 = CRT(moduli)
        res = crt1.remainder(remainders)
        assert res == expected

    testCRT((3,5), (2,3), 8)

    testCRT((3,5,7), (2,3,2), 23)
