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

"""
Chinese remainder for constant moduli.
"""

class CRT(object):
    def __init__(self, moduli):
        self.moduli = moduli
        self.P = reduce(lambda a, b: a*b, moduli)
        self.invppmds = []
        for md in self.moduli:
            pp = self.P // md
            inv = mul_inv(pp, md)
            self.invppmds.append(inv * md)

    def remainder(self, remainders):
        assert len(remainders) == len(self.moduli)
        s = 0
        for rem, invppmd in zip(remainders, self.invppmds):
            s += rem * invppmd
        res = s % self.P
        return res



if __name__ == "__main__":

    def testCRT1():
        moduli = (2,3,7)
        crt1 = CRT(moduli)
        remainders = (2,3,2)
        res = crt1.remainder(remainders)
        assert res == 23

    testCRT1()

