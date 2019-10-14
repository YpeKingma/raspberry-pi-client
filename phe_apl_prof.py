""" Profile phe_apl.py """
import cProfile
from phe_apl import generateKeysPaillierScheme1, mpz

def testAdditionPerf(m, nBitSize):
    print("testAdditionPerf", m , nBitSize)
    # Compute sum(1..m) via encryption and check the decryption for correctness.
    m = mpz(m)
    total = (1 + m) * m // 2 # see https://en.wikipedia.org/wiki/Triangular_number
    print("total", total)
    pub, prv = generateKeysPaillierScheme1(nBitSize=nBitSize, useSecureRandom=False)
    fac = (pub.n - 1) // total # get many non zero bits to be encrypted
    fac = 2 ** 512 # passes, but test fails for 2 ** 1024, why?
    fac = 2 ** 1024
    print("fac", fac)
    print("total * fac", total * fac)
    totalEnc = 1
    for i in range(1,m+1):
        enc = pub.encrypt(i * fac) # runtime is dominated by powmod call in encrypt()
        totalEnc = (totalEnc * enc) % pub.nSquared
        decTotal = prv.decrypt(totalEnc) # outside loop to test encryption only.
        subtotal = (1 + i) * i // 2
        assert decTotal == subtotal * fac
    print("decTotal", decTotal)
    assert decTotal == (total * fac)
    print("testAdditionPerf passed")

#cProfile.run('testAdditionPerf(m=100000, nBitSize=1024)')
testAdditionPerf(m=100, nBitSize=2048)
