""" Profile phe_apl.py """
import cProfile
from phe_apl import generateKeysPaillierScheme1, mpz

def testAdditionPerf(m, nBitSize):
    print("testAdditionPerf", m , nBitSize)
    # Compute sum(1..m) via encryption and check the decryption for correctness.
    m = mpz(m)
    total = (1 + m) * m // 2 # see https://en.wikipedia.org/wiki/Triangular_number
    pub, prv = generateKeysPaillierScheme1(nBitSize=nBitSize, useSecureRandom=False)
    fac = pub.n // total # get many non zero bits to be encrypted
    totalEnc = 1
    for i in range(1,m+1):
        enc = pub.encrypt(i * fac) # runtime is dominated by powmod call in encrypt()
        totalEnc = (totalEnc * enc) % pub.nSquared
        decTotal = prv.decrypt(totalEnc) # outside loop to test encryption only.
    assert decTotal == (total * fac)
    print("testAdditionPerf passed")

#cProfile.run('testAdditionPerf(m=100000, nBitSize=1024)')
testAdditionPerf(m=100, nBitSize=2048)
