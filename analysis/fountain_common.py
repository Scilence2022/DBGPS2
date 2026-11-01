"""DNA Fountain Python-3 companion: exact reimplementation of the parts of
TeamErlich/dna-fountain (python2) needed to verify DBGPS2 output.

- byte <-> DNA: A=00, C=01, G=10, T=11, MSB first (utils.int_to_dna)
- droplet wire format: seed(4B) || payload(32B) || RS(2B), RS = reedsolo RSCodec(2)
- PRNG: robust soliton sampler; python2 random.sample semantics (float multiply,
  not getrandbits) so seeds regenerate the same segment subsets as python2.
"""
import math
import random

DNA_TO_INT = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
INT_TO_DNA = 'ACGT'


def dna_to_bytes(dna):
    bits = ''.join(format(DNA_TO_INT[b], '02b') for b in dna)
    return bytes(int(bits[i:i + 8], 2) for i in range(0, len(bits), 8))


def bytes_to_dna(data):
    bits = ''.join(format(b, '08b') for b in data)
    return ''.join(INT_TO_DNA[int(bits[i:i + 2], 2)] for i in range(0, len(bits), 2))


def seed_of_droplet(dna_block):
    """First 16 nt of the 152-nt data block decode to the 4-byte big-endian seed."""
    return int.from_bytes(dna_to_bytes(dna_block[:16]), 'big')


def gen_tau(S, K, delta):
    pivot = int(math.floor(K / S))
    val1 = [S / K * 1 / d for d in range(1, pivot)]
    val2 = [S / K * math.log(S / delta)]
    val3 = [0.0 for d in range(pivot, K)]
    return val1 + val2 + val3


def gen_rho(K):
    return [1.0 / K] + [1.0 / (d * (d - 1)) for d in range(2, K + 1)]


def gen_rsd_cdf(K, S, delta):
    tau = gen_tau(S, K, delta)
    rho = gen_rho(K)
    Z = sum(rho) + sum(tau)
    mu = [(rho[d] + tau[d]) / Z for d in range(K)]
    cdf = []
    acc = 0.0
    for m in mu:
        acc += m
        cdf.append(acc)
    return cdf, Z


def py2_random_sample(rnd, n, k):
    """python2.7 random.sample(range(n), k) using rnd.random() draws only.
    Both the 'pool' and 'selected set' branches consume randomness the same way
    (one random() per pick, plus none for sorting)."""
    result = [None] * k
    setsize = 21
    if k > 5:
        setsize += 4 ** int(math.ceil(math.log(k * 3, 4)))
    if n <= setsize:
        pool = list(range(n))
        for i in range(k):
            j = int(rnd.random() * (n - i))
            result[i] = pool[j]
            pool[j] = pool[n - i - 1]
    else:
        selected = set()
        selected_add = selected.add
        for i in range(k):
            j = int(rnd.random() * n)
            while j in selected:
                j = int(rnd.random() * n)
            selected_add(j)
            result[i] = j
        result.sort()
    return result


class PRNG:
    """Robust solition PRNG matching robust_solition.pyx (python random path)."""

    def __init__(self, K, delta, c):
        self.K = float(K)
        self.K_int = int(K)
        self.delta = delta
        self.c = c
        S = self.c * math.log(self.K / self.delta) * math.sqrt(self.K)
        self.S = S
        self.cdf, self.Z = gen_rsd_cdf(K, S, delta)
        self.state = 1

    def set_seed(self, seed):
        self.state = seed

    def get_src_blocks(self, seed=None):
        if seed:
            self.state = seed
        blockseed = self.state
        rnd = random.Random()
        rnd.seed(self.state)
        p = rnd.random()
        d = None
        for ix, v in enumerate(self.cdf):
            if v > p:
                d = ix + 1
                break
        if d is None:
            d = ix + 1
        nums = py2_random_sample(rnd, self.K_int, d)
        return blockseed, d, nums


def screen_dna(dna, max_homopolymer=3, gc_dev=0.05):
    """screen_repeat_dna over the 0/1/2/3 representation, but on ACGT here."""
    if len(dna) < 4:
        return False
    run = 1
    for i in range(1, len(dna)):
        if dna[i] == dna[i - 1]:
            run += 1
            if run > max_homopolymer:
                return False
        else:
            run = 1
    gc = sum(1 for b in dna if b in 'GC') / len(dna)
    return 0.5 - gc_dev <= gc <= 0.5 + gc_dev


class ReedSolomon2:
    """reedsolo.RSCodec(2) wrapper; keeps pure detection/correction semantics."""

    def __init__(self):
        from reedsolo import RSCodec
        self.codec = RSCodec(2)

    def check(self, word38):
        """Return corrected 36-byte data if the 38-byte word is valid or
        correctable, else None."""
        try:
            dec = self.codec.decode(bytes(word38))
            return list(dec[0])
        except Exception:
            return None

    def is_codeword(self, word38):
        """True only for an exact RS(38,36) codeword. check() also returns data
        for words in which only the parity needs correction, so comparing its
        result with word[:36] is not a detect-only test."""
        word = bytes(word38)
        return bytes(self.codec.encode(word[:36])) == word

    def encode(self, data36):
        return list(self.codec.encode(bytes(data36)))
