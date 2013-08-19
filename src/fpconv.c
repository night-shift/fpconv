#include <stdbool.h>
#include <string.h>
#include <assert.h>

#include "fpconv.h"
#include "powers.h"

#define MMSK 0x000FFFFFFFFFFFFFULL /* mantissa */
#define EMSK 0x7FF0000000000000ULL /* exponent */
#define HMSK 0x0010000000000000ULL /* hidden bit */
#define SMSK 0x8000000000000000ULL /* sign bit */
#define EBIAS (1023 + 52) /* exponent bias */

#define ONE_LOG_10 0.30102999566398114 /* 1 / lg 10 */

#define LOMSK 0x00000000FFFFFFFFULL

static uint64_t tens[] = {
    10000000000000000000U, 1000000000000000000U, 100000000000000000U,
    10000000000000000U, 1000000000000000U, 100000000000000U,
    10000000000000U, 1000000000000U, 100000000000U,
    10000000000U, 1000000000U, 100000000U,
    10000000U, 1000000U, 100000U,
    10000U, 1000U, 100U,
    10U, 1U 
};

static inline int absv(int n)
{
    return n < 0 ? -n : n;
}

static inline int ceilv(double fp)
{
    return (int)fp + (1 - (int)((int)(fp + 1) - fp));
}

static inline uint64_t get_dbits(double d)
{
    union {
        double   dbl;
        uint64_t i;
    } dbl_bits = { d };

    return dbl_bits.i;
}

static Fp build_fp(double d)
{
    uint64_t bits = get_dbits(d);

    Fp fp;
    fp.frac = bits & MMSK;
    fp.exp = (bits & EMSK) >> 52;

    if(fp.exp) {
        fp.frac += HMSK;
        fp.exp -= EBIAS;

    } else {
        fp.exp = -EBIAS + 1;
    }
    
    return fp;
}

static void normalize(Fp* fp)
{
    while ((fp->frac & HMSK) == 0) {
        fp->frac <<= 1;
        fp->exp--;
    }

    /* exponent + hidden bit */
    int shift = 64 - 52 - 1;
    fp->frac <<= shift;
    fp->exp -= shift;
}

static Fp get_power10(int k)
{
    /* powers of ten exponent offset */
    int index = 343 + k;

    assert(index >= 0 && index < 687);

    return powers10[index];
}

static int k_comp(int exp, int alpha, int q)
{
    return ceilv((alpha - (exp + q) + (q - 1)) * ONE_LOG_10);
}

static void get_normalized_boundaries(Fp* fp, Fp* lower, Fp* upper)
{
    upper->frac = (fp->frac << 1) + 1; 
    upper->exp  = fp->exp - 1;

    while (!(upper->frac & (HMSK << 1))) {
        upper->frac <<= 1;
        upper->exp--;
    }

    int u_shift = 64 - 52 - 2;

    upper->frac <<= u_shift;
    upper->exp = upper->exp - u_shift;

    bool fraction = fp->frac == HMSK;
    int l_shift = fraction ? 2 : 1;

    lower->frac = (fp->frac << l_shift) - 1;
    lower->exp = fp->exp - l_shift;

    lower->frac <<= lower->exp - upper->exp;
    lower->exp = upper->exp;
}

static Fp multiply(Fp* a, Fp* b)
{
    uint64_t ah_bl = (a->frac >> 32)   * (b->frac & LOMSK);
    uint64_t al_bh = (a->frac & LOMSK) * (b->frac >> 32);
    uint64_t al_bl = (a->frac & LOMSK) * (b->frac & LOMSK);
    uint64_t ah_bh = (a->frac >> 32)   * (b->frac >> 32);

    uint64_t tmp = (ah_bl & LOMSK) + (al_bh & LOMSK) + (al_bl >> 32); 
    /* round up */
    tmp += 1U << 31;

    Fp fp = {
        ah_bh + (ah_bl >> 32) + (al_bh >> 32) + (tmp >> 32),
        a->exp + b->exp + 64
    };
    
    return fp;
}

static void round_digit(char* digits, int ndigits, uint64_t delta, uint64_t rem, uint64_t kappa, uint64_t frac)
{
    while (rem < frac && delta - rem >= kappa &&
           (rem + kappa < frac || frac - rem > rem + kappa - frac)) {

        digits[ndigits - 1]--;
        rem += kappa;
    }
}

static int digit_gen(Fp* fp, Fp* upper, Fp* lower, char* digits, int* K)
{
    uint64_t wfrac = upper->frac - fp->frac;
    uint64_t delta = upper->frac - lower->frac;

    Fp one;
    one.frac = 1ULL << -upper->exp;
    one.exp  = upper->exp;
    
    uint64_t part1 = upper->frac >> -one.exp;
    uint64_t part2 = upper->frac & (one.frac - 1);

    int idx = 0, kappa = 10;
    /* 1000000000 */
    for(uint64_t* divp = tens + 9; kappa > 0; divp++) {

        uint64_t div = *divp;
        unsigned digit = part1 / div;

        if (digit || idx) {
            digits[idx++] = digit + '0';
        }

        part1 -= digit * div;
        kappa--;

        uint64_t tmp = (part1 <<-one.exp) + part2;
        if (tmp <= delta) {
            *K += kappa;
            round_digit(digits, idx, delta, tmp, div << -one.exp, wfrac);

            return idx;
        }
    }

    /* 10 */
    uint64_t* unit = tens + 17;

    while(true) {
        part2 *= 10;
        delta *= 10;
        kappa--;

        unsigned digit = part2 >> -one.exp;
        if (digit || idx) {
            digits[idx++] = digit + '0';
        }

        part2 &= one.frac - 1;
        if (part2 < delta) {
            *K += kappa;
            round_digit(digits, idx, delta, part2, one.frac, wfrac * *unit);

            return idx;
        }

        unit--;
    }
}

static int grisu2(double d, char* digits, int* K)
{
    const int q = 64, alpha = -59;

    Fp w = build_fp(d);

    Fp lower, upper;
    get_normalized_boundaries(&w, &lower, &upper);

    normalize(&w);

    int k = k_comp(upper.exp, alpha, q);
    Fp cp = get_power10(k);

    w     = multiply(&w,     &cp);
    upper = multiply(&upper, &cp);
    lower = multiply(&lower, &cp);

    lower.frac++;
    upper.frac--;

    *K = -k;

    return digit_gen(&w, &upper, &lower, digits, K);
}

static int emit_digits(char* digits, int ndigits, char* buf, int K)
{
    int exp = absv(K + ndigits - 1);

    /* write plain integer */
    if(K >= 0 && (exp < (ndigits + 4))) {
        memcpy(buf, digits, ndigits);
        memset(buf + ndigits, '0', K);

        return ndigits + K;
    }

    /* write decimal w/o scientific notation */
    if(K < 0 && (K > -7 || exp < 4)) {
        int offset = ndigits - absv(K);
        /* fp < 1.0 -> write leading zero */
        if(offset <= 0) {
            offset = -offset;
            buf[0] = '0';
            buf[1] = '.';
            memset(buf + 2, '0', offset);
            memcpy(buf + offset + 2, digits, ndigits);

            return ndigits + 2 + offset;

        /* fp > 1.0 */
        } else {
            memcpy(buf, digits, offset);
            buf[offset] = '.';
            memcpy(buf + offset + 1, digits + offset, ndigits - offset);

            return ndigits + 1;
        }
    }

    /* write decimal w/ scientific notation */

    int idx = 0;
    buf[idx++] = digits[0];

    if(ndigits > 1) {
        buf[idx++] = '.';
        memcpy(buf + idx, digits + 1, ndigits - 1);
        idx += ndigits - 1;
    }

    buf[idx++] = 'e';

    char sign = K + ndigits - 1 < 0 ? '-' : '+';
    buf[idx++] = sign;

    int div = 0;

    if(exp > 99) {
        div = exp / 100;
        buf[idx++] = div + '0';
        exp -= div * 100;
    }
    if(exp > 9) {
        int dec = exp / 10;
        buf[idx++] = dec + '0';
        exp -= dec * 10;

    } else if(div) {
        buf[idx++] = '0';
    }

    buf[idx++] = exp % 10 + '0';

    return idx;
}

static int filter_special(double fp, char* buf)
{
    if(fp == 0.0) {
        buf[0] = '0';
        return 1;
    }

    uint64_t bits = get_dbits(fp);

    bool nan = (bits & EMSK) == EMSK;

    if(!nan) {
        return 0;
    }

    if(bits & MMSK) {
        buf[0] = 'n'; buf[1] = 'a'; buf[2] = 'n';

    } else {
        buf[0] = 'i'; buf[1] = 'n'; buf[2] = 'f';
    }
    
    return 3;
}

int fpconv_dtoa(double d, char dest[24])
{
    char digits[18];

    int str_len = 0;

    if(get_dbits(d) & SMSK) {
        dest[0] = '-';
        str_len++;
    }

    int spec = filter_special(d, dest + str_len);

    if(spec) {
        return str_len + spec;
    }

    int K = 0;
    int ndigits = grisu2(d, digits, &K);

    str_len += emit_digits(digits, ndigits, dest + str_len, K);

    return str_len;
}
