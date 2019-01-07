Minimalistic C / D implementation of Fabian Loitsch's Grisu-algorithm [[pdf]](http://florian.loitsch.com/publications/dtoa-pldi2010.pdf).
Grisu converts floating point numbers to an optimal decimal string representation without loss of precision.

### C Api
```c
int fpconv_dtoa(double fp, char dest[24]);
```
* Writes string representation of ```fp``` to ```dest``` and returns the number of written characters
* The emitted string will never exceed 24 characters
* Does not null terminate the string
* Assumes ```fp``` is an IEEE 64-bit floating point number

### Example usage
```c
void print(double d)
{
    char buf[24 + 1]; /* reserve space for null terminator */
    int str_len = fpconv_dtoa(d, buf);

    buf[str_len] = '\0';
    printf("%s", buf);
}
```

### Why not just use `snprintf`?
Convert doubles faster to shortest strings without precision loss.

Average processing time on a mix of "long" and "short" doubles in nanoseconds:
```
             short long
snprintf %g : 515  700
     %0.18g : 989  1171
         %f : 737  3047
fpconv_dtoa : 165  193

snprintf overhead : 71
```
Measured with `gcc-4.8 -O3` and `glibc 2.17`.

