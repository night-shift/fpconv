
Minimalistic C implementation of [Fabian Loitsch's Grisu-algorithim](http://florian.loitsch.com/publications/dtoa-pldi2010.pdf)
Grisu converts floating point numbers to an optimal decimal string representation without precision loss.

#Api

```c
int fpconv_dtoa(double fp, char[24] dest);
```
    *Writes string representation of fp to dest and returns the number of written characters
    *The emitted string will never exceed 24 characters
    *Does not null terminate the string
    
Exemplary usage
```c
void print(double d)
{
    char buf[24 + 1]; /* reserve space for null terminator */
    int str_len = fpconv_dtoa(d, buf);
    
    buf[str_len] = '\0';
    printf("%s", buf);
}
```
    



