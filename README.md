# Durandal

Durandal is a rank-based signature scheme based on the Schnorr-Lyubashevsky approach. See [the ePrint](https://eprint.iacr.org/2018/1192.pdf) for more details.

## Dependencies

This software requires :

* g++
* make
* [NTL](http://shoup.net/ntl/)
* [GMP](https://gmplib.org/)

## Compiling and running

Use `make durandal` to compile the signature executable. The executable is placed in the `bin/` folder. Run `bin/durandal` to the program.

By default, the software runs the first set of parameters given in the paper. To switch between parameter set `I` and `II` simply overwrite `src/parameters.h` with `src/parameters-I.h` or `src/parameters-II.h` respectively.

## Example outputs

Parameter set `I` :

```
B4A244C0E85DBEDC018E40BD7C8F9FD69BA99D71761F84264F57AB94BA1A64A79F38961BF1806BBBB8218AB6B3D8CEE1A331557BE57BBB41E46B566CC23BF89860A04B80E0694D53A632B5EC
B4A244C0E85DBEDC018E40BD7C8F9FD69BA99D71761F84264F57AB94BA1A64A79F38961BF1806BBBB8218AB6B3D8CEE1A331557BE57BBB41E46B566CC23BF89860A04B80E0694D53A632B5EC
Valid signature
Keygen : 0.004128
Offline signature : 0.580611
Online Signature : 0.007313
Verification : 0.006905
```

Parameter set `II` :

```
38695A4A0D4F3CDCF4C0E2D3D0B4B480A19FEC08CADC6EABD18B5BA6C176D6908F48BCA07F241A2187761401D97DD6D1EFF1A7B6A91FF78AA165DA8E5C8F50CB207D96C6B7830BD913E81118AEA38136533C8C7725C2F706D28DB2853BBA86836591B8
38695A4A0D4F3CDCF4C0E2D3D0B4B480A19FEC08CADC6EABD18B5BA6C176D6908F48BCA07F241A2187761401D97DD6D1EFF1A7B6A91FF78AA165DA8E5C8F50CB207D96C6B7830BD913E81118AEA38136533C8C7725C2F706D28DB2853BBA86836591B8
Valid signature
Keygen : 0.00532
Offline signature : 1.36114
Online Signature : 0.011028
Verification : 0.010145
```

Timings are given in ms.

## Performances

Here are tables comparing performance of this software with timings given in the paper :

| I | Keygen | Offline signature | Online signature | Verification |
|:---:|:---:|:---:|:---:|:---:|
| Paper | 4 | 350 | 4 | 5 |
| This software | 3.2 | 637 | 7.6 | 7.2 |

| II | Keygen | Offline signature | Online signature | Verification |
|:---:|:---:|:---:|:---:|:---:|
| Paper | 5 | 700 | 5 | 6 |
| This software | 4.6 | 1412 | 11 | 10.1 |

The timings in the paper were based on the implementation of the most costly part for each algorithm, which explains the gap.