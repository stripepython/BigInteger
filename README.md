# Operation Documentation

BigInteger version 3.0 supports various operations. In the following complexity analysis, $w=8, w'=4$.

## Initialization

- `BigInteger()`: Creates a new `BigInteger` with default value $0$.
- `BigInteger(const BigInteger& x)`: Creates a new `BigInteger` with value $x$, time complexity $O(\dfrac{n}{w})$ where $n$ is the length of $x$.
- `BigInteger(int64_t x)`: Creates a new `BigInteger` with value $x$, time complexity $O(\log x)$.
- `BigInteger(const std::string& s)`: Creates a new `BigInteger` from a string, time complexity $O(n)$ where $n$ is the string length. Valid strings must consist of one or more `-` signs followed by numeric characters.
- `BigInteger(const std::vector<bool>& v)`: Creates a new `BigInteger` from binary representation, time complexity $O(n^2)$ where $n$ is the length of binary representation.
- `BigInteger.from_int128(__int128 x)`: A `static` function that creates a new `BigInteger` from `__int128` type with value $x$, time complexity $O(\log x)$. In environments without `__int128` support, this operation is unavailable.

## I/O

- `std::cin >> x`: Inputs a `BigInteger` value, time complexity $O(n)$ where $n$ is the string length. Valid input requirements are the same as string initialization.
- `std::cout << x`: Outputs a `BigInteger` value, time complexity $O(n)$ where $n$ is the length of the integer.

## Type Conversion

- `a.to_string()`: Returns `std::string` type, the string representation of `a`. Time complexity $O(n)$ where $n$ is the length of the integer.
- `a.to_int64()`: Returns `int64_t` type, the 64-bit integer conversion of `a`, time complexity $O(\dfrac{n}{w})$ where $n$ is the length of the integer. Behavior is undefined if overflow occurs.
- `a.to_binary()`: Returns `std::vector<bool>` type, the binary representation of $a$, time complexity $O(n^2)$ where $n$ is the length of the integer.
- `a.to_int128()`: Returns `__int128` type, the 128-bit integer conversion of `a`, time complexity $O(\dfrac{n}{w})$ where $n$ is the length of the integer. Behavior is undefined if overflow occurs. In environments without `__int128` support, this operation is unavailable.

## Basic Operations

- `a.zero()`: Checks if $a$ equals $0$, time complexity $O(1)$.
- `!a`: Checks if $a$ is not $0$, time complexity $O(1)$.
- `a.positive()`: Checks if $a$ is positive, time complexity $O(1)$. **$0$ is not positive.**
- `a.negative()`: Checks if $a$ is negative, time complexity $O(1)$.

- `a.compare(const BigInteger& b)`: Returns comparison result between $a$ and $b$ as `int` type. Returns $-1$ if $a<b$, $0$ if $a=b$, and $1$ if $a>b$. Time complexity $O(\dfrac{n}{w})$ where $n$ is the maximum length of the two integers.
- `a <=> b, a <= b, a < b, a == b, a != b, a > b, a >= b`: Returns corresponding comparison results, time complexity $O(\dfrac{n}{w})$ where $n$ is the maximum length of the two integers.
- `-a`: Returns $-a$, time complexity $O(\dfrac{n}{w})$ where $n$ is the length of the integer.
- `~a`: Returns $-a-1$, time complexity $O(\dfrac{n}{w})$ where $n$ is the length of the integer.
- `a.abs()`: Returns $|a|$, time complexity $O(\dfrac{n}{w})$ where $n$ is the length of the integer.
- `a + b`: Returns $a+b$, time complexity $O(\dfrac{n}{w})$ where $n$ is the maximum length of the two integers. Supports in-place addition `a += b`. Faster when `b` is `int32_t` type.
- `a - b`: Returns $a-b$, time complexity $O(\dfrac{n}{w})$ where $n$ is the maximum length of the two integers. Supports in-place subtraction `a -= b`. Faster when `b` is `int32_t` type.
- `a * b`: Returns $a \times b$, time complexity $O(\dfrac{n \log n}{w'})$ where $n$ is the maximum length of the two integers. When $n <$ `8 * FFT_LIMIT`, uses $O(\dfrac{n^2}{w^2})$ long multiplication. `FFT_LIMIT` defaults to $8$. When `b` is `int32_t` type, time complexity is $O(\dfrac{n}{w})$ and supports in-place multiplication. Throws `FFTLimitExceededError` when $n > 2^{20}$.
- `a.square()`: Returns $a^2$, time complexity $O(\dfrac{n}{w})$ where $n$ is the length of the integer, faster than `a * a`.

- `a.half()`: Returns $\lfloor \dfrac{a}{2} \rfloor$, time complexity $O(\dfrac{n}{w})$ where $n$ is the length of the integer, faster than `a / 2`.

- `a / b`: Returns $\lfloor \dfrac{a}{b} \rfloor$, time complexity $O(\dfrac{n \log n}{w'})$ where $n$ is the maximum length of the two integers. When $n <$ `8 * NEWTON_DIV_LIMIT`, uses $O(\dfrac{n^2}{w})$ long division. `NEWTON_DIV_LIMIT` defaults to $32$. When `b` is `int64_t` type, time complexity is $O(\dfrac{n}{w})$ and supports in-place division. Throws `ZeroDivisionError` when $b=0$.

- `a % b`: Returns $a \bmod b$, same time complexity as `a / b`. Throws `ZeroDivisionError` when $b=0$.

- `a.divmod(b)`: Returns an `std::pair` of $(\lfloor \dfrac{a}{b} \rfloor, a \bmod b)$, same time complexity as `a / b`, but without optimization when `b` is `int64_t` type. Throws `ZeroDivisionError` when $b=0$.

- `a.mod2()`: Returns $a \bmod 2$, time complexity $O(1)$.

- `a.pow(b)`: Returns $a^b$, time complexity $O(\dfrac{nb \log nb}{w})$ where $n$ is the length of the integer. `b` should be `int64_t` type.

- `a.pow(b, p)`: Returns $a^b \bmod p$, time complexity $O(\dfrac{nb \log nb}{w})$ where $n$ is the length of the integer. `b` should be `int64_t` type, and `p` should be `BigInteger` type.

- `a.sqrt()`: Returns $\lfloor \sqrt{a} \rfloor$, time complexity $O(\dfrac{n \log n}{w})$ where $n$ is the length of the integer. Throws `NegativeRadicandError` when $a<0$.

- `a.root(x)`: Returns $\lfloor \sqrt[x]{a} \rfloor$, time complexity $O(\dfrac{n \log n}{w})$ where $n$ is the length of the integer. Throws `NegativeRadicandError` when $x \le 0$. Throws `NegativeRadicandError` when $2 \mid x$ and $a < 0$.

- `a.gcd(b)`: Returns $\gcd(a,b)$, time complexity $O(\dfrac{n^2}{w})$ where $n$ is the maximum length of the two integers.

- `a.lcm(b)`: Returns $\operatorname{lcm}(a,b)$, time complexity $O(\dfrac{n^2}{w})$ where $n$ is the maximum length of the two integers.

- `a << x`: Returns $a \times 2^x$, time complexity $O(\dfrac{(n+x) \log(n+x)}{w'})$ where $n$ is the length of the integer.

- `a >> x`: Returns $\lfloor \dfrac{a}{2^x} \rfloor$, time complexity $O(\dfrac{(n-x) \log(n-x)}{w'})$ where $n$ is the length of the integer.

- `a & b`: Returns bitwise AND of $a$ and $b$, time complexity $O(n^2)$ where $n$ is the maximum length of the two integers.

- `a | b`: Returns bitwise OR of $a$ and $b$, time complexity $O(n^2)$ where $n$ is the maximum length of the two integers.

- `a ^ b`: Returns bitwise XOR of $a$ and $b$, time complexity $O(n^2)$ where $n$ is the maximum length of the two integers.

## Other Functions

- `factorial(n)`: Returns `BigInteger` type, the value of $n!$, time complexity $O(\dfrac{n^2}{w})$.
- `rand_bigint(n)`: Returns a random `BigInteger` of length $n$, time complexity $O(n)$.

## Internal Functions

These functions are not recommended for use.

- `a._digit_len()`: Returns $\lfloor \dfrac{n}{w} \rfloor$ where $n$ is the length of the integer, time complexity $O(1)$.
- `a._move_l(x)`: Returns $|n \times 10^{wx}|$, time complexity $O(\dfrac{n}{w}+x)$ where $n$ is the length of the integer.
- `a._move_r(x)`: Returns $|\lfloor \dfrac{n}{10^{wx}} \rfloor|$, time complexity $O(\dfrac{n}{w}-x)$ where $n$ is the length of the integer.
- `__FFT::dft(a, n)`: Performs DFT on array `a` of length $n$ with type `__FFT::complex[]`, time complexity $O(n\log n)$. Requires $n$ to be a power of 2 not exceeding $2^{21}$.
- `__FFT::idft(a, n)`: Performs IDFT on array `a` of length $n$ with type `__FFT::complex[]`, time complexity $O(n\log n)$. Requires $n$ to be a power of 2 not exceeding $2^{21}$.
- `__helper(a, b, f)`: Performs bitwise operation `f` on $a$ and $b$, time complexity $O(n^2)$ where $n$ is the maximum length of the two integers. `f` should be a function like `bool f(bool, bool)`.

# Pros and Cons

Pros:

- Highly encapsulated, supports almost all operations of `int` type without manual implementation.
- Fast speed, most operations are optimized to excellent complexity.
- Relatively short code length.

Cons:

- Slow bitwise operations.
- Doesn't support very large integer multiplication/division. Specifically, throws `FFTLimitExceededError` when bit length exceeds $2^{20}=1048576$.
- Occasional bugs.

# Acknowledgments

Thanks to contributors (listed in dictionary order):

- @[bcdmwSjy](https://www.luogu.com.cn/user/514727)
- @[konyakest](https://www.luogu.com.cn/user/482660)
- @[Noiers](https://www.luogu.com.cn/user/1402616)
- @[Xudongning](https://www.luogu.com.cn/user/1636821)
