/// Generic N-form in D dimensions
/// See: https://blog.rust-lang.org/inside-rust/2021/09/06/Splitting-const-generics.html
#[derive(Debug, Clone, Copy)]
pub struct Form<T, const D: usize, const N: usize>
where
    [(); binomial_coeff_recursive(D, N)]:,
{
    pub components: [T; binomial_coeff_recursive(D, N)],
}

impl<T: Default, const D: usize, const N: usize> Default for Form<T, D, N>
where
    [(); binomial_coeff_recursive(D, N)]:,
{
    fn default() -> Self {
        let components = [(); binomial_coeff_recursive(D, N)].map(|_| T::default());

        Self { components }
    }
}

/// Implements the binomial coefficient n over k recursively.
///
/// ```
///  / n \   / n-1 \   / n-1 \      / n \   / n \
///  |   | = |     | + |     |  and |   | = |   | = 1 for all n;
///  \ k /   \ k-1 /   \  k  /      \ 0 /   \ n /
/// ```
///
pub const fn binomial_coeff_recursive(n: usize, k: usize) -> usize {
    match (n, k) {
        (_, 0) => 1,
        (n, k) if n == k => 1,
        _ => binomial_coeff_recursive(n - 1, k - 1) + binomial_coeff_recursive(n - 1, k),
    }
}

#[test]
fn test_binomial_coeff_recursive() {
    assert_eq!(binomial_coeff_recursive(4, 2), 6);
    assert_eq!(binomial_coeff_recursive(6, 3), 20);
    assert_eq!(binomial_coeff_recursive(9, 6), 84);
}
