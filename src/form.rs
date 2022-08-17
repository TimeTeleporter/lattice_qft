/// Generic N-form in D dimensions
/// See: https://blog.rust-lang.org/inside-rust/2021/09/06/Splitting-const-generics.html
#[derive(Debug, Clone, Copy)]
pub struct Form<T, const D: usize, const N: usize>([T; binomial_coeff_recursive(D, N)])
where
    [(); binomial_coeff_recursive(D, N)]:;

impl<T: Default, const D: usize, const N: usize> Default for Form<T, D, N>
where
    [(); binomial_coeff_recursive(D, N)]:,
{
    fn default() -> Self {
        Self([(); binomial_coeff_recursive(D, N)].map(|_| T::default()))
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
    assert_eq!(binomial_coeff_recursive(3, 0), 1);
    assert_eq!(binomial_coeff_recursive(4, 2), 6);
    assert_eq!(binomial_coeff_recursive(6, 3), 20);
    assert_eq!(binomial_coeff_recursive(9, 6), 84);
}

#[test]
fn test_lattice3d_forms() {
    let zero_form = Form::<i32, 3, 0>::default();
    assert_eq!(zero_form.0, [0]);

    let one_form = Form::<i32, 3, 1>::default();
    assert_eq!(one_form.0, [0, 0, 0]);

    let two_form = Form::<i32, 3, 2>::default();
    assert_eq!(two_form.0, [0, 0, 0]);

    let three_form = Form::<i32, 3, 3>::default();
    assert_eq!(three_form.0, [0]);
}
