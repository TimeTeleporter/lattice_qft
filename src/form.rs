//! If we consider a discrete euclidean spacetime, we can identify different objects as living on different parts of the lattice.
//! As an exmaple, consider a 2+1 spacetime lattice. We identify the objects
//!
//! 1. Lattice points - zero-forms.
//! 2. Links between lattice points - one-forms.
//! 3. Plaquettes - two-forms.
//! 4. Cubes - three-forms.
//!
//! For each of these objects we can define fields.
//! Further realize, that there are equally as many cubes as lattice points in this case.
//! We also have three times the amount of links than we have lattice points.
//! We have the same amount of links as plaquettes.
//!
//! Further we can compose the forms to get new ones.

use num::Zero; // This allowes us to use generics that are able to become zero.

/// Generic N-form in D dimensions
#[derive(Debug, Clone, Copy)]
pub struct Form<T, const D: usize, const N: usize>
where
    T: Copy,
    [(); binomial_coeff_recursive(D, N)]:,
{
    components: [T; binomial_coeff_recursive(D, N)],
}

impl<T, const D: usize, const N: usize> Default for Form<T, D, N>
where
    T: Default + Copy,
    [(); binomial_coeff_recursive(D, N)]:,
{
    fn default() -> Self {
        let components = [(); binomial_coeff_recursive(D, N)].map(|_| T::default());

        Self { components }
    }
}

impl<T, const D: usize, const N: usize> Form<T, D, N>
where
    T: Default + Copy,
    [(); binomial_coeff_recursive(D, N)]:,
{
    fn new() -> Self {
        Self::default()
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
