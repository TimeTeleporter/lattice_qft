/// A 3 dimensional spacetime lattice.
#[derive(Debug)]
pub struct Lattice<T, const L: usize, const W: usize, const H: usize>([[[T; H]; W]; L]);

impl<T, const L: usize, const W: usize, const H: usize> Default for Lattice<T, L, W, H>
where
    T: Default,
{
    fn default() -> Lattice<T, L, W, H> {
        let ary = [(); L].map(|_| [(); W].map(|_| [(); H].map(|_| T::default())));
        Lattice(ary)
    }
}

impl<T, const L: usize, const W: usize, const H: usize> Lattice<T, L, W, H>
where
    T: Default + Copy,
{
    pub fn new() -> Lattice<T, L, W, H> {
        Lattice::<T, L, W, H>::default()
    }

    pub fn value(&mut self, x: usize, y: usize, t: usize) -> T {
        self.0[x][y][t]
    }
}
