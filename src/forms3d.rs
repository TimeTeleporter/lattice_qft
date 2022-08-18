use std::ops::Sub;

use crate::{
    indexed_lattice3d::IndexedLattice3d,
    lattice3d::{Direction, Lattice},
};

/// For fields living on the sites of the lattice.
pub type ZeroForm<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> =
    IndexedLattice3d<T, MAX_X, MAX_Y, MAX_T>;

/// For fields on the links of the lattice, saved per site, one per direction.
pub type OneForm<T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> =
    IndexedLattice3d<[T; 3], MAX_X, MAX_Y, MAX_T>;

impl<
        T: Default + Sub<Output = T> + Copy,
        const MAX_X: usize,
        const MAX_Y: usize,
        const MAX_T: usize,
    > Codifferential<OneForm<T, MAX_X, MAX_Y, MAX_T>> for ZeroForm<T, MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    /// Builds a 1-form from a 0-form.
    fn increase_degree(&self) -> OneForm<T, MAX_X, MAX_Y, MAX_T> {
        let mut one_form = OneForm::<T, MAX_X, MAX_Y, MAX_T>::default();

        for (index, value) in self.values.0.iter().enumerate() {
            one_form.values.0[index][0] =
                self.values.0[self.next_neighbour_index(index, Direction::X)] - (*value);
            one_form.values.0[index][1] =
                self.values.0[self.next_neighbour_index(index, Direction::Y)] - (*value);
            one_form.values.0[index][2] =
                self.values.0[self.next_neighbour_index(index, Direction::T)] - (*value);
        }

        one_form
    }
}

pub trait Codifferential<T> {
    fn increase_degree(&self) -> T;
}
