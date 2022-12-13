//! This module implements a methods on fields that are of a HeightVariable type.

use std::{
    fmt::Display,
    ops::{Add, Div, Mul, Sub},
};

use rand::prelude::*;

use crate::field::Field;

/// A collection of traits to be fullfilled by the field type, i.e. i8, i32
pub trait HeightVariable<T>:
    Copy
    + Add<Output = T>
    + Sub<Output = T>
    + Mul<Output = T>
    + Div<Output = T>
    + Default
    + From<i8>
    + Into<f64>
    + Display
    + PartialOrd
{
}

impl<T> HeightVariable<T> for T where
    T: Copy
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Div<Output = T>
        + Default
        + From<i8>
        + Into<f64>
        + Display
        + PartialOrd
{
}

/// For [Field]s of [HeightVariable]s we implement possible methods.
pub trait HeightField<T, const D: usize, const SIZE: usize>
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    /// Shifts all values by the value of a random chosen height value.
    fn normalize_random(&mut self);

    /// Subtracts a shift from all values of the field.
    fn shift_values(&mut self, shift: T);
}

impl<'a, T, const D: usize, const SIZE: usize> HeightField<T, D, SIZE> for Field<'a, T, D, SIZE>
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    fn normalize_random(&mut self) {
        let mut rng = ThreadRng::default();
        let &shift = self.values.choose(&mut rng).unwrap();
        self.shift_values(shift);
    }

    /// Subtracts a shift from all values of the field.
    fn shift_values(&mut self, shift: T) {
        for value in self.values.iter_mut() {
            *value = *value - shift;
        }
    }
}

#[test]
fn test_field_shift() {
    use crate::lattice::Lattice;

    let lattice: Lattice<4, 81> = Lattice::new([3, 3, 3, 3]);

    let mut field: Field<i8, 4, 81> = Field::new(&lattice);

    field.shift_values(3);

    assert_eq!(field.values[80], -3);
}

/// This trait implements methods on fields with height variables to calculate
/// some kind of action.
pub trait Action<T: HeightVariable<T>> {
    fn integer_action(&self) -> T;

    /// Calculates the action of the lattice with the coupling constant.
    fn action_observable(&self, temp: f64) -> f64 {
        temp * <T as Into<f64>>::into(self.integer_action())
    }

    /// Calculates the action of a indexed lattice site.
    fn local_action(&self, index: usize) -> T;

    /// Calculates the action of a indexed lattice site, where the value of the
    /// height variable at the indexed site can be given arbitrarily.
    fn assumed_local_action(&self, index: usize, value: T) -> T;

    /// The action of the bond going from the indexed site in the given
    /// direction.
    fn bond_action(&self, index: usize, direction: usize) -> T;

    /// The action of the bond going from the indexed site in the given
    /// direction, where the value of the height variable at the indexed site
    /// can be given arbitrarily.
    fn assumed_bond_action(&self, index: usize, direction: usize, value: T) -> T;
}

impl<'a, T, const D: usize, const SIZE: usize> Action<T> for Field<'a, T, D, SIZE>
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    fn integer_action(&self) -> T {
        let mut sum: T = T::default();
        for index in 0..SIZE {
            for direction in 0..D {
                sum = sum + self.bond_action(index, direction);
            }
        }
        sum
    }

    fn local_action(&self, index: usize) -> T {
        self.assumed_local_action(index, self.values[index])
    }

    fn assumed_local_action(&self, index: usize, value: T) -> T {
        let mut sum: T = T::default();
        for neighbour in self.lattice.get_neighbours_array(index) {
            let diff: T = value - self.values[neighbour];
            sum = sum + (diff * diff);
        }
        sum
    }

    fn bond_action(&self, index: usize, direction: usize) -> T {
        self.assumed_bond_action(index, direction, self.values[index])
    }

    fn assumed_bond_action(&self, index: usize, direction: usize, value: T) -> T {
        let neighbour_index: usize = self.lattice.values[index][direction];
        let diff: T = value - self.values[neighbour_index];
        diff * diff
    }
}
