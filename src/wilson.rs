use std::fmt::Display;

use serde::{Deserialize, Serialize};

use crate::{
    field::{Field, HeightVariable},
    lattice::Lattice,
};

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Wilsonloop {
    width: usize,
    height: usize,
    right: bool,
}

impl Wilsonloop {
    pub fn new(width: usize, height: usize, right: bool) -> Wilsonloop {
        Wilsonloop {
            width,
            height,
            right,
        }
    }
}

impl Display for Wilsonloop {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.right {
            true => write!(f, "{}x{} right-handed Wilson loop", self.width, self.height),
            false => write!(f, "{}x{} left-handed Wilson loop", self.width, self.height),
        }
    }
}

impl Wilsonloop {
    /// Checks if the given bond goes through the wilson loop.
    ///
    /// This only works for 3-dimensional lattices accuratly, as wilson loops
    /// are undifined in other geometries.
    pub fn is_bond_through_loop<'a, const D: usize, const SIZE: usize>(
        &self,
        lattice: &'a Lattice<D, SIZE>,
        index: usize,
        direction: usize,
    ) -> bool
    where
        [(); D * 2_usize]:,
    {
        let Wilsonloop { width, height, .. } = self;
        let coords: [usize; D] = lattice.calc_coords_from_index(index).into_array();
        let yep: bool = ((direction == 1 && coords[1] == 0)
            || (direction == D + 1 && coords[1] == 1))
            && coords[0] < *width
            && coords[2] < *height;
        yep
    }

    pub fn direction_bond_formula<T: HeightVariable<T>, const D: usize, const SIZE: usize>(
        &self,
        field: &Field<T, D, SIZE>,
        value: T,
        index: usize,
        direction: usize,
    ) -> T
    where
        [(); D * 2_usize]:,
    {
        let modifier: T = if self.is_bond_through_loop(field.lattice, index, direction) {
            Into::<T>::into(1_i8)
        } else {
            Into::<T>::into(0_i8)
        };
        let diff: T =
            value - field.values[field.lattice.get_neighbours_array(index)[direction]] + modifier;
        diff * diff
    }

    pub fn direction_bond_action<T: HeightVariable<T>, const D: usize, const SIZE: usize>(
        &self,
        field: &Field<T, D, SIZE>,
        index: usize,
        direction: usize,
    ) -> T
    where
        [(); D * 2_usize]:,
    {
        self.direction_bond_formula(field, field.values[index], index, direction)
    }

    pub fn calculate_assumed_action<T: HeightVariable<T>, const D: usize, const SIZE: usize>(
        &self,
        field: &Field<T, D, SIZE>,
        index: usize,
        value: T,
    ) -> T
    where
        [(); D * 2_usize]:,
    {
        let mut sum: T = T::default();
        for direction in 0..(D * 2_usize) {
            sum = sum + self.direction_bond_formula(field, value, index, direction);
        }
        sum
    }

    pub fn integer_action<T: HeightVariable<T>, const D: usize, const SIZE: usize>(
        &self,
        field: &Field<T, D, SIZE>,
    ) -> T
    where
        [(); D * 2_usize]:,
    {
        let mut sum: T = T::default();
        for index in 0..SIZE {
            for direction in 0..D {
                sum = sum + self.direction_bond_action(field, index, direction);
            }
        }
        sum
    }

    pub fn action_observable<T: HeightVariable<T>, const D: usize, const SIZE: usize>(
        &self,
        field: &Field<T, D, SIZE>,
        temp: f64,
    ) -> f64
    where
        [(); D * 2_usize]:,
    {
        self.integer_action(field).into() * temp
    }
}
