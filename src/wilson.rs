use std::ops::Deref;

use crate::{
    field::Field,
    heightfield::{Action, HeightVariable},
};

pub struct WilsonField<'a, T, const D: usize, const SIZE: usize>
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    field: Field<'a, T, D, SIZE>,
    width: usize,
    height: usize,
    spin: T,
}

impl<'a, T, const D: usize, const SIZE: usize> Deref for WilsonField<'a, T, D, SIZE>
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    type Target = Field<'a, T, D, SIZE>;

    fn deref(&self) -> &Self::Target {
        &self.field
    }
}

impl<'a, T, const D: usize, const SIZE: usize> Action<T> for WilsonField<'a, T, D, SIZE>
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    fn integer_action(&self) -> T {
        todo!()
    }

    fn assumed_local_action(&self, index: usize, value: T) -> T {
        todo!()
    }

    fn bond_action(&self, index: usize, direction: usize) -> T {
        todo!()
    }

    fn assumed_bond_action(&self, index: usize, direction: usize, value: T) -> T {
        todo!()
    }
}
