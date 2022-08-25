use std::ops::{Add, Mul};

use crate::forms3d::{Form, OneForm};

/// Datatype to calculate the action of the whole lattice.
pub type Action<'a, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> =
    Form<'a, f64, MAX_X, MAX_Y, MAX_T>;

impl<'a, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Action<'a, MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    pub fn get_lattice_action(&self) -> f64 {
        self.values.iter().product()
    }
}

impl<'a, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Action<'a, MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    pub fn from_one_form<T>(one_form: &OneForm<'a, T, MAX_X, MAX_Y, MAX_T>, e: f64) -> Self
    where
        T: Into<f64> + Copy + Add<Output = T> + Mul<Output = T>,
    {
        let mut action = Action::<'a, MAX_X, MAX_Y, MAX_T>::from_lattice(one_form.lattice);

        for (index, value) in action.values.iter_mut().enumerate() {
            *value = one_form.values[index].square().into();

            *value = 0.0 - *value * e * e / 2.0;

            *value = f64::exp(*value)
        }

        action
    }
}
