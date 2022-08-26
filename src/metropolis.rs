use std::ops::{Add, Mul, Sub};

use rand::{random, Rng};

use crate::{field3d::Field3d, lattice3d::Directions};

pub trait Metropolis<const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    type Input: Add<Output = Self::Input>
        + Default
        + Sub<Output = Self::Input>
        + From<i8>
        + Copy
        + Mul<Output = Self::Input>;

    fn metropolis_sweep(&mut self)
    where
        f64: From<<Self as Metropolis<MAX_X, MAX_Y, MAX_T>>::Input>;

    fn calculate_action(&self) -> Self::Input;
}

impl<'a, T, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    Metropolis<MAX_X, MAX_Y, MAX_T> for Field3d<'a, T, MAX_X, MAX_Y, MAX_T>
where
    T: Copy
        + Mul<Output = T>
        + From<i8>
        + Sub<Output = T>
        + Sub<Output = T>
        + Default
        + Add<Output = T>,
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    type Input = T;

    fn metropolis_sweep(&mut self)
    where
        f64: From<<Self as Metropolis<MAX_X, MAX_Y, MAX_T>>::Input>,
    {
        let mut rng = rand::thread_rng();

        for (index, _) in self.lattice.0.iter().enumerate() {
            let value = self.values[index].clone();
            let coin: bool = random();
            let new_value: Self::Input = match coin {
                true => value + 1_i8.into(),
                false => value - 1_i8.into(),
            };

            let mut sum_differenced: Self::Input = Self::Input::default();
            let mut new_sum_differenced: Self::Input = Self::Input::default();

            for direction in Directions::as_array() {
                let next_neighbour = (self.get_next_neighbour_value(index, direction)).clone();
                let prev_neighbour = (self.get_prev_neighbour_value(index, direction)).clone();
                calculate_single_direction_action(
                    &mut sum_differenced,
                    next_neighbour,
                    value,
                    prev_neighbour,
                );
                calculate_single_direction_action(
                    &mut new_sum_differenced,
                    next_neighbour,
                    value,
                    prev_neighbour,
                );
            }

            let difference = sum_differenced - new_sum_differenced;

            let exponent: f64 = f64::exp(difference.into());

            let rng: f64 = rng.gen_range(0.0..1.0);

            if rng <= exponent {
                self.set_value(new_value, index);
            }
        }
    }

    fn calculate_action(&self) -> T {
        let mut action: T = T::default();
        for (index, &value) in self.values.iter().enumerate() {
            for direction in Directions::as_array() {
                let next_neighbour = (self.get_next_neighbour_value(index, direction)).clone();
                let prev_neighbour = (self.get_prev_neighbour_value(index, direction)).clone();
                calculate_single_direction_action(
                    &mut action,
                    next_neighbour,
                    value,
                    prev_neighbour,
                );
            }
        }
        action
    }
}

fn calculate_single_direction_action<T>(
    sum_differenced: &mut T,
    next_neighbour: T,
    value: T,
    prev_neighbour: T,
) where
    T: Copy + Mul<Output = T> + From<i8> + Sub<Output = T> + Sub<Output = T> + Add<Output = T>,
{
    *sum_differenced = *sum_differenced
        + next_neighbour * next_neighbour
        + value * 2_i8.into() * (value - next_neighbour - prev_neighbour)
        + prev_neighbour * prev_neighbour;
}

#[test]
fn test_calculate_single_direction_action() {
    let mut target = 0;
    calculate_single_direction_action::<i32>(&mut target, 100, 3, -200);
    assert_eq!(target, 50618);
}
