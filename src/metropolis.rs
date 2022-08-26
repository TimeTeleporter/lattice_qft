use std::ops::{Add, Mul, Sub};

use rand::{random, Rng};

use crate::{field3d::Field3d, lattice3d::Directions};

#[cfg(test)]
use crate::{field3d::ConvertField, lattice3d::Lattice3d};

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

    const TEMP: f64;

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

    const TEMP: f64 = 1.0;

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
                sum_differenced = sum_differenced
                    + calculate_single_direction_action(next_neighbour, value, prev_neighbour);
                new_sum_differenced = new_sum_differenced
                    + calculate_single_direction_action(next_neighbour, value, prev_neighbour);
            }

            let difference = new_sum_differenced - sum_differenced;

            let action: f64 = difference.into();
            let action = action * Self::TEMP;

            let exponent: f64 = f64::exp(action);

            let rng: f64 = rng.gen_range(0.0..1.0);

            println!(
                "Probability is {}, {} was drawn. Action: {}",
                exponent, rng, action
            );

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
                action = action + value * value + next_neighbour * next_neighbour
                    - value * next_neighbour * 2_i8.into();
            }
        }
        action
    }
}

fn calculate_single_direction_action<T>(next_neighbour: T, value: T, prev_neighbour: T) -> T
where
    T: Copy + Mul<Output = T> + From<i8> + Sub<Output = T> + Sub<Output = T> + Add<Output = T>,
{
    next_neighbour * next_neighbour
        + value * 2_i8.into() * (value - next_neighbour - prev_neighbour)
        + prev_neighbour * prev_neighbour
}

#[test]
fn test_calculate_single_direction_action() {
    let target = calculate_single_direction_action::<i32>(100, 3, -200);
    assert_eq!(target, 50618);
}

#[test]
fn test_action_adding_one_everywhere() {
    const TEST_X: usize = 10;
    const TEST_Y: usize = 10;
    const TEST_T: usize = 10;

    let lattice: Lattice3d<TEST_X, TEST_Y, TEST_T> = Lattice3d::default();

    let field: Field3d<i8, TEST_X, TEST_Y, TEST_T> = Field3d::random(&lattice);
    let mut field: Field3d<i32, TEST_X, TEST_Y, TEST_T> = Field3d::from_field(field);

    field.print_values_formated();

    let action = field.calculate_action().clone();

    field.values = field.values.into_iter().map(|x| x + 1).collect();

    field.print_values_formated();

    assert_eq!(action, field.calculate_action());
}
