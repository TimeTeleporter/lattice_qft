use std::ops::{Add, Mul, Sub};

use rand::{rngs::ThreadRng, Rng};

use crate::{field3d::Field3d, lattice3d::Directions};

#[cfg(test)]
use crate::{field3d::ConvertField, lattice3d::Lattice3d};

type MyInt = i64;

pub trait Metropolis<const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    type FieldType: Add<Output = Self::FieldType>
        + Default
        + Sub<Output = Self::FieldType>
        + From<i8>
        + Copy
        + Mul<Output = Self::FieldType>;

    const TEMP: f64;

    //fn metropolis_single(&mut self, index: usize);

    fn metropolis_sweep(&mut self);

    fn calculate_action(&self) -> Self::FieldType;
}

impl<'a, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Metropolis<MAX_X, MAX_Y, MAX_T>
    for Field3d<'a, i64, MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    type FieldType = i64;

    const TEMP: f64 = 0.01;

    fn metropolis_sweep(&mut self) {
        let mut rng = ThreadRng::default();

        for (index, _) in self.lattice.0.iter().enumerate() {
            let value = self.values[index].clone();
            let coin: bool = rng.gen();
            let new_value: Self::FieldType = match coin {
                true => value + 1_i8 as Self::FieldType,
                false => value - 1_i8 as Self::FieldType,
            };

            let mut sum_differenced: Self::FieldType = Self::FieldType::default();
            let mut new_sum_differenced: Self::FieldType = Self::FieldType::default();

            for direction in Directions::as_array() {
                let next_neighbour = (self.get_next_neighbour_value(index, direction)).clone();
                let prev_neighbour = (self.get_prev_neighbour_value(index, direction)).clone();
                sum_differenced = sum_differenced
                    + calculate_single_direction_action(next_neighbour, value, prev_neighbour);
                new_sum_differenced = new_sum_differenced
                    + calculate_single_direction_action(next_neighbour, new_value, prev_neighbour);
            }

            let difference = new_sum_differenced - sum_differenced;

            let action: f64 = difference as f64;
            let action = action * Self::TEMP;

            let exponent: f64 = f64::exp(action);

            let rng: f64 = rng.gen_range(0.0..1.0);

            //println!("{}, {} drawn. Action {}", exponent, rng, action);

            if rng <= exponent {
                self.set_value(new_value, index);
            } else {
                //println!("Was rejected.");
            }
        }
    }

    fn calculate_action(&self) -> Self::FieldType {
        let mut action: Self::FieldType = Self::FieldType::default();
        for (index, &value) in self.values.iter().enumerate() {
            for direction in Directions::as_array() {
                let next_neighbour = (self.get_next_neighbour_value(index, direction)).clone();
                action = action + value * value + next_neighbour * next_neighbour
                    - value * next_neighbour * 2_i8 as Self::FieldType;
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
/*
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
*/
