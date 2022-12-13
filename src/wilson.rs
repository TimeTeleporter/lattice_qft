use std::ops::{Deref, DerefMut};

use crate::{
    algorithm::bonds::BondsField,
    field::Field,
    heightfield::{Action, HeightVariable},
    lattice::Lattice,
};

pub struct WilsonField<'a, T, const SIZE: usize>
where
    T: HeightVariable<T>,
{
    field: Field<'a, T, 3, SIZE>,
    bonds: BondsField<'a, 3, SIZE>,
    spin: bool,
}

impl<'a, T: HeightVariable<T>, const SIZE: usize> Deref for WilsonField<'a, T, SIZE> {
    type Target = Field<'a, T, 3, SIZE>;

    fn deref(&self) -> &Self::Target {
        &self.field
    }
}

impl<'a, T: HeightVariable<T>, const SIZE: usize> DerefMut for WilsonField<'a, T, SIZE> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.field
    }
}

impl<'a, T, const SIZE: usize> WilsonField<'a, T, SIZE>
where
    T: HeightVariable<T>,
{
    pub fn new(
        lattice: &'a Lattice<3, SIZE>,
        width: usize,
        height: usize,
    ) -> WilsonField<'a, T, SIZE> {
        let field: Field<'a, T, 3, SIZE> = Field::new(lattice);
        Self::from_field(field, width, height)
    }

    pub fn from_field(
        field: Field<'a, T, 3, SIZE>,
        width: usize,
        height: usize,
    ) -> WilsonField<'a, T, SIZE> {
        let mut bonds: BondsField<'a, 3, SIZE> = BondsField::new(field.lattice);
        for index in 0..SIZE {
            let [x, y, t]: [usize; 3] = field.lattice.calc_coords_from_index(index).into_array();
            if y == 0 && x < width && t < height {
                bonds.activate(index, 1);
            }
        }
        bonds.check_coherence();
        WilsonField {
            field,
            bonds,
            spin: true,
        }
    }
}

impl<'a, T, const SIZE: usize> Action<T> for WilsonField<'a, T, SIZE>
where
    T: HeightVariable<T>,
{
    fn integer_action(&self) -> T {
        let mut sum: T = T::default();
        for index in 0..SIZE {
            for direction in 0..3_usize {
                sum = sum + self.bond_action(index, direction);
            }
        }
        sum
    }

    fn local_action(&self, index: usize) -> T {
        self.assumed_local_action(index, self.field.values[index])
    }

    fn assumed_local_action(&self, index: usize, value: T) -> T {
        let mut sum: T = T::default();
        for direction in 0..(3 * 2_usize) {
            sum = sum + self.assumed_bond_action(index, direction, value);
        }
        sum
    }

    fn bond_action(&self, index: usize, direction: usize) -> T {
        self.assumed_bond_action(index, direction, self.field.values[index])
    }

    fn assumed_bond_action(&self, index: usize, direction: usize, value: T) -> T {
        if self.bonds.is_active(index, direction) {
            println!("Modified");
            let spin_modifier: T = ((match direction {
                1 => 1_i8,
                4 => -1_i8,
                _ => 0_i8,
            }) * (match self.spin {
                true => 1_i8,
                false => -1_i8,
            }))
            .into();
            let neighbour_index: usize = self.field.lattice.values[index][direction];
            let diff: T = value - self.field.values[neighbour_index] + spin_modifier;
            diff * diff
        } else {
            self.field.assumed_bond_action(index, direction, value)
        }
    }
}

#[test]
fn test_wilson_field() {
    const SIZE_ARY: [usize; 3] = [2, 2, 2];
    const SIZE: usize = SIZE_ARY[0] * SIZE_ARY[1] * SIZE_ARY[2];
    let lattice: Lattice<3, SIZE> = Lattice::new(SIZE_ARY);
    let field: Field<i32, 3, SIZE> = Field::new(&lattice);
    let field: WilsonField<i32, SIZE> = WilsonField::from_field(field, 1, 1);

    // Checking all field values to be zero
    for index in 0..SIZE {
        assert_eq!(field.field.values[index], 0);
    }

    // Checking only the given bonds to be active
    for index in 0..SIZE {
        for direction in 0..(3 * 2_usize) {
            match (index, direction) {
                (0, 1) | (2, 4) => assert!(field.bonds.values[index][direction]),
                (_, _) => assert!(!field.bonds.values[index][direction]),
            }
        }
    }

    // Checking the modifiers to be correct
    for index in 0..SIZE {
        for direction in 0..(3 * 2_usize) {
            let bond = field.bond_action(index, direction);
            match (index, direction) {
                (0, 1) => assert_eq!(bond, 1),
                (2, 4) => assert_eq!(bond, -1),
                (_, _) => assert_eq!(bond, 0),
            }
        }
    }
}
