use crate::{
    bonds::LinksField,
    field::Field,
    heightfield::{Action, HeightField, HeightVariable},
    lattice::Lattice,
};

pub struct WilsonField<'a, T, const SIZE: usize>
where
    T: HeightVariable<T>,
{
    pub field: Field<'a, T, 3, SIZE>,
    links: LinksField<'a, 3, SIZE>,
    spin: bool,
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
        let mut links: LinksField<'a, 3, SIZE> = LinksField::new(field.lattice);
        for index in 0..SIZE {
            let [x, y, t]: [usize; 3] = field.lattice.calc_coords_from_index(index).into_array();
            if y == 0 && x < width && t < height {
                links.activate(index, 1);
            }
        }
        links.check_coherence();
        WilsonField {
            field,
            links,
            spin: true,
        }
    }
}

impl<'a, T: HeightVariable<T>, const SIZE: usize> HeightField<T, 3, SIZE>
    for WilsonField<'a, T, SIZE>
{
    fn normalize_random(&mut self) {
        self.field.normalize_random()
    }

    fn shift_values(&mut self, shift: T) {
        self.field.shift_values(shift)
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
        if self.links.is_active(index, direction) {
            let direction_modifier: T = T::from(match direction {
                1 => 1_i8,
                4 => (-1_i8) as i8,
                _ => 0_i8,
            });
            let spin_modifier: T = T::from(match self.spin {
                true => 1_i8,
                false => -1_i8,
            });
            let neighbour_index: usize = self.field.lattice.values[index][direction];
            let diff: T =
                value - self.field.values[neighbour_index] + direction_modifier * spin_modifier;
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
    let field: Field<i8, 3, SIZE> = Field::random(&lattice);
    let field: Field<i32, 3, SIZE> = Field::from_field(field);
    let field: WilsonField<i32, SIZE> = WilsonField::from_field(field, 1, 1);

    // Checking only the given links to be active
    for index in 0..SIZE {
        for direction in 0..(3 * 2_usize) {
            match (index, direction) {
                (0, 1) | (2, 4) => assert!(field.links.values[index][direction]),
                (_, _) => assert!(!field.links.values[index][direction]),
            }
        }
    }

    // Checking that each bond only ever evaluates to one value.
    for index in 0..SIZE {
        for direction in 0..(3 * 2_usize) {
            let neighbour_index: usize = field.field.lattice.get_neighbours_array(index)[direction];
            let inverse_direction: usize = (direction + 3) % 6;
            assert_eq!(
                field.bond_action(index, direction),
                field.bond_action(neighbour_index, inverse_direction)
            );
        }
    }
}
