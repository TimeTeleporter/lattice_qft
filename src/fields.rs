//! This module implements fields of values on lattices.

use std::{
    collections::VecDeque,
    fmt::Display,
    ops::{Add, Deref, DerefMut, Div, Mul, Sub},
};

use rand::{distributions::Standard, prelude::*};

use crate::lattice::{Lattice, LatticeCoords};

// - Field --------------------------------------------------------------------

/// A field is a set of values assigned to each lattice site.
#[derive(Debug, Clone)]
pub struct Field<'a, T, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    pub values: Vec<T>,
    pub lattice: &'a Lattice<D, SIZE>,
}

// - Field - Constructors -----------------------------------------------------

impl<'a, T, const D: usize, const SIZE: usize> Field<'a, T, D, SIZE>
where
    [(); D * 2_usize]:,
{
    /// Constructor for a new field on the lattice initialized to be zero everywhere.
    pub fn new(lattice: &'a Lattice<D, SIZE>) -> Self
    where
        T: Default,
    {
        let values: Vec<()> = vec![(); lattice.size.iter().product()];
        let values: Vec<T> = values.iter().map(|_| T::default()).collect();

        Field::<'a, T, D, SIZE> { values, lattice }
    }

    /// Constructor for a new field on the lattice initialized to be random everywhere.
    pub fn random(lattice: &'a Lattice<D, SIZE>) -> Self
    where
        Standard: Distribution<T>,
    {
        let values: Vec<()> = vec![(); lattice.size.iter().product()];
        let values: Vec<T> = values.iter().map(|_| random()).collect();

        Field::<'a, T, D, SIZE> { values, lattice }
    }

    /// Convert a field of one type into a field of a different type losslessly.
    pub fn from_field<U: Into<T>>(field: Field<'a, U, D, SIZE>) -> Self {
        let values = field.values.into_iter().map(|x| x.into()).collect();

        Field::<'a, T, D, SIZE> {
            values,
            lattice: field.lattice,
        }
    }
}

// - Field - Printing ---------------------------------------------------------

// Printing 3-dimensional fields nicely.
impl<T, const SIZE: usize> Field<'_, T, 3, SIZE>
where
    T: std::fmt::Debug,
{
    /// Prints the value in a intuitive way.
    pub fn print_values_formated(&self, size: [usize; 3]) {
        for t in 0..size[2] {
            println!("t = {}", t);
            for y in 0..size[1] {
                print!("[");
                for x in 0..size[0] {
                    if x == size[0] - 1 {
                        println!(
                            "{:?} ]",
                            self.values[self
                                .lattice
                                .calc_index_from_coords(LatticeCoords::new([x, y, t]))]
                        );
                    } else {
                        print!(
                            "{:?}, ",
                            self.values[self
                                .lattice
                                .calc_index_from_coords(LatticeCoords::new([x, y, t]))]
                        );
                    }
                }
            }
        }
    }
}

// Printing two-dimensional fields nicely.
impl<T, const SIZE: usize> Field<'_, T, 2, SIZE>
where
    T: std::fmt::Debug,
{
    /// Prints the value in a intuitive way.
    pub fn print_values_formated(&self, size: [usize; 2]) {
        for y in 0..size[1] {
            print!("[");
            for x in 0..size[0] {
                if x == size[0] - 1 {
                    println!(
                        "{:?} ]",
                        self.values[self
                            .lattice
                            .calc_index_from_coords(LatticeCoords::new([x, y]))]
                    );
                } else {
                    print!(
                        "{:?}, ",
                        self.values[self
                            .lattice
                            .calc_index_from_coords(LatticeCoords::new([x, y]))]
                    );
                }
            }
        }
    }
}

// - Bonds --------------------------------------------------------------------

/// Models a field that has values on the bonds between sites.
pub struct BondsField<'a, T, const D: usize, const SIZE: usize>(
    Field<'a, [T; D * 2_usize], D, SIZE>,
)
where
    [(); D * 2_usize]:;

// - Bonds - Constructors -----------------------------------------------------

impl<'a, T, const D: usize, const SIZE: usize> BondsField<'a, T, D, SIZE>
where
    [(); D * 2_usize]:,
{
    /// Constructor for a new BondsField with the default value.
    pub fn new(lattice: &'a Lattice<D, SIZE>) -> Self
    where
        T: Default,
    {
        let values: Vec<()> = vec![(); SIZE];
        let values: Vec<[T; D * 2_usize]> = values
            .into_iter()
            .map(|_| {
                let ary: [(); D * 2_usize] = [(); D * 2_usize];
                let ary: [T; D * 2_usize] = ary.map(|_| T::default());
                ary
            })
            .collect();

        BondsField(Field::<'a, [T; D * 2_usize], D, SIZE> { values, lattice })
    }

    pub fn random(lattice: &'a Lattice<D, SIZE>) -> Self
    where
        Standard: Distribution<T>,
    {
        let values: Vec<()> = vec![(); SIZE];
        let values: Vec<[T; D * 2_usize]> = values
            .into_iter()
            .map(|_| {
                let ary: [(); D * 2_usize] = [(); D * 2_usize];
                let ary: [T; D * 2_usize] = ary.map(|_| random());
                ary
            })
            .collect();

        BondsField(Field::<'a, [T; D * 2_usize], D, SIZE> { values, lattice })
    }

    /// Convert a field of one type into a field of a different type losslessly.
    pub fn from_field<U: Into<T>>(field: BondsField<'a, U, D, SIZE>) -> Self {
        let lattice: &'a Lattice<D, SIZE> = field.lattice;

        let mut input: Vec<[U; D * 2_usize]> = field.0.values;
        input.reverse();

        let mut values: Vec<[T; D * 2_usize]> = Vec::with_capacity(SIZE);
        while let Some(ary) = input.pop() {
            let ary: [T; D * 2_usize] = ary.map(|x| x.into());
            values.push(ary);
        }
        BondsField(Field::<'a, [T; D * 2_usize], D, SIZE> { values, lattice })
    }
}

// - Bonds - Deref ------------------------------------------------------------

impl<'a, T, const D: usize, const SIZE: usize> Deref for BondsField<'a, T, D, SIZE>
where
    [(); D * 2_usize]:,
{
    type Target = Field<'a, [T; D * 2_usize], D, SIZE>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<'a, T, const D: usize, const SIZE: usize> DerefMut for BondsField<'a, T, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

// - Bonds - Auxilliary -------------------------------------------------------

impl<'a, T, const D: usize, const SIZE: usize> BondsField<'a, T, D, SIZE>
where
    T: Clone,
    [(); D * 2_usize]:,
{
    /// Activate a link of the LinksField
    pub fn set_value(&mut self, index: usize, direction: usize, value: T) {
        let neighbour = self.lattice.get_neighbours_array(index)[direction];
        self.0.values[index][direction] = value.clone();
        self.0.values[neighbour][(direction + D) % (D * 2_usize)] = value;
    }
}

impl<'a, T, const D: usize, const SIZE: usize> BondsField<'a, T, D, SIZE>
where
    T: PartialEq + std::fmt::Debug,
    [(); D * 2_usize]:,
{
    /// A method to check if the links, which are saved on both of the
    /// sites, are correct at both points.
    pub fn check_coherence(&self) {
        for (index, links) in self.values.iter().enumerate() {
            for (direction, link) in links.iter().enumerate() {
                assert_eq!(
                    *link,
                    self.values[self.lattice.get_neighbours_array(index)[direction]]
                        [(direction + D) % (D * 2)]
                );
            }
        }
    }
}

// - Links --------------------------------------------------------------------

/// Models the activation of outgoing links from a lattice site
pub struct LinksField<'a, const D: usize, const SIZE: usize>(BondsField<'a, bool, D, SIZE>)
where
    [(); D * 2_usize]:;

impl<'a, const D: usize, const SIZE: usize> LinksField<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    /// Constructor for a new LinksField on the lattice initialized to be false everywhere.
    pub fn new(lattice: &'a Lattice<D, SIZE>) -> Self {
        assert_eq!(bool::default(), false);
        LinksField(BondsField::<bool, D, SIZE>::new(lattice))
    }
}

impl<'a, const D: usize, const SIZE: usize> Deref for LinksField<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    type Target = BondsField<'a, bool, D, SIZE>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<'a, const D: usize, const SIZE: usize> LinksField<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    /// Activate a link of the LinksField
    pub fn activate(&mut self, index: usize, direction: usize) {
        self.0.set_value(index, direction, true);
    }

    pub fn is_active(&self, index: usize, direction: usize) -> bool {
        self.values[index][direction]
    }

    /// Return collection of all clusters
    pub fn collect_clusters(&self) -> Vec<Vec<usize>> {
        // Remember if a given site has already been considered for a cluster
        let mut unvisited: Vec<bool> = vec![true; SIZE];
        let mut clusters: Vec<Vec<usize>> = Vec::new();
        for index in 0..SIZE {
            if unvisited[index] {
                // Build a new cluster if you come across a site not yet in one.
                clusters.push(self.build_cluster(index, &mut unvisited));
            }
        }
        clusters
    }

    /// Build a cluster from activated links
    fn build_cluster(&self, index: usize, unvisited: &mut Vec<bool>) -> Vec<usize> {
        let mut cluster: Vec<usize> = Vec::with_capacity(SIZE / 2);
        let mut checklist: VecDeque<usize> = VecDeque::with_capacity(SIZE / 2);
        // Add the new site to a checklist for sites to check for
        // activated, unvisited neighbours to be added to the cluster
        checklist.push_back(index);
        unvisited[index] = false;

        // Continuously working through the checklist, for each entry...
        while let Some(check) = checklist.pop_front() {
            // ...add them to the cluster list and...
            cluster.push(check);

            // ...check each neighbour...
            for (direction, neighbour) in self
                .lattice
                .get_neighbours_array(check)
                .into_iter()
                .enumerate()
            {
                // ...if it hasn't been visited and the link is active,
                // add it to the checklist and mark it as visited.
                if unvisited[neighbour] && self.values[check][direction] {
                    checklist.push_back(neighbour);
                    unvisited[neighbour] = false;
                }
            }
        }

        cluster.shrink_to_fit();
        cluster
    }
}

// - WilsonField --------------------------------------------------------------

pub struct WilsonField<'a, T, const SIZE: usize> {
    pub field: Field<'a, T, 3, SIZE>,
    links: LinksField<'a, 3, SIZE>,
    spin: bool,
}

impl<'a, T, const SIZE: usize> WilsonField<'a, T, SIZE> {
    pub fn new(
        lattice: &'a Lattice<3, SIZE>,
        width: usize,
        height: usize,
    ) -> WilsonField<'a, T, SIZE>
    where
        T: Default,
    {
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

// - HeightVariable -----------------------------------------------------------

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

// - HeightField --------------------------------------------------------------

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

// - Action -------------------------------------------------------------------

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

// - Field - Tests ------------------------------------------------------------

#[test]
fn test_field_conversion() {
    let lattice: Lattice<3, 8> = Lattice::new([2, 2, 2]);

    let field8: Field<i8, 3, 8> = Field::new(&lattice);

    let field16: Field<i16, 3, 8> = Field::new(&lattice);

    let field: Field<i16, 3, 8> = Field::from_field(field8);

    assert_eq!(field16.values, field.values);
}

// - Links - Tests ------------------------------------------------------------

#[test]
fn test_links_field() {
    const D: usize = 3;
    const SIZE_ARY: [usize; D] = [4, 5, 3];
    const SIZE: usize = SIZE_ARY[0] * SIZE_ARY[1] * SIZE_ARY[2];

    let lattice: Lattice<D, SIZE> = Lattice::new(SIZE_ARY);
    let mut links_field: LinksField<D, SIZE> = LinksField::new(&lattice);

    // Setting the link at [0, 2, 1] (index = 28) in the third direction
    // (direction = 2).
    assert_eq!(links_field.values[28][2], false);
    assert_eq!(links_field.values[48][5], false);
    links_field.activate(28, 2);
    //links_field.print_values_formated(SIZE_ARY);
    assert_eq!(links_field.values[28][2], true);
    assert_eq!(links_field.values[48][5], true);

    // Setting the link at [2, 1, 0] (index = 6) in the negatice first direction
    // (direction = 3).
    assert_eq!(links_field.values[6][3], false);
    assert_eq!(links_field.values[5][0], false);
    links_field.activate(6, 3);
    assert_eq!(links_field.values[6][3], true);
    assert_eq!(links_field.values[5][0], true);
}

#[test]
fn test_all_links_deactivated() {
    const D: usize = 3;
    const SIZE_ARY: [usize; D] = [4, 5, 3];
    const SIZE: usize = SIZE_ARY[0] * SIZE_ARY[1] * SIZE_ARY[2];

    let lattice: Lattice<D, SIZE> = Lattice::new(SIZE_ARY);
    let links_field: LinksField<D, SIZE> = LinksField::new(&lattice);

    let clusters: Vec<Vec<usize>> = links_field.collect_clusters();

    let mut test: Vec<Vec<usize>> = Vec::new();
    for index in 0..SIZE {
        let mut vect: Vec<usize> = Vec::new();
        vect.push(index);
        test.push(vect);
    }

    assert_eq!(test, clusters);
}

#[test]
fn test_all_links_active() {
    const D: usize = 3;
    const SIZE_ARY: [usize; D] = [4, 5, 3];
    const SIZE: usize = SIZE_ARY[0] * SIZE_ARY[1] * SIZE_ARY[2];

    let lattice: Lattice<D, SIZE> = Lattice::new(SIZE_ARY);
    let mut links_field: LinksField<D, SIZE> = LinksField::new(&lattice);

    for links in links_field.0 .0.values.iter_mut() {
        for link in links.iter_mut() {
            *link = true;
        }
    }

    let clusters: Vec<Vec<usize>> = links_field.collect_clusters();

    let mut test: Vec<Vec<usize>> = Vec::new();
    let mut cluster: Vec<usize> = Vec::with_capacity(SIZE);
    for index in 0..SIZE {
        cluster.push(index);
    }
    test.push(cluster);

    assert_eq!(test.len(), clusters.len());
}

// - HeightField - Tests ------------------------------------------------------

#[test]
fn test_field_shift() {
    use crate::lattice::Lattice;

    let lattice: Lattice<4, 81> = Lattice::new([3, 3, 3, 3]);

    let mut field: Field<i8, 4, 81> = Field::new(&lattice);

    field.shift_values(3);

    assert_eq!(field.values[80], -3);
}

// - Action - Tests -----------------------------------------------------------

#[test]
fn test_wilson_field_action() {
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
