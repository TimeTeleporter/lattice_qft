use std::ops::{Add, Div, Mul, Sub};

use rand::{rngs::ThreadRng, Rng};

use crate::{action::Action, field::Field};

use self::bonds::BondsField;

pub trait Cluster<T>: Action {
    fn cluster_sweep(&mut self, temp: f64) -> usize;

    fn reflect_value(&self, index: usize, plane: T, modifier: T) -> T;
}

impl<'a, T, const D: usize, const SIZE: usize> Cluster<T> for Field<'a, T, D, SIZE>
where
    f64: From<T>,
    T: Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Div<Output = T>
        + Default
        + From<i8>
        + Into<i64>
        + Into<f64>
        + PartialOrd
        + Copy,
    [(); D * 2_usize]:,
{
    /// Implementation of the cluster algortihm. It identifies a mirror plane
    /// and constructs clusters of lattice sites with values on the same side
    /// of the mirror. For each cluster it is randomly decided if all of the
    /// values of its sites are mirrored or not.
    fn cluster_sweep(&mut self, temp: f64) -> usize{
        let mut rng = ThreadRng::default();

        // Set the mirror plane randomly on a height value
        let plane: T = self.values[rng.gen_range(0..SIZE)];
        let modifier: T = match rng.gen::<bool>() {
            true => 0_i8.into(),
            false => match rng.gen::<bool>() {
                true => (-1_i8).into(),
                false => 1_i8.into(),
            },
        };
        //println!("Plane and modifier set");

        // Initialize memory to save the activated bonds and set all to false
        let mut bonds: BondsField<D, SIZE> = BondsField::new(self.lattice);
        //println!("Bonds initialized");

        // Going through the lattice sites...
        for index in 0..SIZE {
            // ...for each neighbour in positive coordinate direction...
            for (direction, neighbour) in self
                .lattice
                .pos_neighbours_array(index)
                .into_iter()
                .enumerate()
            {
                // ...calculate both the normal and the reflected action of the
                // link between them.
                let action: T = <Self as Action>::calculate_link_action(
                    self.values[index],
                    self.values[neighbour],
                );
                let reflected_action: T = <Self as Action>::calculate_link_action(
                    self.reflect_value(index, plane, modifier),
                    self.values[neighbour],
                );

                // Dont activate a bond if both are on the same side.
                if action >= reflected_action {
                    continue;
                }

                // Activate the link with probability 1 - exp(S - S').
                let draw: f64 = rng.gen_range(0.0..=1.0);
                let prob: f64 = 1.0 - (Into::<f64>::into(action - reflected_action) * temp).exp();
                if draw <= prob {
                    bonds.activate(index, direction)
                }
            }
        }
        //println!("Bonds activated");

        // Build the clusters
        let clusters: Vec<Vec<usize>> = bonds.collect_clusters();
        let clusters_amount: usize = clusters.len();

        // For each cluster decide to flip it
        for cluster in clusters {
            let coin: bool = rng.gen();
            if coin {
                for index in cluster {
                    self.values[index] = self.reflect_value(index, plane, modifier);
                }
            }
        }

        clusters_amount
    }

    fn reflect_value(&self, index: usize, plane: T, modifier: T) -> T {
        <i8 as Into<T>>::into(2_i8) * plane + modifier - self.values[index]
    }
}

pub(self) mod bonds {
    use std::{collections::VecDeque, ops::Deref};

    use crate::{field::Field, lattice::Lattice};

    /// Models the activation of outgoing bonds from a lattice site
    pub struct BondsField<'a, const D: usize, const SIZE: usize>(
        Field<'a, [bool; D * 2_usize], D, SIZE>,
    )
    where
        [(); D * 2_usize]:;

    impl<'a, const D: usize, const SIZE: usize> BondsField<'a, D, SIZE>
    where
        [(); D * 2_usize]:,
    {
        /// Constructor for a new BondsField on the lattice initialized to be false everywhere.
        pub fn new(lattice: &'a Lattice<D, SIZE>) -> Self {
            let values: Vec<()> = vec![(); SIZE];
            let values: Vec<[bool; D * 2_usize]> =
                values.into_iter().map(|_| [false; D * 2_usize]).collect();

            BondsField(Field::<'a, [bool; D * 2_usize], D, SIZE> { values, lattice })
        }

        /// Activate a link of the BondsField
        pub fn activate(&mut self, index: usize, direction: usize) {
            let neighbour = self.0.lattice.get_neighbours_array(index)[direction];
            self.0.values[index][direction] = true;
            self.0.values[neighbour][(direction + D) % (D * 2_usize)] = true;
        }

        /// Return collection of all clusters
        pub fn collect_clusters(self) -> Vec<Vec<usize>> {
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

        /// Build a cluster from activated bonds
        fn build_cluster(&self, index: usize, unvisited: &mut Vec<bool>) -> Vec<usize> {
            let mut cluster: Vec<usize> = Vec::with_capacity(SIZE);
            let mut checklist: VecDeque<usize> = VecDeque::with_capacity(SIZE);
            // Add the new site to a checklist for sites to check for
            // activated, unvisited neighbours to be added to the cluster
            checklist.push_back(index);
            unvisited[index] = false;

            // Continuously working through the checklist, for each entry...
            while let Some(index) = checklist.pop_front() {
                // ...add them to the cluster list and...
                cluster.push(index);

                // ...check each neighbour...
                for (direction, neighbour) in self
                    .0
                    .lattice
                    .get_neighbours_array(index)
                    .into_iter()
                    .enumerate()
                {
                    // ...if it hasn't been visited and the link is active,
                    // add it to the checklist and mark it as visited.
                    if self.0.values[neighbour][direction] && unvisited[neighbour] {
                        checklist.push_back(neighbour);
                        unvisited[neighbour] = false;
                    }
                }
            }

            cluster.shrink_to_fit();
            cluster
        }
    }

    impl<'a, const D: usize, const SIZE: usize> Deref for BondsField<'a, D, SIZE>
    where
        [(); D * 2_usize]:,
    {
        type Target = Field<'a, [bool; D * 2_usize], D, SIZE>;

        fn deref(&self) -> &Self::Target {
            &self.0
        }
    }

    #[test]
    fn test_bonds_field() {
        const D: usize = 3;
        const SIZE_ARY: [usize; D] = [4, 5, 3];
        const SIZE: usize = SIZE_ARY[0] * SIZE_ARY[1] * SIZE_ARY[2];

        let lattice: Lattice<D, SIZE> = Lattice::new(SIZE_ARY);
        let mut bonds_field: BondsField<D, SIZE> = BondsField::new(&lattice);

        // Setting the bond at [0, 2, 1] (index = 28) in the third direction
        // (direction = 2).
        assert_eq!(bonds_field.values[28][2], false);
        assert_eq!(bonds_field.values[48][5], false);
        bonds_field.activate(28, 2);
        //bonds_field.print_values_formated(SIZE_ARY);
        assert_eq!(bonds_field.values[28][2], true);
        assert_eq!(bonds_field.values[48][5], true);

        // Setting the bond at [2, 1, 0] (index = 6) in the negatice first direction
        // (direction = 3).
        assert_eq!(bonds_field.values[6][3], false);
        assert_eq!(bonds_field.values[5][0], false);
        bonds_field.activate(6, 3);
        assert_eq!(bonds_field.values[6][3], true);
        assert_eq!(bonds_field.values[5][0], true);
    }

    #[test]
    fn test_all_bonds_deactivated() {
        const D: usize = 3;
        const SIZE_ARY: [usize; D] = [4, 5, 3];
        const SIZE: usize = SIZE_ARY[0] * SIZE_ARY[1] * SIZE_ARY[2];

        let lattice: Lattice<D, SIZE> = Lattice::new(SIZE_ARY);
        let bonds_field: BondsField<D, SIZE> = BondsField::new(&lattice);

        let clusters: Vec<Vec<usize>> = bonds_field.collect_clusters();

        let mut test: Vec<Vec<usize>> = Vec::new();
        for index in 0..SIZE {
            let mut vect: Vec<usize> = Vec::new();
            vect.push(index);
            test.push(vect);
        }

        assert_eq!(test, clusters);
    }

    #[test]
    fn test_all_bonds_active() {
        const D: usize = 3;
        const SIZE_ARY: [usize; D] = [4, 5, 3];
        const SIZE: usize = SIZE_ARY[0] * SIZE_ARY[1] * SIZE_ARY[2];

        let lattice: Lattice<D, SIZE> = Lattice::new(SIZE_ARY);
        let mut bonds_field: BondsField<D, SIZE> = BondsField::new(&lattice);

        for bonds in bonds_field.0.values.iter_mut() {
            for bond in bonds.iter_mut() {
                *bond = true;
            }
        }

        let clusters: Vec<Vec<usize>> = bonds_field.collect_clusters();

        let mut test: Vec<Vec<usize>> = Vec::new();
        let mut cluster: Vec<usize> = Vec::with_capacity(SIZE);
        for index in 0..SIZE {
            cluster.push(index);
        }
        test.push(cluster);

        assert_eq!(test.len(), clusters.len());
    }
}

#[test]
fn test_field_mirror() {
    use crate::action::Action;
    use crate::lattice::Lattice;

    const D: usize = 4;
    const SIZE_ARY: [usize; D] = [3, 3, 3, 3];
    const SIZE: usize = SIZE_ARY[0] * SIZE_ARY[1] * SIZE_ARY[2] * SIZE_ARY[3];

    let lattice: Lattice<D, SIZE> = Lattice::new(SIZE_ARY);
    let field: Field<i8, D, SIZE> = Field::random(&lattice);
    let mut field: Field<i32, D, SIZE> = Field::from_field(field);

    let old_action: i64 = field.sum_link_actions();

    for index in 0..SIZE {
        field.values[index] = field.reflect_value(index, 0, 0);
    }

    assert_eq!(old_action, field.sum_link_actions());
}
