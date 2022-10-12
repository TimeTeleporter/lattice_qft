use std::ops::{Add, Div, Mul, Sub};

use rand::{rngs::ThreadRng, Rng};

use crate::{action::Action, field::Field};

use self::bonds::BondsField;

pub trait Cluster: Action {
    fn cluster_sweep(&mut self, temp: f64);
}

impl<'a, T, const D: usize, const SIZE: usize> Cluster for Field<'a, T, D, SIZE>
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
    fn cluster_sweep(&mut self, temp: f64) {
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

        let mut reflected: Field<T, D, SIZE> = self.clone();
        reflected.mirror_values(plane, modifier);
        let reflected: Field<T, D, SIZE> = reflected; // Make reflected unmutable

        // Activate links if they are on the same side of the plane:
        let mut bonds: BondsField<D, SIZE> = BondsField::new(self.lattice);
        // Going through the lattice sites...
        for index in 0..SIZE {
            // ...for each neighbour in positive coordinate direction...
            for (direction, &neighbour) in
                self.lattice.pos_neighbours_array(index).iter().enumerate()
            {
                // ...calculate both the normal and the reflected action of the
                // link between them.
                let action: T = <Self as Action>::calculate_link_action(
                    self.values[index],
                    self.values[neighbour],
                );
                let reflected_action: T = <Self as Action>::calculate_link_action(
                    reflected.values[index],
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

        // Build the clusters
        let clusters: Vec<Vec<usize>> = bonds.build_clusters();

        // For each cluster decide to flip it
        for cluster in clusters {
            let coin: bool = rng.gen();
            if coin {
                for index in cluster {
                    self.values[index] = reflected.values[index];
                }
            }
        }
    }
}

pub(self) mod bonds {
    use std::ops::Deref;

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
            let values: Vec<()> = vec![(); lattice.size.iter().product()];
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

        /// Build a cluster from activated bonds.
        pub fn build_clusters(&self) -> Vec<Vec<usize>> {
            let mut unvisited: Vec<bool> = vec![true; SIZE];
            let mut clusters: Vec<Vec<usize>> = Vec::new();
            for index in 0..SIZE {
                if unvisited[index] {
                    let mut cluster: Vec<usize> = Vec::new();
                    self.cluster_up(index, &mut unvisited, &mut cluster);
                    clusters.push(cluster);
                }
            }

            clusters
        }

        /// Recursively add sites to a cluster
        fn cluster_up(&self, index: usize, unvisited: &mut Vec<bool>, cluster: &mut Vec<usize>) {
            if unvisited[index] {
                unvisited[index] = false;
                cluster.push(index);
                for neighbour_index in self.0.values[index]
                    .iter()
                    .enumerate()
                    .filter(|(_, &bond)| bond)
                    .map(|(direction, _)| self.0.lattice.get_neighbours_array(index)[direction])
                {
                    self.cluster_up(neighbour_index, unvisited, cluster);
                }
            }
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

        let clusters: Vec<Vec<usize>> = bonds_field.build_clusters();

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

        let clusters: Vec<Vec<usize>> = bonds_field.build_clusters();

        let mut test: Vec<Vec<usize>> = Vec::new();
        let mut cluster: Vec<usize> = Vec::with_capacity(SIZE);
        for index in 0..SIZE {
            cluster.push(index);
        }
        test.push(cluster);

        assert_eq!(test.len(), clusters.len());
    }
}
