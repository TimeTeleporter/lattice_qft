use std::fmt::Display;

use rand::prelude::*;
use serde::{Deserialize, Serialize};

use crate::{
    field::Field,
    heightfield::{Action, HeightField, HeightVariable},
    pause,
};

use self::bonds::BondsField;

const VERBOSE: bool = false;

#[derive(Debug, Serialize, Deserialize, Clone, Copy)]
pub enum AlgorithmType {
    Metropolis,
    Cluster,
}

impl Display for AlgorithmType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AlgorithmType::Metropolis => write!(f, "Metropolis"),
            AlgorithmType::Cluster => write!(f, "Cluster   "),
        }
    }
}

pub trait Algorithm<T, const D: usize, const SIZE: usize>
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    fn field_sweep(&self, field: &mut Field<T, D, SIZE>, temp: f64) -> usize;
}

impl<T, const D: usize, const SIZE: usize> Algorithm<T, D, SIZE> for AlgorithmType
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    fn field_sweep(&self, field: &mut Field<T, D, SIZE>, temp: f64) -> usize {
        match self {
            AlgorithmType::Metropolis => Metropolis::field_sweep(&Metropolis, field, temp),
            AlgorithmType::Cluster => Cluster::field_sweep(&Cluster, field, temp),
        }
    }
}

struct Metropolis;
struct Cluster;

impl<T, const D: usize, const SIZE: usize> Algorithm<T, D, SIZE> for Metropolis
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    fn field_sweep(&self, field: &mut Field<T, D, SIZE>, temp: f64) -> usize {
        let mut rng = ThreadRng::default();
        let mut acceptance: usize = 0;
        for index in 0..SIZE {
            if Metropolis::metropolis_single(field, index, temp, &mut rng) {
                acceptance += 1;
            };
        }
        acceptance
    }
}

impl Metropolis {
    fn metropolis_single<T: HeightVariable<T>, const D: usize, const SIZE: usize>(
        field: &mut Field<T, D, SIZE>,
        index: usize,
        temp: f64,
        rng: &mut ThreadRng,
    ) -> bool
    where
        [(); D * 2_usize]:,
    {
        // Initialize the change to be measured
        let coin: bool = rng.gen();
        let old_value: T = field.values[index];
        let new_value: T = match coin {
            true => field.values[index] + T::from(1_i8),
            false => field.values[index] - T::from(1_i8),
        };

        // Calculate the action of both possibilities
        let old_action = field.assumed_local_action(index, old_value);
        let new_action = field.assumed_local_action(index, new_value);

        // Accept the new action if its lower than the previous.
        // Else accept it with a proportional probability.
        let draw: f64 = rng.gen_range(0.0..=1.0);
        let prob: f64 = (Into::<f64>::into(old_action - new_action) * temp).exp();
        if draw <= prob {
            field.values[index] = new_value;
            return true;
        }
        false
    }
}

impl<T, const D: usize, const SIZE: usize> Algorithm<T, D, SIZE> for Cluster
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    /// Implementation of the cluster algorithm.
    fn field_sweep(&self, field: &mut Field<T, D, SIZE>, temp: f64) -> usize {
        let mut rng = ThreadRng::default();

        // Set the mirror plane randomly on a height value
        let plane: T = field.values[rng.gen_range(0..SIZE)];
        let modifier: T = match rng.gen::<bool>() {
            true => 0_i8.into(),
            false => match rng.gen::<bool>() {
                true => (-1_i8).into(),
                false => 1_i8.into(),
            },
        };
        if VERBOSE {
            println!("Mirror plane set on {plane} and {modifier}");
        }

        // Initialize memory to save the activated bonds and set all to false
        let mut bonds: BondsField<D, SIZE> = BondsField::new(field.lattice);

        // Going through the lattice sites...
        for index in 0..SIZE {
            // ...for each neighbour in positive coordinate direction...
            for (direction, neighbour_index) in field
                .lattice
                .pos_neighbours_array(index)
                .into_iter()
                .enumerate()
            {
                // ...calculate both the normal and the reflected action of the
                // link between them.
                let neighbour: T = field.values[neighbour_index];
                let standard_action: T = field.bond_action(index, direction);
                let reflected_action: T = field.assumed_bond_action(
                    index,
                    direction,
                    <Field<T, D, SIZE> as HeightField<T, D, SIZE>>::reflect_value(
                        neighbour, plane, modifier,
                    ),
                );

                let action_difference: T = standard_action - reflected_action;

                // Dont activate a bond if both are on the same side.
                if action_difference >= Into::<T>::into(0_i8) {
                    continue;
                };

                // Activate the link with probability 1 - exp(S - S').
                let draw: f64 = rng.gen_range(0.0..=1.0);
                let prob: f64 = 1.0 - (Into::<f64>::into(action_difference) * temp).exp();
                if draw <= prob {
                    bonds.activate(index, direction)
                }
            }
        }

        // Build the clusters
        let clusters: Vec<Vec<usize>> = bonds.collect_clusters();

        let clusters_amount: usize = clusters.len();

        // For each cluster decide to flip it
        for cluster in clusters {
            if VERBOSE {
                println!("Looking at the cluster {:?}", cluster);
                pause();
                for &index in cluster.iter() {
                    let neighbours: [usize; D * 2_usize] =
                        field.lattice.get_neighbours_array(index);
                    println!(
                        "{index} has the neighbours {:?} and bonds {:?}",
                        neighbours, bonds.values[index]
                    );
                    for (direction, &bond) in bonds.values[index].iter().enumerate() {
                        let neighbour = neighbours[direction];
                        if bond {
                            assert!(cluster.contains(&neighbour));
                            println!("We found {neighbour} to be also in the cluster");
                        }
                    }
                }
            }

            let coin: bool = rng.gen();
            if coin {
                for index in cluster {
                    field.values[index] =
                        <Field<T, D, SIZE> as HeightField<T, D, SIZE>>::reflect_value(
                            field.values[index],
                            plane,
                            modifier,
                        );
                }
            }
        }

        clusters_amount
    }
}

pub(self) mod bonds {
    use std::{collections::VecDeque, ops::Deref};

    use crate::{field::Field, lattice::Lattice, pause};

    use super::VERBOSE;

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

        /// Build a cluster from activated bonds
        fn build_cluster(&self, index: usize, unvisited: &mut Vec<bool>) -> Vec<usize> {
            let mut cluster: Vec<usize> = Vec::with_capacity(SIZE / 2);
            let mut checklist: VecDeque<usize> = VecDeque::with_capacity(SIZE / 2);
            // Add the new site to a checklist for sites to check for
            // activated, unvisited neighbours to be added to the cluster
            if VERBOSE {
                println!("Starting to build the cluster around {index}");
            }
            checklist.push_back(index);
            unvisited[index] = false;

            // Continuously working through the checklist, for each entry...
            while let Some(check) = checklist.pop_front() {
                // ...add them to the cluster list and...
                cluster.push(check);
                if VERBOSE {
                    pause();
                    println!(
                        "Looking at {check}, with bonds {:?} and neighbours {:?}",
                        self.0.values[check],
                        self.0.lattice.get_neighbours_array(check)
                    );
                }

                // ...check each neighbour...
                for (direction, neighbour) in self
                    .0
                    .lattice
                    .get_neighbours_array(check)
                    .into_iter()
                    .enumerate()
                {
                    // ...if it hasn't been visited and the link is active,
                    // add it to the checklist and mark it as visited.
                    if unvisited[neighbour] && self.0.values[check][direction] {
                        if VERBOSE {
                            println!("{neighbour} gets added to the checklist :)");
                        }
                        checklist.push_back(neighbour);
                        unvisited[neighbour] = false;
                    } else if VERBOSE {
                        println!("{neighbour} is not added to the checklist :(")
                    }
                }
            }

            cluster.shrink_to_fit();
            cluster
        }

        /// A method to check if the bonds, which are saved on both of the
        /// sites, are correct at both points.
        #[allow(dead_code)]
        pub fn check_coherence(&self) {
            for (index, bonds) in self.0.values.iter().enumerate() {
                for (direction, &bond) in bonds.iter().enumerate() {
                    assert_eq!(
                        bond,
                        self.0.values[self.0.lattice.get_neighbours_array(index)[direction]]
                            [(direction + D) % (D * 2)]
                    );
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
