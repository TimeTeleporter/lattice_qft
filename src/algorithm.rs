use std::fmt::Display;

use rand::prelude::*;

use crate::fields::{Action, Field, HeightVariable, LinksField, WilsonField};

// - AlgorithmType ------------------------------------------------------------

/// Lists all availible algorithms
#[derive(Debug, Clone)]
pub enum AlgorithmType {
    Metropolis(Metropolis),
    Cluster(Cluster),
}

// Here we implement Display in order to convert into String.
impl Display for AlgorithmType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AlgorithmType::Metropolis(_) => write!(f, "Metropolis"),
            AlgorithmType::Cluster(_) => write!(f, "Cluster   "),
        }
    }
}

impl AlgorithmType {
    pub fn new_metropolis() -> AlgorithmType {
        AlgorithmType::Metropolis(Metropolis)
    }

    pub fn new_cluster() -> AlgorithmType {
        AlgorithmType::Cluster(Cluster)
    }
}

// - Algorithm ----------------------------------------------------------------

/// Requirements for algorithms
pub trait Algorithm<T, const D: usize, const SIZE: usize>
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    /// This abstract method takes references to self, a field and and a temp
    /// and performes a Markov chain step from one field configuration to the
    /// next. It returns a statistic as a usize.
    fn field_sweep(&self, field: &mut Field<T, D, SIZE>, temp: f64, rng: &mut ThreadRng) -> usize;
}

// Deconstructs the method call to their own struct.
impl<T, const D: usize, const SIZE: usize> Algorithm<T, D, SIZE> for AlgorithmType
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    fn field_sweep(&self, field: &mut Field<T, D, SIZE>, temp: f64, rng: &mut ThreadRng) -> usize {
        match self {
            AlgorithmType::Metropolis(algo) => algo.field_sweep(field, temp, rng),
            AlgorithmType::Cluster(algo) => algo.field_sweep(field, temp, rng),
        }
    }
}

// - Wilson -------------------------------------------------------------------

pub trait WilsonAlgorithm<T, const SIZE: usize>
where
    T: HeightVariable<T>,
{
    /// This abstract method takes references to self, a field and and a temp
    /// and performes a Markov chain step from one field configuration to the
    /// next. It returns a statistic as a usize.
    fn wilson_sweep(&self, field: &mut WilsonField<T, SIZE>, temp: f64) -> usize;
}

// Deconstructs the method call to their own struct.
impl<T, const SIZE: usize> WilsonAlgorithm<T, SIZE> for AlgorithmType
where
    T: HeightVariable<T>,
{
    fn wilson_sweep(&self, field: &mut WilsonField<T, SIZE>, temp: f64) -> usize {
        match self {
            AlgorithmType::Metropolis(algo) => algo.wilson_sweep(field, temp),
            AlgorithmType::Cluster(algo) => algo.wilson_sweep(field, temp),
        }
    }
}

// - Metropolis ---------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct Metropolis;

impl<T, const D: usize, const SIZE: usize> Algorithm<T, D, SIZE> for Metropolis
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    fn field_sweep(&self, field: &mut Field<T, D, SIZE>, temp: f64, rng: &mut ThreadRng) -> usize {
        let mut acceptance: usize = 0;
        for index in 0..SIZE {
            if Metropolis.metropolis_single(field, index, temp, rng) {
                acceptance += 1;
            };
        }
        acceptance
    }
}

impl Metropolis {
    fn metropolis_single<T: HeightVariable<T>, const D: usize, const SIZE: usize>(
        &self,
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

impl<T, const SIZE: usize> WilsonAlgorithm<T, SIZE> for Metropolis
where
    T: HeightVariable<T>,
{
    fn wilson_sweep(&self, field: &mut WilsonField<T, SIZE>, temp: f64) -> usize {
        let mut rng = ThreadRng::default();
        let mut acceptance: usize = 0;
        for index in 0..SIZE {
            if Metropolis.wilson_single(field, index, temp, &mut rng) {
                acceptance += 1;
            };
        }
        acceptance
    }
}

impl Metropolis {
    fn wilson_single<T: HeightVariable<T>, const SIZE: usize>(
        &self,
        field: &mut WilsonField<T, SIZE>,
        index: usize,
        temp: f64,
        rng: &mut ThreadRng,
    ) -> bool {
        // Initialize the change to be measured
        let coin: bool = rng.gen();
        let old_value: T = field.field.values[index];
        let new_value: T = match coin {
            true => field.field.values[index] + T::from(1_i8),
            false => field.field.values[index] - T::from(1_i8),
        };

        // Calculate the action of both possibilities
        let old_action = field.assumed_local_action(index, old_value);
        let new_action = field.assumed_local_action(index, new_value);

        // Accept the new action if its lower than the previous.
        // Else accept it with a proportional probability.
        let draw: f64 = rng.gen_range(0.0..=1.0);
        let prob: f64 = (Into::<f64>::into(old_action - new_action) * temp).exp();
        if draw <= prob {
            field.field.values[index] = new_value;
            return true;
        }
        false
    }
}

// - Cluster ------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct Cluster;

impl<T, const D: usize, const SIZE: usize> Algorithm<T, D, SIZE> for Cluster
where
    T: HeightVariable<T>,
    [(); D * 2_usize]:,
{
    /// Implementation of the cluster algorithm.
    fn field_sweep(&self, field: &mut Field<T, D, SIZE>, temp: f64, rng: &mut ThreadRng) -> usize {
        // Set the mirror plane randomly on a height value
        let plane: T = field.values[rng.gen_range(0..SIZE)];
        let modifier: T = match rng.gen::<bool>() {
            true => 0_i8.into(),
            false => match rng.gen::<bool>() {
                true => (-1_i8).into(),
                false => 1_i8.into(),
            },
        };

        // Initialize memory to save the activated links and set all to false
        let mut links: LinksField<D, SIZE> = LinksField::new(field.lattice);

        // Going through the lattice sites...
        for index in 0..SIZE {
            // ...for each neighbour in positive coordinate direction...
            for direction in 0..D {
                // ...calculate both the normal and the reflected action of the
                // link between them.
                let action: T = field.bond_action(index, direction);
                let reflected_action: T = field.assumed_bond_action(
                    index,
                    direction,
                    reflect_value(field.values[index], plane, modifier),
                );

                let action_difference: T = action - reflected_action;

                // Dont activate a link if both are on the same side.
                if action_difference >= 0_i8.into() {
                    continue;
                };

                // Activate the link with probability 1 - exp(S - S').
                let draw: f64 = rng.gen_range(0.0..=1.0);
                let prob: f64 = 1.0 - (Into::<f64>::into(action_difference) * temp).exp();
                if draw <= prob {
                    links.activate(index, direction)
                }
            }
        }

        // Build the clusters
        let clusters: Vec<Vec<usize>> = links.collect_clusters();

        let clusters_amount: usize = clusters.len();

        // For each cluster decide to flip it
        for cluster in clusters {
            let coin: bool = rng.gen();
            if coin {
                for index in cluster {
                    field.values[index] = reflect_value(field.values[index], plane, modifier);
                }
            }
        }

        clusters_amount
    }
}

/// Reflects a hieght value on a given plane, with the modifier moving the
/// plane by half steps.
fn reflect_value<T: HeightVariable<T>>(value: T, plane: T, modifier: T) -> T {
    Into::<T>::into(2_i8) * plane + modifier - value
}

impl<T: HeightVariable<T>, const SIZE: usize> WilsonAlgorithm<T, SIZE> for Cluster {
    fn wilson_sweep(&self, field: &mut WilsonField<T, SIZE>, temp: f64) -> usize {
        let _ = (field, temp);
        todo!()
    }
}
