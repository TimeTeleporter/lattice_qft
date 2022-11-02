use std::{ops::Deref, time::Instant};

use serde::{Deserialize, Serialize};

use crate::{
    action::Action,
    cluster::Cluster,
    error::ObsChain,
    field::Field,
    lattice::{Lattice, Lattice3d},
    metropolis::Metropolis,
};

/// Datatype to save and read simulation output.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimResult {
    name: String,
    pub temp: f64,
    burnin: usize,
    iterations: usize,
    observable: f64,
    error: Option<f64>,
}

impl SimResult {
    pub fn new(
        name: String,
        temp: f64,
        burnin: usize,
        iterations: usize,
        observable: f64,
        error: Option<f64>,
    ) -> Self {
        SimResult {
            name,
            temp,
            burnin,
            iterations,
            observable,
            error,
        }
    }

    pub fn set_error(&mut self, error_option: Option<f64>) {
        self.error = error_option;
    }
}

#[derive(Debug, Clone, Copy)]
pub enum SimulationType {
    ClusterSim,
    MetropolisSim,
}

impl SimulationType {
    fn single_sweep<const D: usize, const SIZE: usize>(
        &self,
        field: &mut Field<i32, D, SIZE>,
        temp: f64,
    ) where
        [(); D * 2_usize]:,
    {
        match self {
            SimulationType::ClusterSim => field.cluster_sweep(temp),
            SimulationType::MetropolisSim => field.metropolis_sweep(temp),
        }
    }
}

pub struct Simulation<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    name: String,
    sim_type: SimulationType,
    size_normalized: bool,
    lattice: &'a Lattice<D, SIZE>,
    temp: f64,
    burnin: usize,
    iterations: usize,
}

impl<'a, const D: usize, const SIZE: usize> Simulation<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    pub fn new(
        name: String,
        sim_type: SimulationType,
        size_normalized: bool,
        lattice: &'a Lattice<D, SIZE>,
        temp: f64,
        burnin: usize,
        iterations: usize,
    ) -> Self {
        Simulation {
            name,
            sim_type,
            size_normalized,
            lattice,
            temp,
            burnin,
            iterations,
        }
    }

    pub fn run(&self) -> (SimResult, ObsChain) {
        println!("Started {} {}", self.temp, self.name);
        let time = Instant::now();
        let field: Field<i8, D, SIZE> = Field::random(self.lattice);
        let mut field: Field<i32, D, SIZE> = Field::from_field(field);

        // Burnin: compute an amount of sweeps to achieve equilibrium
        for _step in 0..(self.burnin) {
            self.sim_type.single_sweep(&mut field, self.temp);
            field.normalize_random();
        }

        // After having reached equilibrium, for each consecutive field configuration
        // that is generated by the markov chain, we calculate the observable and average it.
        let mut observable_array: Vec<f64> = Vec::with_capacity(self.iterations); // Simulation observable arrray
        for _step in 0..(self.iterations) {
            self.sim_type.single_sweep(&mut field, self.temp);
            observable_array.push(match self.size_normalized {
                true => field.size_normalized_action_observable(self.temp),
                false => field.action_observable(self.temp),
            });
            field.normalize_random();
        }

        // Mean value of the lattice actions generated by the simulation
        let observable: f64 = observable_array.iter().sum::<f64>() / observable_array.len() as f64;

        // Save the field configurations in order to do further calculations
        let data: ObsChain = ObsChain::new(observable_array);

        // Construct the return data type
        let results: SimResult = self.produce_simresult(observable);

        let duration = time.elapsed().as_secs_f32();
        println!("{} {} took {} secs", self.temp, self.name, duration);
        println!("{:?}", results);

        (results, data)
    }

    fn produce_simresult(&self, observable: f64) -> SimResult {
        SimResult {
            name: self.name.clone(),
            temp: self.temp,
            burnin: self.burnin,
            iterations: self.iterations,
            observable,
            error: None,
        }
    }
}

pub struct Simulation3d<'a, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>(
    Simulation<'a, 3, { MAX_X * MAX_Y * MAX_T }>,
)
where
    [(); MAX_X * MAX_Y * MAX_T]:;

impl<'a, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    Simulation3d<'a, MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    pub fn new(
        name: String,
        sim_type: SimulationType,
        size_normalized: bool,
        lattice: &'a Lattice3d<MAX_X, MAX_Y, MAX_T>,
        temp: f64,
        burnin: usize,
        iterations: usize,
    ) -> Self {
        let sim: Simulation<3, { MAX_X * MAX_Y * MAX_T }> =
            Simulation::new(name, sim_type, size_normalized, lattice.deref(), temp, burnin, iterations);
        Simulation3d(sim)
    }
}

impl<'a, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Deref
    for Simulation3d<'a, MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    type Target = Simulation<'a, 3, { MAX_X * MAX_Y * MAX_T }>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
