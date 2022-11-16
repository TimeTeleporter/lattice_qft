use std::{error::Error, time::Instant};

use serde::{Deserialize, Serialize};

use crate::{
    action::Action,
    cluster::Cluster,
    error::ObsChain,
    field::Field,
    lattice::{Lattice, Lattice3d},
    metropolis::Metropolis,
};

#[derive(Debug, Serialize, Deserialize)]
pub enum ComputationType {
    ClusterSimulation,
    MetropolisSimulation,
    Test,
}

#[derive(Debug, Serialize, Deserialize, Clone, Copy)]
#[serde(tag = "type")]
pub enum Observable {
    Action,
    Wilson { width: usize, height: usize },
    SizeNormalized,
}

impl Observable {
    fn observe<const D: usize, const SIZE: usize>(
        &self,
        field: &Field<i32, D, SIZE>,
        temp: f64,
    ) -> f64
    where
        [(); D * 2_usize]:,
    {
        match self {
            Observable::Action => field.action_observable(temp),
            Observable::Wilson { width, height } => {
                field.wilson_loop_observable(temp, *width, *height)
            }
            Observable::SizeNormalized => field.size_normalized_action_observable(temp),
        }
    }
}

impl ToString for Observable {
    fn to_string(&self) -> String {
        match self {
            Observable::Action => "Action observable".to_string(),
            Observable::Wilson { width, height } => format!("{width}x{height} Wilson loop"),
            Observable::SizeNormalized => "Action observable size normalized".to_string(),
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ComputationResults {
    dimensions: usize,
    size: usize,
    x: Option<usize>,
    y: Option<usize>,
    t: Option<usize>,
    computation: ComputationType,
    observable: String,
    temp: f64,
    range: Option<usize>,
    burnin: Option<usize>,
    iterations: Option<usize>,
    output: f64,
    error: Option<f64>,
}

#[derive(Debug)]
pub struct Computation<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    lattice: &'a Lattice<D, SIZE>,
    x: Option<usize>,
    y: Option<usize>,
    t: Option<usize>,
    computation: ComputationType,
    observable: Observable,
    temp: f64,
    range: Option<usize>,
    burnin: Option<usize>,
    iterations: Option<usize>,
}

pub trait Computatable<'a, const D: usize, const SIZE: usize>
where
    Self: Sized,
    [(); D * 2_usize]:,
{
    fn from_computation(computation: Computation<'a, D, SIZE>) -> Result<Self, Box<dyn Error>>;
    fn run(self) -> Result<ComputationResults, Box<dyn Error>>;
    fn into_result(self, observable: f64) -> Result<ComputationResults, Box<dyn Error>>;
}

impl<'a, const D: usize, const SIZE: usize> Computatable<'a, D, SIZE> for Computation<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn from_computation(computation: Computation<'a, D, SIZE>) -> Result<Self, Box<dyn Error>> {
        Ok(computation)
    }

    fn run(self) -> Result<ComputationResults, Box<dyn Error>> {
        match self.computation {
            ComputationType::ClusterSimulation | ComputationType::MetropolisSimulation => {
                Simulation::<D, SIZE>::from_computation(self)?.run()
            }
            ComputationType::Test => Test::<D, SIZE>::from_computation(self)?.run(),
        }
    }

    fn into_result(self, observable: f64) -> Result<ComputationResults, Box<dyn Error>> {
        match self.computation {
            ComputationType::ClusterSimulation | ComputationType::MetropolisSimulation => {
                Simulation::from_computation(self)?.into_result(observable)
            }
            ComputationType::Test => Test::from_computation(self)?.into_result(observable),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum Algorithm {
    Cluster,
    Metropolis,
}

impl Algorithm {
    /// Performs a sweep with the given algorithm. Returns the corresponding
    /// statistics.
    fn field_sweep<const D: usize, const SIZE: usize>(
        &self,
        field: &mut Field<i32, D, SIZE>,
        temp: f64,
    ) -> usize
    where
        [(); D * 2_usize]:,
    {
        // Here we may catch the return values of the sweeps.
        match self {
            Algorithm::Cluster => field.cluster_sweep(temp),
            Algorithm::Metropolis => field.metropolis_sweep(temp),
        }
    }
}

impl Algorithm {
    /// A method to convert Algorithm into ComputationType.
    fn into_computation_type(self) -> ComputationType {
        match self {
            Algorithm::Cluster => ComputationType::ClusterSimulation,
            Algorithm::Metropolis => ComputationType::MetropolisSimulation,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Simulation<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    lattice: &'a Lattice<D, SIZE>,
    x: Option<usize>,
    y: Option<usize>,
    t: Option<usize>,
    algorithm: Algorithm,
    observable: Observable,
    temp: f64,
    burnin: usize,
    iterations: usize,
}

impl<'a, const D: usize, const SIZE: usize> Simulation<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    pub fn new_compuatation(
        lattice: &'a Lattice<D, SIZE>,
        algorithm: Algorithm,
        observable: Observable,
        temp: f64,
        burnin: usize,
        iterations: usize,
    ) -> Computation<'a, D, SIZE> {
        let computation: ComputationType = algorithm.into_computation_type();
        Computation {
            lattice,
            x: None,
            y: None,
            t: None,
            computation,
            observable,
            temp,
            range: None,
            burnin: Some(burnin),
            iterations: Some(iterations),
        }
    }
}

impl<'a, const D: usize, const SIZE: usize> Computatable<'a, D, SIZE> for Simulation<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    /// Builds a simulation from computation data. May fail.
    fn from_computation(computation: Computation<'a, D, SIZE>) -> Result<Self, Box<dyn Error>> {
        Ok(Simulation {
            lattice: computation.lattice,
            algorithm: match computation.computation {
                ComputationType::ClusterSimulation => Algorithm::Cluster,
                ComputationType::MetropolisSimulation => Algorithm::Metropolis,
                ComputationType::Test => Err("ComputationType::Test not allowed for Simulation")?,
            },
            x: computation.x,
            y: computation.y,
            t: computation.t,
            observable: computation.observable,
            temp: computation.temp,
            burnin: computation
                .burnin
                .ok_or("Burnin undefined for Simulation")?,
            iterations: computation
                .iterations
                .ok_or("Iterations undefined for Simulation")?,
        })
    }

    /// This is the main function in which the simulation is run. It
    /// initializes a field, whose configuration is theniteratively altered by
    /// sweep functions depending on the chosen simulation type. Finally a
    /// result is calculated and returned.
    fn run(self) -> Result<ComputationResults, Box<dyn Error>> {
        let time = Instant::now();
        let field: Field<i8, D, SIZE> = Field::random(self.lattice);
        let mut field: Field<i32, D, SIZE> = Field::from_field(field);

        // Burnin: compute an amount of sweeps to achieve equilibrium
        for _step in 0..(self.burnin) {
            //println!("Sweep {_step}");
            self.algorithm.field_sweep(&mut field, self.temp);
            field.normalize_random();
        }

        // After having reached equilibrium, for each consecutive field configuration
        // that is generated by the markov chain, we calculate the observable and average it.
        let mut observable_array: Vec<f64> = Vec::with_capacity(self.iterations); // Simulation observable arrray
        for _step in 0..(self.iterations) {
            //println!("Sweep {_step}");
            self.algorithm.field_sweep(&mut field, self.temp);
            observable_array.push(self.observable.observe(&field, self.temp));
            field.normalize_random();
        }

        let obs: ObsChain = ObsChain::new(observable_array.clone());
        let fehler: Option<f64> = obs.get_error_tanh(self.temp, 2).ok();

        // Mean value of the lattice actions generated by the simulation
        let observable: f64 = observable_array.iter().sum::<f64>() / observable_array.len() as f64;

        // Construct the return data type
        let mut result: ComputationResults = self.clone().into_result(observable)?;
        result.error = fehler;

        let duration: f32 = time.elapsed().as_secs_f32();
        println!(
            "{} {:?} Simulation took {duration} secs",
            self.temp, self.algorithm
        );

        Ok(result)
    }

    fn into_result(self, output: f64) -> Result<ComputationResults, Box<dyn Error>> {
        Ok(ComputationResults {
            dimensions: D,
            size: SIZE,
            x: self.x,
            y: self.y,
            t: self.t,
            computation: match self.algorithm {
                Algorithm::Cluster => ComputationType::ClusterSimulation,
                Algorithm::Metropolis => ComputationType::MetropolisSimulation,
            },
            observable: self.observable.to_string(),
            temp: self.temp,
            range: None,
            burnin: Some(self.burnin),
            iterations: Some(self.iterations),
            output,
            error: None,
        })
    }
}

pub struct Simulation3d<const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
where
    [(); MAX_X * MAX_Y * MAX_T]:;

impl<'a, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    Simulation3d<MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    pub fn new_compuatation(
        lattice: &'a Lattice3d<MAX_X, MAX_Y, MAX_T>,
        algorithm: Algorithm,
        observable: Observable,
        temp: f64,
        burnin: usize,
        iterations: usize,
    ) -> Computation<'a, 3, { MAX_X * MAX_Y * MAX_T }> {
        let computation: ComputationType = algorithm.into_computation_type();
        Computation {
            lattice,
            x: Some(MAX_X),
            y: Some(MAX_Y),
            t: Some(MAX_T),
            computation,
            observable,
            temp,
            range: None,
            burnin: Some(burnin),
            iterations: Some(iterations),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Test<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    lattice: &'a Lattice<D, SIZE>,
    x: Option<usize>,
    y: Option<usize>,
    t: Option<usize>,
    observable: Observable,
    temp: f64,
    range: usize,
}

impl<'a, const D: usize, const SIZE: usize> Test<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    pub fn new_computation(
        lattice: &'a Lattice<D, SIZE>,
        observable: Observable,
        temp: f64,
        range: usize,
    ) -> Computation<'a, D, SIZE> {
        Computation {
            lattice,
            x: None,
            y: None,
            t: None,
            computation: ComputationType::Test,
            observable,
            temp,
            range: Some(range),
            burnin: None,
            iterations: None,
        }
    }
}

impl<'a, const D: usize, const SIZE: usize> Computatable<'a, D, SIZE> for Test<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn from_computation(comp: Computation<'a, D, SIZE>) -> Result<Self, Box<dyn Error>> {
        Ok(Test {
            lattice: comp.lattice,
            x: comp.x,
            y: comp.y,
            t: comp.t,
            observable: comp.observable,
            temp: comp.temp,
            range: comp.range.ok_or("Range undefined for Test")?,
        })
    }

    /// This function calculates the observables over all configurations, which
    /// only works for small lattices.
    fn run(self) -> Result<ComputationResults, Box<dyn Error>> {
        let time = Instant::now();

        // 16 ^ 7 = 268’435’456
        let Some(permutations): Option<u64> = (self.range as u64)
            .checked_pow(SIZE as u32 - 1) else { panic!("Permutations overflow.")};

        let field: Field<i8, D, SIZE> = Field::new(self.lattice);
        let mut field: Field<i32, D, SIZE> = Field::from_field(field);

        assert_eq!(SIZE, field.values.len());

        field.values[SIZE - 1] = (self.range / 2) as i32;

        // The partition function is the sum over all Bolzmann weights
        let mut partfn: f64 = 0.0;

        // The observable is the sum over all weights (Bolzmann times observable),
        // devided by the partition function.
        let mut test: f64 = 0.0;

        let boundary: i32 = self.range as i32 - 1;

        for _ in 0..permutations {
            'updateconfig: for index in 0..SIZE {
                match field.values[index] {
                    x if x < boundary => {
                        field.values[index] = field.values[index] + 1;
                        break 'updateconfig;
                    }
                    x if x == boundary => {
                        field.values[index] = 0;
                    }
                    _ => {
                        panic!("config entry out of bounds.");
                    }
                }
            }
            let bolz: f64 = (-field.action_observable(self.temp)).exp();
            test = test + (self.observable.observe(&field, self.temp) * bolz);
            partfn = partfn + bolz;
        }

        let result: ComputationResults = self.clone().into_result(test / partfn)?;
        let duration = time.elapsed().as_secs_f32();
        println!("{} Test took {duration} secs", self.temp);

        Ok(result)
    }

    fn into_result(self, output: f64) -> Result<ComputationResults, Box<dyn Error>> {
        Ok(ComputationResults {
            dimensions: D,
            size: SIZE,
            x: self.x,
            y: self.y,
            t: self.t,
            computation: ComputationType::Test,
            observable: self.observable.to_string(),
            temp: self.temp,
            range: Some(self.range),
            burnin: None,
            iterations: None,
            output: output,
            error: None,
        })
    }
}

pub struct Test3d<const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
where
    [(); MAX_X * MAX_Y * MAX_T]:;

impl<'a, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Test3d<MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    pub fn new_computation(
        lattice: &'a Lattice3d<MAX_X, MAX_Y, MAX_T>,
        observable: Observable,
        temp: f64,
        range: usize,
    ) -> Computation<'a, 3, { MAX_X * MAX_Y * MAX_T }> {
        Computation {
            lattice,
            x: Some(MAX_X),
            y: Some(MAX_Y),
            t: Some(MAX_T),
            computation: ComputationType::Test,
            observable,
            temp,
            range: Some(range),
            burnin: None,
            iterations: None,
        }
    }
}
