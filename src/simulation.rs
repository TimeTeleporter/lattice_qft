use std::{error::Error, ops::Deref, time::Instant};

use serde::{Deserialize, Serialize};

use crate::{
    action::Action,
    cluster::Cluster,
    error::ObsChain,
    export::{clean_csv, CsvData},
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
pub enum Observable {
    Action,
    /// This defines the Wilson loop observable.
    ///
    /// # Parameters
    ///
    /// 1. width
    /// 2. height
    /// 3. temp
    ///
    Wilson(usize, usize, f64),
    SizeNormalized,
}

impl Observable {
    /// Calculates an observable from the field variables.
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
            Observable::Wilson(width, height, temp) => {
                field.wilson_loop_observable(*temp, *width, *height)
                    * (-1.0 * temp * *height as f64 * *width as f64).exp()
            }
            Observable::SizeNormalized => field.size_normalized_action_observable(temp),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum Algorithm {
    Cluster,
    Metropolis,
}

pub enum SweepReturn {
    ClustersAmount(usize),
    AcceptanceRate(usize),
}

impl SweepReturn {
    fn print_stats(mut ary: Vec<SweepReturn>) -> Result<(), Box<dyn Error>> {
        let entries: f64 = ary.len() as f64;
        let poppy: SweepReturn = ary.pop().ok_or("No sweep stats to analyse.")?;
        let mut sum = match poppy {
            SweepReturn::ClustersAmount(entry) => {
                print!("Average cluster size of ");
                entry
            }
            SweepReturn::AcceptanceRate(entry) => {
                print!("Average acceptance rate of ");
                entry
            }
        } as f64;
        for sweepy in ary {
            sum += match sweepy {
                SweepReturn::ClustersAmount(addend) => addend,
                SweepReturn::AcceptanceRate(addend) => addend,
            } as f64;
        }
        let mean: f64 = sum / entries;
        println!("{mean}");

        Ok(())
    }
}

impl Algorithm {
    /// Performs a single sweep depending on the simulation type. Returns some statistics from the sweep
    fn single_sweep<const D: usize, const SIZE: usize>(
        &self,
        field: &mut Field<i32, D, SIZE>,
        temp: f64,
    ) -> SweepReturn
    where
        [(); D * 2_usize]:,
    {
        match self {
            Algorithm::Cluster => SweepReturn::ClustersAmount(field.cluster_sweep(temp)),
            Algorithm::Metropolis => SweepReturn::AcceptanceRate(field.metropolis_sweep(temp)),
        }
    }
}

pub enum ComputationResult {
    SimResult(SimResult),
    SimResultObs((SimResult, ObsChain)),
    TestResult(TestResult),
}

impl ComputationResult {
    pub fn write_to_csv(
        self,
        sim_path: &str,
        bin_path: &str,
        test_path: &str,
    ) -> Result<(), Box<dyn Error>> {
        match self {
            ComputationResult::SimResultObs((res, obs)) => {
                let temp: f64 = res.temp;
                res.read_write_csv(sim_path)?;
                clean_csv(bin_path)?;
                for bin in obs.calculate_binnings(temp, 2) {
                    bin.read_write_csv(bin_path)?;
                }
                Ok(())
            }
            ComputationResult::TestResult(test) => test.read_write_csv(test_path),
            ComputationResult::SimResult(res) => res.read_write_csv(sim_path),
        }
    }
}

pub enum ComputationType<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    Simulation(Simulation<'a, D, SIZE>),
    Test(TestSim<'a, D, SIZE>),
}

impl<'a, const D: usize, const SIZE: usize> ComputationType<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    pub fn run(&self) -> ComputationResult {
        match self {
            ComputationType::Simulation(sim) => ComputationResult::SimResult({
                let (res, _) = sim.run();
                res
            }),
            ComputationType::Test(test) => ComputationResult::TestResult(test.run()),
        }
    }
}

pub enum ComputationType3d<'a, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    Simulation(Simulation3d<'a, MAX_X, MAX_Y, MAX_T>),
    Test(TestSim3d<'a, MAX_X, MAX_Y, MAX_T>),
}

impl<'a, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    ComputationType3d<'a, MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    pub fn run(&self) -> ComputationResult {
        match self {
            ComputationType3d::Simulation(sim) => ComputationResult::SimResult({
                let (res, _) = sim.run();
                res
            }),
            ComputationType3d::Test(test) => ComputationResult::TestResult(test.run()),
        }
    }
}

pub struct Simulation<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    name: String,
    algorithm: Algorithm,
    observable: Observable,
    lattice: &'a Lattice<D, SIZE>,
    pub temp: f64,
    burnin: usize,
    iterations: usize,
}

impl<'a, const D: usize, const SIZE: usize> Simulation<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    pub fn new(
        name: String,
        algorithm: Algorithm,
        observable: Observable,
        lattice: &'a Lattice<D, SIZE>,
        temp: f64,
        burnin: usize,
        iterations: usize,
    ) -> Self {
        Simulation {
            name,
            algorithm,
            observable,
            lattice,
            temp,
            burnin,
            iterations,
        }
    }

    /// This is the main function in which the simulation is run. It
    /// initializes a field, whose configuration is theniteratively altered by
    /// sweep functions depending on the chosen simulation type. Finally a
    /// result is calculated and returned.
    pub fn run(&self) -> (SimResult, ObsChain) {
        println!("Started {} {}", self.temp, self.name);
        let time = Instant::now();
        let field: Field<i8, D, SIZE> = Field::random(self.lattice);
        let mut field: Field<i32, D, SIZE> = Field::from_field(field);

        // Burnin: compute an amount of sweeps to achieve equilibrium
        for _step in 0..(self.burnin) {
            println!("Sweep {_step}");
            self.algorithm.single_sweep(&mut field, self.temp);
            field.normalize_random();
        }

        // After having reached equilibrium, for each consecutive field configuration
        // that is generated by the markov chain, we calculate the observable and average it.
        let mut observable_array: Vec<f64> = Vec::with_capacity(self.iterations); // Simulation observable arrray
        let mut sweepstats_ary: Vec<SweepReturn> = Vec::with_capacity(self.iterations);
        for _step in 0..(self.iterations) {
            println!("Sweep {_step}");
            sweepstats_ary.push(self.algorithm.single_sweep(&mut field, self.temp));
            observable_array.push(self.observable.observe(&field, self.temp));
            field.normalize_random();
        }

        if let Err(err) = SweepReturn::print_stats(sweepstats_ary) {
            eprintln!("{err}");
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
        algorithm: Algorithm,
        observable: Observable,
        lattice: &'a Lattice3d<MAX_X, MAX_Y, MAX_T>,
        temp: f64,
        burnin: usize,
        iterations: usize,
    ) -> Self {
        let sim: Simulation<3, { MAX_X * MAX_Y * MAX_T }> = Simulation::new(
            name,
            algorithm,
            observable,
            lattice.deref(),
            temp,
            burnin,
            iterations,
        );
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

/// Datatype to save and read simulation output.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestResult {
    name: String,
    pub temp: f64,
    range: usize,
    observable: f64,
    error: Option<f64>,
}

impl TestResult {
    pub fn new(name: String, temp: f64, range: usize, observable: f64, error: Option<f64>) -> Self {
        TestResult {
            name,
            temp,
            range,
            observable,
            error,
        }
    }

    pub fn set_error(&mut self, error_option: Option<f64>) {
        self.error = error_option;
    }
}

/// A struct to initialize a test of a lattice with lattice sites in a given
/// range.
pub struct TestSim<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    name: String,
    observable: Observable,
    lattice: &'a Lattice<D, SIZE>,
    pub temp: f64,
    range: usize,
}

impl<'a, const D: usize, const SIZE: usize> TestSim<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    pub fn new(
        name: String,
        observable: Observable,
        lattice: &'a Lattice<D, SIZE>,
        temp: f64,
        range: usize,
    ) -> Self {
        TestSim {
            name,
            observable,
            lattice,
            temp,
            range,
        }
    }

    /// This is the main function in which the simulation is run. It
    /// initializes a field, whose configuration is theniteratively altered by
    /// sweep functions depending on the chosen simulation type. Finally a
    /// result is calculated and returned.
    pub fn run(&self) -> TestResult {
        let permutations: usize = self.range.pow(SIZE as u32 - 1); // 16 ^ 7 = 268’435’456

        print!("Started to calculate the observable over the range");
        println!(" -{} to +{}.", self.range / 2, self.range / 2);
        println!("Those are {} field configurations.", permutations);

        let field: Field<i8, D, SIZE> = Field::new(self.lattice);
        let mut field: Field<i32, D, SIZE> = Field::from_field(field);

        field.values[7] = (self.range / 2) as i32;

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

        self.produce_simresult(test / partfn)
    }

    fn produce_simresult(&self, observable: f64) -> TestResult {
        TestResult {
            name: self.name.clone(),
            temp: self.temp,
            range: self.range,
            observable,
            error: None,
        }
    }
}

pub struct TestSim3d<'a, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>(
    TestSim<'a, 3, { MAX_X * MAX_Y * MAX_T }>,
)
where
    [(); MAX_X * MAX_Y * MAX_T]:;

impl<'a, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize>
    TestSim3d<'a, MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    pub fn new(
        name: String,
        observable: Observable,
        lattice: &'a Lattice3d<MAX_X, MAX_Y, MAX_T>,
        temp: f64,
        range: usize,
    ) -> Self {
        TestSim3d(TestSim {
            name,
            observable,
            lattice,
            temp,
            range,
        })
    }
}

impl<'a, const MAX_X: usize, const MAX_Y: usize, const MAX_T: usize> Deref
    for TestSim3d<'a, MAX_X, MAX_Y, MAX_T>
where
    [(); MAX_X * MAX_Y * MAX_T]:,
{
    type Target = TestSim<'a, 3, { MAX_X * MAX_Y * MAX_T }>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
