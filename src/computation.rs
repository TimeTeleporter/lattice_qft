use std::{error::Error, fmt::Display, time::Instant};

use rand::rngs::ThreadRng;
use serde::{Deserialize, Serialize};

use crate::{
    algorithm::{Algorithm, AlgorithmType, WilsonAlgorithm},
    export::CsvData,
    fields::{Field, HeightField, WilsonField},
    lattice::Lattice,
    outputdata::{Observe, OutputData, OutputDataType, Plotting, UpdateOutputData},
};

const BURNIN_ANNOUNCEMENT: bool = false;

// - Computation --------------------------------------------------------------

#[derive(Debug, Clone)]
pub enum Computation<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    Simulation(Simulation<'a, D, SIZE>),
    Test(Test<'a, D, SIZE>),
    WilsonSim(WilsonSim<'a, SIZE>),
    WilsonTest(WilsonTest<'a, SIZE>),
}

// - Computation - Formatting -------------------------------------------------

impl<'a, const D: usize, const SIZE: usize> Display for Computation<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Computation::Simulation(Simulation { algorithm, .. }) => {
                write!(f, "{} Simulation", algorithm)
            }
            Computation::Test(Test { .. }) => write!(f, "Test"),
            Computation::WilsonSim(WilsonSim {
                algorithm,
                width,
                height,
                ..
            }) => write!(f, "{}x{} Wilson {} Simulation", width, height, algorithm),
            Computation::WilsonTest(WilsonTest { width, height, .. }) => {
                write!(f, "{}x{} Wilson Test", width, height)
            }
        }
    }
}

// - Computation - Constructors -----------------------------------------------

impl<'a, const D: usize, const SIZE: usize> Computation<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    pub fn new_simulation(
        lattice: &'a Lattice<D, SIZE>,
        temp: f64,
        algorithm: AlgorithmType,
        burnin: u64,
        iterations: u64,
        output: Vec<OutputData<'a, D, SIZE>>,
    ) -> Computation<'a, D, SIZE> {
        Computation::Simulation(Simulation {
            lattice,
            temp,
            algorithm,
            burnin,
            iterations,
            output,
            duration: None,
        })
    }
}

impl<'a, const SIZE: usize> Computation<'a, 3, SIZE> {
    pub fn new_wilson_sim(
        lattice: &'a Lattice<3, SIZE>,
        temp: f64,
        algorithm: AlgorithmType,
        burnin: u64,
        iterations: u64,
        width: usize,
        height: usize,
        output: Vec<OutputData<'a, 3, SIZE>>,
    ) -> Computation<'a, 3, SIZE> {
        Computation::WilsonSim(WilsonSim {
            lattice,
            temp,
            algorithm,
            burnin,
            iterations,
            width,
            height,
            output,
            duration: None,
        })
    }
}

/// Shared behaviour of all [Computation] variants.
pub trait Compute<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    fn run(self) -> Result<Computation<'a, D, SIZE>, Box<dyn Error>>;
    fn into_computation(self) -> Computation<'a, D, SIZE>;
}

impl<'a, const SIZE: usize> Compute<'a, 3, SIZE> for Computation<'a, 3, SIZE> {
    fn run(self) -> Result<Computation<'a, 3, SIZE>, Box<dyn Error>> {
        match self {
            Computation::Simulation(sim) => sim.run(),
            Computation::Test(test) => test.run(),
            Computation::WilsonSim(sim) => sim.run(),
            Computation::WilsonTest(test) => test.run(),
        }
    }

    fn into_computation(self) -> Computation<'a, 3, SIZE> {
        self
    }
}

pub fn parse_simulation_results<const SIZE: usize>(data: Vec<Computation<3, SIZE>>) {
    for (_, computation) in data.into_iter().enumerate() {
        // Fetching the index to append to
        let index: u64 = match ComputationSummary::fetch_csv_data(crate::RESULTS_PATH, true)
            .and_then(|summary| {
                summary
                    .last()
                    .map(|last| last.index + 1)
                    .ok_or("No last element, starting anew.".into())
            }) {
            Ok(index) => index,
            Err(err) => {
                eprint!("{}", err);
                0
            }
        };

        // Building the computation summary and handling outputs.
        let (mut summary, outputs) = ComputationSummary::from_computation(computation, index);

        for output in outputs.into_iter() {
            match output.data {
                OutputDataType::ActionObservable(obs) => {
                    let path: &str = &(crate::ACTION_DATA_PATH_INCOMPLETE.to_owned()
                        + &index.to_string()
                        + &".csv");
                    if let Err(err) = obs.results().read_write_csv(path, false) {
                        eprint!("Writing action data: {}", err);
                    }
                    summary = summary.set_action_data();
                }
                OutputDataType::DifferencePlot(obs) => {
                    for (direction, plot) in obs.plot().into_iter().enumerate() {
                        let path: &str = &(crate::DIFFERENCE_PLOT_PATH_INCOMPLETE.to_owned()
                            + &index.to_string()
                            + &"_"
                            + &direction.to_string()
                            + &".csv");
                        if let Err(err) = plot.read_write_csv(path, true) {
                            eprint!("Writing difference plot data: {}", err);
                        };
                    }
                    summary = summary.set_bonds_data();
                }
                OutputDataType::CorrelationData(obs) => {
                    let path: &str = &(crate::CORR_FN_PATH_INCOMPLETE.to_owned()
                        + &"correlation_"
                        + &index.to_string()
                        + &".csv");
                    if let Err(err) = obs.plot().overwrite_csv(path) {
                        eprint!("Writing correlation data: {}", err);
                    };
                    summary = summary.set_correlation_data();
                }
                OutputDataType::EnergyObservable(obs) => {
                    let path: &str = &(crate::ENERGY_DATA_PATH_INCOMPLETE.to_owned()
                        + &index.to_string()
                        + &".csv");
                    if let Err(err) = obs.results().read_write_csv(path, false) {
                        eprint!("Writing energy data: {}", err);
                    }
                    summary = summary.set_energy_data();
                }
            }
        }

        if let Err(err) = summary.read_write_csv(crate::RESULTS_PATH, true) {
            eprint!("{}", err);
        };
    }
}

// - Simulation ---------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct Simulation<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    lattice: &'a Lattice<D, SIZE>,
    temp: f64,
    algorithm: AlgorithmType,
    burnin: u64,
    iterations: u64,
    output: Vec<OutputData<'a, D, SIZE>>,
    duration: Option<u64>,
}

impl<'a, const D: usize, const SIZE: usize> Compute<'a, D, SIZE> for Simulation<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn run(mut self) -> Result<Computation<'a, D, SIZE>, Box<dyn Error>> {
        let mut rng = ThreadRng::default();

        let time = Instant::now();
        let field: Field<i8, D, SIZE> = Field::random(self.lattice);
        let mut field: Field<i32, D, SIZE> = Field::from_field(field);

        // Burnin: compute an amount of sweeps to achieve equilibrium
        for step in 0..(self.burnin) {
            self.algorithm.field_sweep(&mut field, self.temp, &mut rng);
            if BURNIN_ANNOUNCEMENT {
                if step == 0 {
                    println!("{:.2} {: <12} Started", self.temp, "Burnin:");
                } else if step == self.burnin / 4 {
                    println!("{:.2} {: <12} 25%", self.temp, "Burnin:");
                } else if step == self.burnin / 2 {
                    println!("{:.2} {: <12} 50%", self.temp, "Burnin:");
                } else if step == self.burnin / 4 * 3 {
                    println!("{:.2} {: <12} 75%", self.temp, "Burnin:");
                }
            }
            if step % 10 == 0 {
                field.normalize_random();
            }
        }

        // After having reached equilibrium, for each consecutive field configuration
        // that is generated by the markov chain, we calculate the observable and average it.
        for step in 0..(self.iterations) {
            if step == 0 {
                println!("{:.2} {: <12} Started", self.temp, "Simulation:");
            } else if step == self.iterations / 4 {
                println!("{:.2} {: <12} 25%", self.temp, "Simulation:");
            } else if step == self.iterations / 2 {
                println!("{:.2} {: <12} 50%", self.temp, "Simulation:");
            } else if step == self.iterations / 4 * 3 {
                println!("{:.2} {: <12} 75%", self.temp, "Simulation:");
            }
            //println!("Sweep {_step}");
            self.algorithm.field_sweep(&mut field, self.temp, &mut rng);

            self.output
                .iter_mut()
                .filter(|data| data.is_step_for_next_update(step))
                .for_each(|data| data.update(&field, &mut rng));

            if step % 10 == 0 {
                field.normalize_random();
            }
        }

        let duration: u64 = time.elapsed().as_secs();
        println!(
            "{} {:?} Simulation took {duration} secs",
            self.temp, self.algorithm
        );

        self.duration = Some(duration);

        Ok(self.into_computation())
    }

    fn into_computation(self) -> Computation<'a, D, SIZE> {
        Computation::Simulation(self)
    }
}

// - WilsonSim ----------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct WilsonSim<'a, const SIZE: usize> {
    lattice: &'a Lattice<3, SIZE>,
    temp: f64,
    algorithm: AlgorithmType,
    burnin: u64,
    iterations: u64,
    width: usize,
    height: usize,
    output: Vec<OutputData<'a, 3, SIZE>>,
    duration: Option<u64>,
}

impl<'a, const SIZE: usize> Compute<'a, 3, SIZE> for WilsonSim<'a, SIZE> {
    fn run(mut self) -> Result<Computation<'a, 3, SIZE>, Box<dyn Error>> {
        let mut rng = ThreadRng::default();
        let time = Instant::now();
        let field: Field<i8, 3, SIZE> = Field::random(self.lattice);
        let field: Field<i32, 3, SIZE> = Field::from_field(field);
        let mut wilson: WilsonField<i32, SIZE> =
            WilsonField::from_field(field, self.width, self.height);

        // Burnin: compute an amount of sweeps to achieve equilibrium
        for step in 0..(self.burnin) {
            //println!("Sweep {_step}");
            self.algorithm.wilson_sweep(&mut wilson, self.temp);
            if step % 10 == 0 {
                wilson.normalize_random();
            }
        }

        // After having reached equilibrium, for each consecutive field configuration
        // that is generated by the markov chain, we calculate the observable and average it.
        for step in 0..(self.iterations) {
            if step == 0 {
                println!("{:.2} {: <12} Started", self.temp, "Wilson:");
            } else if step == self.iterations / 4 {
                println!("{:.2} {: <12} 25%", self.temp, "Wilson:");
            } else if step == self.iterations / 2 {
                println!("{:.2} {: <12} 50%", self.temp, "Wilson:");
            } else if step == self.iterations / 4 * 3 {
                println!("{:.2} {: <12} 75%", self.temp, "Wilson:");
            }
            self.algorithm.wilson_sweep(&mut wilson, self.temp);
            for data in self.output.iter_mut() {
                data.update(&wilson.field, &mut rng)
            }
            if step % 10 == 0 {
                wilson.normalize_random();
            }
        }

        let duration: u64 = time.elapsed().as_secs();
        println!(
            "{} {:?} Simulation took {duration} secs",
            self.temp, self.algorithm
        );
        self.duration = Some(duration);

        Ok(self.into_computation())
    }

    fn into_computation(self) -> Computation<'a, 3, SIZE> {
        Computation::WilsonSim(self)
    }
}

// - Test ---------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct Test<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    lattice: &'a Lattice<D, SIZE>,
    temp: f64,
    range: i32,
    permutations: u64,
    output: Vec<OutputData<'a, D, SIZE>>,
    duration: Option<u64>,
}

impl<'a, const D: usize, const SIZE: usize> Compute<'a, D, SIZE> for Test<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    /// This function calculates the observables over all configurations, which
    /// only works for small lattices.
    fn run(mut self) -> Result<Computation<'a, D, SIZE>, Box<dyn Error>> {
        let mut rng = ThreadRng::default();
        let time = Instant::now();

        let field: Field<i8, D, SIZE> = Field::new(self.lattice);
        let mut field: Field<i32, D, SIZE> = Field::from_field(field);

        assert_eq!(SIZE, field.values.len());

        field.values[SIZE - 1] = (self.range / 2) as i32;

        let boundary: i32 = self.range as i32 - 1;

        println!("Test: {} permutations", self.permutations);

        for step in 0..self.permutations {
            if step == 0 {
                println!("{} {: <12} Started", self.temp, "Test:");
            } else if step == self.permutations / 4 {
                println!("{} {: <12} 25%", self.temp, "Test:");
            } else if step == self.permutations / 2 {
                println!("{} {: <12} 50%", self.temp, "Test:");
            } else if step == self.permutations / 4 * 3 {
                println!("{} {: <12} 75%", self.temp, "Test:");
            }
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
            for data in self.output.iter_mut() {
                data.update(&field, &mut rng);
            }
        }

        let duration: u64 = time.elapsed().as_secs();
        println!("{} Test took {duration} secs", self.temp);

        self.duration = Some(duration);

        Ok(self.into_computation())
    }

    fn into_computation(self) -> Computation<'a, D, SIZE> {
        Computation::Test(self)
    }
}

// - WilsonTest ---------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct WilsonTest<'a, const SIZE: usize> {
    lattice: &'a Lattice<3, SIZE>,
    temp: f64,
    range: usize,
    permutations: u64,
    width: usize,
    height: usize,
    output: Vec<OutputData<'a, 3, SIZE>>,
    duration: Option<u64>,
}

impl<'a, const SIZE: usize> Compute<'a, 3, SIZE> for WilsonTest<'a, SIZE> {
    /// This function calculates the observables over all configurations, which
    /// only works for small lattices.
    fn run(mut self) -> Result<Computation<'a, 3, SIZE>, Box<dyn Error>> {
        let mut rng = ThreadRng::default();
        let time = Instant::now();

        let field: Field<i8, 3, SIZE> = Field::new(self.lattice);
        let field: Field<i32, 3, SIZE> = Field::from_field(field);
        let mut wilson: WilsonField<'a, i32, SIZE> =
            WilsonField::from_field(field, self.width, self.height);

        assert_eq!(SIZE, wilson.field.values.len());

        wilson.field.values[SIZE - 1] = (self.range / 2) as i32;

        let boundary: i32 = self.range as i32 - 1;

        println!("WilsonTest: {} permutations", self.permutations);

        for step in 0..self.permutations {
            if step == 0 {
                println!("{} {: <12} Started", self.temp, "WilsonTest:");
            } else if step == self.permutations / 4 {
                println!("{} {: <12} 25%", self.temp, "WilsonTest:");
            } else if step == self.permutations / 2 {
                println!("{} {: <12} 50%", self.temp, "WilsonTest:");
            } else if step == self.permutations / 4 * 3 {
                println!("{} {: <12} 75%", self.temp, "WilsonTest:");
            }
            'updateconfig: for index in 0..SIZE {
                match wilson.field.values[index] {
                    x if x < boundary => {
                        wilson.field.values[index] = wilson.field.values[index] + 1;
                        break 'updateconfig;
                    }
                    x if x == boundary => {
                        wilson.field.values[index] = 0;
                    }
                    _ => {
                        panic!("config entry out of bounds.");
                    }
                }
            }
            for data in self.output.iter_mut() {
                data.update(&wilson, &mut rng);
            }
        }

        let duration: u64 = time.elapsed().as_secs();
        println!("{} Test took {duration} secs", self.temp);

        self.duration = Some(duration);

        Ok(self.into_computation())
    }

    fn into_computation(self) -> Computation<'a, 3, SIZE> {
        Computation::WilsonTest(self)
    }
}

// - Computation - Summary ----------------------------------------------------

#[derive(Debug, Default, Serialize, Deserialize, Clone)]
/// We export the [`ComputationResult`] in a csv file with this formating
pub struct ComputationSummary {
    pub index: u64,
    pub d: usize,
    pub size: usize,
    pub x: usize,
    pub y: usize,
    pub t: usize,
    pub temp: f64,
    pub comptype: String,
    pub burnin: Option<u64>,
    pub iterations: u64,
    pub duration: Option<u64>,
    pub algorithm_statistic: Option<u64>,
    pub action_data: bool,
    pub energy_data: bool,
    pub difference_data: bool,
    pub correlation_data: bool,
}

impl ComputationSummary {
    pub fn from_computation<const SIZE: usize>(
        computation: Computation<3, SIZE>,
        index: u64,
    ) -> (Self, Vec<OutputData<3, SIZE>>) {
        let comptype: String = computation.to_string();
        match computation {
            Computation::Simulation(comp) => {
                let (d, size) = (3, SIZE);
                let [x, y, t] = comp.lattice.size;
                let temp: f64 = comp.temp;
                let (burnin, iterations) = (Some(comp.burnin), comp.iterations);
                let new: ComputationSummary = ComputationSummary {
                    index,
                    d,
                    size,
                    x,
                    y,
                    t,
                    temp,
                    comptype,
                    burnin,
                    iterations,
                    duration: comp.duration,
                    ..Default::default()
                };
                let outputs: Vec<OutputData<3, SIZE>> = comp.output;
                (new, outputs)
            }
            Computation::Test(comp) => {
                let (d, size) = (3, SIZE);
                let [x, y, t] = comp.lattice.size;
                let temp: f64 = comp.temp;
                let iterations = comp.permutations;
                let new: ComputationSummary = ComputationSummary {
                    index,
                    d,
                    size,
                    x,
                    y,
                    t,
                    temp,
                    comptype,
                    iterations,
                    ..Default::default()
                };
                let outputs: Vec<OutputData<3, SIZE>> = comp.output;
                (new, outputs)
            }
            Computation::WilsonSim(comp) => {
                let (d, size) = (3, SIZE);
                let [x, y, t] = comp.lattice.size;
                let temp: f64 = comp.temp;
                let (burnin, iterations) = (Some(comp.burnin), comp.iterations);
                let new: ComputationSummary = ComputationSummary {
                    index,
                    d,
                    size,
                    x,
                    y,
                    t,
                    temp,
                    comptype,
                    burnin,
                    iterations,
                    ..Default::default()
                };
                let outputs: Vec<OutputData<3, SIZE>> = comp.output;
                (new, outputs)
            }
            Computation::WilsonTest(comp) => {
                let (d, size) = (3, SIZE);
                let [x, y, t] = comp.lattice.size;
                let temp: f64 = comp.temp;
                let iterations = comp.permutations;
                let new: ComputationSummary = ComputationSummary {
                    index,
                    d,
                    size,
                    x,
                    y,
                    t,
                    temp,
                    comptype,
                    iterations,
                    ..Default::default()
                };
                let outputs: Vec<OutputData<3, SIZE>> = comp.output;
                (new, outputs)
            }
        }
    }

    pub fn set_action_data(mut self) -> Self {
        self.action_data = true;
        self
    }

    pub fn set_energy_data(mut self) -> Self {
        self.energy_data = true;
        self
    }

    pub fn set_bonds_data(mut self) -> Self {
        self.difference_data = true;
        self
    }

    pub fn set_correlation_data(mut self) -> Self {
        self.correlation_data = true;
        self
    }
}
