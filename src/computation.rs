use std::{error::Error, fmt::Display, time::Instant};

use serde::{Deserialize, Serialize};

use crate::{
    algorithm::{Algorithm, AlgorithmType},
    error::ObsChain,
    field::Field,
    heightfield::{Action, HeightField},
    lattice::Lattice,
    observable::{Observable, ObservableType},
};

#[derive(Debug, Serialize, Deserialize, Clone)]
pub enum ComputationType {
    Simulation {
        algorithm: AlgorithmType,
        burnin: usize,
        iterations: usize,
    },
    Test {
        range: usize,
        permutations: u64,
    },
}

impl Display for ComputationType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ComputationType::Simulation {
                algorithm,
                iterations,
                ..
            } => {
                write!(f, "{} Simulation ({})", algorithm, iterations)
            }
            ComputationType::Test { range, .. } => write!(f, "Test (0 - {})", range),
        }
    }
}

#[derive(Debug, Clone)]
struct Simulation<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    lattice: &'a Lattice<D, SIZE>,
    temp: f64,
    algorithm: AlgorithmType,
    burnin: usize,
    iterations: usize,
    observable: ObservableType,
}

#[derive(Debug, Clone)]
struct Test<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    lattice: &'a Lattice<D, SIZE>,
    temp: f64,
    range: usize,
    permutations: u64,
    observable: ObservableType,
}

#[derive(Debug, Clone)]
pub struct Computation<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    lattice: &'a Lattice<D, SIZE>,
    temp: f64,
    comptype: ComputationType,
    observable: ObservableType,
}

impl<'a, const D: usize, const SIZE: usize> Computation<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    pub fn new_simulation(
        lattice: &'a Lattice<D, SIZE>,
        temp: f64,
        algorithm: AlgorithmType,
        burnin: usize,
        iterations: usize,
        observable: ObservableType,
    ) -> Computation<'a, D, SIZE> {
        Computation {
            lattice,
            temp,
            comptype: ComputationType::Simulation {
                algorithm,
                burnin,
                iterations,
            },
            observable,
        }
    }

    pub fn new_test(
        lattice: &'a Lattice<D, SIZE>,
        temp: f64,
        range: usize,
        observable: ObservableType,
    ) -> Computation<'a, D, SIZE> {
        // 16 ^ 7 = 268’435’456
        let Some(permutations): Option<u64> = (range as u64)
            .checked_pow(SIZE as u32 - 1) else { panic!("Permutations overflow.")};
        Computation {
            lattice,
            temp,
            comptype: ComputationType::Test {
                range,
                permutations,
            },
            observable,
        }
    }

    fn into_simulation(self) -> Result<Simulation<'a, D, SIZE>, Box<dyn Error>> {
        if let ComputationType::Simulation {
            algorithm,
            burnin,
            iterations,
        } = self.comptype
        {
            Ok(Simulation {
                lattice: self.lattice,
                temp: self.temp,
                algorithm,
                burnin,
                iterations,
                observable: self.observable,
            })
        } else {
            Err("Unable to convert Computation into Simulation.".into())
        }
    }

    fn into_test(self) -> Result<Test<'a, D, SIZE>, Box<dyn Error>> {
        if let ComputationType::Test {
            range,
            permutations,
        } = self.comptype
        {
            Ok(Test {
                lattice: self.lattice,
                temp: self.temp,
                range,
                permutations,
                observable: self.observable,
            })
        } else {
            Err("Unable to convert Computation into Test.".into())
        }
    }
}

pub trait Compute {
    fn run(self) -> Result<ComputationResult, Box<dyn Error>>;
    fn into_comp_result(self, result: f64) -> ComputationResult;
}

impl<'a, const D: usize, const SIZE: usize> Compute for Computation<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn run(self) -> Result<ComputationResult, Box<dyn Error>> {
        match self.comptype {
            ComputationType::Simulation { .. } => self.into_simulation()?.run(),
            ComputationType::Test { .. } => self.into_test()?.run(),
        }
    }

    fn into_comp_result(self, result: f64) -> ComputationResult {
        let [x, y, t]: [Option<usize>; 3] = if D == 3 {
            [
                Some(self.lattice.size[0]),
                Some(self.lattice.size[1]),
                Some(self.lattice.size[2]),
            ]
        } else {
            [None, None, None]
        };
        ComputationResult {
            d: D,
            size: SIZE,
            x,
            y,
            t,
            temp: self.temp,
            comptype: self.comptype,
            observable: self.observable,
            result,
            error: None,
        }
    }
}

impl<'a, const D: usize, const SIZE: usize> Compute for Simulation<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn run(self) -> Result<ComputationResult, Box<dyn Error>> {
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
        let fehler: f64 = obs.get_error_tanh(self.temp, 2)?;

        // Mean value of the lattice actions generated by the simulation
        let result: f64 = observable_array.iter().sum::<f64>() / observable_array.len() as f64;

        // Construct the return data type
        let mut result: ComputationResult = self.clone().into_comp_result(result);
        result.add_error(fehler);

        let duration: f32 = time.elapsed().as_secs_f32();
        println!(
            "{} {:?} Simulation took {duration} secs",
            self.temp, self.algorithm
        );

        Ok(result)
    }

    fn into_comp_result(self, result: f64) -> ComputationResult {
        let [x, y, t]: [Option<usize>; 3] = if D == 3 {
            [
                Some(self.lattice.size[0]),
                Some(self.lattice.size[1]),
                Some(self.lattice.size[2]),
            ]
        } else {
            [None, None, None]
        };
        ComputationResult {
            d: D,
            size: SIZE,
            x,
            y,
            t,
            temp: self.temp,
            comptype: ComputationType::Simulation {
                algorithm: self.algorithm,
                burnin: self.burnin,
                iterations: self.iterations,
            },
            observable: self.observable,
            result,
            error: None,
        }
    }
}

impl<'a, const D: usize, const SIZE: usize> Compute for Test<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    /// This function calculates the observables over all configurations, which
    /// only works for small lattices.
    fn run(self) -> Result<ComputationResult, Box<dyn Error>> {
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

        let result: ComputationResult = self.clone().into_comp_result(test / partfn);
        let duration = time.elapsed().as_secs_f32();
        println!("{} Test took {duration} secs", self.temp);

        Ok(result)
    }

    fn into_comp_result(self, result: f64) -> ComputationResult {
        let [x, y, t]: [Option<usize>; 3] = if D == 3 {
            [
                Some(self.lattice.size[0]),
                Some(self.lattice.size[1]),
                Some(self.lattice.size[2]),
            ]
        } else {
            [None, None, None]
        };
        ComputationResult {
            d: D,
            size: SIZE,
            x,
            y,
            t,
            temp: self.temp,
            comptype: ComputationType::Test {
                range: self.range,
                permutations: self.permutations,
            },
            observable: self.observable,
            result,
            error: None,
        }
    }
}

// ----------------------------------------------------------------------------

#[derive(Debug, Serialize, Deserialize)]
/// The result of a [`Computation`] that has finished.
pub struct ComputationResult {
    d: usize,
    size: usize,
    x: Option<usize>,
    y: Option<usize>,
    t: Option<usize>,
    temp: f64,
    comptype: ComputationType,
    observable: ObservableType,
    result: f64,
    error: Option<f64>,
}

impl ComputationResult {
    /// Adds a error to the [`ComputationResult`] data.
    fn add_error(&mut self, error: f64) {
        self.error = Some(error);
    }

    pub fn into_export(self) -> ComputationExport {
        ComputationExport {
            d: self.d,
            size: self.size,
            x: self.x,
            y: self.y,
            t: self.t,
            temp: self.temp,
            comptype: self.comptype.to_string(),
            observable: self.observable.to_string(),
            result: self.result,
            error: self.error,
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
/// We export the [`ComputationResult`] in a csv file with this formating
pub struct ComputationExport {
    d: usize,
    size: usize,
    x: Option<usize>,
    y: Option<usize>,
    t: Option<usize>,
    temp: f64,
    comptype: String,
    observable: String,
    result: f64,
    error: Option<f64>,
}
