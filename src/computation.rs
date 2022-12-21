use std::{error::Error, fmt::Display, time::Instant};

use serde::{Deserialize, Serialize};

use crate::{
    algorithm::{Algorithm, AlgorithmType, WilsonAlgorithm},
    fields::{Action, Field, HeightField, WilsonField},
    lattice::Lattice,
    observable::{Observable, ObservableType, ObservableValue},
};

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
            Computation::Simulation(Simulation {
                algorithm,
                iterations,
                ..
            }) => write!(f, "{} Simulation ({})", algorithm, iterations),
            Computation::Test(Test { range, .. }) => write!(f, "Test (0 - {})", range),
            Computation::WilsonSim(WilsonSim {
                algorithm,
                iterations,
                width,
                height,
                ..
            }) => write!(
                f,
                "{}x{} Wilson {} Simulation ({})",
                width, height, algorithm, iterations
            ),
            Computation::WilsonTest(WilsonTest {
                range,
                width,
                height,
                ..
            }) => write!(f, "{}x{} Wilson Test (0 - {})", width, height, range),
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
        burnin: usize,
        iterations: usize,
        observable: ObservableType,
    ) -> Computation<'a, D, SIZE> {
        Computation::Simulation(Simulation {
            lattice,
            temp,
            algorithm,
            burnin,
            iterations,
            observable,
        })
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
        Computation::Test(Test {
            lattice,
            temp,
            range,
            permutations,
            observable,
        })
    }
}

impl<'a, const SIZE: usize> Computation<'a, 3, SIZE> {
    pub fn new_wilson_sim(
        lattice: &'a Lattice<3, SIZE>,
        temp: f64,
        algorithm: AlgorithmType,
        burnin: usize,
        iterations: usize,
        width: usize,
        height: usize,
        observable: ObservableType,
    ) -> Computation<'a, 3, SIZE> {
        Computation::WilsonSim(WilsonSim {
            lattice,
            temp,
            algorithm,
            burnin,
            iterations,
            width,
            height,
            observable,
        })
    }
    pub fn new_wilson_test(
        lattice: &'a Lattice<3, SIZE>,
        temp: f64,
        range: usize,
        width: usize,
        height: usize,
        observable: ObservableType,
    ) -> Computation<'a, 3, SIZE> {
        // 16 ^ 7 = 268’435’456
        let Some(permutations): Option<u64> = (range as u64)
            .checked_pow(SIZE as u32 - 1) else { panic!("Permutations overflow.")};
        Computation::WilsonTest(WilsonTest {
            lattice,
            temp,
            range,
            permutations,
            width,
            height,
            observable,
        })
    }
}

/// Shared behaviour of all [Computation] variants.
pub trait Compute<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    fn run(self) -> Result<ComputationResult<'a, D, SIZE>, Box<dyn Error>>;
    fn into_comp_result(self, result: f64) -> ComputationResult<'a, D, SIZE>;
}

impl<'a, const D: usize, const SIZE: usize> Compute<'a, D, SIZE> for Computation<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn run(self) -> Result<ComputationResult<'a, D, SIZE>, Box<dyn Error>> {
        match self {
            Computation::Simulation(sim) => sim.run(),
            Computation::Test(test) => test.run(),
            Computation::WilsonSim(sim) => sim.run(),
            Computation::WilsonTest(test) => test.run(),
        }
    }

    fn into_comp_result(self, result: f64) -> ComputationResult<'a, D, SIZE> {
        match self {
            Computation::Simulation(sim) => sim.into_comp_result(result),
            Computation::Test(test) => test.into_comp_result(result),
            Computation::WilsonSim(sim) => sim.into_comp_result(result),
            Computation::WilsonTest(test) => test.into_comp_result(result),
        }
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
    burnin: usize,
    iterations: usize,
    observable: ObservableType,
}

impl<'a, const D: usize, const SIZE: usize> Compute<'a, D, SIZE> for Simulation<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    fn run(self) -> Result<ComputationResult<'a, D, SIZE>, Box<dyn Error>> {
        let time = Instant::now();
        let field: Field<i8, D, SIZE> = Field::random(self.lattice);
        let mut field: Field<i32, D, SIZE> = Field::from_field(field);

        let mut obs: Observable<i32, D, SIZE> =
            Observable::new(self.lattice, self.observable.clone(), self.temp);

        // Burnin: compute an amount of sweeps to achieve equilibrium
        for _step in 0..(self.burnin) {
            //println!("Sweep {_step}");
            self.algorithm.field_sweep(&mut field, self.temp);
            field.normalize_random();
        }

        // After having reached equilibrium, for each consecutive field configuration
        // that is generated by the markov chain, we calculate the observable and average it.
        for _step in 0..(self.iterations) {
            //println!("Sweep {_step}");
            self.algorithm.field_sweep(&mut field, self.temp);
            obs.update(&field);
            field.normalize_random();
        }

        let result: f64 = obs.result();

        let duration: f32 = time.elapsed().as_secs_f32();
        println!(
            "{} {:?} Simulation took {duration} secs",
            self.temp, self.algorithm
        );

        Ok(self.into_comp_result(result))
    }

    fn into_comp_result(self, result: f64) -> ComputationResult<'a, D, SIZE> {
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
            comptype: Computation::Simulation(self.clone()),
            observable: self.observable,
            result,
            error: None,
        }
    }
}

// - WilsonSim ----------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct WilsonSim<'a, const SIZE: usize> {
    lattice: &'a Lattice<3, SIZE>,
    temp: f64,
    algorithm: AlgorithmType,
    burnin: usize,
    iterations: usize,
    width: usize,
    height: usize,
    observable: ObservableType,
}

impl<'a, const D: usize, const SIZE: usize> Compute<'a, D, SIZE> for WilsonSim<'a, SIZE>
where
    [(); D * 2_usize]:,
{
    fn run(self) -> Result<ComputationResult<'a, D, SIZE>, Box<dyn Error>> {
        let time = Instant::now();
        let field: Field<i8, 3, SIZE> = Field::random(self.lattice);
        let field: Field<i32, 3, SIZE> = Field::from_field(field);
        let mut wilson: WilsonField<i32, SIZE> =
            WilsonField::from_field(field, self.width, self.height);

        let mut obs: Observable<i32, 3, SIZE> =
            Observable::new(self.lattice, self.observable.clone(), self.temp);

        // Burnin: compute an amount of sweeps to achieve equilibrium
        for _step in 0..(self.burnin) {
            //println!("Sweep {_step}");
            self.algorithm.wilson_sweep(&mut wilson, self.temp);
            wilson.normalize_random();
        }

        // After having reached equilibrium, for each consecutive field configuration
        // that is generated by the markov chain, we calculate the observable and average it.
        for _step in 0..(self.iterations) {
            //println!("Sweep {_step}");
            self.algorithm.wilson_sweep(&mut wilson, self.temp);
            obs.update(&wilson.field);
            wilson.normalize_random();
        }

        let result: f64 = obs.result();

        let duration: f32 = time.elapsed().as_secs_f32();
        println!(
            "{} {:?} Simulation took {duration} secs",
            self.temp, self.algorithm
        );

        Ok(self.into_comp_result(result))
    }

    fn into_comp_result(self, result: f64) -> ComputationResult<'a, D, SIZE> {
        let [x, y, t]: [Option<usize>; 3] = [
            Some(self.lattice.size[0]),
            Some(self.lattice.size[1]),
            Some(self.lattice.size[2]),
        ];
        ComputationResult {
            d: 3,
            size: SIZE,
            x,
            y,
            t,
            temp: self.temp,
            comptype: Computation::WilsonSim(self.clone()),
            observable: self.observable,
            result,
            error: None,
        }
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
    range: usize,
    permutations: u64,
    observable: ObservableType,
}

impl<'a, const D: usize, const SIZE: usize> Compute<'a, D, SIZE> for Test<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    /// This function calculates the observables over all configurations, which
    /// only works for small lattices.
    fn run(self) -> Result<ComputationResult<'a, D, SIZE>, Box<dyn Error>> {
        let time = Instant::now();

        let field: Field<i8, D, SIZE> = Field::new(self.lattice);
        let mut field: Field<i32, D, SIZE> = Field::from_field(field);

        assert_eq!(SIZE, field.values.len());

        field.values[SIZE - 1] = (self.range / 2) as i32;

        // The partition function is the sum over all Bolzmann weights
        let mut partfn: ObservableValue<f64> = ObservableValue::default();

        // The observable is the sum over all weights (Bolzmann times observable),
        // devided by the partition function.
        let mut test: ObservableValue<f64> = ObservableValue::default();

        let boundary: i32 = self.range as i32 - 1;

        for _ in 0..self.permutations {
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
            test.update(self.observable.observe(&field, self.temp) * bolz);
            partfn.update(bolz);
        }
        let result: f64 = test.result() / partfn.result();

        let duration = time.elapsed().as_secs_f32();
        println!("{} Test took {duration} secs", self.temp);

        Ok(self.into_comp_result(result))
    }

    fn into_comp_result(self, result: f64) -> ComputationResult<'a, D, SIZE> {
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
            comptype: Computation::Test(self.clone()),
            observable: self.observable,
            result,
            error: None,
        }
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
    observable: ObservableType,
}

impl<'a, const D: usize, const SIZE: usize> Compute<'a, D, SIZE> for WilsonTest<'a, SIZE>
where
    [(); D * 2_usize]:,
{
    /// This function calculates the observables over all configurations, which
    /// only works for small lattices.
    fn run(self) -> Result<ComputationResult<'a, D, SIZE>, Box<dyn Error>> {
        let time = Instant::now();

        let field: Field<i8, 3, SIZE> = Field::new(self.lattice);
        let field: Field<i32, 3, SIZE> = Field::from_field(field);
        let mut wilson: WilsonField<'a, i32, SIZE> =
            WilsonField::from_field(field, self.width, self.height);

        assert_eq!(SIZE, wilson.field.values.len());

        wilson.field.values[SIZE - 1] = (self.range / 2) as i32;

        // The partition function is the sum over all Bolzmann weights
        let mut partfn: ObservableValue<f64> = ObservableValue::default();

        // The observable is the sum over all weights (Bolzmann times observable),
        // devided by the partition function.
        let mut test: ObservableValue<f64> = ObservableValue::default();

        let boundary: i32 = self.range as i32 - 1;

        for _ in 0..self.permutations {
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
            let bolz: f64 = (-wilson.action_observable(self.temp)).exp();
            test.update(self.observable.observe(&wilson.field, self.temp) * bolz);
            partfn.update(bolz);
        }

        let result: f64 = test.result() / partfn.result();

        let duration = time.elapsed().as_secs_f32();
        println!("{} Test took {duration} secs", self.temp);

        Ok(self.into_comp_result(result))
    }

    fn into_comp_result(self, result: f64) -> ComputationResult<'a, D, SIZE> {
        let [x, y, t]: [Option<usize>; 3] = [
            Some(self.lattice.size[0]),
            Some(self.lattice.size[1]),
            Some(self.lattice.size[2]),
        ];
        ComputationResult {
            d: D,
            size: SIZE,
            x,
            y,
            t,
            temp: self.temp,
            comptype: Computation::WilsonTest(self.clone()),
            observable: self.observable,
            result,
            error: None,
        }
    }
}

// - Results ------------------------------------------------------------------

#[derive(Debug)]
/// The result of a [`Computation`] that has finished.
pub struct ComputationResult<'a, const D: usize, const SIZE: usize>
where
    [(); D * 2_usize]:,
{
    d: usize,
    size: usize,
    x: Option<usize>,
    y: Option<usize>,
    t: Option<usize>,
    temp: f64,
    comptype: Computation<'a, D, SIZE>,
    observable: ObservableType,
    result: f64,
    error: Option<f64>,
}

impl<'a, const D: usize, const SIZE: usize> ComputationResult<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
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
