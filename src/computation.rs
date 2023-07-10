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

    pub fn new_test(
        lattice: &'a Lattice<D, SIZE>,
        temp: f64,
        range: usize,
    ) -> Computation<'a, D, SIZE> {
        // 16 ^ 7 = 268’435’456
        let Some(permutations): Option<u64> = (range as u64)
            .checked_pow(SIZE as u32 - 1) else { panic!("Permutations overflow.")};
        let mut output: Vec<OutputData<D, SIZE>> = Vec::new();
        output.push(OutputData::new_test_action_observable(temp));
        Computation::Test(Test {
            lattice,
            temp,
            range,
            permutations,
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
        burnin: usize,
        iterations: usize,
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
    pub fn new_wilson_test(
        lattice: &'a Lattice<3, SIZE>,
        temp: f64,
        range: usize,
        width: usize,
        height: usize,
    ) -> Computation<'a, 3, SIZE> {
        // 16 ^ 7 = 268’435’456
        let Some(permutations): Option<u64> = (range as u64)
            .checked_pow(SIZE as u32 - 1) else { panic!("Permutations overflow.")};
        let mut output: Vec<OutputData<'a, 3, SIZE>> = Vec::new();
        output.push(OutputData::new_test_action_observable(temp));
        Computation::WilsonTest(WilsonTest {
            lattice,
            temp,
            range,
            permutations,
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
        let index: usize = match ComputationSummary::fetch_csv_data(crate::RESULTS_PATH, true)
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
                    summary = summary.set_action(obs.result());
                }
                OutputDataType::TestActionObservable(obs) => {
                    summary = summary.set_action(obs.result());
                }
                OutputDataType::DifferencePlot(obs) => {
                    for (direction, plot) in obs.plot().into_iter().enumerate() {
                        let path: &str = &(crate::DIFFERENCE_PLOT_PATH_INCOMPLETE.to_owned()
                            + &"difference_"
                            + &index.to_string()
                            + &"_"
                            + &direction.to_string()
                            + &".csv");
                        if let Err(err) = plot.read_write_csv(path, true) {
                            eprint!("{}", err);
                        };
                    }
                    summary = summary.set_bonds_data();
                }
                OutputDataType::EnergyPlot(obs) => {
                    for (direction, plot) in obs.plot().into_iter().enumerate() {
                        let path: &str = &(crate::ENERGY_PLOT_PATH_INCOMPLETE.to_owned()
                            + &"energy_"
                            + &index.to_string()
                            + &"_"
                            + &direction.to_string()
                            + &".csv");
                        if let Err(err) = plot.read_write_csv(path, true) {
                            eprint!("{}", err);
                        };
                    }
                    summary = summary.set_energy_data();
                }
                OutputDataType::CorrelationData(obs) => {
                    if let Some(plot) = obs.plot().into_iter().next() {
                        let path: &str = &(crate::CORR_FN_PATH_INCOMPLETE.to_owned()
                            + &"correlation_"
                            + &index.to_string()
                            + &".csv");
                        if let Err(err) = plot.overwrite_csv(path) {
                            eprint!("{}", err);
                        };
                        summary = summary.set_correlation_data();
                    }
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
    burnin: usize,
    iterations: usize,
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
            if step == 0 {
                println!("{:.2} {: <12} Started", self.temp, "Burnin:");
            } else if step == self.burnin / 4 {
                println!("{:.2} {: <12} 25%", self.temp, "Burnin:");
            } else if step == self.burnin / 2 {
                println!("{:.2} {: <12} 50%", self.temp, "Burnin:");
            } else if step == self.burnin / 4 * 3 {
                println!("{:.2} {: <12} 75%", self.temp, "Burnin:");
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
    burnin: usize,
    iterations: usize,
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
    range: usize,
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

#[derive(Debug, Serialize, Deserialize, Clone)]
/// We export the [`ComputationResult`] in a csv file with this formating
pub struct ComputationSummary {
    pub index: usize,
    pub d: Option<usize>,
    pub size: Option<usize>,
    pub x: Option<usize>,
    pub y: Option<usize>,
    pub t: Option<usize>,
    pub temp: Option<f64>,
    pub comptype: Option<String>,
    pub comptime: Option<u64>,
    pub action: Option<f64>,
    pub energy_data: bool,
    pub difference_data: bool,
    pub correlation_data: bool,
    pub corr12: Option<f64>,
    pub corr12_err: Option<f64>,
}

impl ComputationSummary {
    pub fn new(index: usize) -> Self {
        ComputationSummary {
            index,
            d: None,
            size: None,
            x: None,
            y: None,
            t: None,
            temp: None,
            comptype: None,
            comptime: None,
            action: None,
            energy_data: false,
            difference_data: false,
            correlation_data: false,
            corr12: None,
            corr12_err: None,
        }
    }

    pub fn from_computation<const SIZE: usize>(
        computation: Computation<3, SIZE>,
        index: usize,
    ) -> (Self, Vec<OutputData<3, SIZE>>) {
        let comptype: String = computation.to_string();
        match computation {
            Computation::Simulation(comp) => {
                let new: ComputationSummary = ComputationSummary::new(index)
                    .set_size(comp.lattice.size)
                    .set_computation(comp.temp, comptype, comp.duration);
                let outputs: Vec<OutputData<3, SIZE>> = comp.output;
                (new, outputs)
            }
            Computation::Test(comp) => {
                let new: ComputationSummary = ComputationSummary::new(index)
                    .set_size(comp.lattice.size)
                    .set_computation(comp.temp, comptype, comp.duration);
                let outputs: Vec<OutputData<3, SIZE>> = comp.output;
                (new, outputs)
            }
            Computation::WilsonSim(comp) => {
                let new: ComputationSummary = ComputationSummary::new(index)
                    .set_size(comp.lattice.size)
                    .set_computation(comp.temp, comptype, comp.duration);
                let outputs: Vec<OutputData<3, SIZE>> = comp.output;
                (new, outputs)
            }
            Computation::WilsonTest(comp) => {
                let new: ComputationSummary = ComputationSummary::new(index)
                    .set_size(comp.lattice.size)
                    .set_computation(comp.temp, comptype, comp.duration);
                let outputs: Vec<OutputData<3, SIZE>> = comp.output;
                (new, outputs)
            }
        }
    }

    pub fn set_size<const D: usize>(mut self, size_ary: [usize; D]) -> Self
    where
        [(); D]:,
    {
        self.d = Some(D);
        self.size = Some(size_ary.iter().product());
        match D {
            1 => self.x = Some(size_ary[0]),
            2 => {
                self.x = Some(size_ary[0]);
                self.y = Some(size_ary[1]);
            }
            3 => {
                self.x = Some(size_ary[0]);
                self.y = Some(size_ary[1]);
                self.t = Some(size_ary[2]);
            }
            _ => {}
        }
        self
    }

    pub fn set_computation(mut self, temp: f64, comptype: String, duration: Option<u64>) -> Self {
        self.temp = Some(temp);
        self.comptype = Some(comptype);
        self.comptime = duration;
        self
    }

    pub fn set_action(mut self, result: f64) -> Self {
        self.action = Some(result);
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

    pub fn set_correlation_length(mut self, corr12: f64) -> Self {
        self.corr12 = Some(corr12);
        self
    }

    pub fn set_correlation_length_error(mut self, corr12_err: f64) -> Self {
        self.corr12_err = Some(corr12_err);
        self
    }
}

// - Field - Exporting --------------------------------------------------------

#[derive(Serialize, Deserialize)]
pub struct FieldExport3d<T>
where
    T: Serialize,
{
    x: usize,
    y: usize,
    t: usize,
    val: T,
}

impl<'a, T, const SIZE: usize> Field<'a, T, 3, SIZE>
where
    T: Serialize,
{
    pub fn into_export(self) -> Vec<FieldExport3d<T>> {
        let mut ary: Vec<FieldExport3d<T>> = Vec::new();
        for (index, val) in self.values.into_iter().enumerate() {
            let [x, y, t] = self.lattice.calc_coords_from_index(index).into_array();
            ary.push(FieldExport3d { x, y, t, val })
        }
        ary
    }
}
