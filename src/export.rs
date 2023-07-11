use std::{error::Error, fs::File};

use serde::{Deserialize, Serialize};

use crate::computation::{ComputationSummary, FieldExport3d};

pub trait CsvData
where
    Self: Serialize + for<'a> Deserialize<'a>,
{
    fn read_write_csv(self, path: &str, has_headers: bool) -> Result<(), Box<dyn Error>> {
        // Initialize data storage
        let mut storage: Vec<Self> = Self::fetch_csv_data(path, has_headers)?;

        // Append the new entry
        storage.push(self);

        // Write to the file
        write_to_csv(path, storage)?;

        Ok(())
    }

    fn overwrite_csv(self, path: &str) -> Result<(), Box<dyn Error>> {
        clean_csv(path)?;
        self.read_write_csv(path, false)
    }

    fn fetch_csv_data(path: &str, has_headers: bool) -> Result<Vec<Self>, Box<dyn Error>> {
        let mut storage: Vec<Self> = Vec::new();

        read_from_csv(path, &mut storage, has_headers)?;

        Ok(storage)
    }
}

fn write_to_csv<T>(path: &str, storage: Vec<T>) -> Result<(), Box<dyn Error>>
where
    T: Serialize + for<'a> Deserialize<'a>,
{
    let mut wtr = csv::Writer::from_path(path)?;
    for data in storage {
        wtr.serialize(data)?;
    }
    wtr.flush()?;
    Ok(())
}

fn read_from_csv<T>(
    path: &str,
    storage: &mut Vec<T>,
    has_headers: bool,
) -> Result<(), Box<dyn Error>>
where
    T: Serialize + for<'a> Deserialize<'a>,
{
    let file: File = match File::open(path) {
        Ok(file) => file,
        Err(err) => match err.kind() {
            std::io::ErrorKind::NotFound => {
                File::create(path)?;
                File::open(path)?
            }
            _ => Err(err)?,
        },
    };
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(has_headers)
        .from_reader(file);
    Ok(for result in rdr.deserialize() {
        let test_data: T = result?;
        storage.push(test_data);
    })
}

pub fn clean_csv(path: &str) -> Result<(), Box<dyn Error>> {
    let storage: Vec<()> = Vec::new();
    write_to_csv(path, storage)?;
    Ok(())
}

pub fn get_correlation_fn(index: usize, corr_fn_path: &str) -> Result<Vec<f64>, Box<dyn Error>> {
    let path: &str = &(corr_fn_path.to_owned() + &"correlation_" + &index.to_string() + &".csv");
    let corr_fn: Result<Vec<f64>, Box<dyn Error>> = f64::fetch_csv_data(path, false)
        .map_err(|err| format!("Fetching {}: {}", path, err).into());
    corr_fn
}

pub fn get_correlation_fn_with_err(
    index: usize,
    corr_fn_path: &str,
) -> Result<(Vec<f64>, Option<Vec<f64>>), Box<dyn Error>> {
    let corr_fn: Vec<f64> = get_correlation_fn(index, corr_fn_path)?;
    let path: &str =
        &(corr_fn_path.to_owned() + &"correlation_" + &index.to_string() + &"_err.csv");
    let corr_fn_err: Option<Vec<f64>> = f64::fetch_csv_data(path, false).ok();
    Ok((corr_fn, corr_fn_err))
}

impl CsvData for ComputationSummary {}
impl CsvData for FieldExport3d<f64> {}
impl CsvData for f64 {}

impl<T: CsvData> CsvData for Vec<T> {
    fn read_write_csv(mut self, path: &str, has_headers: bool) -> Result<(), Box<dyn Error>> {
        // Initialize data storage
        let mut storage: Vec<T> = <T as CsvData>::fetch_csv_data(path, has_headers)?;

        // Append the new entry
        storage.append(&mut self);

        // Write to the file
        write_to_csv(path, storage)?;

        Ok(())
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FitResult {
    index: usize,
    pub m: f64,
    pub n: f64,
    pub a: f64,
    pub b: f64,
    pub corr: f64,
}

impl FitResult {
    pub fn new(index: usize, m: f64, n: f64, a: f64, b: f64) -> FitResult {
        FitResult {
            index,
            m,
            n,
            a,
            b,
            corr: 1.0 / m,
        }
    }
}

impl CsvData for FitResult {}

#[derive(Debug, Serialize, Deserialize, Clone, Default)]
pub struct CorrelationLengths {
    index: usize,
    m12: f64,
    m23: f64,
    m34: f64,
    m13: f64,
    m24: f64,
    m14: f64,
    m12_err: Option<f64>,
    m23_err: Option<f64>,
    m34_err: Option<f64>,
    m13_err: Option<f64>,
    m24_err: Option<f64>,
    m14_err: Option<f64>,
    pub corr12: f64,
    corr23: f64,
    corr34: f64,
    corr13: f64,
    corr24: f64,
    corr14: f64,
    pub corr12_err: Option<f64>,
    corr23_err: Option<f64>,
    corr34_err: Option<f64>,
    corr13_err: Option<f64>,
    corr24_err: Option<f64>,
    corr14_err: Option<f64>,
}

impl CorrelationLengths {
    pub fn new(index: usize, values: [f64; 6]) -> CorrelationLengths {
        let [m12, m23, m34, m13, m24, m14] = values;
        let corrs: [f64; 6] = values.map(|x| 1.0 / x);
        let [corr12, corr23, corr34, corr13, corr24, corr14] = corrs;
        CorrelationLengths {
            index,
            m12,
            m23,
            m34,
            m13,
            m24,
            m14,
            corr12,
            corr23,
            corr34,
            corr13,
            corr24,
            corr14,
            ..CorrelationLengths::default()
        }
    }

    pub fn set_errors(&mut self, errors: [Option<f64>; 6]) {
        [
            self.m12_err,
            self.m23_err,
            self.m34_err,
            self.m13_err,
            self.m24_err,
            self.m14_err,
        ] = errors;
        let values = [self.m12, self.m23, self.m34, self.m13, self.m24, self.m14];
        let corr_errors: [Option<f64>; 6] = errors
            .zip(values)
            .map(|(delta_x, x)| delta_x.map(|delta| delta / (x * x)));
        [
            self.corr12_err,
            self.corr23_err,
            self.corr34_err,
            self.corr13_err,
            self.corr24_err,
            self.corr14_err,
        ] = corr_errors;
    }
}

impl CsvData for CorrelationLengths {}

#[derive(Debug, Serialize, Deserialize)]
pub struct OldComputationSummary {
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

impl CsvData for OldComputationSummary {}
