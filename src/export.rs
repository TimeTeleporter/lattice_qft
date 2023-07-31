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
impl CsvData for ComputationSummary {}
impl CsvData for FieldExport3d<f64> {}
impl CsvData for f64 {}
impl CsvData for (u64, Vec<f64>) {}

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
    index: u64,
    bin_size: u64,
    pub m: f64,
    pub n: f64,
    pub a: f64,
    pub b: f64,
    pub corr: f64,
}

impl FitResult {
    pub fn new(index: u64, bin_size: u64, m: f64, n: f64, a: f64, b: f64) -> FitResult {
        FitResult {
            index,
            bin_size,
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
    pub index: u64,
    pub bin_size: u64,
    m12: Option<f64>,
    m23: Option<f64>,
    m34: Option<f64>,
    m13: Option<f64>,
    m24: Option<f64>,
    m14: Option<f64>,
    m12_err: Option<f64>,
    m23_err: Option<f64>,
    m34_err: Option<f64>,
    m13_err: Option<f64>,
    m24_err: Option<f64>,
    m14_err: Option<f64>,
    m12_sigma: Option<f64>,
    pub corr12: Option<f64>,
    corr23: Option<f64>,
    corr34: Option<f64>,
    corr13: Option<f64>,
    corr24: Option<f64>,
    corr14: Option<f64>,
    pub corr12_err: Option<f64>,
    corr23_err: Option<f64>,
    corr34_err: Option<f64>,
    corr13_err: Option<f64>,
    corr24_err: Option<f64>,
    corr14_err: Option<f64>,
    corr12_sigma: Option<f64>,
}

impl CorrelationLengths {
    pub fn new(index: u64, bin_size: u64, m12: f64) -> Self {
        CorrelationLengths {
            index,
            bin_size,
            m12: Some(m12),
            corr12: Some(1.0 / m12),
            ..Default::default()
        }
    }

    pub fn set_corr12_error(mut self, m12_err: f64) -> Self {
        self.m12_err = Some(m12_err);
        self.corr12_err = self.m12.map(|x| m12_err / (x * x));
        self
    }

    pub fn set_corr12_sigma(mut self, m12_sigma: f64) -> Self {
        self.m12_sigma = Some(m12_sigma);
        self.corr12_sigma = self.m12.map(|x| m12_sigma / (x * x));
        self
    }

    pub fn set_other_values(mut self, values: [Option<f64>; 5]) -> Self {
        [self.m23, self.m34, self.m13, self.m24, self.m14] = values;
        [
            self.corr23,
            self.corr34,
            self.corr13,
            self.corr24,
            self.corr14,
        ] = values.map(|opt| opt.map(|x| 1.0 / x));
        self
    }
}

impl CsvData for CorrelationLengths {}
