use std::{error::Error, fs::File};

use serde::{Deserialize, Serialize};

pub trait CsvData {
    fn write_to_csv(self, path: &str) -> Result<(), Box<dyn Error>>;
}

impl<T> CsvData for T
where
    T: Serialize + for<'a> Deserialize<'a>,
{
    fn write_to_csv(self, path: &str) -> Result<(), Box<dyn Error>> {
        // Initialize data storage
        let mut storage: Vec<T> = Vec::new();

        // Open the file and append it to the storage
        let file = File::open(path)?;
        let mut rdr = csv::Reader::from_reader(file);
        for result in rdr.deserialize() {
            let test_data: T = result?;
            storage.push(test_data);
        }

        // Append the new entry
        storage.push(self);

        // Write to the file
        let mut wtr = csv::Writer::from_path(path)?;
        for data in storage {
            wtr.serialize(data)?;
        }

        wtr.flush()?;

        Ok(())
    }
}
