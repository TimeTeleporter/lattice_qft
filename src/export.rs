use std::{error::Error, fs::File};

use serde::{Deserialize, Serialize};

pub trait CsvData {
    fn read_write_csv(self, path: &str) -> Result<(), Box<dyn Error>>;
    fn overwrite_csv(self, path: &str) -> Result<(), Box<dyn Error>>;
}

impl<T> CsvData for T
where
    T: Serialize + for<'a> Deserialize<'a>,
{
    fn read_write_csv(self, path: &str) -> Result<(), Box<dyn Error>> {
        // Initialize data storage
        let mut storage: Vec<T> = Vec::new();

        read_from_csv(path, &mut storage)?;

        // Append the new entry
        storage.push(self);

        // Write to the file
        write_to_csv(path, storage)?;

        Ok(())
    }

    fn overwrite_csv(self, path: &str) -> Result<(), Box<dyn Error>> {
        // Initialize data storage
        let mut storage: Vec<T> = Vec::new();

        // Append the new entry
        storage.push(self);

        write_to_csv(path, storage)?;

        Ok(())
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

fn read_from_csv<T>(path: &str, storage: &mut Vec<T>) -> Result<(), Box<dyn Error>>
where
    T: Serialize + for<'a> Deserialize<'a>,
{
    let file = File::open(path)?;
    let mut rdr = csv::Reader::from_reader(file);
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
