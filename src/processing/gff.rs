use bio::io::gff;
use bio_types::strand::Strand;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use crate::error::{Result, VisiogenError};

pub fn coords_from_gene_name(
    gff_path: &String,
    gene: &String,
) -> Result<Option<(u64, u64, Strand)>> {
    let path = Path::new(gff_path);
    let file = File::open(path).map_err(|e| VisiogenError::IoError(e))?;
    let reader = BufReader::new(file);
    let mut gff_reader = gff::Reader::new(reader, gff::GffType::GFF3);

    for record in gff_reader.records() {
        let rec = record.map_err(|e| VisiogenError::IoError(e.into()))?;
        if let Some(attributes) = rec.attributes().get("Name") {
            if attributes == gene {
                return Ok(Some((
                    *rec.start(),
                    *rec.end(),
                    rec.strand().unwrap_or(Strand::Forward),
                )));
            }
        }
    }
    Ok(None)
}
