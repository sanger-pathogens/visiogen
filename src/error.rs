use thiserror::Error;

#[derive(Error, Debug)]
pub enum VisiogenError {
    #[error("Failed to parse GFF3 file: {0}")]
    GffParseError(String),

    #[error("Failed to parse GFA file: {0}")]
    GfaParseError(String),

    #[error("Failed to build off-target indexes: {0}")]
    IndexBuildError(String),

    #[error("Failed to query off-target indexes: {0}")]
    IndexQueryError(String),

    #[error("Missing required argument: {0}")]
    MissingArgument(String),

    #[error("Gene processing error: {0}")]
    GeneProcessingError(String),

    #[error("I/O error: {0}")]
    IoError(#[from] std::io::Error),

    #[error("no unique k‑mers found in '{gene}' ({start}:{end})")]
    NoUniqueKmers { gene: String, start: u64, end: u64 },

    #[error("Other error: {0}")]
    Other(String),
}

pub type Result<T> = std::result::Result<T, VisiogenError>;
