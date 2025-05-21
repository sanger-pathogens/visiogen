use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct FilteredKmers {
    pub gene: String,
    pub start: u64,
    pub end: u64,
    pub kmers: HashMap<String, Vec<usize>>,
    pub strand: String,
}
