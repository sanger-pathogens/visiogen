use chrono::Local;
use log::*;
use simplelog::*;
use std::fs::File;

pub fn set_up_logging() {
    let current_time = Local::now().format("%m-%d_%H-%M-%S").to_string();
    let log_filename = format!("visiogen_{}.log", current_time);

    CombinedLogger::init(vec![
        TermLogger::new(
            LevelFilter::Warn,
            Config::default(),
            TerminalMode::Mixed,
            ColorChoice::Auto,
        ),
        WriteLogger::new(
            LevelFilter::Info,
            Config::default(),
            File::create(log_filename).unwrap(),
        ),
    ])
    .unwrap();
}
