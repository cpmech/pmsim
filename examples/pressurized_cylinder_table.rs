use pmsim::prelude::*;
use russell_lab::*;
use std::collections::HashMap;
use std::env;
use std::fmt::Write;
use std::fs;
use std::fs::File;
use std::io::Write as io_write;
use std::path::Path;
use std::process::Command;

const OUT_DIR: &str = "/tmp/pmsim/latex/";

fn main() -> Result<(), StrError> {
    // run(
    //     "pressurized_cylinder2d_elast",
    //     "_tri",
    //     &["tri3", "tri6", "tri10", "tri15"],
    // )?;
    run(
        "pressurized_cylinder2d_elast",
        "_qua",
        &["qua4", "qua8", "qua9", "qua12", "qua16", "qua17"],
    )?;
    // run(
    //     "pressurized_cylinder3d_elast",
    //     "",
    //     &["tet4", "tet10", "tet20", "hex8", "hex20", "hex32"],
    // )?;
    Ok(())
}

fn run(name: &str, suffix: &str, str_kinds: &[&str]) -> Result<(), StrError> {
    // genies
    let str_genies = &["mumps", "umfpack", "inteldss"];

    // aux strings
    let row_genies =
        "NDOF & \\multicolumn{3}{c|}{{MUMPS}} & \\multicolumn{3}{c|}{{UMFPACK}} & \\multicolumn{3}{c}{{Intel DSS}} \\\\";

    // table env
    let mut buf = String::new();
    writeln!(
        &mut buf,
        "\\begin{{table}}\n\
         \\begin{{tabular}}{{@{{}}r|rrr|rrr|rrr@{{}}}}\\toprule\n\
         {}",
        row_genies
    )
    .unwrap();

    // run for each kind
    for (counter, str_kind) in str_kinds.iter().enumerate() {
        let mut ndof: Vec<usize> = Vec::new();
        let mut mumps = HashMap::new();
        let mut umfpack = HashMap::new();
        let mut inteldss = HashMap::new();

        // read data
        for str_genie in str_genies {
            // check if results are available
            let json = format!("/tmp/pmsim/{}_{}_{}.json", name, str_genie, str_kind);
            let path_json = Path::new(&json);
            if path_json.exists() {
                // load results
                let results = ConvergenceResults::from(&path_json)?;
                assert_eq!(results.name, *str_kind);
                match *str_genie {
                    "mumps" => {
                        ndof = results.ndof.clone();
                        for i in 0..results.ndof.len() {
                            mumps.insert(results.ndof[i], (results.time[i], results.error[i]));
                        }
                    }
                    "umfpack" => {
                        for i in 0..results.ndof.len() {
                            umfpack.insert(results.ndof[i], (results.time[i], results.error[i]));
                        }
                    }
                    "inteldss" => {
                        for i in 0..results.ndof.len() {
                            inteldss.insert(results.ndof[i], (results.time[i], results.error[i]));
                        }
                    }
                    _ => panic!("genie not available"),
                }
            }
        }

        if ndof.len() == 0 {
            continue;
        } else {
            let row_labels = "\\bf{Time} & \\bf{Error} & \\bf{RT} & \\bf{Time} & \\bf{Error} & \\bf{RT} & \\bf{Time} & \\bf{Error} & \\bf{RT} \\\\\\midrule";
            writeln!(&mut buf, "\\bf{{{}}} & {}", up_first(str_kind), row_labels).unwrap();
        }

        // generate the contents
        for n in &ndof {
            let e;
            write!(&mut buf, "{:>6}", n).unwrap();
            let ref_time;
            if let Some(data) = mumps.get(n) {
                ref_time = data.0;
                e = format!("{:>8.2e}", data.1);
                write!(&mut buf, " & {:>8} & {:>9} & {:>4.2}", ts(data.0), e, 1.0).unwrap();
            } else {
                panic!("reference data (MUMPS) is not available");
            }
            if let Some(data) = umfpack.get(n) {
                let ee = format!("{:>8.2e}", data.1);
                let es = if ee == e { ee } else { format!("*{}", ee) };
                let rel_time = (data.0 as f64) / (ref_time as f64);
                write!(&mut buf, " & {:>8} & {:>9} & {:>4.2}", ts(data.0), es, rel_time).unwrap();
            }
            if let Some(data) = inteldss.get(n) {
                let ee = format!("{:>8.2e}", data.1);
                let es = if ee == e { ee } else { format!("*{}", ee) };
                let rel_time = (data.0 as f64) / (ref_time as f64);
                write!(&mut buf, " & {:>8} & {:>9} & {:>4.2}", ts(data.0), es, rel_time).unwrap();
            }
            writeln!(&mut buf, "\\\\").unwrap();
        }

        // gap
        if counter != str_kinds.len() - 1 {
            writeln!(&mut buf, "[0.5em]").unwrap();
        }
    }

    write!(
        &mut buf,
        "\\bottomrule\n\
         \\end{{tabular}}\n\
         \\end{{table}}"
    )
    .unwrap();

    let key = format!("{}{}", name, suffix);
    call_latex(&key, &buf)
}

fn call_latex(key: &String, buffer: &String) -> Result<(), StrError> {
    // create directory
    let name = format!("{}.tex", key);
    let filepath = format!("{}/{}", OUT_DIR, name);
    let path = Path::new(&filepath);
    if let Some(p) = path.parent() {
        fs::create_dir_all(p).map_err(|_| "cannot create directory")?;
    }

    // write file
    let contents = format!(
        "\\documentclass[10pt,a4paper]{{article}}\n\
         \\usepackage[a4paper,margin=1in]{{geometry}}\n\
         \\usepackage{{booktabs}}\n\
         \\usepackage{{graphicx}}\n\
         \\begin{{document}}\n\
         {}\n\
         \\end{{document}}",
        buffer
    );
    let mut file = File::create(path).map_err(|_| "cannot create file")?;
    file.write_all(contents.as_bytes()).map_err(|_| "cannot write file")?;

    // force sync
    file.sync_all().map_err(|_| "cannot sync file")?;

    // change working dir
    let out_dir = format!("{}", OUT_DIR);
    let working_dir = Path::new(&out_dir);
    env::set_current_dir(&working_dir).unwrap();

    // execute file
    let output = Command::new("latexmk")
        .args([
            "-pdf",
            "-shell-escape",
            "-halt-on-error",
            "-interaction=batchmode",
            &filepath,
        ])
        .output()
        .map_err(|_| "cannot run latexmk")?;

    // output from stderr
    let err = String::from_utf8(output.stderr).unwrap();
    if err != "" {
        println!("{}", err);
        return Err("latexmk failed");
    }
    Ok(())
}

const NS_PER_NANOSECOND: u128 = 1;
const NS_PER_MICROSECOND: u128 = 1000 * NS_PER_NANOSECOND;
const NS_PER_MILLISECOND: u128 = 1000 * NS_PER_MICROSECOND;
const NS_PER_SECOND: u128 = 1000 * NS_PER_MILLISECOND;
const NS_PER_MINUTE: u128 = 60 * NS_PER_SECOND;
const NS_PER_HOUR: u128 = 60 * NS_PER_MINUTE;

/// Formats the nanoseconds in 1 second; value < 1 second
fn format_nanoseconds_in_seconds(buf: &mut String, value: u128) {
    if value < NS_PER_MICROSECOND {
        write!(buf, "{:.2}ns", value).unwrap();
    } else if value < NS_PER_MILLISECOND {
        write!(buf, "{:.2}µs", (value as f64) / (NS_PER_MICROSECOND as f64)).unwrap();
    } else {
        write!(buf, "{:.2}ms", (value as f64) / (NS_PER_MILLISECOND as f64)).unwrap();
    }
}

/// Returns a pretty string representing the value in nanoseconds
pub fn ts(nanoseconds: u128) -> String {
    if nanoseconds == 0 {
        return "0ns".to_string();
    }

    let mut value = nanoseconds;
    let mut buf = String::new();
    if value < NS_PER_SECOND {
        // nanoseconds is smaller than a second => use small units such as 2.5ms
        format_nanoseconds_in_seconds(&mut buf, value);
    } else {
        // nanoseconds is greater than a second => use large units such as 3m2.5s
        if value >= NS_PER_HOUR {
            let hours = value / NS_PER_HOUR;
            value -= hours * NS_PER_HOUR;
            write!(&mut buf, "{}h", hours).unwrap();
        }
        if value >= NS_PER_MINUTE {
            let minutes = value / NS_PER_MINUTE;
            value -= minutes * NS_PER_MINUTE;
            write!(&mut buf, "{}m", minutes).unwrap();
        }
        if value > 0 {
            if value < NS_PER_SECOND {
                format_nanoseconds_in_seconds(&mut buf, value);
            } else {
                let seconds = (value as f64) / (NS_PER_SECOND as f64);
                write!(&mut buf, "{:.2}s", &seconds).unwrap();
            }
        }
    }

    buf
}

fn up_first(s: &str) -> String {
    let mut c = s.chars();
    match c.next() {
        None => String::new(),
        Some(f) => f.to_uppercase().collect::<String>() + c.as_str(),
    }
}
