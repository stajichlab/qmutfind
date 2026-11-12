//! End-to-end test of `qmutfind fetch` against a local mock of the FungiDB web service.

use std::fs;
use std::io::{Read, Write};
use std::net::TcpListener;
use std::path::Path;
use std::process::Command;
use std::sync::{Arc, Mutex};
use std::thread;

const KEY: &str = "secret-test-key-123";
const FASTA: &str = ">chr1\nACGTACGT\n";
const GFF: &str = "##gff-version 3\nchr1\tsrc\tgene\t1\t8\t.\t+\t.\tID=g1\n";

struct Mock {
    base: String,
    auth_headers: Arc<Mutex<Vec<String>>>,
}

/// A server that answers three paths. `raw-files` needs `Authorization: Bearer KEY` (or the bare key when `bare`).
fn mock(bare: bool) -> Mock {
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let base = format!("http://{}", listener.local_addr().unwrap());
    let auth_headers = Arc::new(Mutex::new(Vec::new()));
    let seen = auth_headers.clone();
    let base2 = base.clone();
    thread::spawn(move || {
        for stream in listener.incoming() {
            let mut stream = match stream { Ok(s) => s, Err(_) => break };
            let mut buf = [0u8; 8192];
            let n = stream.read(&mut buf).unwrap_or(0);
            let req = String::from_utf8_lossy(&buf[..n]).to_string();
            let path = req.split_whitespace().nth(1).unwrap_or("").to_string();
            let auth = req.lines().find_map(|l| l.strip_prefix("Authorization: ").or_else(|| l.strip_prefix("authorization: "))).map(str::to_string);
            if path.contains("raw-files") || path.contains("downloads") { seen.lock().unwrap().push(auth.clone().unwrap_or_default()); }
            let ok_auth = auth.as_deref() == Some(&format!("Bearer {KEY}")) && !bare || auth.as_deref() == Some(KEY) && bare;
            let (status, body): (&str, String) = if path.contains("GenomeDataTypes") {
                let json = format!(
                    r#"{{"records":[{{"id":[],"attributes":{{"primary_key":"Aspergillus fumigatus Af293","URLGenomeFasta":"{base2}/a/service/raw-files/Current_Release/Af/fasta/data/FungiDB-71_Af_Genome.fasta","URLgff":"{base2}/common/downloads/Current_Release/Af/gff/data/FungiDB-71_Af.gff"}}}},{{"id":[],"attributes":{{"primary_key":"Aspergillus fumigatus A1163","URLGenomeFasta":"x","URLgff":"y"}}}}]}}"#
                );
                ("200 OK", json)
            } else if path.contains("raw-files") {
                if ok_auth { ("200 OK", if path.contains("/gff/") { GFF } else { FASTA }.to_string()) } else { ("401 Unauthorized", "Valid API Key required for this endpoint.".to_string()) }
            } else if path.contains("downloads") {
                // The real site lists /common/downloads URLs that now answer 404; the file is under raw-files.
                ("404 Not Found", String::new())
            } else {
                ("404 Not Found", String::new())
            };
            let resp = format!("HTTP/1.1 {status}\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{body}", body.len());
            let _ = stream.write_all(resp.as_bytes());
        }
    });
    Mock { base, auth_headers }
}

fn run(dir: &Path, base: &str, key: Option<&str>, organism: &str) -> std::process::Output {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_qmutfind"));
    cmd.args(["fetch", "--organism", organism, "--out-dir", "out", "--base-url", base]).current_dir(dir).env_remove("EUPATHDB_KEY");
    if let Some(k) = key { cmd.env("EUPATHDB_KEY", k); }
    cmd.output().unwrap()
}

fn gunzip(path: &Path) -> String {
    let mut s = String::new();
    flate2::read::MultiGzDecoder::new(fs::File::open(path).unwrap()).read_to_string(&mut s).unwrap();
    s
}

#[test]
fn fetch_downloads_the_fasta_and_gff_as_gzip_and_sends_the_key_as_a_bearer_token() {
    let m = mock(false);
    let dir = tempfile::tempdir().unwrap();

    let out = run(dir.path(), &m.base, Some(KEY), "Aspergillus fumigatus Af293");

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    assert_eq!(gunzip(&dir.path().join("out/FungiDB-71_Af_Genome.fasta.gz")), FASTA);
    assert_eq!(gunzip(&dir.path().join("out/FungiDB-71_Af.gff.gz")), GFF);
    assert!(m.auth_headers.lock().unwrap().iter().any(|h| h == &format!("Bearer {KEY}")));
}

#[test]
fn the_key_never_appears_in_the_output() {
    let m = mock(false);
    let dir = tempfile::tempdir().unwrap();

    let out = run(dir.path(), &m.base, Some(KEY), "Aspergillus fumigatus Af293");

    let text = format!("{}{}", String::from_utf8_lossy(&out.stdout), String::from_utf8_lossy(&out.stderr));
    assert!(!text.contains(KEY));
}

#[test]
fn a_rejected_bearer_token_is_retried_with_the_bare_key() {
    let m = mock(true);
    let dir = tempfile::tempdir().unwrap();

    let out = run(dir.path(), &m.base, Some(KEY), "Aspergillus fumigatus Af293");

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    assert_eq!(gunzip(&dir.path().join("out/FungiDB-71_Af_Genome.fasta.gz")), FASTA);
}

#[test]
fn without_a_key_the_error_names_the_variable_and_no_file_is_left_behind() {
    let m = mock(false);
    let dir = tempfile::tempdir().unwrap();

    let out = run(dir.path(), &m.base, None, "Aspergillus fumigatus Af293");

    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("EUPATHDB_KEY"), "{err}");
    assert!(!dir.path().join("out/FungiDB-71_Af_Genome.fasta.gz").exists());
    assert!(!dir.path().join("out/FungiDB-71_Af_Genome.fasta.gz.part").exists());
}

#[test]
fn an_unknown_organism_lists_close_matches() {
    let m = mock(false);
    let dir = tempfile::tempdir().unwrap();

    let out = run(dir.path(), &m.base, Some(KEY), "fumigatus");

    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("Aspergillus fumigatus Af293") && err.contains("Aspergillus fumigatus A1163"), "{err}");
}

#[test]
fn an_organism_name_broken_across_lines_by_a_paste_still_matches() {
    let m = mock(false);
    let dir = tempfile::tempdir().unwrap();

    let out = run(dir.path(), &m.base, Some(KEY), "  Aspergillus   fumigatus\n  Af293 ");

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
}
