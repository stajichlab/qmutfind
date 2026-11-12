//! `qmutfind fetch`: download the genome FASTA and GFF of one organism from FungiDB.
//!
//! The organism list (with the file URLs) is public. The FASTA download needs an API key, which is read from the
//! environment variable `EUPATHDB_KEY`, or from `--eupathdb-key-file` (which wins when both are given). The key goes only to the host of `--base-url`, and it is never printed.

use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::Path;

const DEFAULT_BASE: &str = "https://fungidb.org";
const KEY_VARIABLE: &str = "EUPATHDB_KEY";
const MAX_REDIRECTS: usize = 5;

/// Percent-encode everything except the unreserved characters.
fn percent_encode(text: &str) -> String {
    let mut out = String::new();
    for b in text.bytes() {
        if b.is_ascii_alphanumeric() || matches!(b, b'-' | b'_' | b'.' | b'~') {
            out.push(b as char);
        } else {
            out.push_str(&format!("%{b:02X}"));
        }
    }
    out
}

/// The host part of an http(s) URL, lower case. `None` when the URL has another form.
fn host_of(url: &str) -> Option<String> {
    let rest = url.strip_prefix("https://").or_else(|| url.strip_prefix("http://"))?;
    Some(rest.split(['/', '?', '#']).next()?.to_ascii_lowercase())
}

struct Organism {
    name: String,
    fasta: String,
    gff: String,
}

/// Ask the site for all organisms with their file URLs. This request needs no key.
fn list_organisms(agent: &ureq::Agent, base: &str) -> Result<Vec<Organism>, String> {
    let config = r#"{"attributes":["primary_key","URLGenomeFasta","URLgff"],"tables":[],"attributeFormat":"text"}"#;
    let url = format!(
        "{base}/fungidb/service/record-types/organism/searches/GenomeDataTypes/reports/standard?reportConfig={}",
        percent_encode(config)
    );
    let body = agent.get(&url).call().map_err(|e| format!("cannot list organisms from {base}: {e}"))?.into_string().map_err(|e| e.to_string())?;
    let json: serde_json::Value = serde_json::from_str(&body).map_err(|e| format!("the organism list is not JSON: {e}"))?;
    let records = json["records"].as_array().ok_or("the organism list has no records")?;
    Ok(records
        .iter()
        .filter_map(|r| {
            let a = &r["attributes"];
            Some(Organism { name: a["primary_key"].as_str()?.to_string(), fasta: a["URLGenomeFasta"].as_str()?.to_string(), gff: a["URLgff"].as_str()?.to_string() })
        })
        .collect())
}

/// Pick the organism whose name equals the request (ignoring case). Otherwise list the names that contain it.
fn choose<'a>(organisms: &'a [Organism], wanted: &str) -> Result<&'a Organism, String> {
    let w = wanted.split_whitespace().collect::<Vec<_>>().join(" ").to_lowercase();
    if let Some(o) = organisms.iter().find(|o| o.name.to_lowercase() == w) {
        return Ok(o);
    }
    let close: Vec<&str> = organisms.iter().filter(|o| o.name.to_lowercase().contains(&w)).map(|o| o.name.as_str()).collect();
    if close.is_empty() {
        Err(format!("no organism named \"{wanted}\" at this site"))
    } else {
        Err(format!("no organism named exactly \"{wanted}\". Close matches:\n  {}", close.join("\n  ")))
    }
}

/// One GET with manual redirects. The key is sent only while the host stays the same as `key_host`.
fn get(agent: &ureq::Agent, url: &str, auth: Option<&str>, key_host: &str) -> Result<ureq::Response, ureq::Error> {
    let mut url = url.to_string();
    for _ in 0..=MAX_REDIRECTS {
        let mut request = agent.get(&url);
        if let (Some(value), Some(host)) = (auth, host_of(&url)) {
            if host == key_host {
                request = request.set("Authorization", value);
            }
        }
        let response = request.call()?;
        if matches!(response.status(), 301 | 302 | 303 | 307 | 308) {
            match response.header("Location") {
                Some(next) => url = next.to_string(),
                None => return Ok(response),
            }
        } else {
            return Ok(response);
        }
    }
    Err(ureq::Error::Status(310, ureq::Response::new(310, "Too many redirects", "").unwrap()))
}

/// Download `url` into `out_dir`, compressing with gzip unless the name already ends in `.gz`.
fn download(agent: &ureq::Agent, url: &str, key: Option<&str>, base: &str, out_dir: &Path) -> Result<String, String> {
    // Never send the key over plain http when the site is https.
    // The organism list gives `/common/downloads/` URLs. They answer 404 on the real site (2026-11-12). The same
    // file is served, with the key, under `/a/service/raw-files/`.
    let url = url.replace("/common/downloads/", "/a/service/raw-files/");
    let url = url.as_str();
    let url = if base.starts_with("https://") && url.starts_with("http://") { format!("https://{}", &url["http://".len()..]) } else { url.to_string() };
    let key_host = host_of(base).ok_or_else(|| format!("--base-url must start with http:// or https://, not {base}"))?;
    let name = url.rsplit('/').next().filter(|n| !n.is_empty()).ok_or_else(|| format!("cannot find a file name in {url}"))?;
    let already_gz = name.ends_with(".gz");
    let final_name = if already_gz { name.to_string() } else { format!("{name}.gz") };

    let bearer = key.map(|k| format!("Bearer {k}"));
    let schemes: Vec<Option<String>> = match key {
        Some(k) => vec![bearer, Some(k.to_string())],
        None => vec![None],
    };
    let mut response = None;
    for (i, auth) in schemes.iter().enumerate() {
        match get(agent, &url, auth.as_deref(), &key_host) {
            Ok(r) => {
                response = Some(r);
                break;
            }
            Err(ureq::Error::Status(401, _)) if key.is_none() => {
                return Err(format!("{url} needs an API key. Set the environment variable {KEY_VARIABLE} (get the key from your FungiDB account page, Web Services Access)"));
            }
            Err(ureq::Error::Status(401, _)) if i + 1 < schemes.len() => continue,
            Err(ureq::Error::Status(401, _)) => return Err(format!("the server rejected the key in {KEY_VARIABLE} (HTTP 401)")),
            Err(ureq::Error::Status(code, _)) => return Err(format!("{url}: HTTP {code}")),
            Err(e) => return Err(format!("{url}: {e}")),
        }
    }
    let response = response.ok_or("no response")?;

    fs::create_dir_all(out_dir).map_err(|e| format!("cannot create {}: {e}", out_dir.display()))?;
    let final_path = out_dir.join(&final_name);
    let part_path = out_dir.join(format!("{final_name}.part"));
    let result = (|| -> io::Result<()> {
        let file = BufWriter::new(File::create(&part_path)?);
        let mut reader = response.into_reader();
        if already_gz {
            let mut file = file;
            io::copy(&mut reader, &mut file)?;
            file.flush()?;
        } else {
            let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
            io::copy(&mut reader, &mut encoder)?;
            encoder.finish()?.flush()?;
        }
        fs::rename(&part_path, &final_path)
    })();
    if let Err(e) = result {
        let _ = fs::remove_file(&part_path);
        return Err(format!("cannot save {final_name}: {e}"));
    }
    Ok(final_path.display().to_string())
}

/// Read the key for `host` from a key file. Each line is `KEY` (used for any site) or `HOST KEY` (used for that
/// site only). A line for the host wins over a line without a host. Blank lines and lines that start with `#` are
/// ignored. Keys are never printed.
fn key_from_file(path: &str, host: &str) -> Result<String, String> {
    use std::os::unix::fs::PermissionsExt;
    let meta = fs::metadata(path).map_err(|e| format!("cannot read key file {path}: {e}"))?;
    if meta.permissions().mode() & 0o077 != 0 {
        eprintln!("qmutfind: warning: key file {path} is readable by other users. Run: chmod 600 {path}");
    }
    let text = fs::read_to_string(path).map_err(|e| format!("cannot read key file {path}: {e}"))?;
    let mut default = None;
    for line in text.lines().map(str::trim).filter(|l| !l.is_empty() && !l.starts_with('#')) {
        let mut words = line.split_whitespace();
        match (words.next(), words.next(), words.next()) {
            (Some(key), None, _) => default = default.or(Some(key)),
            (Some(h), Some(key), None) if h.eq_ignore_ascii_case(host) => return Ok(key.to_string()),
            _ => {}
        }
    }
    default.map(str::to_string).ok_or_else(|| format!("no key for {host} in key file {path}"))
}

pub fn fetch(args: &[String]) -> Result<(), String> {
    let flags = crate::parse_flags(args)?;
    let wanted = crate::required(&flags, "organism")?;
    let out_dir = Path::new(crate::required(&flags, "out-dir")?);
    let base = flags.get("base-url").map(|s| s.trim_end_matches('/')).unwrap_or(DEFAULT_BASE);
    let key = match flags.get("eupathdb-key-file") {
        Some(path) => Some(key_from_file(path, &host_of(base).ok_or_else(|| format!("--base-url must start with http:// or https://, not {base}"))?)?),
        None => std::env::var(KEY_VARIABLE).ok().map(|k| k.trim().to_string()).filter(|k| !k.is_empty()),
    };

    let agent = ureq::AgentBuilder::new().redirects(0).timeout_connect(std::time::Duration::from_secs(30)).build();
    let organisms = list_organisms(&agent, base)?;
    let organism = choose(&organisms, wanted)?;
    for url in [&organism.fasta, &organism.gff] {
        let saved = download(&agent, url, key.as_deref(), base, out_dir)?;
        eprintln!("qmutfind: wrote {saved}");
    }
    Ok(())
}
