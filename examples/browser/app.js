// The page: collect the files, start the worker, show the progress and the results.
const $ = (id) => document.getElementById(id);
const worker = new Worker('worker.js', { type: 'module' });

function show(id, text) {
  $(id).textContent = text;
}

function download(name, text, type) {
  const url = URL.createObjectURL(new Blob([text], { type }));
  const link = $('download-' + name);
  link.href = url;
  link.download = `${$('sample').value || 'sample'}.${name === 'vcf' ? 'vcf' : 'stats.json'}`;
  link.hidden = false;
}

$('run').addEventListener('click', () => {
  const reference = $('reference').files[0];
  const reads = [...$('reads').files];
  if (!reference || reads.length === 0) {
    show('status', 'Choose a reference FASTA and at least one reads file.');
    return;
  }
  $('run').disabled = true;
  show('status', 'Starting');
  worker.postMessage({
    reference,
    reads,
    gff: $('gff').files[0] || null,
    sample: $('sample').value || 'sample',
    kind: $('kind').value,
    region: $('region').value.trim() || null,
  });
});

worker.onmessage = (event) => {
  const message = event.data;
  if (message.type === 'status') show('status', message.text);
  if (message.type === 'progress') show('progress', `${message.reads.toLocaleString()} reads mapped`);
  if (message.type === 'error') {
    show('status', `Error: ${message.message}`);
    $('run').disabled = false;
  }
  if (message.type === 'done') {
    show('status', `Done: ${message.reads.toLocaleString()} reads in ${message.seconds.toFixed(1)} s`);
    const summary = JSON.parse(message.stats);
    show('summary', `Calls: ${summary.total}. Ts/Tv: ${summary.ts_tv === null ? 'not defined' : summary.ts_tv.toFixed(2)}. ` +
      `Calls with a fraction from 0.2 to 0.8: ${summary.mid_fraction_calls}.`);
    download('vcf', message.vcf, 'text/plain');
    download('stats', message.stats, 'application/json');
    if (message.report) $('report').srcdoc = message.report;
    $('run').disabled = false;
  }
};
