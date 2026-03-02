nextflow.enable.dsl=2

def timestamp = new java.util.Date().format( 'yyyy-MM-dd_HH-mm' )

workflow {
  // 1. Feste Pfade für die Multi-Locus Analyse
  def db_cox_path = '/mnt/d/Epi2Me_Datenbanken/oomycetes_cox_ref.fasta'
  def db_its_path = '/mnt/d/Epi2Me_Datenbanken/oomycetes_its_ref.fasta'

  if( !params.reads ) {
    log.info "Bitte Reads angeben."
    return
  }

  // 2. Pre-Processing (Läuft nur einmal)
  ch_samples = Channel
    .fromPath("${params.reads}/*/*.fastq.gz")
    .map { f -> tuple(f.parent.name, f) }
    .groupTuple()
    .map { sample, files -> tuple(sample, files.sort()) }

  merged    = MERGE_FASTQ(ch_samples)
  fasta     = FASTQ_TO_FASTA(merged)
  clustered = CLUSTER_VSEARCH(fasta)
  kept      = FILTER_CLUSTERS(clustered)

  // 3. Datenbanken aufbauen (Parallel)
  db_cox_ch = MAKEBLASTDB_COX(Channel.fromPath(db_cox_path).first())
  db_its_ch = MAKEBLASTDB_ITS(Channel.fromPath(db_its_path).first())

  // 4. BLAST parallel ausführen
  kept_cox = kept.combine(db_cox_ch)
  kept_its = kept.combine(db_its_ch)

  blast_cox = BLASTN_COX(kept_cox)
  blast_its = BLASTN_ITS(kept_its)

  // 5. Taxonomy-Tabellen pro Locus erstellen
  tax_cox = JOIN_COUNTS_BLAST_COX(blast_cox)
  tax_its = JOIN_COUNTS_BLAST_ITS(blast_its)

  // 6. Sammeln und Consensus bilden
  tax_cox_list = tax_cox.map { sample, taxfile -> taxfile }.collect()
  tax_its_list = tax_its.map { sample, taxfile -> taxfile }.collect()

  summary = CONSENSUS_RESULTS(tax_cox_list, tax_its_list)
  REPORT_HTML(summary)
}

process MERGE_FASTQ {
  tag "$sample"
  container = null 
  publishDir "${params.out_dir}/per_sample/${sample}/reads", mode: 'copy', overwrite: true
  input: tuple val(sample), path(reads)
  output: tuple val(sample), path("${sample}.fastq.gz")
  script: "cat ${reads.join(' ')} > ${sample}.fastq.gz"
}

process FASTQ_TO_FASTA {
  tag "$sample"
  input: tuple val(sample), path(fq)
  output: tuple val(sample), path("${sample}.fasta")
  script: "seqkit fq2fa ${fq} > ${sample}.fasta"
}

process CLUSTER_VSEARCH {
  tag "$sample"
  input: tuple val(sample), path(fa)
  output: tuple val(sample), path("${sample}.centroids.fasta"), path("${sample}.clusters.uc")
  script: "vsearch --cluster_fast ${fa} --id 0.85 --minseqlength 10 --strand both --threads ${task.cpus} --uc ${sample}.clusters.uc --centroids ${sample}.centroids.fasta"
}

process FILTER_CLUSTERS {
  tag "$sample"
  container = null 
  input: tuple val(sample), path(centroids), path(uc)
  output: tuple val(sample), path("${sample}.cluster_counts.tsv"), path("${sample}.centroids.kept.fasta")
  script:
  """
  awk -F'\\t' 'BEGIN{OFS="\\t"} \$1=="S"{cl=\$2; id=\$9; c[cl]=id; n[id]=1} \$1=="H"{n[c[\$2]]++} END{for(i in n) print i,n[i]}' ${uc} | sort -k2,2nr > ${sample}.cluster_counts.tsv
  awk '\$2>=1{print \$1}' ${sample}.cluster_counts.tsv > keep.txt
  if [ -s keep.txt ]; then cp ${centroids} ${sample}.centroids.kept.fasta; else touch ${sample}.centroids.kept.fasta; fi
  """
}

// ==========================================
// PARALLELE PROZESSE FÜR COX UND ITS
// ==========================================

process MAKEBLASTDB_COX {
  input: path(db_fasta)
  output: path("blastdb_cox")
  script: "mkdir -p blastdb_cox && cp ${db_fasta} blastdb_cox/db.fasta && makeblastdb -in blastdb_cox/db.fasta -dbtype nucl -out blastdb_cox/cox_db"
}

process MAKEBLASTDB_ITS {
  input: path(db_fasta)
  output: path("blastdb_its")
  script: "mkdir -p blastdb_its && cp ${db_fasta} blastdb_its/db.fasta && makeblastdb -in blastdb_its/db.fasta -dbtype nucl -out blastdb_its/its_db"
}

process BLASTN_COX {
  tag "$sample"
  input: tuple val(sample), path(counts), path(centroids_fa), path(dbdir)
  output: tuple val(sample), path(counts), path("${sample}.cox.blast.tsv")
  script: "if [ ! -s ${centroids_fa} ]; then touch ${sample}.cox.blast.tsv; exit 0; fi; blastn -query ${centroids_fa} -db ${dbdir}/cox_db -max_target_seqs 5 -num_threads ${task.cpus} -outfmt '6 qseqid sseqid pident length qlen bitscore evalue qcovs stitle' > ${sample}.cox.blast.tsv"
}

process BLASTN_ITS {
  tag "$sample"
  input: tuple val(sample), path(counts), path(centroids_fa), path(dbdir)
  output: tuple val(sample), path(counts), path("${sample}.its.blast.tsv")
  script: "if [ ! -s ${centroids_fa} ]; then touch ${sample}.its.blast.tsv; exit 0; fi; blastn -query ${centroids_fa} -db ${dbdir}/its_db -max_target_seqs 5 -num_threads ${task.cpus} -outfmt '6 qseqid sseqid pident length qlen bitscore evalue qcovs stitle' > ${sample}.its.blast.tsv"
}

process JOIN_COUNTS_BLAST_COX {
  tag "$sample"
  container = null
  input: tuple val(sample), path(counts_tsv), path(blast_tsv)
  output: tuple val(sample), path("${sample}.cox.taxonomy.tsv")
  script:
  """
  python - << 'PY'
import csv
min_id, min_cov = float("${params.min_identity}"), float("${params.min_coverage}")
hits = {}
try:
    with open("${blast_tsv}") as f:
        for row in csv.reader(f, delimiter='\\t'):
            if not row or row[0] in hits: continue
            if float(row[2]) < min_id or float(row[7]) < min_cov: continue
            hits[row[0]] = {'stitle': row[8] if len(row) > 8 else 'NA'}
except Exception: pass

with open("${sample}.cox.taxonomy.tsv", 'w') as out:
    out.write("cluster_id\\tread_count\\tbest_hit_title\\n")
    with open("${counts_tsv}") as f:
        for row in csv.reader(f, delimiter='\\t'):
            if row: out.write(f"{row[0]}\\t{row[1]}\\t{hits.get(row[0], {'stitle':'Unclassified'})['stitle']}\\n")
PY
  """
}

process JOIN_COUNTS_BLAST_ITS {
  tag "$sample"
  container = null
  input: tuple val(sample), path(counts_tsv), path(blast_tsv)
  output: tuple val(sample), path("${sample}.its.taxonomy.tsv")
  script:
  """
  python - << 'PY'
import csv
min_id, min_cov = float("${params.min_identity}"), float("${params.min_coverage}")
hits = {}
try:
    with open("${blast_tsv}") as f:
        for row in csv.reader(f, delimiter='\\t'):
            if not row or row[0] in hits: continue
            if float(row[2]) < min_id or float(row[7]) < min_cov: continue
            hits[row[0]] = {'stitle': row[8] if len(row) > 8 else 'NA'}
except Exception: pass

with open("${sample}.its.taxonomy.tsv", 'w') as out:
    out.write("cluster_id\\tread_count\\tbest_hit_title\\n")
    with open("${counts_tsv}") as f:
        for row in csv.reader(f, delimiter='\\t'):
            if row: out.write(f"{row[0]}\\t{row[1]}\\t{hits.get(row[0], {'stitle':'Unclassified'})['stitle']}\\n")
PY
  """
}

// ==========================================
// CONSENSUS LOGIK
// ==========================================

process CONSENSUS_RESULTS {
  container = null
  publishDir "${params.out_dir}/summary", mode: 'copy', overwrite: true
  input: 
  path(cox_tables)
  path(its_tables)
  
  output: 
  tuple path("multilocus_species_counts.csv"), path("multilocus_abundance_matrix.csv")

  script:
  """
  python - << 'PY'
import csv
import re
import glob

data = {}
all_samples = set()
all_consensus_names = set()

# Der robuste Namens-Reiniger aus v1.1.2
def clean_name(raw_title):
    if "Unclassified" in raw_title: return "Unclassified"
    raw_title = raw_title.replace(',', ' ').replace(';', ' ').replace('UNVERIFIED:', '').strip()
    parts = raw_title.split()
    while parts and (any(char.isdigit() for char in parts[0]) or parts[0].count('.') > 0):
        parts.pop(0)
    if not parts: return "Unknown"
    cleaned_title = " ".join(parts)
    match = re.search(r'^([A-Z][a-z]+\\s+[a-z0-9\\-]+(\\s+(f\\.\\s*sp\\.|f\\.sp\\.|var\\.|subsp\\.|f\\.)\\s+[a-z0-9\\-]+)?)', cleaned_title)
    if match: return match.group(1).strip()
    return f"{parts[0]} {parts[1]}" if len(parts) >= 2 else parts[0]

# Datenstruktur: dict[sample][cluster_id] = {'count': x, 'cox': name, 'its': name}
parsed_data = {}

# Lese COX Tabellen
for f in glob.glob("*.cox.taxonomy.tsv"):
    sample = f.replace(".cox.taxonomy.tsv", "")
    all_samples.add(sample)
    if sample not in parsed_data: parsed_data[sample] = {}
    with open(f, 'r') as file:
        reader = csv.DictReader(file, delimiter='\\t')
        for row in reader:
            cid = row['cluster_id']
            if cid not in parsed_data[sample]: parsed_data[sample][cid] = {'count': int(row['read_count']), 'cox': 'Unclassified', 'its': 'Unclassified'}
            parsed_data[sample][cid]['cox'] = clean_name(row['best_hit_title'])

# Lese ITS Tabellen
for f in glob.glob("*.its.taxonomy.tsv"):
    sample = f.replace(".its.taxonomy.tsv", "")
    all_samples.add(sample)
    if sample not in parsed_data: parsed_data[sample] = {}
    with open(f, 'r') as file:
        reader = csv.DictReader(file, delimiter='\\t')
        for row in reader:
            cid = row['cluster_id']
            # Falls der Cluster nur in ITS existiert (sollte durch Pipeline-Design nicht vorkommen, aber zur Sicherheit)
            if cid not in parsed_data[sample]: parsed_data[sample][cid] = {'count': int(row['read_count']), 'cox': 'Unclassified', 'its': 'Unclassified'}
            parsed_data[sample][cid]['its'] = clean_name(row['best_hit_title'])

# Consensus bilden
for sample, clusters in parsed_data.items():
    for cid, info in clusters.items():
        cox_name = info['cox']
        its_name = info['its']
        count = info['count']
        
        consensus = ""
        if cox_name == "Unclassified" and its_name == "Unclassified":
            consensus = "Unclassified"
        elif cox_name != "Unclassified" and its_name == "Unclassified":
            consensus = f"{cox_name} (COX only)"
        elif cox_name == "Unclassified" and its_name != "Unclassified":
            consensus = f"{its_name} (ITS only)"
        elif cox_name == its_name:
            consensus = f"{cox_name} (Confirmed)"
        else:
            consensus = f"Conflict: COX={cox_name} | ITS={its_name}"
            
        all_consensus_names.add(consensus)
        if consensus not in data: data[consensus] = {}
        data[consensus][sample] = data[consensus].get(sample, 0) + count

sorted_samples = sorted(list(all_samples))
sorted_consensus = sorted(list(all_consensus_names))

# Tabellen schreiben
with open("multilocus_species_counts.csv", 'w') as out:
    header = ['consensus_species'] + sorted_samples + ['total']
    out.write(','.join(header) + '\\n')
    for sp in sorted_consensus:
        total = sum(data[sp].values())
        row = [f'"{sp}"'] + [str(data[sp].get(sa, 0)) for sa in sorted_samples] + [str(total)]
        out.write(','.join(row) + '\\n')

with open("multilocus_abundance_matrix.csv", 'w') as out:
    out.write(','.join(['consensus_species'] + sorted_samples) + '\\n')
    for sp in sorted_consensus:
        out.write(','.join([f'"{sp}"'] + [str(data[sp].get(sa, 0)) for sa in sorted_samples]) + '\\n')
PY
  """
}

process REPORT_HTML {
  container = null
  publishDir "${params.out_dir}", mode: 'copy', overwrite: true
  input: tuple path(metagenomics_csv), path(matrix_csv)
  output: path("wf-metabarcoding-report.html")
  script:
  """
  python - << 'PY'
from pathlib import Path
html = f'''<!doctype html><html><head><meta charset="utf-8"/><title>Multi-Locus Report</title>
<style>body{{font-family:sans-serif;margin:30px;}}.box{{border:1px solid #ddd;padding:20px;background:#fafafa;}}a{{color:#0366d6;text-decoration:none;font-weight:bold;}}</style>
</head><body><h1>Multi-Locus Metabarcoding Report (COX + ITS)</h1><div class="box"><h3>Ergebnisse</h3><ul>
<li><a href="summary/{Path("${metagenomics_csv}").name}" target="_blank">Consensus Counts</a></li>
<li><a href="summary/{Path("${matrix_csv}").name}" target="_blank">Abundance Matrix</a></li>
</ul></div></body></html>'''
with open("wf-metabarcoding-report.html", "w", encoding="utf-8") as f: f.write(html)
PY
  """
}
