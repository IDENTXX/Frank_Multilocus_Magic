nextflow.enable.dsl=2

def timestamp = new java.util.Date().format( 'yyyy-MM-dd_HH-mm' )

workflow {
  // Der dynamische Katalog für Multi-Locus Kombinationen
  def db_catalog = [
    'oomycetes_cox_its' : [
        path1: '/mnt/d/Epi2Me_Datenbanken/oomycetes_cox_ref.fasta', name1: 'COX',
        path2: '/mnt/d/Epi2Me_Datenbanken/oomycetes_its_ref.fasta', name2: 'ITS'
    ],
    'fusarium_tef_its' : [
        path1: '/mnt/d/Epi2Me_Datenbanken/fusarium_TEF_ref.fasta',  name1: 'TEF1',
        path2: '/mnt/d/Epi2Me_Datenbanken/oomycetes_its_ref.fasta', name2: 'ITS' 
    ]
  ]

  def choice = db_catalog.containsKey(params.db_choice) ? params.db_choice : 'oomycetes_cox_its'
  def selected_db = db_catalog[choice]

  if( !params.reads ) {
    log.info "Bitte Reads angeben."
    return
  }

  // 1. Pre-Processing
  ch_samples = Channel
    .fromPath("${params.reads}/*/*.fastq.gz")
    .map { f -> tuple(f.parent.name, f) }
    .groupTuple()
    .map { sample, files -> tuple(sample, files.sort()) }

  merged    = MERGE_FASTQ(ch_samples)
  fasta     = FASTQ_TO_FASTA(merged)
  clustered = CLUSTER_VSEARCH(fasta)
  kept      = FILTER_CLUSTERS(clustered)

  // 2. Datenbanken aufbauen 
  db_m1_ch = MAKEBLASTDB_M1(Channel.fromPath(selected_db.path1).first())
  db_m2_ch = MAKEBLASTDB_M2(Channel.fromPath(selected_db.path2).first())

  // 3. BLAST parallel ausführen
  kept_m1 = kept.combine(db_m1_ch)
  kept_m2 = kept.combine(db_m2_ch)

  blast_m1 = BLASTN_M1(kept_m1)
  blast_m2 = BLASTN_M2(kept_m2)

  // 4. Taxonomy-Tabellen pro Locus erstellen (jetzt mit pident)
  tax_m1 = JOIN_COUNTS_BLAST_M1(blast_m1)
  tax_m2 = JOIN_COUNTS_BLAST_M2(blast_m2)

  // 5. Sammeln und Consensus bilden
  tax_m1_list = tax_m1.map { sample, taxfile -> taxfile }.collect()
  tax_m2_list = tax_m2.map { sample, taxfile -> taxfile }.collect()

  summary = CONSENSUS_RESULTS(tax_m1_list, tax_m2_list, selected_db.name1, selected_db.name2)
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

process MAKEBLASTDB_M1 {
  input: path(db_fasta)
  output: path("blastdb_m1")
  script: "mkdir -p blastdb_m1 && cp ${db_fasta} blastdb_m1/db.fasta && makeblastdb -in blastdb_m1/db.fasta -dbtype nucl -out blastdb_m1/m1_db"
}

process MAKEBLASTDB_M2 {
  input: path(db_fasta)
  output: path("blastdb_m2")
  script: "mkdir -p blastdb_m2 && cp ${db_fasta} blastdb_m2/db.fasta && makeblastdb -in blastdb_m2/db.fasta -dbtype nucl -out blastdb_m2/m2_db"
}

process BLASTN_M1 {
  tag "$sample"
  input: tuple val(sample), path(counts), path(centroids_fa), path(dbdir)
  output: tuple val(sample), path(counts), path("${sample}.m1.blast.tsv")
  script: "if [ ! -s ${centroids_fa} ]; then touch ${sample}.m1.blast.tsv; exit 0; fi; blastn -query ${centroids_fa} -db ${dbdir}/m1_db -max_target_seqs 5 -num_threads ${task.cpus} -outfmt '6 qseqid sseqid pident length qlen bitscore evalue qcovs stitle' > ${sample}.m1.blast.tsv"
}

process BLASTN_M2 {
  tag "$sample"
  input: tuple val(sample), path(counts), path(centroids_fa), path(dbdir)
  output: tuple val(sample), path(counts), path("${sample}.m2.blast.tsv")
  script: "if [ ! -s ${centroids_fa} ]; then touch ${sample}.m2.blast.tsv; exit 0; fi; blastn -query ${centroids_fa} -db ${dbdir}/m2_db -max_target_seqs 5 -num_threads ${task.cpus} -outfmt '6 qseqid sseqid pident length qlen bitscore evalue qcovs stitle' > ${sample}.m2.blast.tsv"
}

process JOIN_COUNTS_BLAST_M1 {
  tag "$sample"
  container = null
  input: tuple val(sample), path(counts_tsv), path(blast_tsv)
  output: tuple val(sample), path("${sample}.m1.taxonomy.tsv")
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
            # Speichere stitle UND pident (Homologie in %)
            hits[row[0]] = {'stitle': row[8] if len(row) > 8 else 'NA', 'pident': float(row[2])}
except Exception: pass

with open("${sample}.m1.taxonomy.tsv", 'w') as out:
    out.write("cluster_id\\tread_count\\tbest_hit_title\\tpident\\n")
    with open("${counts_tsv}") as f:
        for row in csv.reader(f, delimiter='\\t'):
            if row: 
                hit = hits.get(row[0], {'stitle':'Unclassified', 'pident': 0.0})
                out.write(f"{row[0]}\\t{row[1]}\\t{hit['stitle']}\\t{hit['pident']}\\n")
PY
  """
}

process JOIN_COUNTS_BLAST_M2 {
  tag "$sample"
  container = null
  input: tuple val(sample), path(counts_tsv), path(blast_tsv)
  output: tuple val(sample), path("${sample}.m2.taxonomy.tsv")
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
            hits[row[0]] = {'stitle': row[8] if len(row) > 8 else 'NA', 'pident': float(row[2])}
except Exception: pass

with open("${sample}.m2.taxonomy.tsv", 'w') as out:
    out.write("cluster_id\\tread_count\\tbest_hit_title\\tpident\\n")
    with open("${counts_tsv}") as f:
        for row in csv.reader(f, delimiter='\\t'):
            if row: 
                hit = hits.get(row[0], {'stitle':'Unclassified', 'pident': 0.0})
                out.write(f"{row[0]}\\t{row[1]}\\t{hit['stitle']}\\t{hit['pident']}\\n")
PY
  """
}

process CONSENSUS_RESULTS {
  container = null
  publishDir "${params.out_dir}/summary", mode: 'copy', overwrite: true
  input: 
  path(m1_tables)
  path(m2_tables)
  val(name1)
  val(name2)
  
  output: 
  tuple path("multilocus_species_counts.csv"), path("multilocus_abundance_matrix.csv"), path("homology_qc_report.csv")

  script:
  """
  python - << 'PY'
import csv
import re
import glob

data = {}
all_samples = set()
all_consensus_names = set()

def clean_name(raw_title):
    if "Unclassified" in raw_title: return "Unclassified"
    raw_title = raw_title.replace(',', ' ').replace(';', ' ').replace('UNVERIFIED:', '').replace('UNVERIFIED_ORG:', '').strip()
    parts = raw_title.split()
    while parts and (any(char.isdigit() for char in parts[0]) or parts[0].count('.') > 0):
        parts.pop(0)
    if not parts: return "Unknown"
    cleaned_title = " ".join(parts)
    match = re.search(r'^([A-Z][a-z]+\\s+[a-z0-9\\-]+(\\s+(f\\.\\s*sp\\.|f\\.sp\\.|var\\.|subsp\\.|f\\.)\\s+[a-z0-9\\-]+)?)', cleaned_title)
    if match: return match.group(1).strip()
    return f"{parts[0]} {parts[1]}" if len(parts) >= 2 else parts[0]

parsed_data = {}

# Lese Marker 1 Tabellen
for f in glob.glob("*.m1.taxonomy.tsv"):
    sample = f.replace(".m1.taxonomy.tsv", "")
    all_samples.add(sample)
    if sample not in parsed_data: parsed_data[sample] = {}
    with open(f, 'r') as file:
        reader = csv.DictReader(file, delimiter='\\t')
        for row in reader:
            cid = row['cluster_id']
            if cid not in parsed_data[sample]: parsed_data[sample][cid] = {'count': int(row['read_count']), 'm1': 'Unclassified', 'm2': 'Unclassified', 'm1_pid': 0.0, 'm2_pid': 0.0}
            parsed_data[sample][cid]['m1'] = clean_name(row['best_hit_title'])
            parsed_data[sample][cid]['m1_pid'] = float(row['pident'])

# Lese Marker 2 Tabellen
for f in glob.glob("*.m2.taxonomy.tsv"):
    sample = f.replace(".m2.taxonomy.tsv", "")
    all_samples.add(sample)
    if sample not in parsed_data: parsed_data[sample] = {}
    with open(f, 'r') as file:
        reader = csv.DictReader(file, delimiter='\\t')
        for row in reader:
            cid = row['cluster_id']
            if cid not in parsed_data[sample]: parsed_data[sample][cid] = {'count': int(row['read_count']), 'm1': 'Unclassified', 'm2': 'Unclassified', 'm1_pid': 0.0, 'm2_pid': 0.0}
            parsed_data[sample][cid]['m2'] = clean_name(row['best_hit_title'])
            parsed_data[sample][cid]['m2_pid'] = float(row['pident'])

homology_stats = {}

# Consensus mit dynamischen Namen bilden & Homologie sammeln
for sample, clusters in parsed_data.items():
    for cid, info in clusters.items():
        m1_name = info['m1']
        m2_name = info['m2']
        count = info['count']
        
        consensus = ""
        if m1_name == "Unclassified" and m2_name == "Unclassified":
            consensus = "Unclassified"
        elif m1_name != "Unclassified" and m2_name == "Unclassified":
            consensus = f"{m1_name} (${name1} only)"
        elif m1_name == "Unclassified" and m2_name != "Unclassified":
            consensus = f"{m2_name} (${name2} only)"
        elif m1_name == m2_name:
            consensus = f"{m1_name} (Confirmed)"
        else:
            consensus = f"Conflict: ${name1}={m1_name} | ${name2}={m2_name}"
            
        all_consensus_names.add(consensus)
        if consensus not in data: data[consensus] = {}
        data[consensus][sample] = data[consensus].get(sample, 0) + count

        # Homologie Werte speichern
        if consensus not in homology_stats:
            homology_stats[consensus] = {'m1_pids': [], 'm2_pids': []}
        
        if m1_name != "Unclassified":
            homology_stats[consensus]['m1_pids'].append((info['m1_pid'], count))
        if m2_name != "Unclassified":
            homology_stats[consensus]['m2_pids'].append((info['m2_pid'], count))

sorted_samples = sorted(list(all_samples))
sorted_consensus = sorted(list(all_consensus_names))

# Output 1: Counts
with open("multilocus_species_counts.csv", 'w') as out:
    header = ['consensus_species'] + sorted_samples + ['total']
    out.write(','.join(header) + '\\n')
    for sp in sorted_consensus:
        total = sum(data[sp].values())
        row = [f'"{sp}"'] + [str(data[sp].get(sa, 0)) for sa in sorted_samples] + [str(total)]
        out.write(','.join(row) + '\\n')

# Output 2: Matrix
with open("multilocus_abundance_matrix.csv", 'w') as out:
    out.write(','.join(['consensus_species'] + sorted_samples) + '\\n')
    for sp in sorted_consensus:
        out.write(','.join([f'"{sp}"'] + [str(data[sp].get(sa, 0)) for sa in sorted_samples]) + '\\n')

# Output 3: NEU - Homologie QC Report
with open("homology_qc_report.csv", 'w') as out_qc:
    out_qc.write(f"consensus_species,${name1}_mean_id_%,${name1}_min_%,${name1}_max_%,${name2}_mean_id_%,${name2}_min_%,${name2}_max_%\\n")
    for sp in sorted_consensus:
        m1_data = homology_stats.get(sp, {}).get('m1_pids', [])
        m2_data = homology_stats.get(sp, {}).get('m2_pids', [])

        def get_stats(data_list):
            if not data_list: return "NA", "NA", "NA"
            total_reads = sum(c for p, c in data_list)
            if total_reads == 0: return "NA", "NA", "NA"
            # Gewichteter Durchschnitt basierend auf Reads pro Cluster
            mean_val = sum(p * c for p, c in data_list) / total_reads
            min_val = min(p for p, c in data_list)
            max_val = max(p for p, c in data_list)
            return f"{mean_val:.2f}", f"{min_val:.2f}", f"{max_val:.2f}"

        m1_mean, m1_min, m1_max = get_stats(m1_data)
        m2_mean, m2_min, m2_max = get_stats(m2_data)

        out_qc.write(f'"{sp}",{m1_mean},{m1_min},{m1_max},{m2_mean},{m2_min},{m2_max}\\n')
PY
  """
}

process REPORT_HTML {
  container = null
  publishDir "${params.out_dir}", mode: 'copy', overwrite: true
  input: tuple path(metagenomics_csv), path(matrix_csv), path(qc_csv)
  output: path("wf-metabarcoding-report.html")
  script:
  """
  python - << 'PY'
from pathlib import Path
html = f'''<!doctype html><html><head><meta charset="utf-8"/><title>Multi-Locus Report</title>
<style>body{{font-family:sans-serif;margin:30px;}}.box{{border:1px solid #ddd;padding:20px;background:#fafafa;}}a{{color:#0366d6;text-decoration:none;font-weight:bold;}}</style>
</head><body><h1>Multi-Locus Metabarcoding Report</h1><div class="box"><h3>Ergebnisse</h3><ul>
<li><a href="summary/{Path("${metagenomics_csv}").name}" target="_blank">Consensus Counts</a></li>
<li><a href="summary/{Path("${matrix_csv}").name}" target="_blank">Abundance Matrix</a></li>
<li><a href="summary/{Path("${qc_csv}").name}" target="_blank">Homology QC Report</a></li>
</ul></div></body></html>'''
with open("wf-metabarcoding-report.html", "w", encoding="utf-8") as f: f.write(html)
PY
  """
}