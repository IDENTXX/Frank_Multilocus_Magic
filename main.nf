nextflow.enable.dsl=2

def timestamp = new java.util.Date().format( 'yyyy-MM-dd_HH-mm' )

workflow {
  def db_catalog = [
    'oomycetes_cox' : '/mnt/d/Epi2Me_Datenbanken/oomycetes_cox_ref.fasta',
    'fusarium_tef'  : '/mnt/d/Epi2Me_Datenbanken/fusarium_TEF_ref.fasta',
    'oomycetes_its' : '/mnt/d/Epi2Me_Datenbanken/oomycetes_its_ref.fasta'
  ]

  def resolved_db_path = params.custom_db_path
  if (params.db_choice != 'custom_fasta') {
      resolved_db_path = db_catalog[params.db_choice]
  }

  if( !params.reads ) {
    log.info "Bitte Reads angeben."
    return
  }

  ch_samples = Channel
    .fromPath("${params.reads}/*/*.fastq.gz")
    .map { f -> tuple(f.parent.name, f) }
    .groupTuple()
    .map { sample, files -> tuple(sample, files.sort()) }

  merged    = MERGE_FASTQ(ch_samples)
  fasta     = FASTQ_TO_FASTA(merged)
  clustered = CLUSTER_VSEARCH(fasta)
  kept      = FILTER_CLUSTERS(clustered)

  ref_fa    = Channel.fromPath(resolved_db_path).first()
  dbdir_ch  = MAKEBLASTDB(ref_fa)

  kept_with_db = kept.combine(dbdir_ch)
  blasted      = BLASTN(kept_with_db)
  tax          = JOIN_COUNTS_BLAST(blasted)

  tax_tables_list = tax.map { sample, taxfile -> taxfile }.collect()
  summary = AGGREGATE_RESULTS(tax_tables_list)
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

process MAKEBLASTDB {
  input: path(db_fasta)
  output: path("blastdb")
  script:
  """
  mkdir -p blastdb
  cp ${db_fasta} blastdb/db.fasta
  makeblastdb -in blastdb/db.fasta -dbtype nucl -out blastdb/generic_db
  """
}

process BLASTN {
  tag "$sample"
  input: tuple val(sample), path(counts), path(centroids_fa), path(dbdir)
  output: tuple val(sample), path(counts), path("${sample}.blast.tsv")
  script:
  """
  if [ ! -s ${centroids_fa} ]; then touch ${sample}.blast.tsv; exit 0; fi
  blastn -query ${centroids_fa} -db ${dbdir}/generic_db -max_target_seqs 5 -num_threads ${task.cpus} -outfmt "6 qseqid sseqid pident length qlen bitscore evalue qcovs stitle" > ${sample}.blast.tsv
  """
}

process JOIN_COUNTS_BLAST {
  tag "$sample"
  container = null
  input: tuple val(sample), path(counts_tsv), path(blast_tsv)
  output: tuple val(sample), path("${sample}.taxonomy.tsv")
  script:
  """
  python - << 'PY'
import csv
min_id = float("${params.min_identity}")
min_cov = float("${params.min_coverage}")
hits = {}
try:
    with open("${blast_tsv}") as f:
        reader = csv.reader(f, delimiter='\\t')
        for row in reader:
            if not row: continue
            qid = row[0]
            if qid in hits: continue
            if float(row[2]) < min_id or float(row[7]) < min_cov: continue
            stitle = row[8] if len(row) > 8 else 'NA'
            hits[qid] = {'sseqid': row[1], 'pident': row[2], 'qcovs': row[7], 'stitle': stitle}
except Exception: pass

with open("${sample}.taxonomy.tsv", 'w') as out:
    out.write("cluster_id\\tread_count\\tbest_sseqid\\tpident\\tqcovs\\tbest_hit_title\\n")
    with open("${counts_tsv}") as f:
        reader = csv.reader(f, delimiter='\\t')
        for row in reader:
            if not row: continue
            h = hits.get(row[0], {'sseqid':'NA', 'pident':'0', 'qcovs':'0', 'stitle':'Unclassified'})
            out.write(f"{row[0]}\\t{row[1]}\\t{h['sseqid']}\\t{h['pident']}\\t{h['qcovs']}\\t{h['stitle']}\\n")
PY
  """
}

process AGGREGATE_RESULTS {
  container = null
  publishDir "${params.out_dir}/summary", mode: 'copy', overwrite: true
  input: path(tables)
  output: tuple path("wf-metagenomics-counts-species.csv"), path("abundance_matrix.csv")
  script:
  """
  printf "sample\\tcluster_id\\tread_count\\tbest_sseqid\\tpident\\tqcovs\\tbest_hit_title\\n" > raw_combined.tsv
  for f in ${tables}; do
    fname=\$(basename \$f)
    s=\${fname%%.taxonomy.tsv}
    awk -v s=\$s 'NR>1{print s"\\t"\$0}' \$f >> raw_combined.tsv
  done

  python - << 'PY'
import csv
import re

data, all_samples, all_species = {}, set(), set()

def clean_name(raw_title):
    if "Unclassified" in raw_title: return "Unclassified"
    raw_title = raw_title.replace(',', ' ').replace(';', ' ').replace('UNVERIFIED:', '').strip()
    
    # Accession am Anfang aggressiver entfernen (JQ429332.1 etc.)
    parts = raw_title.split()
    while parts and (any(char.isdigit() for char in parts[0]) or parts[0].count('.') > 0):
        parts.pop(0)
    
    if not parts: return "Unknown"
    cleaned_title = " ".join(parts)
    
    # Finaler Regex für lycopersici und Varietäten
    match = re.search(r'^([A-Z][a-z]+\\s+[a-z0-9\\-]+(\\s+(f\\.\\s*sp\\.|f\\.sp\\.|var\\.|subsp\\.|f\\.)\\s+[a-z0-9\\-]+)?)', cleaned_title)
    
    if match:
        return match.group(1).strip()
    return f"{parts[0]} {parts[1]}" if len(parts) >= 2 else parts[0]

with open('raw_combined.tsv', 'r') as f:
    for row in csv.DictReader(f, delimiter='\\t'):
        sp = clean_name(row['best_hit_title'])
        sa = row['sample']
        all_samples.add(sa)
        all_species.add(sp)
        if sp not in data: data[sp] = {}
        data[sp][sa] = data[sp].get(sa, 0) + int(row['read_count'] if row['read_count'].isdigit() else 0)

sorted_samples = sorted(list(all_samples))
sorted_species = sorted(list(all_species))

with open('wf-metagenomics-counts-species.csv', 'w') as out:
    header = ['species'] + sorted_samples + ['total', 'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'tax']
    out.write(','.join(header) + '\\n')
    for sp in sorted_species:
        total = sum(data[sp].values())
        genus = sp.split()[0] if sp != "Unclassified" else "NA"
        row = [sp] + [str(data[sp].get(sa, 0)) for sa in sorted_samples] + [str(total), "NA","NA","NA","NA","NA","NA", genus, f"NA;NA;NA;NA;NA;NA;{genus};{sp}"]
        out.write(','.join(row) + '\\n')

with open('abundance_matrix.csv', 'w') as out:
    out.write(','.join(['species'] + sorted_samples) + '\\n')
    for sp in sorted_species:
        out.write(','.join([sp] + [str(data[sp].get(sa, 0)) for sa in sorted_samples]) + '\\n')
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
html = f'''<!doctype html><html><head><meta charset="utf-8"/><title>Metabarcoding Report</title>
<style>body{{font-family:sans-serif;margin:30px;}}.box{{border:1px solid #ddd;padding:20px;background:#fafafa;}}a{{color:#0366d6;text-decoration:none;font-weight:bold;}}</style>
</head><body><h1>Generic Metabarcoding Report</h1><div class="box"><h3>Ergebnisse</h3><ul>
<li><a href="summary/{Path("${metagenomics_csv}").name}" target="_blank">Metagenomics Counts</a></li>
<li><a href="summary/{Path("${matrix_csv}").name}" target="_blank">Abundance Matrix</a></li>
</ul></div></body></html>'''
with open("wf-metabarcoding-report.html", "w", encoding="utf-8") as f: f.write(html)
PY
  """
}